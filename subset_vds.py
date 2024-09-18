"""Create a subset VDS from a larger VDS

This script subsets a given VDS by either a number of samples, a region or set of regions, or both.

Arguments:
    --vds-path: str. Path to the VDS in GCS. Does not need the project. e.g cpg-bioheart-test/vds/bioheart1-0.vds should be entered as vds/bioheart1-0.vds
    --n-samples: int, optional. The number of samples to subset the VDS down to.
    --intervals: str, optional. A (comma separated) string in the format 'chr:start-end' of the interval to subset to, or the path to a file in gcs with one interval per line.
    --output-format: str. A comma separated string of output formats to generate. Only formats in [vcf, bed, vds] can be chosen.

Raises:
    ValueError: only supports values of n-samples that are less than the total number of sampels in the input VDS.
    FileExistsError: Will not overwrite existing files.
    AttributeError: Unrecognised arguments.
    AttributeError: The user must provide one of n-samples or intervals, optherwise the script will not do any subsetting.
"""

from argparse import ArgumentParser
from io import StringIO
from random import sample
from re import fullmatch, split
from typing import Any

from cpg_utils.config import dataset_path, get_access_level, output_path, to_path
from cpg_utils.hail_batch import init_batch
from hail import (
    IntervalExpression,
    MatrixTable,
    Table,
    export_plink,
    export_vcf,
    get_reference,
    parse_locus_interval,
    split_multi_hts,
)
from hail.vds import filter_intervals, filter_samples, read_vds, to_dense_mt
from hail.vds.variant_dataset import VariantDataset


def check_output_already_exists(output_format: list[str], infile_name: str) -> None:
    """Check for existence of output files

    Args:
        output_format (list[str]): list of output formats to generate files for
        infile_name (str): base name of the dataset

    Returns:
        Tuple[bool, str]: A tuple to indicate if any output files already exist, and which ones they are.
    """
    output_errors: str = ""
    files_exist_errors: bool = False
    for format in output_format:
        if format == "vds" and to_path(output_path(f"{infile_name}_subset")).exists():
            output_errors += f"The output VDS {infile_name}_subset already exists. Refusing to overwrite it.\n"
            files_exist_errors = True
        if (
            format == "bed"
            and to_path(output_path(f"{infile_name}_subset.bed")).exists()
        ):
            output_errors += f"The output {to_path(output_path(infile_name))}_subset.bed fileset exists. Refusing to overwrite it.\n"
            files_exist_errors = True
        if (
            format == "vcf"
            and to_path(output_path(f"{infile_name}_subset.vcf.bgz")).exists()
        ):
            output_errors += f"The output file {to_path(output_path(infile_name))}_subset.vcf.bgz exists. Refusing to overwrite it.\n"
            files_exist_errors = True
    if files_exist_errors:
        raise FileExistsError(f"{output_errors}")


def parse_intervals(intervals: str) -> list[str]:
    """read intervals from either the command line or a file

    Args:
        intervals (str): a string of either comma separated intervals, or a path to a file of intervals (one per line) in gcs

    Raises:
        ValueError: the intervals need to match the following pattern: ^[0-9cxymt]+ (does it start with some kind of contig identifier?)

    Returns:
        list[str]: a list of strings representing intervals
    """
    interval_list: list[str] = []
    problem_intervals: list[str] = []
    if to_path(intervals).exists():
        with open(to_path(intervals), "rt") as interval_file:
            interval_list = interval_file.readlines()
    else:
        interval_list = intervals.split(",")
    for interval in interval_list:
        if not fullmatch("^[0-9chrxymt]+", interval.split(":")[0].lower()):
            problem_intervals.append(interval)
    if len(problem_intervals) > 0:
        raise ValueError(
            f"The list of supplied intervals contains the following incorrectly formatted contigs: {[problem for problem in problem_intervals]}"
        )
    return interval_list


def convert_intervals_to_locus(interval_list: list[str]) -> list[IntervalExpression]:
    """convert a list of string intervals to interval expressions that hail understands for locus filtering

    Args:
        interval_list (list[str]): a list of intervals to convert to locuses

    Raises:
        ValueError: chromosome start positions have to be greter than 0
        ValueError: intervals need to be in the format of either 'chrom', 'chrom:pos' or 'chrom:start-end'.

    Returns:
        list[IntervalExpression]: the list of locuses in a format that hail can understand
    """
    locus_list: list[IntervalExpression] = []
    contig: str = ""
    start: str | int = ""
    end: str | int = ""
    for interval in interval_list:
        split_interval: list[str] = split("[:-]", interval)
        if len(split_interval) == 1:
            contig = split_interval[0]
            start = "start"
            end = "end"
        elif len(split_interval) == 2:
            contig = split_interval[0]
            assert int(
                start := split_interval[1]
            ), f"If only one position is specified, it must be an int, not {start}"
            start = int(split_interval[1])
            end = int(split_interval[1]) + 1
        elif len(split_interval) == 3:
            contig = split_interval[0]
            if start != "start":
                assert int(
                    start := split_interval[1]
                ), f"The start value: {start} couldn't be converted to an int"
                if int(start) < 1:
                    raise ValueError(
                        "Chromosome start positions must be greater than 0."
                    )
            start = int(split_interval[1])
            if end != "end":
                assert int(
                    end := split_interval[2]
                ), f"The end value: {end} couldn't be converted to an int."
                if int(end) > get_reference("GRCh38").lengths[contig]:
                    end = get_reference("GRCh38").lengths[contig]
            end = int(split_interval[2])
        else:
            raise ValueError(
                f"The interval: {interval} does not conform to the required format of either 'chrom', 'chrom:pos' or 'chrom:start-end'."
            )
        if start != "start" and end != "end":
            assert int(start) < int(end)
        locus_list.append(
            parse_locus_interval(f"{contig}:{start}-{end}", reference_genome="GRCh38")
        )
    return locus_list


def get_subset_sample_list(input_vds: VariantDataset, n_samples: int) -> list[str]:
    """generates a list of samples from the input vds to use for subsetting

    Args:
        input_vds (VariantDataset): the input vds
        n_samples (int): the number of samples to subset to

    Raises:
        ValueError: if n_samples is 0, throw an error as this would remove all data from the vds
        ValueError: if n_samples is >= to the number of samples in the input vds throw an error, as this would not do any filtering

    Returns:
        list[str]: a random list of samples of length n_samples from the input vds, to use in subsetting
    """
    starting_samples: list[str] = input_vds.variant_data.s.collect()
    n_total_samples: int = len(starting_samples)
    if n_samples == 0:
        raise ValueError(
            "Subsetting to 0 samples will result in an empty MatrixTable, so not doing that."
        )
    if n_samples >= n_total_samples:
        raise ValueError(
            f"The subset sample size {n_samples} is at least the same size as the total number of samples in the input {n_total_samples}. Please choose a smaller number."
        )
    subset_sample_list = sample(starting_samples, n_samples)
    return subset_sample_list


def subset_by_samples(
    input_vds: VariantDataset, subset_sample_list: list[str]
) -> VariantDataset:
    """subset a vds by samples

    Args:
        input_vds (VariantDataset): the input vds to filter on
        subset_sample_list (list[str]): the list of sampels to retain in the vds

    Returns:
        VariantDataset: the subset vds
    """
    subset_vds: VariantDataset
    subset_vds = filter_samples(input_vds, subset_sample_list)
    return subset_vds


def subset_by_locus(
    parsed_locus: list[IntervalExpression], input_vds: VariantDataset
) -> VariantDataset:
    """subset a vds by locus

    Args:
        parsed_locus (list[IntervalExpression]): a list of locuses to filter the vds down to
        input_vds (VariantDataset): the input vds to filter

    Raises:
        ValueError: throws an error if all rows are removed by filtering

    Returns:
        VariantDataset: the subset vds
    """
    subset_vds: VariantDataset
    subset_vds = filter_intervals(input_vds, parsed_locus)
    if subset_vds.variant_data.count_rows() == 0:
        raise ValueError(
            f"No rows remain after applying the following locus filters: {parsed_locus}"
        )
    return subset_vds


def write_outputs(
    output_formats: list[str],
    subset_vds: VariantDataset | None,
    subset_sample_list: list[str] | None,
    infile_name: str,
    is_test: bool,
) -> None:
    """Writes the vds out in the specified formats

    Args:
        output_formats (list[str]): list of output formats to write
        subset_vds (VariantDataset | None): the resulting vds after subsetting by samples and / or intervals
        subset_sample_list (list[str] | None): the samples left in the output data, when subsetting by samples
        infile_name (str): the name of the input vds, to append '_subset' to
        is_test (bool): controls whether to write to test or not

    Raises:
        FileExistsError: throws an error if any of the proposed output paths exist, as it will not overwrite them
    """
    if "vds" in output_formats:
        if to_path(output_path(f"{infile_name}_subset")).exists():
            raise FileExistsError(
                f'The file {to_path(output_path("{infile_name}_subset"))} exists. Refusing to overwrite it.'
            )
        subset_vds.write(output_path(f"{infile_name}_subset", test=is_test))

    subset_dense_mt: MatrixTable | Table | Any = to_dense_mt(subset_vds)
    subset_dense_mt = split_multi_hts(subset_dense_mt)

    for format in output_formats:
        if format == "bed":
            export_plink(
                subset_dense_mt,
                output_path(f"{infile_name}_subset", test=is_test),
                call=subset_dense_mt.LGT,
                ind_id=subset_dense_mt.s,
            )
        if format == "vcf":
            if "gvcf_info" in subset_dense_mt.entry:
                subset_dense_mt = subset_dense_mt.drop("gvcf_info")
            export_vcf(
                subset_dense_mt,
                output_path(f"{infile_name}_subset.vcf.bgz", test=is_test),
                tabix=True,
            )

    if subset_sample_list:
        if to_path(output_path("subset_samples_file.txt")).exists():
            raise FileExistsError(
                f"The file {to_path(output_path('subset_samples_file.txt'))} exists. Refusing to overwrite it."
            )
        outdata = StringIO()
        for single_sample in subset_sample_list:
            outdata.write(f"{single_sample}\n")
        with to_path(output_path("subset_samples_file.txt", test=is_test)).open(
            "wt"
        ) as outfile:
            outfile.write(outdata.getvalue())
            outdata.close()


def main(
    vds_path: str,
    n_samples: int | None,
    intervals: str | None,
    output_formats: list[str],
) -> None:
    allowed_output_formats: list[str] = ["vcf", "bed", "vds"]
    if len(output_formats) == 0:
        raise ValueError(
            f"A value for 'output_formats' needs to be supplied, one of: {allowed_output_formats}"
        )
    if not all([outformat in allowed_output_formats for outformat in output_formats]):
        raise NotImplementedError(
            f"This script only supports {allowed_output_formats} as output formats."
        )

    init_batch()
    input_vds: VariantDataset = read_vds(dataset_path(vds_path))
    infile_name: str = vds_path.split("/")[-1].split(".")[0]
    subset_vds: VariantDataset | None = None
    subset_sample_list: list[str] | None = None
    is_test: bool = True
    access: str = get_access_level()
    parsed_intervals: list[str]
    parsed_locus: list[IntervalExpression]

    check_output_already_exists(output_formats, infile_name)

    if access != "test":
        is_test = False

    if n_samples:
        subset_sample_list = get_subset_sample_list(input_vds, n_samples)
        if intervals:
            parsed_intervals = parse_intervals(intervals)
            parsed_locus = convert_intervals_to_locus(parsed_intervals)
            subset_vds = subset_by_locus(
                parsed_locus,
                subset_by_samples(input_vds, subset_sample_list),
            )
        else:
            subset_vds = subset_by_samples(input_vds, subset_sample_list)
    elif intervals and not n_samples:
        parsed_intervals = parse_intervals(intervals)
        parsed_locus = convert_intervals_to_locus(parsed_intervals)
        subset_vds = subset_by_locus(parsed_locus, input_vds)

    write_outputs(output_formats, subset_vds, subset_sample_list, infile_name, is_test)


if __name__ == "__main__":
    parser = ArgumentParser(description="VDS subsetting script")
    parser.add_argument("--vds-path", help="Path to VDS in GCP.", required=True)
    parser.add_argument(
        "--n-samples",
        help="Number of samples to subset the VDS down to.",
        required=False,
        type=int,
    )
    parser.add_argument(
        "--intervals",
        help="Interval(s) provided either on the command line in the format chr:start-end or in a text file, one per line to keep in the VDS.",
        required=False,
    )
    parser.add_argument(
        "--output-formats",
        help="Comma separated string indicating what output formats you would like. One of [vcf, bed, vds]",
        required=True,
    )
    args, failures = parser.parse_known_args()
    if failures:
        parser.print_help()
        raise AttributeError(f"Invalid arguments {failures}")
    if not args.n_samples and not args.interval:
        raise AttributeError(
            "If neither a number of samples or intervals are provided, this script will not do any subsetting!"
        )
    main(
        vds_path=args.vds_path,
        n_samples=args.n_samples,
        intervals=args.intervals,
        output_formats=args.output_formats.lower().split(","),
    )
