"""Create a subset VDS from a larger VDS

This script subsets a given VDS by either a number of samples, a region or set of regions, or both.

Arguments:
    --vds-path: str. Path to the VDS in GCS. Does not need the project. e.g cpg-bioheart-test/vds/bioheart1-0.vds should be entered as vds/bioheart1-0.vds.
    --n-samples: int, optional. The number of samples to subset the VDS down to.
    --intervals: str, optional. A (comma separated) string in the format 'chr:start-end' of the interval to subset to, or the path to a file in gcs with one interval per line.
    --output-formats: str. A space separated string of output formats to generate. Only formats in [vcf, bed, vds] can be chosen.
    --random-seed: int, optional. An int to control the random number generator (default: 19700101).

Raises:
    ValueError: only supports values of n-samples that are less than the total number of samples in the input VDS.
    FileExistsError: Will not overwrite existing files.
    AttributeError: Unrecognised arguments.
    AttributeError: The user must provide one of n-samples or intervals, otherwise the script will not do any subsetting.

Example usage:

analysis-runner --dataset bioheart \
    --access-level test \
    --output-dir test-subset \
    --description 'Test VDS subsetting script' \
    python3 subset_vds.py --vds-path vds/bioheart1-0.vds \
    --output-formats vcf bed \
    --intervals chr22,chr1:50000-100000 \
    --n-samples 2
"""

from argparse import ArgumentParser, Namespace
from random import sample, seed
from re import split
from typing import Any

from cpg_utils.config import dataset_path, output_path, to_path
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
from hail.utils.java import FatalError
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
        if format == "vds":
            if to_path(
                output_path(f"{infile_name}_subset", category="analysis")
            ).exists():
                output_errors += f"The output VDS {infile_name}_subset.vds already exists. Refusing to overwrite it.\n"
                files_exist_errors = True
            if to_path(
                output_path("subset_samples_file.txt", category="analysis")
            ).exists():
                output_errors += f"The file {to_path(output_path('subset_samples_file.txt'))} exists. Refusing to overwrite it."
                files_exist_errors = True
        if (
            format == "bed"
            and to_path(
                output_path(f"{infile_name}_subset.bed", category="analysis")
            ).exists()
        ):
            output_errors += f"The output {to_path(output_path(infile_name))}_subset.bed fileset exists. Refusing to overwrite it.\n"
            files_exist_errors = True
        if (
            format == "vcf"
            and to_path(
                output_path(f"{infile_name}_subset.vcf.bgz", category="analysis")
            ).exists()
        ):
            output_errors += f"The output file {to_path(output_path(infile_name))}_subset.vcf.bgz exists. Refusing to overwrite it.\n"
            files_exist_errors = True
    if files_exist_errors:
        raise FileExistsError(f"{output_errors}")


def parse_intervals(intervals: str) -> list[str]:
    """read intervals from either the command line or a file

    Args:
        intervals (str): a string of either comma separated intervals, or a path to a file of intervals (one per line) in gcs

    Returns:
        list[str]: a list of strings representing intervals
    """
    interval_list: list[str] = []
    if to_path(intervals).exists():
        with open(to_path(intervals), "rt") as interval_file:
            interval_list = interval_file.readlines()
    else:
        interval_list = intervals.split(",")
    return interval_list


def convert_intervals_to_locus(interval_list: list[str]) -> list[IntervalExpression]:
    """convert a list of string intervals to interval expressions that hail understands for locus filtering

    Args:
        interval_list (list[str]): a list of intervals to convert to locuses

    Raises:
        ValueError: intervals need to be in the format of either 'chrom', 'chrom:pos', 'chrom:start-end' or chrom:pos-chrom:pos. Each contig also needs to be in the list of reference contigs for GRCh38

    Returns:
        list[IntervalExpression]: the list of locuses in a format that hail can understand
    """
    locus_list: list[IntervalExpression] = []
    hail_reference_contigs: list[Any] = get_reference("GRCh38").contigs
    invalid_locuses: list[str] = []
    locus_errors: str = ""
    for interval in interval_list:
        split_interval: list[str] = split("[:-]", interval)
        n_interval_elements: int = len(split_interval)
        if n_interval_elements < 4 and split_interval[0] not in hail_reference_contigs:
            invalid_locuses.append(interval)
            locus_errors += f"The contig: {split_interval[0]} is not in the list of reference contigs: {hail_reference_contigs}."
        elif n_interval_elements == 4 and not all(
            [
                contig in hail_reference_contigs
                for contig in [split_interval[0], split_interval[2]]
            ]
        ):
            invalid_locuses.append(interval)
            locus_errors += f"Not all input contigs {[split_interval[0], split_interval[2]]} are in the list of reference contigs {hail_reference_contigs}.\n"
        elif n_interval_elements > 4:
            invalid_locuses.append(interval)
            locus_errors += f"The interval {interval} does not conform to any acceptable input format.\n See https://hail.is/docs/0.2/functions/genetics.html#hail.expr.functions.parse_locus_interval for the acceptable formats.\n"
        try:
            locus_list.append(parse_locus_interval(interval, reference_genome="GRCh38"))
        except FatalError as e:
            invalid_locuses.append(interval)
            locus_errors += f"{e}\n"

    if len(invalid_locuses) > 0:
        raise ValueError(
            f"The following locus(es) contain errors: {invalid_locuses}. For more information, read the error message(s) \n{locus_errors}"
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
        subset_sample_list (list[str]): the list of samples to retain in the vds

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
        subset_vds.write(output_path(f"{infile_name}_subset.vds", category="analysis"))

    subset_dense_mt: MatrixTable | Table | Any = to_dense_mt(subset_vds)
    subset_dense_mt = split_multi_hts(subset_dense_mt)

    for format in output_formats:
        if format == "bed":
            export_plink(
                subset_dense_mt,
                output_path(f"{infile_name}_subset", category="analysis"),
                call=subset_dense_mt.LGT,
                ind_id=subset_dense_mt.s,
            )
        if format == "vcf":
            if "gvcf_info" in subset_dense_mt.entry:
                subset_dense_mt = subset_dense_mt.drop("gvcf_info")
            export_vcf(
                subset_dense_mt,
                output_path(f"{infile_name}_subset.vcf.bgz", category="analysis"),
                tabix=True,
            )

    if subset_sample_list:
        with to_path(output_path("subset_samples_file.txt", category="analysis")).open(
            "wt"
        ) as outfile:
            outfile.write("\n".join(subset_sample_list))


def main(
    vds_path: str,
    n_samples: int | None,
    intervals: str | None,
    output_formats: list[str],
    random_seed: int,
) -> None:
    seed(random_seed)
    init_batch()
    input_vds: VariantDataset = read_vds(dataset_path(vds_path))
    infile_name: str = vds_path.split("/")[-1].split(".")[0]
    subset_vds: VariantDataset | None = None
    subset_sample_list: list[str] | None = None
    parsed_intervals: list[str]
    parsed_locus: list[IntervalExpression]

    check_output_already_exists(output_formats, infile_name)

    # Always subset by interval first, if possible
    # https://discuss.hail.is/t/filtering-samples-from-vds-in-google-cloud/3718/6
    if n_samples:
        subset_sample_list = get_subset_sample_list(input_vds, n_samples)
        if intervals:
            parsed_intervals = parse_intervals(intervals)
            parsed_locus = convert_intervals_to_locus(parsed_intervals)
            subset_vds = subset_by_samples(
                subset_by_locus(parsed_locus, input_vds), subset_sample_list
            )
        else:
            subset_vds = subset_by_samples(input_vds, subset_sample_list)
    elif intervals:
        parsed_intervals = parse_intervals(intervals)
        parsed_locus = convert_intervals_to_locus(parsed_intervals)
        subset_vds = subset_by_locus(parsed_locus, input_vds)

    write_outputs(output_formats, subset_vds, subset_sample_list, infile_name)


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
        help="Space separated string indicating what output formats you would like. One of [vcf, bed, vds]",
        required=True,
        nargs="+",
        choices=["vcf", "bed", "vds"],
    )
    parser.add_argument(
        "--random-seed",
        help="Seed for the random number generator when subsetting by samples",
        required=False,
        default=19700101,
        type=int,
    )
    args: Namespace = parser.parse_args()
    if not args.n_samples and not args.intervals:
        raise AttributeError(
            "If neither a number of samples or intervals are provided, this script will not do any subsetting!"
        )
    main(
        vds_path=args.vds_path,
        n_samples=args.n_samples,
        intervals=args.intervals,
        output_formats=args.output_formats,
        random_seed=args.random_seed,
    )
