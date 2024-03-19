import argparse
import os
import sys
import logging
import time
import tracegroomer.utils as ut
from tracegroomer.tidy import perform_type_prep


def prep_args() -> argparse.ArgumentParser:
    show_defaults = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(prog="python -m tracegroomer",
                                     formatter_class=show_defaults)

    parser.add_argument('-cf', '--config_file', type=str,
                        help="Configuration file given as absolute path")

    parser.add_argument('-lm', '--labeled_metabo_file', type=str, default=None,
                        help="Labeled metabolomics input file, absolute path")

    parser.add_argument('-tf', '--type_of_file', type=str, default=None,
                        help="One of the following: \
                        IsoCor_out_tsv|rule_tsv|VIBMEC_xlsx|generic_xlsx")

    parser.add_argument('--amountMaterial_path', type=str, default=None,
                        help="absolute path to the file having the amount \
                           of material (number of cells, tissue weight, etc) \
                           by sample, for the normalization")

    parser.add_argument("--alternative_div_amount_material",
                        action=argparse.BooleanOptionalAction, default=True,
                        help="When dividing values by the amount of  \
                        material, also multiplies by 'mean(amountMaterial)' \
                        to stay in abundance units")

    parser.add_argument("--div_isotopologues_by_amount_material",
                        action=argparse.BooleanOptionalAction,
                        default=True,
                        help="Apply normalization by the amount of material, \
                        at the level of isotopologue absolute values. \
                        After this, re-computes all derived metrics. \
                        If False, only total abundances are normalized")

    parser.add_argument("--use_internal_standard", default=None, type=str,
                        help='Internal Standard for performing the division: \
                        total_abundances/internal_standard, \
                        example: --use_internal_standard Myristic_acid_d27. \
                        By default is not performed')

    parser.add_argument("--remove_these_metabolites", default=None, type=str,
                        help="Absolute path to the file with columns:  \
                        compartment, metabolite; listing the metabolites to \
                        be completely excluded (you know what you are doing)")

    # for isotopologues and meanenrich_or_fracfontrib only:

    parser.add_argument("--isotopologues_preview",
                        action=argparse.BooleanOptionalAction, default=False,
                        help="Plot isotopologue values, as given")

    parser.add_argument("--isosprop_min_admitted", default=-0.5, type=float,
                        help="Metabolites whose isotopologues proportions are\
                        less or equal to this cutoff, are removed")

    parser.add_argument("--fractions_stomp_values",
                        action=argparse.BooleanOptionalAction, default=True,
                        help="Stomps fractional contributions (synonym: \
                        mean enrichment), and isotopologue proportions, \
                        to max 1.0 and min 0.0")  #meanenrich_or_fracfontrib

    # for total abundance only if VIB data
    parser.add_argument("--under_detection_limit_set_nan",
                        action=argparse.BooleanOptionalAction,
                        default=True, help="On VIB results. Any abundance \
                           < LOD (Limit Of Detection), is set as NaN")

    parser.add_argument("--subtract_blankavg",
                        action=argparse.BooleanOptionalAction, default=True,
                        help="On VIB results. From samples' abundances,  \
                          subtracts the average of the blanks")

    # output files extension
    parser.add_argument("-ox", "--output_files_extension", type=str,
                        default="csv", help="Extension for the output files, \
                        must be one of the following: csv|tsv|txt")

    return parser


def main() -> int:
    logger = logging.getLogger(__name__)
    logger = ut.reset_log_config(logger)
    logger.info("\nTime: {}".format(time.strftime('%Y%m%d-%H.%M.%S')))

    parser = prep_args()
    args = parser.parse_args()
    logger.info(
        f"Running TraceGroomer with the following parameters:")
    for x in vars(args).keys():
        logger.info(f"{x} = {vars(args)[x]} ")

    supported_types = ["IsoCor_out_tsv", "rule_tsv",
                       "generic_xlsx", "VIBMEC_xlsx"]
    assert args.type_of_file in supported_types, logger.critical(
        f"Error: type_of_file {args.type_of_file} unsupported or misspelled. "
        f"Supported types are: {supported_types}. ")

    confdict = ut.open_config_file(os.path.expanduser(args.config_file))
    expected_keys_confdict = ['metadata', 'abundances', 'mean_enrichment',
                              'isotopologue_proportions', 'isotopologues']

    for k in expected_keys_confdict:
        if k not in confdict.keys():
            logger.warning(f"{k} : missing in configuration file! ")
            confdict[k] = None

    out_path = os.path.expanduser(confdict['groom_out_path'])
    assert os.path.exists(out_path), logger.critical(
        f"ERROR. Directory {out_path} does not exist. Abort.")

    targetedMetabo_path = os.path.expanduser(args.labeled_metabo_file)
    assert os.path.isfile(targetedMetabo_path), logger.critical(
        f"ERROR. Labeled metabolomics file {targetedMetabo_path} does not "
        f"exist. Abort.")

    # metadata is expected in the out path, which is the data sub-folder:
    metadata_used_extension = None
    for extension in [".tsv", ".csv"]:
        if os.path.isfile(os.path.join(
                out_path, f"{confdict['metadata']}{extension}")):
            metadata_used_extension = extension
    assert metadata_used_extension is not None, logger.critical(
            f"ERROR. Metadata file '{confdict['metadata']}' must exist in the"
            f" 'out_path' (data subfolder) location: '{out_path}'. Abort.")

    perform_type_prep(args, confdict, metadata_used_extension,
                      targetedMetabo_path,  out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
