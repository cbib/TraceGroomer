import argparse
import os
import sys
import logging
import tracegroomer.utils as ut
from tracegroomer.tidy import perform_type_prep

def prep_args() -> argparse.ArgumentParser:
    show_defaults = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(prog="python -m tracegroomer",
                                     formatter_class=show_defaults)

    parser.add_argument('-cf', '--config_file', type=str,
                        help="Configuration file given as absolute path.")

    parser.add_argument('-lm', '--labeled_metabo_file', type=str, default=None,
                        help="Labeled metabolomics input file, absolute path.") # TODO absolute? data subfolder ?

    parser.add_argument('-tf', '--type_of_file', type=str, default=None,
                        help="One of the following: \
                        IsoCor_out_tsv|VIBMEC_xlsx|generic_xlsx")

    parser.add_argument("-o", "--out_path", type=str, default=None,
                        help="Output directory (the data sub-folder) in absolute path.")

    # for abundance
    parser.add_argument("--under_detection_limit_set_nan",
                        action=argparse.BooleanOptionalAction,
                        default=True, help="On VIB results. Any abundance inferior \
                        to LOD (Limit Of Detection) is set as NaN.")  # np.nan

    parser.add_argument("--auto_drop_metabolite_LOD_based",
                        action=argparse.BooleanOptionalAction,
                        default=True,
                        help="On VIB results.  By compartment, a metabolite is \
                        automatically rejected if all abundances are under \
                        the detection limit. Has effect in all tables.")

    parser.add_argument("--subtract_blankavg",
                        action=argparse.BooleanOptionalAction, default=True,
                        help="On VIB results. From samples' abundances, subtracts  \
                        the average of the blanks.")

    parser.add_argument('--amountMaterial_path', type=str, default=None,
                        help="absolute path to the file .csv having the amount \
                           of material (number of cells, tissue weight, etc) \
                           by sample, for the total abundances normalization.")

    parser.add_argument("--alternative_div_amount_material",
                        action=argparse.BooleanOptionalAction, default=False,
                        help="When dividing abundances by the amount of  \
                        material, also multiplies by 'mean(amountMaterial)' \
                        to stay in abundance units.")

    parser.add_argument("--use_internal_standard", default=None, type=str,
                        help='Internal Standard for performing the division: \
                        abundances/internal_standard, \
                        example: --use_internal_standard Myristic_acid_d27. \
                        By default is not performed.')

    # for isotopologues and meanenrich_or_fracfontrib:

    parser.add_argument("--isotopologues_preview",
                        action=argparse.BooleanOptionalAction, default=False,
                        help="Plot isotopologue values, as given")

    parser.add_argument("--isosprop_min_admitted", default=-0.5, type=float,
                        help="Metabolites whose isotopologues are less or equal \
                        this cutoff are removed.")

    parser.add_argument("--isosprop_stomp_values",
                        action=argparse.BooleanOptionalAction,
                        default=True,
                        help="Stomps isotopologues' proportions \
                        to max 1.0 and min 0.0")

    parser.add_argument("--meanenrich_or_fracfontrib_stomp_values",
                        action=argparse.BooleanOptionalAction, default=True,
                        help="Stomps fractional contributions (synonym: \
                        mean enrichment) to max 1.0 and min 0.0")

    # for all
    parser.add_argument("--remove_these_metabolites", default=None, type=str,
                        help="Absolute path to the .csv file with columns:  \
                        compartment, metabolite; listing the metabolites to \
                        be completely excluded (you know what you are doing)")
    # all tables affected

    return parser


def main() -> int:
    logging.basicConfig(filename=None, encoding='utf-8', level=logging.DEBUG)
    logger = logging.getLogger()
    parser = prep_args()
    args = parser.parse_args()

    out_path = os.path.expanduser(args.out_path)

    file_handler = logging.FileHandler('logs.log')  # TODO: set a good log dir. e.g out_path  ??
    logger.addHandler(file_handler)

    logger.info(
        f"Running TraceGroomer with the following parameters:\n {args}\n")

    confidic = ut.open_config_file(os.path.expanduser(args.config_file))

    expected_keys_makeready = [
        'metadata',
                     'abundances',
                     'mean_enrichment',
                     'isotopologue_proportions',
                     'isotopologues']
    for k in expected_keys_makeready:
        assert k in confidic.keys(), logger.warning(
            f"{k} : missing in configuration file! ")  # TODO challenge this

    targetedMetabo_path = os.path.expanduser(args.labeled_metabo_file)
    # TODO verify it exists ??? located in data subfolder or not ?

    meta_path = os.path.join(out_path, f"{confidic['metadata']}.csv")
    # metadata is expected in the out path, which is the data sub-folder:
    assert os.path.isfile(meta_path), \
        logger.critical(
                f"ERROR. Metadata file must exist and be located in the "
                f"'out_path' (data subfolder) location: '{out_path}'. Abort!"
        )

    amount_mater_path = args.amountMaterial_path
    if args.amountMaterial_path is not None:
        amount_mater_path = os.path.expanduser(args.amountMaterial_path)

    # TODO: reactivate and continue logger in tidy
    #perform_type_prep(args, confidic, meta_path, targetedMetabo_path,
    #                  amount_mater_path, out_path, logger)
    return 0

if __name__ == "__main__":
    sys.exit(main())