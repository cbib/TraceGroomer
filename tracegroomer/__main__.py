import argparse
import os
import sys

import tracegroomer.utils as fg
from tracegroomer.tidy import perform_type_prep

def prep_args() -> argparse.ArgumentParser:
    show_defaults = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(prog="python -m Tracegroomer.tidy",
                                     formatter_class=show_defaults)

    parser.add_argument('config', type=str,
                        help="Configuration file in absolute path")

    parser.add_argument('--targetedMetabo_path', type=str, default=None,
                        help="the absolute path to your file")

    parser.add_argument('--type_of_file', type=str, default=None,
                        help="IsoCor_out_tsv|VIBMEC_xlsx|generic_xlsx")

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
                           of material (number of cells, tissue weight, etc) by sample")

    parser.add_argument("--alternative_div_amount_material",
                        action=argparse.BooleanOptionalAction, default=False,
                        help="On VIB results, when dividing abundances \
                        by the amount of material, \
                        also multiplies by mean(amountMaterial) \
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
                        compartment, metabolite. This file contains \
                        metabolites to be completely excluded from all  \
                        analysis (you know what you are doing).")
    # all tables affected

    return parser


def main() -> int:
    parser = prep_args()
    args = parser.parse_args()
    configfile = os.path.expanduser(args.config)
    confidic = fg.open_config_file(configfile)

    expected_keys_makeready = [
        'metadata',
                     'abundances',
                     'mean_enrichment',
                     'isotopologue_proportions',
                     'isotopologues']
    for k in expected_keys_makeready:
        assert k in confidic.keys(), f"{k} : missing in configuration file! "

   #meta_path = os.path.expanduser(confidic['metadata'])
    targetedMetabo_path = os.path.expanduser(args.targetedMetabo_path)
    groom_out_path = os.path.expanduser(confidic['groom_out_path'])

    meta_path = os.path.join(groom_out_path, f"{confidic['metadata']}.csv")
    amount_mater_path = args.amountMaterial_path
    if args.amountMaterial_path is not None:
        amount_mater_path = os.path.expanduser(
            args.amountMaterial_path)

    perform_type_prep(args, confidic, meta_path, targetedMetabo_path,
                      amount_mater_path, groom_out_path)
    return 0

if __name__ == "__main__":
    sys.exit(main())