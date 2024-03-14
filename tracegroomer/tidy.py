#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Prepare tables,
user defined config file and args

@author: johanna
"""
import os
import logging
import pandas as pd
import tracegroomer.utils as ut
from typing import Dict, Union


logger = logging.getLogger(__name__)
logger = ut.reset_log_config(logger)


class ComposedData:  # refactor this name
    def __init__(self, type_of_file):
        self.metadata: pd.DataFrame = None
        self.type_of_file: str = type_of_file
        self.material_df: pd.DataFrame = None
        self.instandard_abun_df: pd.DataFrame = None
        self.metabolites_to_drop_df: pd.DataFrame = None
        self.expected_keys_confdict = ['mean_enrichment',
                                       'isotopologue_proportions',
                                       'isotopologues',
                                       'abundances']
        self.frames_dict: Dict[str, pd.DataFrame] = dict()
        self.all_metrics_normalized_by_material: bool = False

    def load_metadata(self, metadata_path):
        self.metadata = ut.open_metadata(metadata_path)
        ut.verify_metadata_sample_not_duplicated(self.metadata)

    def load_amount_material_df(self, amount_material_path):
        if amount_material_path is not None:
            self.material_df = ut.open_amount_material(amount_material_path)

    def isocor_data_load(self, targetedMetabo_path, args, confdict):
        assert self.type_of_file == "IsoCor_out_tsv", \
            "function 'isocor_data_load' called for wrong type of file"
        isocor__df = pd.read_csv(targetedMetabo_path, sep="\t")
        frames_dict = ut.isocor_2_frames_dict(
            isocor__df, confdict, args.use_internal_standard)

        self.frames_dict = frames_dict

    def rule_tsv_data_load(self, targetedMetabo_path, args, confdict):
        assert self.type_of_file == "rule_tsv", "wrong type of file"
        data_df = pd.read_csv(targetedMetabo_path, sep="\t")
        data_df.set_index([confdict[
                               'columns_variable_description']['identifier']],
                          inplace=True)
        variable_df = pd.read_csv(
            os.path.join(confdict['groom_out_path'],
                         f"{confdict['variable_description']}.tsv"), sep="\t"
        )

        isonames = variable_df[confdict[
            'columns_variable_description']['compound']].str.cat(
            variable_df[
                confdict['columns_variable_description']['isotopologue_number']
                        ].astype(str),  sep="_m+")
        variable_df = variable_df.assign(isotopologue_name=isonames)
        variable_df.set_index([confdict[
                                   'columns_variable_description'
                               ]['identifier']], inplace=True)
        # secure names matching:
        tmp = variable_df.T[list(data_df.index)]  # same order as data matrix
        variable_df = tmp.T
        data_df.index = variable_df['isotopologue_name']
        self.frames_dict[confdict["isotopologues"]] = data_df

    def generic_xlsx_load(self, targetedMetabo_path, args, confdict):
        assert self.type_of_file == "generic_xlsx", \
            "function called for wrong type of file, must be generic xlsx"
        frames_dict = ut.excelsheets2frames_dict(targetedMetabo_path, confdict)
        tabs_isotopologues = [
            x for x in [confdict['isotopologue_proportions'],
                        confdict['isotopologues'],
                        ] if x is not None]
        assert len(tabs_isotopologues) >= 1, logger.warning(
              "Bad or no isotopologues input")
        for tab in tabs_isotopologues:  # tabs  not split by compartment yet
            tmp = frames_dict[tab]
            new_col = ut.transformmyisotopologues(tmp.columns, "generic")
            tmp.columns = new_col
            frames_dict[tab] = tmp.T

        self.frames_dict = frames_dict

    def vib_data_load(self, targetedMetabo_path, args, confdict):
        assert self.type_of_file == "VIBMEC_xlsx", \
            "function called for wrong type of file"
        frames_dict, internal_standards_df = reshape_vib_data(
            targetedMetabo_path, args, confdict)
        try:
            tmp = frames_dict[confdict['isotopologue_proportions']]
            new_col = ut.transformmyisotopologues(tmp.columns, "vib")
            tmp.columns = new_col
            frames_dict[confdict['isotopologue_proportions']] = tmp
        except Exception as e:
            logger.info(e)
            logger.error("isotopologue_proportions missing "
                         "or misspelled in config")

        self.instandard_abun_df = internal_standards_df
        self.frames_dict = frames_dict

    def transpose_frames(self):
        for x in self.frames_dict.keys():
            tmp = self.frames_dict[x].T
            self.frames_dict[x] = tmp

    def true_key_value_available_frames(self, confdict):
        reverse_dict = dict()
        avail_dict = dict()
        for m in self.expected_keys_confdict:
            reverse_dict[confdict[m]] = m
        for h in self.frames_dict.keys():
            if h is not None:
                avail_dict[reverse_dict[h]] = h

        self.available_frames = avail_dict

    def load_metabolite_to_isotopologue_df(self, confdict):
        """df of correspondences between isotopologues and metabolites
        proper to the given data"""
        try:
            isotopologues_full = list(self.frames_dict[confdict[
                                     "isotopologue_proportions"]].index)
        except TypeError:
            isotopologues_full = list(self.frames_dict[confdict[
                                      "isotopologues"]].index)
        except KeyError:
            isotopologues_full = list(self.frames_dict[confdict[
                                    "isotopologues"]].index)

        self.metabolites_isos_df = ut.isotopologues_meaning_df(
            isotopologues_full)

    def fill_missing_data(self, confdict) -> Dict[str, str]:
        tmp, confdict_new = ut.complete_missing_frames(
            confdict, self.frames_dict,
            self.metadata, self.metabolites_isos_df)
        self.frames_dict = tmp

        return confdict_new

    def save_isotopologues_preview(self, args, confdict, groom_out_path):
        compartmentalized_dict = ut.df_to__dic_bycomp(
            self.frames_dict[confdict['isotopologue_proportions']],
            self.metadata)
        output_plots_dir = os.path.join(groom_out_path, "preview_plots")
        if args.isotopologues_preview:
            if not os.path.exists(output_plots_dir):
                os.makedirs(output_plots_dir)
        ut.save_isos_preview(
            compartmentalized_dict, self.metadata,
            output_plots_dir, args.isotopologues_preview)

    def pull_internal_standard(self, confdict, args):
        if (self.instandard_abun_df is None) and (
             args.use_internal_standard is not None
        ):
            try:
                x = self.frames_dict[confdict['abundances']].columns.tolist()
                y = self.frames_dict[confdict['abundances']
                                     ].loc[args.use_internal_standard, :].tolist()
                instandard_abun_df = pd.DataFrame(
                    {"sample": x,
                     args.use_internal_standard: y
                     })
                instandard_abun_df.index = instandard_abun_df["sample"]
                self.instandard_abun_df = instandard_abun_df
            except Exception as e:
                logger.info(e)
                logger.warning("internal standard not found, continuing.")

    def normalize_isotopologues_by_material(self, args, confdict):
        """when running normalization on the isotopologues absolute values,
        the other metrics are re-computed"""
        assert args.div_isotopologues_by_amount_material, "Erroneous call"
        self.frames_dict = ut.divide_by_amount_material(
            self.frames_dict, confdict,  self.material_df,
            args.alternative_div_amount_material, metric="isotopologues")
        for k in ["isotopologue_proportions", "mean_enrichment", "abundances"]:
            # reset the metrics (except isotopologue absolute) for recomputing
            del self.frames_dict[confdict[k]]
            confdict[k] = None
        # recompute:
        new_frames_dict, new_confdict = ut.complete_missing_frames(
            confdict, self.frames_dict,
            self.metadata, self.metabolites_isos_df
        )
        self.frames_dict = new_frames_dict
        self.all_metrics_normalized_by_material = True

        return new_confdict

    def normalize_total_abundance_by_material(self, args, confdict):
        """the abundances undergo normalization iif it was not done for the
        isotopologues absolute values  (see boolean attribute)"""
        if not self.all_metrics_normalized_by_material:  # bool
            newframes_dict = ut.divide_by_amount_material(
                self.frames_dict, confdict, self.material_df,
                args.alternative_div_amount_material, metric="abundances")
            self.frames_dict = newframes_dict

    def normalize_by_internal_standard(self, args, confdict):
        """Only the total abundances are divided by the internal standard"""
        if args.use_internal_standard is not None:
            logger.info("computing normalization by internal standard")
            frames_dict = ut.abund_divideby_internalStandard(
                self.frames_dict, confdict, self.instandard_abun_df,
                args.use_internal_standard)
            self.frames_dict = frames_dict

    def compartmentalize_frames_dict(self, confdict):
        for k in self.frames_dict.keys():
            tmp = ut.df_to__dic_bycomp(
                self.frames_dict[k], self.metadata
            )
            self.frames_dict[k] = tmp

    def drop_metabolites_infile(self, exclude_list_file: Union[str, None]):
        if exclude_list_file is not None:
            tmp = self.frames_dict.copy()
            logger.info("removing metabolites as specified by user in file:")
            logger.info(f"{exclude_list_file}")
            exclude_df = pd.read_csv(exclude_list_file, sep="\t", header=0)
            try:
                unwanted_metabolites = dict()
                for co in exclude_df["compartment"].unique():
                    mets_l = exclude_df.loc[
                        exclude_df["compartment"] == co, 'metabolite'].tolist()
                    unwanted_metabolites[co] = mets_l
                tmp = ut.drop__metabolites_by_compart(tmp, unwanted_metabolites) #TODO!!
            except FileNotFoundError as err_file:
                logger.info(err_file)
            except Exception as e:
                logger.info(e)
            self.frames_dict = tmp

    def frames_filterby_min_admited_isotopol_proportions(
            self, confdict, isosprop_min_admitted: float
    ):
        isos_propor_dic = self.frames_dict[confdict['isotopologue_proportions']]
        bad_mets = dict()
        for co in isos_propor_dic.keys():
            tmp = isos_propor_dic[co]
            series_bool = tmp.le(isosprop_min_admitted).any()
            isos_bad = series_bool[series_bool]
            set_mets = set([i.split("_m+")[0] for i in isos_bad.index])
            bad_mets[co] = list(set_mets)

        tmp = ut.drop__metabolites_by_compart(self.frames_dict, bad_mets)  # TODO !!
        self.frames_dict = tmp

    def stomp_fraction_values(self, args,confdict):
        tmp = self.frames_dict.copy()
        tmp = isosprop_stomp_values(tmp, confdict,
                                            args.isosprop_stomp_values)
        tmp = meanenrich_or_fracfontrib_stomp_values(
            tmp, confdict,
            args.meanenrich_or_fracfontrib_stomp_values)
        self.frames_dict = tmp

    def transfer__abund_nan__to_all_tables(self, confdict):
        tmp = ut.transfer__abund_nan__to_all_tables(
            confdict, self.frames_dict, self.metadata)
        self.frames_dict = tmp

# end class


def reshape_vib_data(targetedMetabo_path: str, args, confdict):
    frames_dict = ut.excelsheets2frames_dict(targetedMetabo_path, confdict)
    lod_values, blanks_df, internal_standards_df, bad_x_y = ut.pull_LOD_blanks_IS(
        frames_dict[confdict['abundances']])

    frames_dict = ut.reshape_frames_dict_elems(frames_dict, bad_x_y)
    frames_dict = ut.abund_under_lod_set_nan(
        confdict, frames_dict, lod_values, args.under_detection_limit_set_nan)

    frames_dict = ut.abund_subtract_blankavg(frames_dict, confdict, blanks_df,
                                             args.subtract_blankavg)

    return frames_dict, internal_standards_df


def check_confdict_completeness(confdict):
    """Callable when input follows 'rule_tsv' type of file (3 files) """
    logger.info("further checking configuration, must include variables "
                "description")
    expected_level1 = ['variable_description', 'columns_variable_description']
    assert set(expected_level1).issubset(set(list(confdict.keys()))), \
        logger.critical(
        f"missing one or more of '{expected_level1}' in configuration")
    expected_level2 = ['identifier', 'compound', 'isotopologue_number']
    assert set(expected_level2).issubset(
        set(list(confdict['columns_variable_description'].keys()))), \
        logger.critical(f"missing one or more of '{expected_level2}'"
                        f" in columns_variable_description")


def isosprop_stomp_values(frames_dict, confdict, isosprop_stomp_vals: bool):
    if isosprop_stomp_vals:
        isos_propor_dic = frames_dict[confdict['isotopologue_proportions']]
        for co in isos_propor_dic.keys():
            df = isos_propor_dic[co]
            df[df < 0] = 0
            df[df > 1] = 1
            isos_propor_dic[co] = df
        frames_dict[confdict['isotopologue_proportions']] = isos_propor_dic
    return frames_dict


def meanenrich_or_fracfontrib_stomp_values(
        frames_dict, confdict, meanenri_or_fraccontrib_stomp_vals: bool):
    if meanenri_or_fraccontrib_stomp_vals:
        meorfc_dic = frames_dict[confdict['mean_enrichment']]
        for co in meorfc_dic.keys():
            df = meorfc_dic[co]
            df[df < 0] = 0
            df[df > 1] = 1
            meorfc_dic[co] = df
        frames_dict[confdict['mean_enrichment']] = meorfc_dic
    return frames_dict


def save_tables(frames_dict, groom_out_path) -> None:
    for k in frames_dict.keys():
        # reunify the compartments
        tmpli = list()
        for compartment in frames_dict[k].keys():
            tmp = frames_dict[k][compartment]
            tmpli.append(tmp)

        final_df = tmpli[0]    # reunify the compartments
        for i in range(1, len(tmpli)):
            final_df = pd.merge(final_df, tmpli[i], how='outer',
                                left_index=True, right_index=True)

        final_df.index.name = "ID"
        final_df = final_df.reset_index()
        final_df = final_df.drop_duplicates()
        final_df.to_csv(os.path.join(groom_out_path, f"{k}.csv"),
                        sep='\t', header=True, index=False)

        logger.info(f"File saved to: {os.path.join(groom_out_path, k)}.csv")
        print("File saved to:",  os.path.join(groom_out_path, f"{k}.csv"))
        # note : do not clear zero rows, as gives problem pd.merge


def wrapper_common_steps(combo_data: ComposedData,
                         args, confdict, groom_out_path: str) -> None:
    combo_data.load_metabolite_to_isotopologue_df(confdict)
    confdict = combo_data.fill_missing_data(confdict)
    combo_data.true_key_value_available_frames(confdict)

    combo_data.save_isotopologues_preview(args, confdict, groom_out_path)

    combo_data.pull_internal_standard(confdict, args)

    if combo_data.material_df is not None:
        logger.info("computing normalization by amount of material")
        if args.div_isotopologues_by_amount_material and (
                "isotopologues" in list(combo_data.available_frames.keys())):
            confdict = combo_data.normalize_isotopologues_by_material(
                args, confdict)
        else:
            combo_data.normalize_total_abundance_by_material(args, confdict)
    combo_data.normalize_by_internal_standard(args, confdict)
    combo_data.compartmentalize_frames_dict(confdict)
    # last steps use compartmentalized frames
    combo_data.drop_metabolites_infile(args.remove_these_metabolites)
    combo_data.stomp_fraction_values(args, confdict)
    combo_data.transfer__abund_nan__to_all_tables(confdict)
    save_tables(combo_data.frames_dict, groom_out_path)  # TODO: too many decimals, set to 6 places for all dfs


def perform_type_prep(args, confdict, metadata_used_extension: str,
                      targetedMetabo_path: str, groom_out_path
) -> None:
    combo_data = ComposedData(args.type_of_file)
    logger.info("\nLoading data")
    combo_data.load_metadata(
        os.path.join(groom_out_path,
                     f"{confdict['metadata']}{metadata_used_extension}"))

    combo_data.load_amount_material_df(args.amountMaterial_path)

    if args.type_of_file == 'IsoCor_out_tsv':
        combo_data.isocor_data_load(targetedMetabo_path, args,  confdict)

    elif args.type_of_file == 'rule_tsv':
        check_confdict_completeness(confdict)
        combo_data.rule_tsv_data_load(targetedMetabo_path, args, confdict)

    elif args.type_of_file == 'generic_xlsx':
        combo_data.generic_xlsx_load(targetedMetabo_path, args, confdict)

    elif args.type_of_file == 'VIBMEC_xlsx':
        combo_data.vib_data_load(targetedMetabo_path, args, confdict)
        combo_data.transpose_frames()
    # endif
    wrapper_common_steps(combo_data, args, confdict, groom_out_path)







