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


class CompositeData:  # refactor this name
    def __init__(self, type_of_file):
        self.metadata: pd.DataFrame = None
        self.type_of_file: str = type_of_file
        self.material_df: pd.DataFrame = None
        self.internal_standard_abun_df: pd.DataFrame = None
        self.metabolites_to_drop_df: pd.DataFrame = None
        self.frames_dict: Dict[str, pd.DataFrame] = dict()
        self.all_metrics_normalized_by_material: bool = False
        self.expected_keys_confdict = ['mean_enrichment',
                                       'isotopologue_proportions',
                                       'isotopologues',
                                       'abundances']
        self.metabolites2isotopologues_df = None
        self.user_given_names = dict()

    def load_metadata(self, metadata_path):
        self.metadata = ut.open_metadata(metadata_path)
        ut.verify_metadata_sample_not_duplicated(self.metadata)

    def load_amount_material_df(self, amount_material_path):
        if amount_material_path is not None:
            logger.info(f"loading amount of material: {amount_material_path}")
            self.material_df = ut.open_amount_material(amount_material_path)

    def load_metabolites_to_drop_df(self, exclude_list_file: Union[str, None]):
        if exclude_list_file is not None:
            logger.info(f"loading metabolites to drop: {exclude_list_file}")
            self.metabolites_to_drop_df = ut.open_metabolites_to_drop(
                exclude_list_file)

    def isocor_data_load(self, targetedMetabo_path, confdict):
        assert self.type_of_file == "IsoCor_out_tsv", \
            "function 'isocor_data_load' called for wrong type of file"
        isocor__df = pd.read_csv(targetedMetabo_path, sep="\t")
        frames_dict = ut.isocor_2_frames_dict(isocor__df, confdict)

        self.frames_dict = frames_dict

    def rule_tsv_data_load(self, targetedMetabo_path, confdict):
        assert self.type_of_file == "rule_tsv", "wrong type of file"
        data_df = pd.read_csv(targetedMetabo_path, sep="\t")
        data_df.set_index([confdict[
                               'columns_variable_description']['identifier']],
                          inplace=True)
        variable_df = pd.read_csv(
            os.path.join(confdict['groom_out_path'],
                         f"{confdict['variable_description']}.tsv"), sep="\t")

        isonames = variable_df[confdict[
            'columns_variable_description']['compound']].str.cat(
            variable_df[confdict[
                'columns_variable_description']['isotopologue_number']
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

    def generic_xlsx_load(self, targetedMetabo_path, confdict):
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
            frames_dict[tab] = tmp

        self.frames_dict = frames_dict

    def vib_data_load(self, targetedMetabo_path, args, confdict):
        assert self.type_of_file == "VIBMEC_xlsx", \
            "function called for wrong type of file"
        frames_dict, internal_standards_df = ut.reshape_vib_data(
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

        self.internal_standard_abun_df = internal_standards_df
        self.frames_dict = frames_dict

    def transpose_frames(self):
        for x in self.frames_dict.keys():
            tmp = self.frames_dict[x].T
            self.frames_dict[x] = tmp

    def load_metabolite_to_isotopologue_df(self, confdict):
        """df of correspondences between isotopologues and metabolites
        proper to the given data"""
        if (confdict['isotopologue_proportions'] is not None) or (
                confdict['isotopologues'] is not None):
            try:
                isotopologues_full = list(self.frames_dict[confdict[
                                         "isotopologue_proportions"]].index)
            # ok apply same proposed solution whether Key or Type error:
            except TypeError:
                isotopologues_full = list(self.frames_dict[confdict[
                                          "isotopologues"]].index)
            except KeyError:
                isotopologues_full = list(self.frames_dict[confdict[
                                        "isotopologues"]].index)

            self.metabolites2isotopologues_df = ut.isotopologues_meaning_df(
                isotopologues_full)

    def set_user_given_names(self, confdict):
        """
        Set user given names of quantifications only.
        This will be useful for final_files_names property"""
        user_given_names = dict()
        for keyname_quantif in self.expected_keys_confdict:
            if confdict[keyname_quantif] is not None:
                user_given_names[keyname_quantif] = confdict[
                    keyname_quantif]
        self.user_given_names = user_given_names

    def fill_missing_data(self, confdict) -> Dict[str, str]:
        """
        Computes the quantification that the user set None in confdict.
        This is done using the absolute isotopic values, and if not provided,
        will compute mean enrichment from isotopic proportions.
        """
        tmp, confdict_new = ut.complete_missing_frames(
            confdict, self.frames_dict, self.metabolites2isotopologues_df)
        self.frames_dict = tmp

        return confdict_new

    def update_truly_available_frames(self, confdict):
        """
        Sets the properties .available_frames and .reverse_available_frames,
        which are dictionaries of quantifications names that truly exist
          in the object (in the frames_dict)
        """
        reverse_dict = dict()
        avail_dict = dict()
        true_reverse_dict = dict()
        for m in self.expected_keys_confdict:
            try:
                reverse_dict[confdict[m]] = m
            except KeyError:
                continue
        for h in self.frames_dict.keys():
            # if the quantification content and key exists
            if (self.frames_dict[h] is not None) and (h is not None):
                try:
                    avail_dict[reverse_dict[h]] = h
                    true_reverse_dict[h] = reverse_dict[h]
                except Exception as e:
                    print(e)
                    continue

        self.available_frames = avail_dict
        self.reverse_available_frames = true_reverse_dict

    def save_isotopologues_preview(self, args, confdict, groom_out_path):
        if args.isotopologues_preview:
            compartmentalized_dict = ut.df_to__dic_bycomp(
                self.frames_dict[confdict['isotopologue_proportions']],
                self.metadata)
            output_plots_dir = os.path.join(groom_out_path, "preview_plots")
            logger.info("prepare isotopologue proportions overview figures")
            if not os.path.exists(output_plots_dir):
                os.makedirs(output_plots_dir)
            ut.save_isos_preview(
                compartmentalized_dict, self.metadata, output_plots_dir,
                args.isotopologues_preview, args.output_files_extension)
            logger.info(f"saved figures to {output_plots_dir}")

    def pull_internal_standard(self, confdict, args):
        if (self.internal_standard_abun_df is None) and (
             args.use_internal_standard is not None
        ):
            try:
                x = self.frames_dict[confdict['abundances']].columns.tolist()
                y = self.frames_dict[confdict['abundances']].loc[
                    args.use_internal_standard, :].tolist()
                instandard_abun_df = pd.DataFrame(
                    {"sample": x,
                     args.use_internal_standard: y
                     })
                instandard_abun_df.index = instandard_abun_df["sample"]
                self.internal_standard_abun_df = instandard_abun_df
            except Exception as e:
                logger.info(e)
                logger.warning("internal standard not found. Continue.")

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
            confdict, self.frames_dict, self.metabolites2isotopologues_df
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
                self.frames_dict, confdict, self.internal_standard_abun_df,
                args.use_internal_standard)
            self.frames_dict = frames_dict

    def set_final_files_names(self):
        """
        Set final names of output files:
          1. set final names dictionary and add as attribute, and
          2. assign the values of 1. to the keys of the object.frames_dict
        """
        not_user_defined_dict = ut.retrieve_dict_not_user_defined()
        final_files_names_d = dict()  # set final names dictionary
        for valuename in self.reverse_available_frames.keys():
            keyname = self.reverse_available_frames[valuename]
            if keyname in list(self.user_given_names.keys()):
                final_files_names_d[keyname] = self.user_given_names[keyname]
            else:
                if self.available_frames[keyname] is not None:
                    final_files_names_d[keyname] = not_user_defined_dict[
                        keyname]
        # end for
        frames_names_list = list(self.frames_dict.keys())
        for frame_name in frames_names_list:  # assign final names to frames
            if frame_name == "abundances_computed":
                self.frames_dict[final_files_names_d[
                    "abundances"]] = self.frames_dict[frame_name]
                del self.frames_dict[frame_name]
            if frame_name == "mean_enrichment_computed":
                self.frames_dict[final_files_names_d[
                    "mean_enrichment"]] = self.frames_dict[frame_name]
                del self.frames_dict[frame_name]
            if frame_name == "isotopologue_props_computed":
                self.frames_dict[final_files_names_d[
                    "isotopologue_proportions"]] = self.frames_dict[
                    frame_name]
                del self.frames_dict[frame_name]
        # end for
        self.final_files_names = final_files_names_d

    def compartmentalize_frames_dict(self):
        for k in self.frames_dict.keys():
            tmp = ut.df_to__dic_bycomp(
                self.frames_dict[k], self.metadata
            )
            self.frames_dict[k] = tmp

    def drop_metabolites(self):
        if self.metabolites_to_drop_df is not None:
            tmp = self.frames_dict.copy()
            logger.info("removing user specified metabolites, by compartment")
            try:
                unwanted_metabolites = dict()
                for co in self.metabolites_to_drop_df["compartment"].unique():
                    mets_l = self.metabolites_to_drop_df.loc[
                        self.metabolites_to_drop_df["compartment"] == co,
                        'metabolite'].tolist()
                    unwanted_metabolites[co] = mets_l
                tmp = ut.drop__metabolites_by_compart(
                    tmp, self.reverse_available_frames, unwanted_metabolites)
            except Exception as e:
                logger.warning(e)
                logger.warning("Metabolites removal not possible. Continue")

            self.frames_dict = tmp

    def frames_filterby_min_admited_isotopol_proportions(
            self, confdict, isosprop_min_admitted: float
    ):
        isos_propor_dic = self.frames_dict[
            confdict['isotopologue_proportions']]
        bad_mets = dict()
        for co in isos_propor_dic.keys():
            tmp = isos_propor_dic[co]
            series_bool = tmp.le(isosprop_min_admitted).any()
            isos_bad = series_bool[series_bool]
            set_mets = set([i.split("_m+")[0] for i in isos_bad.index])
            bad_mets[co] = list(set_mets)

        tmp = ut.drop__metabolites_by_compart(
            self.frames_dict, self.reverse_available_frames, bad_mets)
        self.frames_dict = tmp

    def stomp_fraction_values(self, args, final_files_names: dict):
        if args.fractions_stomp_values:
            for frac_type in ["mean_enrichment", "isotopologue_proportions"]:
                try:
                    curr_dict = self.frames_dict[final_files_names[frac_type]]
                    for co in curr_dict.keys():
                        df = curr_dict[co]
                        df[df < 0] = 0
                        df[df > 1] = 1
                        curr_dict[co] = df
                    self.frames_dict[final_files_names[frac_type]] = curr_dict
                except KeyError as e:
                    logger.info(f"{e}: unavailable (stomp fractions)")

    def transfer__abund_nan__to_all_tables(self, final_files_names: dict):
        try:
            tmp = ut.transfer__abund_nan__to_all_tables(
                final_files_names, self.frames_dict, self.metadata)
            self.frames_dict = tmp
        except KeyError as e:
            logger.info(f"{e}: unavailable (propagate NaN values)")

# end class


def check_confdict_completeness(confdict) -> None:
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


def save_tables(frames_dict, groom_out_path, output_extension) -> None:
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
        o_file_path = os.path.join(groom_out_path, f"{k}.{output_extension}")
        final_df.to_csv(o_file_path, sep='\t', header=True, index=False)

        logger.info(f"File saved to: {o_file_path}")
        # note : do not clear zero rows, as gives problem pd.merge


def wrapper_common_steps(combo_data: CompositeData,
                         args, confdict, groom_out_path: str) -> None:
    combo_data.load_metabolite_to_isotopologue_df(confdict)
    combo_data.set_user_given_names(confdict)  # only user def names
    confdict = combo_data.fill_missing_data(confdict)  # critical completion
    combo_data.update_truly_available_frames(confdict)  # update 1
    combo_data.save_isotopologues_preview(args, confdict, groom_out_path)
    combo_data.pull_internal_standard(confdict, args)  # before normalisation

    if combo_data.material_df is not None:
        logger.info("computing normalization by amount of material")
        if args.div_isotopologues_by_amount_material and (
                "isotopologues" in list(combo_data.available_frames.keys())):
            confdict = combo_data.normalize_isotopologues_by_material(
                args, confdict)
        else:
            combo_data.normalize_total_abundance_by_material(args, confdict)

    combo_data.normalize_by_internal_standard(args, confdict)
    combo_data.set_final_files_names()
    combo_data.update_truly_available_frames(combo_data.final_files_names)  # update 2
    combo_data.compartmentalize_frames_dict()
    # last steps use compartmentalized frames
    combo_data.drop_metabolites()
    combo_data.stomp_fraction_values(args, combo_data.final_files_names)
    combo_data.transfer__abund_nan__to_all_tables(
            combo_data.final_files_names)
    save_tables(combo_data.frames_dict, groom_out_path,
                args.output_files_extension)


def perform_type_prep(args, confdict, metadata_used_extension: str,
                      targetedMetabo_path: str, groom_out_path) -> None:
    combo_data = CompositeData(args.type_of_file)
    logger.info("\nLoading data")
    combo_data.load_metadata(
        os.path.join(groom_out_path,
                     f"{confdict['metadata']}{metadata_used_extension}"))

    combo_data.load_amount_material_df(args.amountMaterial_path)
    combo_data.load_metabolites_to_drop_df(args.remove_these_metabolites)

    if args.type_of_file == 'IsoCor_out_tsv':
        combo_data.isocor_data_load(targetedMetabo_path, confdict)

    elif args.type_of_file == 'rule_tsv':
        check_confdict_completeness(confdict)
        combo_data.rule_tsv_data_load(targetedMetabo_path, confdict)

    elif args.type_of_file == 'generic_xlsx':
        combo_data.generic_xlsx_load(targetedMetabo_path, confdict)
        combo_data.transpose_frames()

    elif args.type_of_file == 'VIBMEC_xlsx':
        combo_data.vib_data_load(targetedMetabo_path, args, confdict)
        combo_data.transpose_frames()
    # endif
    wrapper_common_steps(combo_data, args, confdict, groom_out_path)
