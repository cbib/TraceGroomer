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
import numpy as np
import re
import tracegroomer.utils as ut
from typing import Dict


logger = logging.getLogger(__name__)
logger = ut.reset_log_config(logger)


class Myclass:  # refactor this name
    def __init__(self, type_of_file):
        self.metadata: pd.DataFrame = None
        self.type_of_file: str = type_of_file
        self.all_metrics_normalized_by_material = False
        self.material_df: pd.DataFrame = None
        self.instandard_abun_df: pd.DataFrame = None
        self.expected_keys_confdict = ['mean_enrichment',
                                       'isotopologue_proportions',
                                      'isotopologues',
                                       'abundances']

    def initialize_frames_dict(self, confdict):
        my_dict = dict()
        my_dict[confdict["isotopologues"]]: Dict[str, pd.DataFrame]= None,
        my_dict[confdict["isotopologue_proportions"]]: Dict[str, pd.DataFrame]= None,
        my_dict[confdict["abundances"]]: Dict[str, pd.DataFrame]= None,
        my_dict[confdict["mean_enrichment"]]: Dict[str, pd.DataFrame]= None,
        self.frames_dict = my_dict

    def load_metadata(self, metadata_path):
        self.metadata = ut.open_metadata(metadata_path)
        ut.verify_metadata_sample_not_duplicated(self.metadata)

    def load_material(self, amount_material_path):
        if amount_material_path is not None:
            self.material_df = ut.open_amount_material(amount_material_path)

    def isocor_data_load(self, targetedMetabo_path, args, confdict):
        assert self.type_of_file == "IsoCor_out_tsv", \
            "function 'load_isocor_data' called for wrong type of file"
        isocor__df = pd.read_csv(targetedMetabo_path, sep="\t")
        frames_dict, instandard_abun_df = ut.isocor_2_frames_dict(
            isocor__df, self.metadata,
            confdict, args.use_internal_standard)

        self.frames_dict = frames_dict
        self.instandard_abun_df = instandard_abun_df

    def load_metabolite_to_isotopologue_df(self, confdict):
        """df of correspondences between isotopologues and metabolites
        proper to the given data"""
        isotopologues_full = self.frames_dict[confdict[
            "isotopologue_proportions"]].columns
        self.metabolites_isos_df = ut.isotopologues_meaning_df(isotopologues_full)

    def compartmentalize_dict(self, confdict):
        for k in self.expected_keys_confdict:
            tmp = ut.df_to__dic_bycomp(
                self.frames_dict[confdict[k]], self.metadata
            )
            self.frames_dict[confdict[k]] = tmp

    def normalize_isotopologues_by_material(self, args, confdict):
        """when running normalization in the isotopologues absolute values,
        the other metrics are re-computed: the abundances do not undergo
        further normalization by material (see boolean attribute)"""
        assert args.div_isotopologues_by_amount_material, "Error function call"
        self.frames_dict =  ut.isosAbsol_divideby_amount_material(
            self.frames_dict, confdict,  amount_material,
            args.alternative_div_amount_material)

        for k in self.expected_keys_confdict:
            confdict[k] = None  # reset the metrics for recomputing
            self.frames_dict[confdict[k]] = None

        new_frames_dict, new_confdict = ut.complete_missing_frames(
            confdict, self.frames_dict,
            self.metadata, self.metabolites_isos_df
        )
        self.frames_dict = new_frames_dict
        self.all_metrics_normalized_by_material = True

        return new_confdict

    def save_isotopologues_preview(self, args, confdict, groom_out_path):
        output_plots_dir = os.path.join(groom_out_path, "preview_plots")
        if args.isotopologues_preview:
            if not os.path.exists(output_plots_dir):
                os.makedirs(output_plots_dir)
        ut.save_isos_preview(self.frames_dict[confdict['isotopologue_proportions']],
                          self.metadata,
                          output_plots_dir, args.isotopologues_preview)

    def normalize_total_abundance_material(self, args, confdict):
        if not self.all_metrics_normalized_by_material:
            newframes_dict = ut.abund_divideby_amount_material(  # total abundance
                self.frames_dict, confdict, self.material_df,
                args.alternative_div_amount_material)
            self.frames_dict = newframes_dict

    def normalize_by_internal_standard(self, args, confdict):
        """Only the total abundances are divided by the internal standard"""
        frames_dict = ut.abund_divideby_internalStandard(
            self.frames_dict, confdict, self.instandard_abun_df,
            args.use_internal_standard)
        self.frames_dict = frames_dict




def drop__metabolites_by_compart(frames_dict_orig: dict,
                                 bad_metabolites_dic: dict) -> dict:
    frames_dict = frames_dict_orig.copy()
    for tab_name in frames_dict.keys():
        for co in bad_metabolites_dic.keys():

            if "isotopolog" in tab_name.lower():
                tmpdf = frames_dict[tab_name][co]
                to_drop_now_isos = list()
                for i in tmpdf.columns:
                    for j in bad_metabolites_dic[co]:
                        if i.startswith(j):
                            to_drop_now_isos.append(i)
                tmpdf = tmpdf.drop(columns=to_drop_now_isos)
                frames_dict[tab_name][co] = tmpdf

            elif "isotopolog" not in tab_name.lower():
                tmpdf = frames_dict[tab_name][co]
                to_drop_now = bad_metabolites_dic[co]
                tmpdf = tmpdf.drop(columns=to_drop_now)
                frames_dict[tab_name][co] = tmpdf

    return frames_dict


def do_generic_prep(meta_path, targetedMetabo_path, args, confdict,
                    amount_mater_path, output_plots_dir):
    metadata = ut.open_metadata(meta_path)
    #ut.verify_metadata_sample_not_duplicated(metadata) # TODO del
    frames_dict = excelsheets2frames_dict(targetedMetabo_path, confdict)
    tabs_isotopologues = [s for s in frames_dict.keys() if
                          "isotopol" in s.lower()]
    assert len(
        tabs_isotopologues) >= 1, "\nError, bad or no isotopologues input"
    for tab in tabs_isotopologues:  # tabs are not split by compartment here
        tmp = frames_dict[tab]
        new_col = transformmyisotopologues(tmp.columns, "generic")
        tmp.columns = new_col
        frames_dict[tab] = tmp
    # end for

    isotopologues_full = frames_dict[tabs_isotopologues[0]].columns
    metabolites_isos_df = ut.isotopologues_meaning_df(isotopologues_full)

    for k in frames_dict.keys():
        tmp_co_dic = df_to__dic_bycomp(frames_dict[k],
                                       metadata)  # split by compartment
        frames_dict[k] = tmp_co_dic

    frames_dict, confdict_new = ut.complete_missing_frames(confdict, frames_dict,
                                                       metadata,
                                                       metabolites_isos_df)

    if args.use_internal_standard is not None:
        instandard_abun_l = list()
        abu = confdict_new['abundances']
        for co in frames_dict[abu].keys():
            tmp_co = frames_dict[abu][co][args.use_internal_standard]
            instandard_abun_l.append(
                pd.DataFrame({args.use_internal_standard: tmp_co}))
        instandard_abun_df = pd.concat(instandard_abun_l)

        frames_dict = abund_divideby_internalStandard(
            frames_dict, confdict_new,
            instandard_abun_df,
            args.use_internal_standard)
    # end if

    frames_dict = abund_divideby_amount_material(    # generic, only abund
        frames_dict, confdict_new,
        amount_mater_path,
        args.alternative_div_amount_material)



    save_isos_preview(frames_dict[confdict_new['isotopologue_proportions']],
                      metadata,
                      output_plots_dir, args.isotopologues_preview)

    return frames_dict, confdict_new


def drop_metabolites_infile(frames_dict, exclude_list_file: str|None):
    if exclude_list_file is not None:
        logger.info("removing metabolites as specified by user in file:")
        logger.info(exclude_list_file, "\n")
        exclude_df = pd.read_csv(exclude_list_file, sep="\t", header=0)
        try:
            unwanted_metabolites = dict()
            for co in exclude_df["compartment"].unique():
                mets_l = exclude_df.loc[
                    exclude_df["compartment"] == co, 'metabolite'].tolist()
                unwanted_metabolites[co] = mets_l
            frames_dict = drop__metabolites_by_compart(frames_dict,
                                                      unwanted_metabolites)
        except FileNotFoundError as err_file:
            logger.info(err_file)
        except Exception as e:
            logger.info(e)
    return frames_dict


def frames_filterby_min_admited_isosprop(frames_dict, confdict,
                                         isosprop_min_admitted: float):
    isos_propor_dic = frames_dict[confdict['isotopologue_proportions']]
    bad_mets = dict()
    for co in isos_propor_dic.keys():
        tmp = isos_propor_dic[co]
        series_bool = tmp.le(isosprop_min_admitted).any()
        isos_bad = series_bool[series_bool]
        set_mets = set([i.split("_m+")[0] for i in isos_bad.index])
        bad_mets[co] = list(set_mets)

    frames_dict = drop__metabolites_by_compart(frames_dict, bad_mets)

    return frames_dict


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


#transfer_nan_all_tables


def perform_type_prep(args, confdict, metadata_used_extension: str,
        targetedMetabo_path: str, groom_out_path
) -> None:
    myobj = Myclass(args.type_of_file)
    myobj.initialize_frames_dict(confdict)
    myobj.load_metadata(os.path.join(groom_out_path,
                     f"{confdict['metadata']}{metadata_used_extension}"))
    myobj.load_material(args.amountMaterial_path)

    if args.type_of_file == 'IsoCor_out_tsv':
        myobj.isocor_data_load(targetedMetabo_path, args,
                               confdict)
        myobj.load_metabolite_to_isotopologue_df(confdict)
        myobj.compartmentalize_dict(confdict)
        if args.div_isotopologues_by_amount_material and (
                myobj.material_df is not None):
            confdict = myobj.normalize_isotopologues_by_material(
                args, confdict)

        myobj.save_isotopologues_preview(args, confdict, groom_out_path)
        myobj.normalize_total_abundance_material(args, confdict)
        myobj. normalize_by_internal_standard(args, confdict)

   # TODO: elif args.type_of_file... rule_tsv :
    elif args.type_of_file == 'VIBMEC_xlsx':
        frames_dict = do_vib_prep(myobj, targetedMetabo_path, args,
                                 confdict,
                                 amount_mater_path, output_plots_dir)
    elif args.type_of_file == 'generic_xlsx':  # generates also confdict !!!
        frames_dict, confdict = do_generic_prep(myobj, targetedMetabo_path,
                                               args, confdict,
                                               amount_mater_path,
                                               output_plots_dir)
    logger.info("i am hereeeee")
    # common steps to any preparation type:
    frames_dict = drop_metabolites_infile(frames_dict,
                                         args.remove_these_metabolites)
    frames_dict = frames_filterby_min_admited_isosprop(
        frames_dict, confdict,
        args.isosprop_min_admitted)
    frames_dict = isosprop_stomp_values(frames_dict, confdict,
                                       args.isosprop_stomp_values)
    frames_dict = meanenrich_or_fracfontrib_stomp_values(
        frames_dict, confdict,
        args.meanenrich_or_fracfontrib_stomp_values)
    frames_dict = transfer__abund_nan__to_all_tables(confdict, frames_dict,
                                                    meta_path)

    for k in frames_dict.keys():
        # reunify the compartments
        tmpli = list()
        for compartment in frames_dict[k].keys():
            tmp = frames_dict[k][compartment].T
            tmpli.append(tmp)

        finalk = tmpli[0]
        for i in range(1, len(tmpli)):
            finalk = pd.merge(finalk, tmpli[i], how='outer',
                              left_index=True, right_index=True)
        # reunified the compartments
        # note : do not clear zero rows, as gives problem pd.merge
        finalk.index.name = "ID"
        finalk = finalk.reset_index()
        finalk = finalk.drop_duplicates()
        finalk.to_csv(
            os.path.join(groom_out_path, f"{k}.csv"),
            sep='\t', header=True, index=False)
        print("File saved to:",  os.path.join(groom_out_path, f"{k}.csv"))