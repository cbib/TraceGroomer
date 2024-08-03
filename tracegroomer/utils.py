#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 12:24:01 2022

@author: johanna
"""
import sys
import os
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import re
from typing import Dict, List, Tuple, Union


def reset_log_config(logger):
    """parameters for the log file"""
    out_log_file_name = "groom.log"
    logging.basicConfig(encoding='utf-8', stream=sys.stdout)
    file_handler = logging.FileHandler(out_log_file_name)
    logger.addHandler(file_handler)
    logger.setLevel(logging.DEBUG)
    return logger


def open_config_file(config_file_path: str) -> Dict[str, Union[str, dict]]:
    try:
        with open(config_file_path, "r") as f:
            confdict = yaml.load(f, Loader=yaml.Loader)
    except yaml.YAMLError as yam_err:
        print(yam_err)
        confdict = None
    except Exception as e:
        print(e)
        confdict = None

    if confdict is None:
        raise ValueError("\nimpossible to read configuration file")

    return confdict


def open_metadata(file_path: str) -> pd.DataFrame:
    try:
        metadata = pd.read_csv(file_path, sep='\t')
        return metadata
    except Exception as e:
        print(f'{e}\nproblem with opening metadata file')
        metadata = None
    if metadata is None:
        raise ValueError("\nproblem opening metadata file")


def open_amount_material(amount_material_path: str) -> pd.DataFrame:
    if amount_material_path is not None:
        try:
            file = amount_material_path
            material_df = pd.read_csv(file, sep='\t', index_col=0)

            assert material_df.shape[1] == 1, \
                "amountMaterial_path: table must have only 2 columns"

            assert (material_df.iloc[:, 0] <= 0).sum() == 0, "amountMaterial\
                 _path: table must not contain zeros nor negative numbers"

        except FileNotFoundError as err_file:
            print(err_file)
        except UnboundLocalError as uerr:
            print(uerr, "amountMaterial_path:  check spelling")
        except Exception as e:
            print(e, "unexpected error, could not open amountMaterial_path")

    return material_df


def retrieve_dict_not_user_defined() -> Dict[str, str]:
    """
    Yields a dictionary of possible replacement names, usable when the
    user set None for given type(s) of quantification in the configuration.
    Note: If None for a type of quantif, TraceGrommer computed it internally
    """
    not_user_defined_dict = {
        'abundances': 'TotalMetaboliteAbundances',
        'mean_enrichment': 'MeanEnrichmentData',
        'isotopologues': 'IsotopologueAbsoluteValues',
        'isotopologue_proportions': 'IsotopologueProportions'}
    return not_user_defined_dict


def open_metabolites_to_drop(exclude_list_file: str) -> pd.DataFrame:
    try:
        exclude_df = pd.read_csv(exclude_list_file, sep="\t", header=0)
        expected_columns = ["compartment", "metabolite"]
        if not set(list(exclude_df.columns)).issubset(set(expected_columns)):
            print(f"columns: {expected_columns} missing in the file\n")
    except FileNotFoundError as err_file:
        print(err_file)
    except Exception as e:
        print(e)

    return exclude_df


def compute_abund_from_absolute_isotopol(df, metabos_isos_df) -> pd.DataFrame:
    """
    input:
       df: input isotopologues in absolute values
       metabos_isos_df: df created internally (see 'isotopologues_meaning_df')
    output:
       dataframe (samples as columns) of total abundances
    """
    metabos_uniq = metabos_isos_df['metabolite'].unique()
    abundance = pd.DataFrame(index=metabos_uniq, columns=df.columns)
    for m in metabos_uniq:
        isos_here = metabos_isos_df.loc[
            metabos_isos_df['metabolite'] == m, 'isotopologue_name']
        sub_df = df.loc[isos_here, :]
        sub_df_sum = sub_df.sum(axis=0,
                                skipna=False)
        # False makes sure that, if all values nan, sum = nan
        abundance.loc[m, :] = sub_df_sum
    return abundance


def compute_isotopologues_proportions_from_absolute(df, metabos_isos_df):
    """
    input:
       df : input isotopologues in absolute values
       metabos_isos_df : df created internally (see 'isotopologues_meaning_df')
    output:
       dataframe (samples as columns) of isotopologue proportions
    """
    metabos_uniq = metabos_isos_df['metabolite'].unique()
    isos_prop = pd.DataFrame(index=df.index, columns=df.columns)
    for m in metabos_uniq:
        isos_here = metabos_isos_df.loc[
            metabos_isos_df['metabolite'] == m, 'isotopologue_name']
        sub_df = df.loc[isos_here, :]
        sub_df_sum = sub_df.sum(axis=0, skipna=False)
        proportions_m = sub_df.div(sub_df_sum.T, axis=1)
        isos_prop.loc[isos_here.tolist(), :] = proportions_m
    isos_prop = isos_prop.round(decimals=9)
    return isos_prop


def compute_MEorFC_from_isotopologues_proportions(df, metabos_isos_df):
    """
    computes mean enrichment (a.k.a fractional contributions)
    input:
      df : isotopologue proportions (whether computed here, or from input)
      metabos_isos_df : df created internally (see 'isotopologues_meaning_df')
    output:
      dataframe (samples as columns) of mean enrichment
    Note:
    mean enrichment (i) = ( (i_m+0 * 0) + (i_m+1 * 1) +...+ (i_m+n * n)  ) / n
          where
           i is one metabolite,
           n the max nb of carbons that can be marked for that metabolite,
           each i_m+x is a value comprised between 0 and 1.
    """
    isos_prop = df
    metabos_uniq = metabos_isos_df['metabolite'].unique()
    # empty df
    meanenrich_or_fraccontrib = pd.DataFrame(index=metabos_uniq,
                                             columns=isos_prop.columns)
    for m in metabos_uniq:
        isos_here = metabos_isos_df.loc[
            metabos_isos_df['metabolite'] == m, 'isotopologue_name']
        coefs = [int(i.split("_m+")[1]) for i in isos_here.tolist()]
        sub_df = isos_prop.loc[isos_here, :]
        sub_df['coefs'] = coefs  # 0,1,3..,n integers are the coefficients
        # compute the factors, produce another pandas df coefs_fracs_prod
        coefs_fracs_prod = sub_df.multiply(sub_df['coefs'], axis=0)
        # line above makes coefs column be multiplied by itself,
        # In a future improve this, for now just dropping the coefs col :
        coefs_fracs_prod.drop(columns=['coefs'], inplace=True)
        # sum the factors
        numerator_val = coefs_fracs_prod.sum(axis=0, skipna=False)
        # divide by n, place that scalar in the empty df
        me_fc_this_metabolite = numerator_val / max(coefs)
        me_fc_this_metabolite.name = m
        meanenrich_or_fraccontrib.loc[m, :] = me_fc_this_metabolite
    meanenrich_or_fraccontrib = meanenrich_or_fraccontrib.round(decimals=9)
    return meanenrich_or_fraccontrib


def complete_missing_frames(confdict, frames_dict, metabolites_isos_df
                            ) -> Tuple[Dict[str, pd.DataFrame], dict]:
    """can apply to any type of inputs, compartmentalized"""
    confdict_new = confdict.copy()

    if confdict['abundances'] is None:
        if confdict['isotopologues'] is not None:
            frames_dict["abundances_computed"] = dict()
            tmp = compute_abund_from_absolute_isotopol(
                frames_dict[confdict['isotopologues']], metabolites_isos_df)
            frames_dict["abundances_computed"] = tmp.astype(float)
            confdict_new['abundances'] = "abundances_computed"
        elif confdict['isotopologues'] is None:
            print(" isotopologues' absolute values not available,\
                 impossible to get abundance")

    if confdict['isotopologue_proportions'] is None:
        if confdict['isotopologues'] is not None:
            frames_dict['isotopologue_props_computed'] = dict()
            tmp = compute_isotopologues_proportions_from_absolute(
                frames_dict[confdict['isotopologues']], metabolites_isos_df
            )
            frames_dict["isotopologue_props_computed"] = tmp.astype(float)
            confdict_new[
                'isotopologue_proportions'] = "isotopologue_props_computed"
        elif confdict['isotopologues'] is None:
            print(" isotopologues' absolute values not available, \
                impossible to get proportions")

    if confdict['mean_enrichment'] is None:
        try:
            frames_dict["mean_enrichment_computed"] = dict()
            tmp = compute_MEorFC_from_isotopologues_proportions(
                frames_dict[confdict_new['isotopologue_proportions']],
                metabolites_isos_df
            )
            frames_dict["mean_enrichment_computed"] = tmp.astype(float)
            confdict_new['mean_enrichment'] = "mean_enrichment_computed"
        except Exception as e:
            print(f"{e}\n impossible to calculate: mean enrichment or fract"
                  f" contribution. Isotopologue proportions not found")

    return frames_dict, confdict_new


def df_to__dic_bycomp(df: pd.DataFrame, metadata: pd.DataFrame) -> dict:
    # splits df into dictionary of dataframes, each for one compartment:
    out_dic = dict()
    for co in metadata['compartment'].unique():
        metada_co = metadata.loc[metadata['compartment'] == co, :]
        df_co = df[list(metada_co['original_name'])]
        out_dic[co] = df_co
    return out_dic


def isocor_2_frames_dict(isocor_input_df,  confdict) -> Dict[str, pd.DataFrame]:
    """ function exclusive for IsoCor type of file
        all the 4 types of quantification are transformed in independent
        data frames, and kept into a dictionary"""
    df = isocor_input_df[
        ['sample', 'metabolite', 'isotopologue', 'corrected_area',
         'isotopologue_fraction', 'mean_enrichment']]
    # converting to the "m+x" style, the isotopologues names :
    isonames = df.metabolite.str.cat(df.isotopologue.astype(str),
                                     sep="_m+")
    df = df.assign(isotopologue_name=isonames)

    metabos_isos_df = df[['metabolite', 'isotopologue_name']]
    metabos_isos_df = metabos_isos_df.drop_duplicates()

    me_or_fc_melted = df[['sample', 'metabolite', 'mean_enrichment']]
    me_or_fc_melted = me_or_fc_melted.drop_duplicates()
    me_or_fc = me_or_fc_melted.pivot(index='metabolite', columns='sample')
    me_or_fc = me_or_fc['mean_enrichment'].reset_index()
    me_or_fc = me_or_fc.set_index('metabolite')

    isos_prop_melted = df[
        ['sample', 'isotopologue_name', 'isotopologue_fraction']]
    isos_prop_melted = isos_prop_melted.drop_duplicates()
    isos_prop = isos_prop_melted.pivot(index='isotopologue_name',
                                       columns='sample')
    isos_prop = isos_prop['isotopologue_fraction'].reset_index()
    isos_prop = isos_prop.set_index('isotopologue_name')

    isos_absolute_melted = df[
        ['sample', 'isotopologue_name', 'corrected_area']]
    isos_absolute_melted = isos_absolute_melted.drop_duplicates()
    isos_absolute = isos_absolute_melted.pivot(index='isotopologue_name',
                                               columns='sample')
    isos_absolute = isos_absolute['corrected_area'].reset_index()
    isos_absolute = isos_absolute.set_index('isotopologue_name')

    abundance = compute_abund_from_absolute_isotopol(isos_absolute,
                                                     metabos_isos_df)
    frames_dict = dict()
    frames_dict[confdict['mean_enrichment']] = me_or_fc
    frames_dict[confdict['isotopologue_proportions']] = isos_prop
    frames_dict[confdict['isotopologues']] = isos_absolute
    frames_dict[confdict['abundances']] = abundance

    return frames_dict


def check_config_and_sheets_match(sheetsnames: List[str],
                                  list_config_tabs: List[str]) -> None:
    name_notfound = set(list_config_tabs) - set(sheetsnames)
    message = f"One or more name_ arguments in config file not matching \
    \nthe excel sheets names:  {name_notfound}. Check spelling!"
    if len(list(name_notfound)) > 0:
        print(message)
    assert len(list(name_notfound)) == 0, message


def fullynumeric(mystring: str) -> bool:
    """
    tests if a string can be converted into float or not.
    e.g. "44" = True, "7A" = False
    """
    try:
        float(mystring)
        return True
    except ValueError:
        return False
    except Exception as e:
        print(f'{e}\n[utils > fullynumeric] {mystring} is not only numeric!')
        return False


def verify_metadata_sample_not_duplicated(metadata_df) -> None:
    def yield_repeated_elems(mylist):
        occur_dic = dict(map(lambda x: (x, list(mylist).count(x)),
                             mylist))  # credits: w3resource.com
        repeated_elems = list()
        for k in occur_dic.keys():
            if occur_dic[k] > 1:
                repeated_elems.append(k)
        return repeated_elems

    sample_duplicated = yield_repeated_elems(list(metadata_df['name_to_plot']))
    if len(sample_duplicated) > 0:
        txt_errors = f"-> duplicated sample names: {sample_duplicated}\n"
        raise ValueError(
            f"Error, found these conflicts in your metadata:\n{txt_errors}")


def isotopologues_meaning_df(isotopologues_full_list: List[str]):
    """
    input: list of isotopologues ['cit_m+0', 'cit_m+1', ...]
       note: extracted from the colnames of the input isotopologues
       (auto-detected any table of isotopologues)
    output: a dataframe in this style:
        metabolite   m+x    isotopologue_name
        cit          m+0    cit_m+0
        cit          m+1    cit_m+1
        ...
        cit          m+6    cit_m+6
        PEP          m+0    PEP_m+0
        ...
    """
    xu = {"metabolite": [], "m+x": [], "isotopologue_name": []}
    for ch in isotopologues_full_list:
        elems = ch.split("_m+")
        xu["metabolite"].append(elems[0])
        xu["m+x"].append("m+{}".format(elems[-1].split("-")[-1]))
        xu["isotopologue_name"].append(ch)
    df = pd.DataFrame.from_dict(xu)
    return df

# from here, functions for isotopologue preview


def save_heatmap_sums_isos(thesums, figuretitle, outputfigure) -> None:
    fig, ax = plt.subplots(figsize=(9, 10))
    sns.heatmap(thesums,
                annot=True, fmt=".1f", cmap="crest",
                square=True,
                annot_kws={
                    'fontsize': 6
                },
                ax=ax)
    plt.xticks(rotation=90)
    plt.title(figuretitle)
    plt.savefig(outputfigure)
    plt.close()


def save_isos_preview(dict_isos_prop, metadata, output_plots_dir,
                      the_boolean_arg, output_extension) -> None:
    if the_boolean_arg:
        for k in metadata['compartment'].unique().tolist():
            df = dict_isos_prop[k]
            samples_co = metadata.loc[
                metadata["compartment"] == k, "original_name"]
            df = df[samples_co]
            df = df.astype(float)
            df = df.assign(metabolite=df.index.str.split("_m+").str[0])
            df = df.assign(isotopologue_type=df.index.str.split("_m+").str[1])
            df["isotopologue_type"] = df['isotopologue_type'].astype(int)

            thesums = compute_sums_isotopol_props(df)

            thesums = thesums.drop(
                columns=['isotopologue_type', 'metabolite'])

            thesums = thesums.astype(float).round(3)
            ff = os.path.join(output_plots_dir, f"sums_Iso_{k}.pdf")
            figuretitle = f"Sums of isotopologue proportions ({k}) "
            save_heatmap_sums_isos(thesums, figuretitle, ff)

            dfmelt = pd.melt(df, id_vars=['metabolite', 'isotopologue_type'])
            dfmelt = impute_custom_levels_to_df(dfmelt)
            table_minimalbymet(dfmelt, os.path.join(output_plots_dir,
                               f"minextremesIso_{k}.{output_extension}"))
            outputfigure = os.path.join(output_plots_dir, f"allsampleIsos_{k}.pdf")
            figtitle = f"{k} compartment, Isotopologues (proportions) \
            across all samples"
            save_rawisos_plot(dfmelt, figuretitle=figtitle,
                              outputfigure=outputfigure)


def impute_custom_levels_to_df(melted: pd.DataFrame):
    """if plotting true, re-level metabolites based on smallest values seen"""
    another = melted.copy()
    another = another.groupby('metabolite').min()
    another = another.sort_values(by='value', ascending=False)
    levelsmetabolites = another.index
    tmp = melted['metabolite']
    melted['metabolite'] = pd.Categorical(tmp, categories=levelsmetabolites)
    return melted


def table_minimalbymet(melted, fileout) -> None:
    try:
        another = melted.copy()
        another = another.groupby('metabolite', observed=False).min()
        another = another.sort_values(by='value', ascending=False)
        another.to_csv(fileout, sep='\t', header=True)
    except Exception as e:
        print(f'[table_minimalbymet > isotopologues preview] {e}. Continue.')
        pass


def save_rawisos_plot(dfmelt, figuretitle, outputfigure) -> None:
    fig, ax = plt.subplots(1, 1, figsize=(16, 10))
    sns.stripplot(ax=ax, data=dfmelt, x="value", y="metabolite", jitter=False,
                  hue="isotopologue_type", size=4, palette="tab20")
    plt.axvline(x=0,
                ymin=0,
                ymax=1,
                linestyle="--", color="gray")
    plt.axvline(x=1,
                ymin=0,
                ymax=1,
                linestyle="--", color="gray")
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.title(figuretitle)
    plt.xlabel("fraction")
    plt.savefig(outputfigure)
    plt.close()

# end functions for isotopologue preview


def abund_divideby_internalStandard(frames_dict, confdict,
                                    internal_standard_df,
                                    use_internal_standard: Union[str, None]):
    # compulsory to subset strict one column as some formats can give > 1
    try:
        internal_standard_df = internal_standard_df[use_internal_standard]
        min_not_zero = internal_standard_df[internal_standard_df > 0].min()
        internal_standard_df.replace(0, min_not_zero, inplace=True)
        # internal_standard_df[internal_standard_df == 0] = \
        #     internal_standard_df[internal_standard_df > 0].min()

        abund_df = frames_dict[confdict['abundances']].copy()
        tmp = abund_df.div(internal_standard_df, axis=1)
        frames_dict[confdict['abundances']] = tmp
    except KeyError:
        print("Error in internal standard normalization, not performed")
        pass

    return frames_dict


def divide_by_amount_material(frames_dict: dict, confdict: dict,
                              material_df: pd.DataFrame,
                              alternative_method: bool,
                              metric: str):
    if material_df is not None:
        abund_df = frames_dict[confdict[metric]].copy()
        if alternative_method:
            material_avg = material_df.iloc[:, 0].mean()
            material_avg_ser = pd.Series([float(material_avg) for i in
                                          range(material_df.shape[0])],
                                         index=material_df.index)
            tmp = abund_df.div(material_df.iloc[:, 0], axis=1)
            tmp = tmp.mul(material_avg_ser, axis=1)
        else:
            tmp = abund_df.div(material_df.iloc[:, 0],  axis=1)

        frames_dict[confdict[metric]] = tmp

    return frames_dict


def drop__metabolites_by_compart(frames_dict: Dict[str, pd.DataFrame],
                                 reverse_available_frames: dict,
                                 bad_metabolites_dict: dict) -> dict:
    for tab_name in frames_dict.keys():
        for co in bad_metabolites_dict.keys():
            if reverse_available_frames[tab_name] in \
                    ["isotopologues", "isotopologue_proportions"]:
                tmpdf = frames_dict[tab_name][co]
                to_drop__isotopologues_names = list()
                for i in list(tmpdf.index):
                    for j in bad_metabolites_dict[co]:
                        if i.startswith(j):
                            to_drop__isotopologues_names.append(i)
                tmpdf = tmpdf.drop(index=to_drop__isotopologues_names)
            else:
                tmpdf = frames_dict[tab_name][co]
                to_drop_now = bad_metabolites_dict[co]
                tmpdf = tmpdf.drop(index=to_drop_now)

            frames_dict[tab_name][co] = tmpdf

    return frames_dict


def transfer__abund_nan__to_all_tables(final_files_names: dict,
                                       frames_dict, metadata):
    """propagates nan from abundance
    # to isotopologues and fractional contributions"""
    # isos_tables = [x for x in [final_files_names['isotopologues'],
    #                            final_files_names['isotopologue_proportions']
    #                            ] if x is not None]
    isos_tables = list()
    for tmp_key in ['isotopologues', 'isotopologue_proportions']:
        try:
            isos_tables.append(final_files_names[tmp_key])
        except KeyError:
            continue
    for co in metadata['compartment'].unique().tolist():
        abu_co = frames_dict[final_files_names['abundances']][co]
        frac_co = frames_dict[final_files_names['mean_enrichment']][co]
        tt = frac_co.mask(abu_co.isnull())

        frames_dict[final_files_names['mean_enrichment']][co] = tt
        # propagation to isotopologues, both prop and absolutes:
        if len(isos_tables) > 0:
            for isoname in isos_tables:
                isoname_df_co = frames_dict[isoname][co]
                tmpfill = list()  # list of dataframes to concatenate
                for metabolite in list(abu_co.index):
                    isoshere = [k for k in list(isoname_df_co.index) if
                                k.startswith(metabolite)]
                    sub_iso_df_co = isoname_df_co.loc[isoshere, :].T
                    sub_iso_df_co = sub_iso_df_co.assign(
                        abu_val=abu_co.loc[metabolite, :].tolist())
                    sub_iso_df_co.loc[sub_iso_df_co['abu_val'].isna(),
                                      :] = np.nan
                    sub_iso_df_co = sub_iso_df_co.drop(columns=['abu_val'])
                    tmpfill.append(sub_iso_df_co.T)
                frames_dict[isoname][co] = pd.concat(tmpfill, axis=0)
    return frames_dict


def compute_sums_isotopol_props(metabolite_isotopologue_names_df):
    sums_df = pd.DataFrame(
        index=metabolite_isotopologue_names_df['metabolite'].unique(),
        columns=metabolite_isotopologue_names_df.columns)
    for metabolite in metabolite_isotopologue_names_df['metabolite'].unique():
        df_sub = metabolite_isotopologue_names_df.loc[
                 metabolite_isotopologue_names_df['metabolite'] == metabolite, :]
        summa = df_sub.sum(axis=0, skipna=False)
        sums_df.loc[metabolite, :] = summa
    return sums_df


def excelsheets2frames_dict(excel_file: str, confdict: dict) -> dict:
    """Extracts data from VIB or generic xlsx files"""
    frames_dict = dict()
    xl = pd.ExcelFile(excel_file)
    sheetsnames = xl.sheet_names
    list_config_tabs = [confdict['abundances'],
                        confdict['mean_enrichment'],
                        confdict['isotopologue_proportions'],
                        confdict['isotopologues']]
    list_config_tabs = [i for i in list_config_tabs if i is not None]

    check_config_and_sheets_match(sheetsnames, list_config_tabs)

    for i in list_config_tabs:
        tmp = pd.read_excel(excel_file, sheet_name=i, engine='openpyxl',
                            header=0, index_col=0)

        badcols = [i for i in list(tmp.columns) if i.startswith("Unnamed")]
        tmp = tmp.loc[:, ~tmp.columns.isin(badcols)]
        tmp.columns = tmp.columns.str.replace(" ", "_")
        tmp.index = tmp.index.str.replace(" ", "_")
        tmp = tmp.replace(" ", "_", regex=False)
        tmp = tmp.dropna(axis=0, how="all")
        frames_dict[i] = tmp

    return frames_dict


# ############ VIB dedicated:


def reshape_vib_data(targetedMetabo_path: str, args, confdict):
    frames_dict = excelsheets2frames_dict(targetedMetabo_path, confdict)
    lod_values, blanks_df, internal_standards_df, bad_x_y = pull_LOD_blanks_IS(
        frames_dict[confdict['abundances']])

    frames_dict = reshape_frames_dict_elems(frames_dict, bad_x_y)
    frames_dict = abund_under_lod_set_nan(
        confdict, frames_dict, lod_values, args.under_detection_limit_set_nan)

    frames_dict = abund_subtract_blankavg(frames_dict, confdict, blanks_df,
                                          args.subtract_blankavg)

    return frames_dict, internal_standards_df


def abund_subtract_blankavg(frames_dict: dict, confdict: dict,
                            blanks_df: pd.Series, subtract_blankavg: bool):
    """on VIB data"""
    abund_df = frames_dict[confdict['abundances']].copy()
    if subtract_blankavg:
        blank_avg = blanks_df.mean(axis=0)
        tmp = abund_df.T.subtract(blank_avg, axis='index')
        tmp[tmp < 0] = 0
        abund_df = tmp.T

        frames_dict[confdict['abundances']] = abund_df

    return frames_dict


def pull_LOD_blanks_IS(abund_df) -> tuple[pd.Series, pd.DataFrame,
                                          pd.DataFrame, dict]:
    """
    Extracts data parts, from VIB total abundance data:
     - Limit of Detection, across samples, into pd Series
     - 'Blanks' values, across variables, into data frame
     - Internal Standard (IS), across samples, , into data frame
    """
    internal_st_precol = tuple()
    pathways_by_vib = list()
    for i in range(len(abund_df.columns)):
        # 00_Internal_Standard, ...., 01_Glycolysis ....
        if "internal_standard" in str(abund_df.columns[i].lower()):
            internal_st_precol = (i, abund_df.columns[i])
        elif re.search(".._", abund_df.columns[i]) and \
                fullynumeric(abund_df.columns[i].split("_")[0]):
            # often '01_Glycolysis' and so on
            pathways_by_vib.append((i, abund_df.columns[i]))

    icolIS = range(internal_st_precol[0] + 1, pathways_by_vib[0][0])
    colIS = [abund_df.columns[i] for i in icolIS]
    internal_standards_df = abund_df[colIS]

    blanks_rows = [i for i in abund_df.index if
                   (i.lower().startswith("blank") or
                    i.lower().startswith("mock"))]
    # synonyms: blank, mock
    blanks_df = abund_df.loc[blanks_rows, :]

    # refine dfs
    elems_x_todrop = [internal_st_precol[1]]
    elems_x_todrop.extend([i[1] for i in pathways_by_vib])
    elems_x_todrop.extend(list(internal_standards_df.columns))
    elems_y_todrop = ['LOD'] + blanks_rows
    # lod_values = lod_values.loc[~lod_values.index.isin(elems_x_todrop)]
    internal_standards_df = internal_standards_df.loc[
        ~internal_standards_df.index.isin(elems_y_todrop)]
    blanks_df = blanks_df.loc[:, ~blanks_df.columns.isin(elems_x_todrop)]

    todrop_x_y = {'x': elems_x_todrop,
                  'y': elems_y_todrop}
    # these x and y just as vib originals (not transposed)
    # * new
    # re-calculate lod_values (ok equal as in VIB excel, verified)
    std_blanks = blanks_df.std(axis=0, ddof=1).multiply(3)
    lod_values = blanks_df.mean() + std_blanks

    return lod_values, blanks_df, internal_standards_df, todrop_x_y


def reshape_frames_dict_elems(frames_dict: dict,
                              todrop_x_y: dict):
    """
    Give proper format to VIB data:
    exclude from each dataframe the rows and columns
    that are specified in the todrop_x_y dictionary
    """
    for k in frames_dict.keys():
        df = frames_dict[k]
        df = df.loc[:, ~df.columns.isin(todrop_x_y["x"])]
        df = df.loc[~df.index.isin(todrop_x_y["y"]), :]
        frames_dict[k] = df

    return frames_dict


def abund_under_lod_set_nan(confdict, frames_dict, lod_values: pd.Series,
                            under_detection_limit_set_nan: bool) -> dict:
    """
    on VIB total abundances, set NaN any value that is below the
    limit of detection (LOD).
    """
    if under_detection_limit_set_nan:
        abund_T = frames_dict[confdict['abundances']].T
        metabolites_below_limit_of_detection = list()
        for i, r in abund_T.iterrows():
            tmp = abund_T.loc[i, :].copy()
            tmp.loc[tmp < lod_values[i]] = np.nan
            if np.isnan(np.array(tmp)).any():
                metabolites_below_limit_of_detection.append(i)
            abund_T.loc[i, :] = tmp

        frames_dict[confdict['abundances']] = abund_T.T

    return frames_dict


def transformmyisotopologues(isos_list, style) -> List[str]:
    """only applies to variables in the VIB and generic formats"""
    if "vib" in style.lower():
        outli = list()
        for ch in isos_list:
            if "_C13-label-" in ch:
                elems = ch.split("_C13-label-")
                metabolite = elems[0]
                species = "m+{}".format(elems[-1].split("-")[-1])
            elif "_PARENT" in ch:
                elems = ch.split("_PARENT")
                metabolite = elems[0]
                species = "m+0"
            else:
                metabolite = ch
                species = "m+?"
            outli.append(metabolite + "_" + species)
    elif "generic" in style.lower():
        try:
            outli = [i.replace("label", "m+") for i in isos_list]
        except Exception as e:
            print(f"{e}\nnot possible to change the isotopologues name style")
            outli = isos_list
    else:
        outli = isos_list
        raise ValueError("isotopologues style not vib nor generic")
    return outli


