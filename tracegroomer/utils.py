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
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import time
from typing import Union


def reset_log_config(logger):
    """parameters for the log file"""
    out_log_file_name = "groom.log"
    logging.basicConfig(encoding='utf-8',
                        stream=sys.stdout,  level=logging.DEBUG)
    file_handler = logging.FileHandler(out_log_file_name)
    logger.addHandler(file_handler)
    return logger


def open_config_file(confifile):
    try:
        with open(confifile, "r") as f:
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
        print(e)
        print('problem with opening metadata file')
        metadata = None
    if metadata is None:
        raise ValueError("\nproblem opening metadata file")


def open_amount_material(amount_material_path: str):
    if amount_material_path is not None:
        try:
            file = amount_material_path
            material_df = pd.read_csv(file, sep='\t', index_col=0)

            assert material_df.shape[1] == 1,\
                "amountMaterial table must have only 2 columns"

            assert (material_df.iloc[:, 0] <= 0).sum() == 0, "amountMaterial table\
                 must not contain zeros nor negative numbers"

        except FileNotFoundError as err_file:
            print(err_file)
        except UnboundLocalError as uerr:
            print(uerr, "config amountMaterial_path:  check spelling")
        except Exception as e:
            print(e)

    return material_df


def compute_abund_from_absolute_isotopol(df, metabos_isos_df):
    """
    input:
       df : input isotopologues in absolute values
       metabos_isos_df : df created internally (see 'isotopologues_meaning_df')
    output:
       transposed dataframe (samples as columns) of total abundances
    """
    df = df.T
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
    return abundance.T


def compute_isotopologues_proportions_from_absolute(df, metabos_isos_df):
    """
    input:
       df : input isotopologues in absolute values
       metabos_isos_df : df created internally (see 'isotopologues_meaning_df')
    output:
       transposed dataframe (samples as columns) of isotopologue proportions
    """
    df = df.T
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
    return isos_prop.T


def compute_MEorFC_from_isotopologues_proportions(df, metabos_isos_df):
    """
    computes mean enrichment (a.k.a fractional contributions)
    input:
      df : isotopologue proportions (whether computed here, or from input)
      metabos_isos_df : df created internally (see 'isotopologues_meaning_df')
    output:
      transposed dataframe (samples as columns) of mean enrichment
    note:
       mean enrichment (i) = ( (i_m+0 * 0) + (i_m+1 * 1) +...+ (i_m+n * n)  ) / n
          where:
           i is one metabolite,
           n the max nb of carbons that can be marked for that metabolite,
           each i_m+x is a value comprised between 0 and 1.
    """
    isos_prop = df.T
    metabos_uniq = metabos_isos_df['metabolite'].unique()
    # empty df
    meanenrich_or_fraccontrib = pd.DataFrame(index=metabos_uniq,
                                             columns=isos_prop.columns)
    for m in metabos_uniq:
        isos_here = metabos_isos_df.loc[
            metabos_isos_df['metabolite'] == m, 'isotopologue_name']
        coefs = [int(i.split("_m+")[1]) for i in isos_here.tolist()] # 0, 1, 2, etc
        sub_df = isos_prop.loc[isos_here, :]
        sub_df['coefs'] = coefs
        # compute the factors, produce another pandas df coefs_fracs_prod
        coefs_fracs_prod = sub_df.multiply(sub_df['coefs'], axis=0)
        # line above makes coefs column be multiplied by itself,
        # TODO : In a future fix this, for now just dropping the coefs col :
        coefs_fracs_prod.drop(columns=['coefs'], inplace=True)
        # sum the factors
        numerator_val = coefs_fracs_prod.sum(axis=0, skipna=False)
        # divide by n, place that scalar in the empty df
        me_fc_this_metabolite = numerator_val / max(coefs)
        me_fc_this_metabolite.name = m
        meanenrich_or_fraccontrib.loc[m, :] = me_fc_this_metabolite
    meanenrich_or_fraccontrib = meanenrich_or_fraccontrib.round(decimals=9)
    return meanenrich_or_fraccontrib.T


def complete_missing_frames(confdict, frames_dict, metadata,
                            metabolites_isos_df) -> dict:
    """can apply to any type of inputs"""
    confdict_new = confdict.copy()

    compartments = metadata['compartment'].unique().tolist()
    if confdict['abundances'] is None:
        if confdict['isotopologues'] is not None:
            frames_dict["abundance_computed"] = dict()
            for co in compartments:
                df_co = frames_dict[confdict['isotopologues']][co]
                tmp_co = compute_abund_from_absolute_isotopol(
                    df_co,
                    metabolites_isos_df)
                frames_dict["abundance_computed"][co] = tmp_co.astype(float)
            confdict_new['abundances'] = "abundance_computed"
        elif confdict['isotopologues'] is None:
            print(" isotopologues' absolute values not available,\
                 impossible to get abundance")
    if confdict['isotopologue_proportions'] is None:
        if confdict['isotopologues'] is not None:
            frames_dict['isotopologues_props_computed'] = dict()
            for co in compartments:
                df_co = frames_dict[confdict['isotopologues']][co]
                tmp_co = compute_isotopologues_proportions_from_absolute(
                    df_co, metabolites_isos_df)
                frames_dict["isotopologues_props_computed"][
                    co] = tmp_co.astype(float)
            confdict_new[
                'isotopologue_proportions'] = "isotopologues_props_computed"
        elif confdict['isotopologues'] is None:
            print(" isotopologues' absolute values not available, \
                impossible to get proportions")
    if confdict['mean_enrichment'] is None:
        try:
            frames_dict["meanEnr_or_FracC_computed"] = dict()
            for co in compartments:
                df_co = frames_dict[confdict_new['isotopologue_proportions']][co]
                tmp_co = compute_MEorFC_from_isotopologues_proportions(
                    df_co,
                    metabolites_isos_df)
                frames_dict["meanEnr_or_FracC_computed"][co] = tmp_co.astype(
                    float)
            confdict_new[
                'mean_enrichment'] = "meanEnr_or_FracC_computed"
        except Exception as e:
            print("impossible to calculate: mean enrichment or fractional \
                  contribution. Isotopologues proportions not found")
            print(e)

    return frames_dict, confdict_new


def df_to__dic_bycomp(df: pd.DataFrame, metadata: pd.DataFrame) -> dict:
    # splits df into dictionary of dataframes, each for one compartment:
    out_dic = dict()
    for co in metadata['compartment'].unique():
        metada_co = metadata.loc[metadata['compartment'] == co, :]
        df_co = df.loc[metada_co['original_name'], :]
        out_dic[co] = df_co
    return out_dic




def isocor_2_frames_dict(isocor_input_df, metadata, confdict,
                        internal_standard: Union[str,None]):
    """ function exclusive for IsoCor type of file"""
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
    me_or_fc = me_or_fc_melted.pivot(index='sample', columns='metabolite')
    me_or_fc = me_or_fc['mean_enrichment'].reset_index()
    me_or_fc = me_or_fc.set_index('sample')

    isos_prop_melted = df[
        ['sample', 'isotopologue_name', 'isotopologue_fraction']]
    isos_prop_melted = isos_prop_melted.drop_duplicates()
    isos_prop = isos_prop_melted.pivot(index='sample',
                                       columns='isotopologue_name')
    isos_prop = isos_prop['isotopologue_fraction'].reset_index()
    isos_prop = isos_prop.set_index('sample')

    isos_absolute_melted = df[
        ['sample', 'isotopologue_name', 'corrected_area']]
    isos_absolute_melted = isos_absolute_melted.drop_duplicates()
    isos_absolute = isos_absolute_melted.pivot(index='sample',
                                               columns='isotopologue_name')
    isos_absolute = isos_absolute['corrected_area'].reset_index()
    isos_absolute = isos_absolute.set_index('sample')

    abundance = compute_abund_from_absolute_isotopol(isos_absolute,
                                                     metabos_isos_df)
    if internal_standard is not None:
        instandard_abun_df = abundance[[internal_standard]]
    else:
        instandard_abun_df = None

    frames_dict = dict()
    frames_dict[confdict['mean_enrichment']] = me_or_fc
    frames_dict[confdict['isotopologue_proportions']] = isos_prop
    frames_dict[confdict['isotopologues']] = isos_absolute
    frames_dict[confdict['abundances']] = abundance

    return frames_dict, instandard_abun_df


def check_config_and_sheets_match(sheetsnames, list_config_tabs):
    name_notfound = set(list_config_tabs) - set(sheetsnames)
    message = f"One or more name_ arguments in config file not matching \
    \nthe excel sheets names:  {name_notfound}. Check spelling!"
    if len(list(name_notfound)) > 0:
        print(message)
    assert len(list(name_notfound)) == 0, message


def fullynumeric(mystring):
    try:
        float(mystring)
        return True
    except ValueError:
        return False
    except Exception as e:
        print(e)
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


def isotopologues_meaning_df(isotopologues_full_list):
    """
    input: list of isotopologues ['cit_m+0', 'cit_m+1', ...]
       note: extracted from the colnames of the input isotopologues (auto-detected any table of isotopologues)
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

def add_metabolite_column(df):
    theindex = df.index
    themetabolites = [i.split("_m+")[0] for i in theindex]
    df = df.assign(metabolite=themetabolites)

    return df


def add_isotopologue_type_column(df):
    theindex = df.index
    preisotopologue_type = [i.split("_m+")[1] for i in theindex]
    theisotopologue_type = [int(i) for i in preisotopologue_type]
    df = df.assign(isotopologue_type=theisotopologue_type)

    return df


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


def save_isos_preview(dic_isos_prop, metadata, output_plots_dir,
                      the_boolean_arg):
    if the_boolean_arg:
        for k in metadata['compartment'].unique().tolist():
            df = dic_isos_prop[k]
            sples_co = metadata.loc[
                metadata["compartment"] == k, "original_name"]
            df = df.loc[sples_co, :]
            df = df.T  # transpose
            df = df.astype(float)
            df = ut.add_metabolite_column(df)
            df = ut.add_isotopologue_type_column(df)

            thesums = compute_sums_isotopol_props(df)

            thesums = thesums.drop(
                columns=['isotopologue_type', 'metabolite'])

            thesums = thesums.astype(float).round(3)
            ff = os.path.join(output_plots_dir, f"sums_Iso_{k}.pdf")
            figuretitle = f"Sums of isotopologue proportions ({k}) "
            ut.save_heatmap_sums_isos(thesums, figuretitle, ff)

            dfmelt = pd.melt(df, id_vars=['metabolite', 'isotopologue_type'])
            dfmelt = ut.givelevels(dfmelt)
            ut.table_minimalbymet(dfmelt,
                                  f"{output_plots_dir}minextremesIso_{k}.csv")
            outputfigure = os.path.join(output_plots_dir ,f"allsampleIsos_{k}.pdf")
            figtitle = f"{k} compartment, Isotopologues (proportions) \
            across all samples"
            ut.save_rawisos_plot(dfmelt, figuretitle=figtitle,
                                 outputfigure=outputfigure)



def givelevels(melted):
    another = melted.copy()
    another = another.groupby('metabolite').min()
    another = another.sort_values(by='value', ascending=False)
    levelsmetabolites = another.index
    tmp = melted['metabolite']
    melted['metabolite'] = pd.Categorical(tmp, categories=levelsmetabolites)

    return melted


def table_minimalbymet(melted, fileout) -> None:
    another = melted.copy()
    another = another.groupby('metabolite').min()
    another = another.sort_values(by='value', ascending=False)
    another.to_csv(fileout, sep='\t', header=True)


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
                                    internal_standards_df,
                                    use_internal_standard: str|None):
    if use_internal_standard is None:
        return frames_dict
    else:
        picked_internal_standard = use_internal_standard
        assert picked_internal_standard in internal_standards_df.columns, "\
               Error, that internal standard is not present in the data"
        abund_dic = frames_dict[confdict['abundances']].copy()
        for compartment in abund_dic.keys():
            inte_sta_co = internal_standards_df.loc[
                          abund_dic[compartment].index, :]
            inte_sta_co_is = inte_sta_co[picked_internal_standard]
            # replace zeros to avoid zero division, uses min of the pd series :
            inte_sta_co_is[inte_sta_co_is == 0] = inte_sta_co_is[
                inte_sta_co_is > 0].min()
            tmp = abund_dic[compartment].div(inte_sta_co_is, axis=0)
            frames_dict[confdict['abundances']][compartment] = tmp
        return frames_dict


def abund_divideby_amount_material(frames_dict: dict, confdict: dict,
                                   amount_material: str,
                                   alternative_method: bool):
    if amount_material is not None:
        try:
            file = amount_material
            material_df = pd.read_csv(file, sep='\t', index_col=0)

            assert material_df.shape[1] == 1,\
                "amountMaterial table must have only 2 columns"

            assert (material_df.iloc[:, 0] <= 0).sum() == 0, "amountMaterial table\
                 must not contain zeros nor negative numbers"

            abund_dic = frames_dict[confdict['abundances']].copy()
            for compartment in abund_dic.keys():
                material_df_s = material_df.loc[
                                list(abund_dic[compartment].index), :]
                if alternative_method:
                    material_avg = material_df_s.iloc[:, 0].mean()
                    material_avg_ser = pd.Series([float(material_avg) for i in
                                                  range(material_df_s.shape[
                                                            0])],
                                                 index=material_df_s.index)
                    tmp = abund_dic[compartment].div(material_df_s.iloc[:, 0],
                                                     axis=0)
                    tmp = tmp.mul(material_avg_ser, axis=0)
                else:
                    tmp = abund_dic[compartment].div(material_df_s.iloc[:, 0],
                                                     axis=0)

                frames_dict[confdict['abundances']][compartment] = tmp

        except FileNotFoundError as err_file:
            print(err_file)
        except UnboundLocalError as uerr:
            print(uerr, "config amountMaterial_path:  check spelling")
        except Exception as e:
            print(e)

    return frames_dict


def isosAbsol_divideby_amount_material(frames_dict: dict, confdict: dict,
                                       material_df: pd.DataFrame,
                                       alternative_method: bool):
    if material_df is not None:
        abund_dic = frames_dict[confdict['isotopologues']].copy()
        for compartment in abund_dic.keys():
            material_df_s = material_df.loc[
                            list(abund_dic[compartment].index), :]
            if alternative_method:
                material_avg = material_df_s.iloc[:, 0].mean()
                material_avg_ser = pd.Series([float(material_avg) for i in
                                              range(material_df_s.shape[
                                                        0])],
                                             index=material_df_s.index)
                tmp = abund_dic[compartment].div(material_df_s.iloc[:, 0],
                                                 axis=0)
                tmp = tmp.mul(material_avg_ser, axis=0)
            else:
                tmp = abund_dic[compartment].div(material_df_s.iloc[:, 0],
                                                 axis=0)

            frames_dict[confdict['isotopologues']][compartment] = tmp
    return frames_dict


def transfer__abund_nan__to_all_tables(confdict, frames_dict, meta_path):
    metadata = ut.open_metadata(meta_path)
    # propagates nan from abundance
    # to isotopologues and fractional contributions
    isos_tables = [s for s in frames_dict.keys() if "isotopol" in s.lower()]
    for co in metadata['compartment'].unique().tolist():
        abu_co = frames_dict[confdict['abundances']][co]
        frac_co = frames_dict[confdict['mean_enrichment']][co]
        tt = frac_co.mask(abu_co.isnull())
        frames_dict[confdict['mean_enrichment']][co] = tt
        # propagation to isotopologues, both prop and absolutes:
        for isoname in isos_tables:
            isoname_df_co = frames_dict[isoname][co]
            tmpfill = list()
            for metabolite in abu_co.columns:
                isoshere = [k for k in isoname_df_co if
                            k.startswith(metabolite)]
                sub_iso_df_co = isoname_df_co[isoshere]
                sub_iso_df_co = sub_iso_df_co.assign(
                    abu_val=abu_co[metabolite])
                sub_iso_df_co.loc[sub_iso_df_co['abu_val'].isna(), :] = np.nan
                sub_iso_df_co = sub_iso_df_co.drop(columns=['abu_val'])
                tmpfill.append(sub_iso_df_co)
            frames_dict[isoname][co] = pd.concat(tmpfill, axis=1)
    return frames_dict


def compute_sums_isotopol_props(dfT):
    sums_df = pd.DataFrame(index=dfT['metabolite'].unique(),
                           columns=dfT.columns)
    for metabolite in dfT['metabolite'].unique():
        df_sub = dfT.loc[dfT['metabolite'] == metabolite, :]
        summa = df_sub.sum(axis=0, skipna=False)
        sums_df.loc[metabolite, :] = summa
    return sums_df


# ############ VIB dedicated:

def abund_subtract_blankavg(frames_dict: dict, confdict: dict,
                            blanks_df: pd.Series, subtract_blankavg: bool):
    """on VIB data"""
    abund_dic = frames_dict[confdict['abundances']].copy()
    if subtract_blankavg:
        for compartment in abund_dic.keys():
            blanks_df_s = blanks_df[list(abund_dic[compartment].columns)]
            blanksAvg_s = blanks_df_s.mean(axis=0)
            abu_df_T = abund_dic[compartment].T
            tmp = abu_df_T.subtract(blanksAvg_s, axis='index')
            tmp[tmp < 0] = 0
            abund_dic[compartment] = tmp.T

        frames_dict[confdict['abundances']] = abund_dic

    return frames_dict

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


def auto_drop_metabolites_uLOD(confdict, frames_dict, metadata, lod_values,
                               auto_drop_metabolite_LOD_based: bool) -> dict:
    """Applied on VIB data as it provides Limit of Detection"""
    # affects all the datasets in frames_dict
    auto_bad_metabolites = dict()
    compartments = metadata['compartment'].unique().tolist()
    for k in compartments:
        auto_bad_metabolites[k] = list()

    if auto_drop_metabolite_LOD_based:
        # drop metabolite if all its values are under LOD
        for co in compartments:
            abund_co = frames_dict[confdict['abundances']][co]
            abund_coT = abund_co.T
            for i, r in abund_coT.iterrows():
                nb_nan = abund_coT.loc[i, :].isna().sum()
                nb_under_LOD = (abund_coT.loc[i, :] < lod_values[i]).sum()
                if (nb_under_LOD == r.size) or (nb_nan == r.size):
                    auto_bad_metabolites[co].append(i)
        frames_dict = drop__metabolites_by_compart(frames_dict,
                                                  auto_bad_metabolites)

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
                ut.fullynumeric(abund_df.columns[i].split("_")[0]):
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
    # lod_values = abund_df.loc['LOD', :]

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


def reshape_frames_dict_elems(frames_dict: dict, metadata: pd.DataFrame,
                             todrop_x_y: dict):
    """
    Give proper format to VIB data:
    exclude from each dataframe the rows and columns
    that are specified in the todrop_x_y dictionary
    """
    trans_dic = dict()
    for k in frames_dict.keys():
        df = frames_dict[k]
        df = df.loc[:, ~df.columns.isin(todrop_x_y["x"])]
        df = df.loc[~df.index.isin(todrop_x_y["y"]), :]
        compartments = metadata['compartment'].unique().tolist()
        trans_dic[k] = dict()
        for co in compartments:
            metada_co = metadata.loc[metadata['compartment'] == co, :]
            df_co = df.loc[metada_co['original_name'], :]
            trans_dic[k][co] = df_co

    frames_dict = trans_dic.copy()
    return frames_dict



def abund_under_lod_set_nan(confdict, frames_dict, metadata,
                            lod_values,
                            under_detection_limit_set_nan) -> dict:
    """on VIB data """
    if under_detection_limit_set_nan:
        for co in metadata['compartment'].unique().tolist():
            abund_co = frames_dict[confdict['abundances']][co]
            abund_coT = abund_co.T
            for i, r in abund_coT.iterrows():
                # avoid future error "FutureWarning: ChainedAssignmentError":
                tmp = abund_coT.loc[i, :].copy()
                tmp.loc[tmp < lod_values[i]] = np.nan
                abund_coT.loc[i, :] = tmp
            frames_dict[confdict['abundances']][co] = abund_coT.T

    return frames_dict


def do_vib_prep(meta_path, targetedMetabo_path, args, confdict,
                amount_mater_path, output_plots_dir):
    # the order of the steps is the one recommended by VIB
    frames_dict = excelsheets2frames_dict(targetedMetabo_path, confdict)
    metadata = ut.open_metadata(meta_path)
    # ut.verify_metadata_sample_not_duplicated(metadata)  # TODO del
    abundance_df = frames_dict[confdict['abundances']]
    lod_values, blanks_df, internal_standards_df, bad_x_y = pull_LOD_blanks_IS(
        abundance_df)

    frames_dict = ut.reshape_frames_dict_elems(frames_dict, metadata, bad_x_y)

    frames_dict = ut.abund_under_lod_set_nan(confdict, frames_dict, metadata,
                                         lod_values,
                                         args.under_detection_limit_set_nan)

    frames_dict = ut.auto_drop_metabolites_uLOD(confdict, frames_dict, metadata,
                                            lod_values, args.
                                            auto_drop_metabolite_LOD_based)
    frames_dict = ut.abund_subtract_blankavg(frames_dict, confdict,
                                         blanks_df, args.subtract_blankavg)

    arg_alt_div_amount_material = args.alternative_div_amount_material
    frames_dict = abund_divideby_amount_material(frames_dict, confdict,  # VIB
                                                amount_mater_path,
                                                arg_alt_div_amount_material)

    frames_dict = abund_divideby_internalStandard(frames_dict, confdict,
                                                 internal_standards_df,
                                                 args.use_internal_standard)

    # transform isotopologues names to the easier "m+x" style:
    for tab in frames_dict.keys():
        if "isotopol" in tab.lower():
            for co in frames_dict[tab]:
                tmp = frames_dict[tab][co]
                new_col = transformmyisotopologues(tmp.columns, "vib")
                tmp.columns = new_col
                frames_dict[tab][co] = tmp
    # end for
    save_isos_preview(frames_dict[confdict['isotopologue_proportions']],
                      metadata,
                      output_plots_dir, args.isotopologues_preview)

    return frames_dict


def transformmyisotopologues(isos_list, style):
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
            print(e)
            print("not possible to change the isotopologues name style")
            outli = isos_list
    else:
        outli = isos_list
        raise ValueError("isotopologues style not vib nor generic")
    return outli



# useful resources:
# count nb of occurrences:
# https://www.w3resource.com/python-exercises/lambda/python-lambda-exercise-49.php
