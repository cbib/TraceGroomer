#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Prepare tables,
user defined config file and args

@author: johanna
"""
import os
import pandas as pd
import numpy as np
import re
import tracegroomer.utils as fg


def excelsheets2frames_dic(excel_file: str, confidic: dict) -> dict:
    frames_dic = dict()
    xl = pd.ExcelFile(excel_file)
    sheetsnames = xl.sheet_names
    list_config_tabs = [confidic['abundances'],
                        confidic['mean_enrichment'],
                        confidic['isotopologue_proportions'],
                        confidic['isotopologues']]
    list_config_tabs = [i for i in list_config_tabs if i is not None]

    def check_config_and_sheets_match(sheetsnames, list_config_tabs):
        name_notfound = set(list_config_tabs) - set(sheetsnames)
        message = f"One or more name_ arguments in config file not matching \
        \nthe excel sheets names:  {name_notfound}. Check spelling!"
        if len(list(name_notfound)) > 0:
            print(message)
        assert len(list(name_notfound)) == 0, message

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
        frames_dic[i] = tmp

    return frames_dic


def pull_LOD_blanks_IS(abund_df) -> tuple[pd.Series, pd.DataFrame,
                                          pd.DataFrame, dict]:
    internal_st_precol = tuple()
    pathways_by_vib = list()
    for i in range(len(abund_df.columns)):
        # 00_Internal_Standard, ...., 01_Glycolysis ....
        if "internal_standard" in str(abund_df.columns[i].lower()):
            internal_st_precol = (i, abund_df.columns[i])
        elif re.search(".._", abund_df.columns[i]) and \
                fg.fullynumeric(abund_df.columns[i].split("_")[0]):
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


def reshape_frames_dic_elems(frames_dic: dict, metadata: pd.DataFrame,
                             todrop_x_y: dict):
    trans_dic = dict()
    for k in frames_dic.keys():
        df = frames_dic[k]
        df = df.loc[:, ~df.columns.isin(todrop_x_y["x"])]
        df = df.loc[~df.index.isin(todrop_x_y["y"]), :]
        compartments = metadata['compartment'].unique().tolist()
        trans_dic[k] = dict()
        for co in compartments:
            metada_co = metadata.loc[metadata['compartment'] == co, :]
            df_co = df.loc[metada_co['original_name'], :]
            trans_dic[k][co] = df_co

    frames_dic = trans_dic.copy()
    return frames_dic


def abund_under_lod_set_nan(confidic, frames_dic, metadata,
                            lod_values,
                            under_detection_limit_set_nan) -> dict:
    if under_detection_limit_set_nan:
        for co in metadata['compartment'].unique().tolist():
            abund_co = frames_dic[confidic['abundances']][co]
            abund_coT = abund_co.T
            for i, r in abund_coT.iterrows():
                # avoid future error "FutureWarning: ChainedAssignmentError":
                tmp = abund_coT.loc[i, :].copy()
                tmp.loc[tmp < lod_values[i]] = np.nan
                abund_coT.loc[i, :] = tmp
            frames_dic[confidic['abundances']][co] = abund_coT.T

    return frames_dic


def drop__metabolites_by_compart(frames_dic_orig: dict,
                                 bad_metabolites_dic: dict) -> dict:
    frames_dic = frames_dic_orig.copy()
    for tab_name in frames_dic.keys():
        for co in bad_metabolites_dic.keys():

            if "isotopolog" in tab_name.lower():
                tmpdf = frames_dic[tab_name][co]
                to_drop_now_isos = list()
                for i in tmpdf.columns:
                    for j in bad_metabolites_dic[co]:
                        if i.startswith(j):
                            to_drop_now_isos.append(i)
                tmpdf = tmpdf.drop(columns=to_drop_now_isos)
                frames_dic[tab_name][co] = tmpdf

            elif "isotopolog" not in tab_name.lower():
                tmpdf = frames_dic[tab_name][co]
                to_drop_now = bad_metabolites_dic[co]
                tmpdf = tmpdf.drop(columns=to_drop_now)
                frames_dic[tab_name][co] = tmpdf

    return frames_dic


def auto_drop_metabolites_uLOD(confidic, frames_dic, metadata, lod_values,
                               auto_drop_metabolite_LOD_based: bool) -> dict:
    # affects all the datasets in frames_dic
    auto_bad_metabolites = dict()
    compartments = metadata['compartment'].unique().tolist()
    for k in compartments:
        auto_bad_metabolites[k] = list()

    if auto_drop_metabolite_LOD_based:
        # drop metabolite if all its values are under LOD
        for co in compartments:
            abund_co = frames_dic[confidic['abundances']][co]
            abund_coT = abund_co.T
            for i, r in abund_coT.iterrows():
                nb_nan = abund_coT.loc[i, :].isna().sum()
                nb_under_LOD = (abund_coT.loc[i, :] < lod_values[i]).sum()
                if (nb_under_LOD == r.size) or (nb_nan == r.size):
                    auto_bad_metabolites[co].append(i)
        frames_dic = drop__metabolites_by_compart(frames_dic,
                                                  auto_bad_metabolites)

    return frames_dic


def abund_subtract_blankavg(frames_dic: dict, confidic: dict,
                            blanks_df: pd.Series, subtract_blankavg: bool):
    abund_dic = frames_dic[confidic['abundances']].copy()
    if subtract_blankavg:
        for compartment in abund_dic.keys():
            blanks_df_s = blanks_df[list(abund_dic[compartment].columns)]
            blanksAvg_s = blanks_df_s.mean(axis=0)
            abu_df_T = abund_dic[compartment].T
            tmp = abu_df_T.subtract(blanksAvg_s, axis='index')
            tmp[tmp < 0] = 0
            abund_dic[compartment] = tmp.T

        frames_dic[confidic['abundances']] = abund_dic

    return frames_dic


def abund_divideby_amount_material(frames_dic: dict, confidic: dict,
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

            abund_dic = frames_dic[confidic['abundances']].copy()
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

                frames_dic[confidic['abundances']][compartment] = tmp

        except FileNotFoundError as err_file:
            print(err_file)
        except UnboundLocalError as uerr:
            print(uerr, "config amountMaterial_path:  check spelling")
        except Exception as e:
            print(e)

    return frames_dic


def isosAbsol_divideby_amount_material(frames_dic: dict, confidic: dict,
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

            abund_dic = frames_dic[confidic['isotopologues']].copy()
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

                frames_dic[confidic['isotopologues']][compartment] = tmp

        except FileNotFoundError as err_file:
            print(err_file)
        except UnboundLocalError as uerr:
            print(uerr, "config amountMaterial_path:  check spelling")
        except Exception as e:
            print(e)

    return frames_dic


def abund_divideby_internalStandard(frames_dic, confidic,
                                    internal_standards_df,
                                    use_internal_standard: str|None):
    if use_internal_standard is None:
        return frames_dic
    else:
        picked_internal_standard = use_internal_standard
        assert picked_internal_standard in internal_standards_df.columns, "\
               Error, that internal standard is not present in the data"
        abund_dic = frames_dic[confidic['abundances']].copy()
        for compartment in abund_dic.keys():
            inte_sta_co = internal_standards_df.loc[
                          abund_dic[compartment].index, :]
            inte_sta_co_is = inte_sta_co[picked_internal_standard]
            # replace zeros to avoid zero division, uses min of the pd series :
            inte_sta_co_is[inte_sta_co_is == 0] = inte_sta_co_is[
                inte_sta_co_is > 0].min()
            tmp = abund_dic[compartment].div(inte_sta_co_is, axis=0)
            frames_dic[confidic['abundances']][compartment] = tmp
        return frames_dic


def transformmyisotopologues(isos_list, style):
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


def compute_sums_isotopol_props(dfT):
    sums_df = pd.DataFrame(index=dfT['metabolite'].unique(),
                           columns=dfT.columns)
    for metabolite in dfT['metabolite'].unique():
        df_sub = dfT.loc[dfT['metabolite'] == metabolite, :]
        summa = df_sub.sum(axis=0, skipna=False)
        sums_df.loc[metabolite, :] = summa
    return sums_df


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
            df = fg.add_metabolite_column(df)
            df = fg.add_isotopologue_type_column(df)

            thesums = compute_sums_isotopol_props(df)

            thesums = thesums.drop(
                columns=['isotopologue_type', 'metabolite'])

            thesums = thesums.astype(float).round(3)
            ff = os.path.join(output_plots_dir, f"sums_Iso_{k}.pdf")
            figuretitle = f"Sums of isotopologue proportions ({k}) "
            fg.save_heatmap_sums_isos(thesums, figuretitle, ff)

            dfmelt = pd.melt(df, id_vars=['metabolite', 'isotopologue_type'])
            dfmelt = fg.givelevels(dfmelt)
            fg.table_minimalbymet(dfmelt,
                                  f"{output_plots_dir}minextremesIso_{k}.csv")
            outputfigure = os.path.join(output_plots_dir ,f"allsampleIsos_{k}.pdf")
            figtitle = f"{k} compartment, Isotopologues (proportions) \
            across all samples"
            fg.save_rawisos_plot(dfmelt, figuretitle=figtitle,
                                 outputfigure=outputfigure)


def do_vib_prep(meta_path, targetedMetabo_path, args, confidic,
                amount_mater_path, output_plots_dir):
    # the order of the steps is the one recommended by VIB
    frames_dic = excelsheets2frames_dic(targetedMetabo_path, confidic)
    metadata = fg.open_metadata(meta_path)
    fg.verify_metadata_sample_not_duplicated(metadata)
    abundance_df = frames_dic[confidic['abundances']]
    lod_values, blanks_df, internal_standards_df, bad_x_y = pull_LOD_blanks_IS(
        abundance_df)

    frames_dic = reshape_frames_dic_elems(frames_dic, metadata, bad_x_y)

    frames_dic = abund_under_lod_set_nan(confidic, frames_dic, metadata,
                                         lod_values,
                                         args.under_detection_limit_set_nan)

    frames_dic = auto_drop_metabolites_uLOD(confidic, frames_dic, metadata,
                                            lod_values, args.
                                            auto_drop_metabolite_LOD_based)
    frames_dic = abund_subtract_blankavg(frames_dic, confidic,
                                         blanks_df, args.subtract_blankavg)

    arg_alt_div_amount_material = args.alternative_div_amount_material
    frames_dic = abund_divideby_amount_material(frames_dic, confidic,
                                                amount_mater_path,
                                                arg_alt_div_amount_material)

    frames_dic = abund_divideby_internalStandard(frames_dic, confidic,
                                                 internal_standards_df,
                                                 args.use_internal_standard)

    # transform isotopologues names to the easier "m+x" style:
    for tab in frames_dic.keys():
        if "isotopol" in tab.lower():
            for co in frames_dic[tab]:
                tmp = frames_dic[tab][co]
                new_col = transformmyisotopologues(tmp.columns, "vib")
                tmp.columns = new_col
                frames_dic[tab][co] = tmp
    # end for
    save_isos_preview(frames_dic[confidic['isotopologue_proportions']],
                      metadata,
                      output_plots_dir, args.isotopologues_preview)

    return frames_dic


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


def complete_missing_frames(confidic, frames_dic, metadata,
                            metabolites_isos_df) -> dict:
    confidic_new = confidic.copy()

    compartments = metadata['compartment'].unique().tolist()
    if confidic['abundances'] is None:
        if confidic['isotopologues'] is not None:
            frames_dic["abundance_computed"] = dict()
            for co in compartments:
                df_co = frames_dic[confidic['isotopologues']][co]
                tmp_co = compute_abund_from_absolute_isotopol(
                    df_co,
                    metabolites_isos_df)
                frames_dic["abundance_computed"][co] = tmp_co.astype(float)
            confidic_new['abundances'] = "abundance_computed"
        elif confidic['isotopologues'] is None:
            print(" isotopologues' absolute values not available,\
                 impossible to get abundance")
    if confidic['isotopologue_proportions'] is None:
        if confidic['isotopologues'] is not None:
            frames_dic['isotopologues_props_computed'] = dict()
            for co in compartments:
                df_co = frames_dic[confidic['isotopologues']][co]
                tmp_co = compute_isotopologues_proportions_from_absolute(
                    df_co, metabolites_isos_df)
                frames_dic["isotopologues_props_computed"][
                    co] = tmp_co.astype(float)
            confidic_new[
                'isotopologue_proportions'] = "isotopologues_props_computed"
        elif confidic['isotopologues'] is None:
            print(" isotopologues' absolute values not available, \
                impossible to get proportions")
    if confidic['mean_enrichment'] is None:
        try:
            frames_dic["meanEnr_or_FracC_computed"] = dict()
            for co in compartments:
                df_co = frames_dic[confidic_new['isotopologue_proportions']][co]
                tmp_co = compute_MEorFC_from_isotopologues_proportions(
                    df_co,
                    metabolites_isos_df)
                frames_dic["meanEnr_or_FracC_computed"][co] = tmp_co.astype(
                    float)
            confidic_new[
                'mean_enrichment'] = "meanEnr_or_FracC_computed"
        except Exception as e:
            print("impossible to calculate: mean enrichment or fractional \
                  contribution. Isotopologues proportions not found")
            print(e)

    return frames_dic, confidic_new


def df_to__dic_bycomp(df: pd.DataFrame, metadata: pd.DataFrame) -> dict:
    # splits df into dictionary of dataframes, each for one compartment:
    out_dic = dict()
    for co in metadata['compartment'].unique():
        metada_co = metadata.loc[metadata['compartment'] == co, :]
        df_co = df.loc[metada_co['original_name'], :]
        out_dic[co] = df_co
    return out_dic


def do_isocorOutput_prep(meta_path, targetedMetabo_path, args, confidic,
                         amount_mater_path, output_plots_dir):
    def isocor_2_frames_dic(isocor_out_df, metadata, confidic,
                            internal_standard: str|None):
        df = isocor_out_df[
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

        frames_dic = dict()
        frames_dic[confidic['mean_enrichment']] = df_to__dic_bycomp(
            me_or_fc, metadata)
        frames_dic[confidic['isotopologue_proportions']] = df_to__dic_bycomp(
            isos_prop, metadata)
        frames_dic[confidic['isotopologues']] = df_to__dic_bycomp(
            isos_absolute, metadata)
        frames_dic[confidic['abundances']] = df_to__dic_bycomp(abundance,
                                                                   metadata)

        return frames_dic, instandard_abun_df

    isocor_output_df = pd.read_csv(targetedMetabo_path, sep="\t")
    metadata = fg.open_metadata(meta_path)
    fg.verify_metadata_sample_not_duplicated(metadata)
    frames_dic, instandard_abun_df = isocor_2_frames_dic(
        isocor_output_df,
        metadata, confidic,
        args.use_internal_standard)
    # at this point isotopologues are in "m+x" format
    frames_dic = abund_divideby_internalStandard(frames_dic, confidic,
                                                 instandard_abun_df,
                                                 args.use_internal_standard)
    frames_dic = abund_divideby_amount_material(
        frames_dic, confidic,
        amount_mater_path,
        args.alternative_div_amount_material)

    # # I tested divide absol isotopols by amount mater => values too small !
    # frames_dic = isosAbsol_divideby_amount_material(
    #     frames_dic, confidic,
    #     amount_mater_path,
    #     args.alternative_div_amount_material)

    save_isos_preview(frames_dic[confidic['isotopologue_proportions']],
                      metadata,
                      output_plots_dir, args.isotopologues_preview)

    return frames_dic


def do_generic_prep(meta_path, targetedMetabo_path, args, confidic,
                    amount_mater_path, output_plots_dir):
    metadata = fg.open_metadata(meta_path)
    fg.verify_metadata_sample_not_duplicated(metadata)
    frames_dic = excelsheets2frames_dic(targetedMetabo_path, confidic)
    tabs_isotopologues = [s for s in frames_dic.keys() if
                          "isotopol" in s.lower()]
    assert len(
        tabs_isotopologues) >= 1, "\nError, bad or no isotopologues input"
    for tab in tabs_isotopologues:  # tabs are not split by compartment here
        tmp = frames_dic[tab]
        new_col = transformmyisotopologues(tmp.columns, "generic")
        tmp.columns = new_col
        frames_dic[tab] = tmp

    # end for

    isotopologues_full = frames_dic[tabs_isotopologues[0]].columns
    metabolites_isos_df = fg.isotopologues_meaning_df(isotopologues_full)

    for k in frames_dic.keys():
        tmp_co_dic = df_to__dic_bycomp(frames_dic[k],
                                       metadata)  # split by compartment
        frames_dic[k] = tmp_co_dic

    frames_dic, confidic_new = complete_missing_frames(confidic, frames_dic,
                                                       metadata,
                                                       metabolites_isos_df)

    if args.use_internal_standard is not None:
        instandard_abun_l = list()
        abu = confidic_new['abundances']
        for co in frames_dic[abu].keys():
            tmp_co = frames_dic[abu][co][args.use_internal_standard]
            instandard_abun_l.append(
                pd.DataFrame({args.use_internal_standard: tmp_co}))
        instandard_abun_df = pd.concat(instandard_abun_l)

        frames_dic = abund_divideby_internalStandard(
            frames_dic, confidic_new,
            instandard_abun_df,
            args.use_internal_standard)
    # end if

    frames_dic = abund_divideby_amount_material(
        frames_dic, confidic_new,
        amount_mater_path,
        args.alternative_div_amount_material)

    # # TODO : think divide or not divide absolute isotopologues
    # frames_dic = isosAbsol_divideby_amount_material(
    #     frames_dic, confidic,
    #     amount_mater_path,
    #     args.alternative_div_amount_material)

    save_isos_preview(frames_dic[confidic_new['isotopologue_proportions']],
                      metadata,
                      output_plots_dir, args.isotopologues_preview)

    return frames_dic, confidic_new


def drop_metabolites_infile(frames_dic, exclude_list_file: str|None):
    if exclude_list_file is not None:
        print("removing metabolites as specified by user in file:")
        print(exclude_list_file, "\n")
        exclude_df = pd.read_csv(exclude_list_file, sep="\t", header=0)
        try:
            unwanted_metabolites = dict()
            for co in exclude_df["compartment"].unique():
                mets_l = exclude_df.loc[
                    exclude_df["compartment"] == co, 'metabolite'].tolist()
                unwanted_metabolites[co] = mets_l
            frames_dic = drop__metabolites_by_compart(frames_dic,
                                                      unwanted_metabolites)
        except FileNotFoundError as err_file:
            print(err_file, "ooo")
        except Exception as e:
            print(e, "iii")
    return frames_dic


def frames_filterby_min_admited_isosprop(frames_dic, confidic,
                                         isosprop_min_admitted: float):
    isos_propor_dic = frames_dic[confidic['isotopologue_proportions']]
    bad_mets = dict()
    for co in isos_propor_dic.keys():
        tmp = isos_propor_dic[co]
        series_bool = tmp.le(isosprop_min_admitted).any()
        isos_bad = series_bool[series_bool]
        set_mets = set([i.split("_m+")[0] for i in isos_bad.index])
        bad_mets[co] = list(set_mets)

    frames_dic = drop__metabolites_by_compart(frames_dic, bad_mets)

    return frames_dic


def isosprop_stomp_values(frames_dic, confidic, isosprop_stomp_vals: bool):
    if isosprop_stomp_vals:
        isos_propor_dic = frames_dic[confidic['isotopologue_proportions']]
        for co in isos_propor_dic.keys():
            df = isos_propor_dic[co]
            df[df < 0] = 0
            df[df > 1] = 1
            isos_propor_dic[co] = df
        frames_dic[confidic['isotopologue_proportions']] = isos_propor_dic
    return frames_dic


def meanenrich_or_fracfontrib_stomp_values(
        frames_dic, confidic, meanenri_or_fraccontrib_stomp_vals: bool):
    if meanenri_or_fraccontrib_stomp_vals:
        meorfc_dic = frames_dic[confidic['mean_enrichment']]
        for co in meorfc_dic.keys():
            df = meorfc_dic[co]
            df[df < 0] = 0
            df[df > 1] = 1
            meorfc_dic[co] = df
        frames_dic[confidic['mean_enrichment']] = meorfc_dic
    return frames_dic


def transfer__abund_nan__to_all_tables(confidic, frames_dic, meta_path):
    metadata = fg.open_metadata(meta_path)
    # propagates nan from abundance
    # to isotopologues and fractional contributions
    isos_tables = [s for s in frames_dic.keys() if "isotopol" in s.lower()]
    for co in metadata['compartment'].unique().tolist():
        abu_co = frames_dic[confidic['abundances']][co]
        frac_co = frames_dic[confidic['mean_enrichment']][co]
        tt = frac_co.mask(abu_co.isnull())
        frames_dic[confidic['mean_enrichment']][co] = tt
        # propagation to isotopologues, both prop and absolutes:
        for isoname in isos_tables:
            isoname_df_co = frames_dic[isoname][co]
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
            frames_dic[isoname][co] = pd.concat(tmpfill, axis=1)
    return frames_dic


def perform_type_prep(args, confidic, meta_path, targetedMetabo_path,
                      amount_mater_path, groom_out_path) -> None:
    output_plots_dir = os.path.join(groom_out_path, "preview_plots")
    if args.isotopologues_preview:
        fg.detect_and_create_dir(output_plots_dir)

    fg.detect_and_create_dir(groom_out_path) # remove in future version

    if args.type_of_file == 'IsoCor_out_tsv':
        frames_dic = do_isocorOutput_prep(meta_path, targetedMetabo_path,
                                          args, confidic,
                                          amount_mater_path, output_plots_dir)
    elif args.type_of_file == 'VIBMEC_xlsx':
        frames_dic = do_vib_prep(meta_path, targetedMetabo_path, args,
                                 confidic,
                                 amount_mater_path, output_plots_dir)
    elif args.type_of_file == 'generic_xlsx':
        frames_dic, confidic = do_generic_prep(meta_path, targetedMetabo_path,
                                               args, confidic,
                                               amount_mater_path,
                                               output_plots_dir)

    # common steps to any preparation type:
    frames_dic = drop_metabolites_infile(frames_dic,
                                         args.remove_these_metabolites)
    frames_dic = frames_filterby_min_admited_isosprop(
        frames_dic, confidic,
        args.isosprop_min_admitted)
    frames_dic = isosprop_stomp_values(frames_dic, confidic,
                                       args.isosprop_stomp_values)
    frames_dic = meanenrich_or_fracfontrib_stomp_values(
        frames_dic, confidic,
        args.meanenrich_or_fracfontrib_stomp_values)
    frames_dic = transfer__abund_nan__to_all_tables(confidic, frames_dic,
                                                    meta_path)

    for k in frames_dic.keys():
        tmpli = list()
        for compartment in frames_dic[k].keys():
            tmp = frames_dic[k][compartment].T
            tmpli.append(tmp)

        finalk = tmpli[0]
        for i in range(1, len(tmpli)):
            finalk = pd.merge(finalk, tmpli[i], how='outer',
                              left_index=True, right_index=True)
        # note : do not clear zero rows, as gives problem pd.merge
        finalk.index.name = "ID"
        finalk = finalk.reset_index()
        finalk = finalk.drop_duplicates()
        finalk.to_csv(
            os.path.join(groom_out_path, f"{k}.csv"),
            sep='\t', header=True, index=False)
        print("File saved to:",  os.path.join(groom_out_path, f"{k}.csv"))


