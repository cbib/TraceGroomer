#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 12:24:01 2022

@author: johanna
"""

import os
import yaml
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import locale


def open_config_file(confifile):
    try:
        with open(confifile, "r") as f:
            confidic = yaml.load(f, Loader=yaml.Loader)
    except yaml.YAMLError as yam_err:
        print(yam_err)
        confidic = None
    except Exception as e:
        print(e)
        confidic = None

    if confidic is None:
        raise ValueError("\nimpossible to read configuration file")

    return confidic


def detect_and_create_dir(namenesteddir):
    if not os.path.exists(namenesteddir):
        os.makedirs(namenesteddir)


def fullynumeric(mystring):
    try:
        float(mystring)
        return True
    except ValueError:
        return False
    except Exception as e:
        print(e)
        return False


def open_metadata(file_path):
    try:
        metadata = pd.read_csv(file_path, sep='\t')
        return metadata
    except Exception as e:
        print(e)
        print('problem with opening metadata file')
        metadata = None
    if metadata is None:
        raise ValueError("\nproblem opening configuration file")


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

# END

# useful resources:
# count nb of occurrences:
# https://www.w3resource.com/python-exercises/lambda/python-lambda-exercise-49.php
