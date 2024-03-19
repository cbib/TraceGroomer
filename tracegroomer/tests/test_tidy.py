from unittest import TestCase

import pandas as pd
import numpy as np

from tracegroomer.tidy import CompositeData


class MimicArg:
    def __init__(self):
        self.fractions_stomp_values = True
        self.use_internal_standard = "W_acid"


class TestCompositeData(TestCase):

    def test_transpose_frames(self):
        df = pd.DataFrame({
            "M0458_label0": [5769.87, 6224.543, 7787.676, 86457.5],
            "M0458_label1": [2879.3, 3676.4, 4023, 4684.4],
            "M0458_label2": [1956, 1845.4, 1005.2, 1309],
            "M0777_label0": [2879.3, 3676.4, 4023, 4684.4],
            "M0777_label1": [79.3, 76.4, 23.56, 44.34],
        })
        df.index = ["sampleA1", "sampleA2", "sampleB1", "sampleB2"]
        myobj = CompositeData("generic_xlsx")
        myobj.frames_dict = {"utopik-info": df}
        myobj.transpose_frames()
        self.assertListEqual(
            list(myobj.frames_dict["utopik-info"].columns),
            ["sampleA1", "sampleA2", "sampleB1", "sampleB2"])
        self.assertTrue(np.all(
            np.array(
                myobj.frames_dict["utopik-info"]["sampleB2"]) == np.array(
                [86457.5, 4684.4, 1309, 4684.4, 44.34]
            )
        ))

    def test_fill_missing_data(self):
        df = pd.DataFrame({
            "sampleA1": [5769.87, 2879.3, 1956.0, 2879.3, 79.3],
            "sampleA2": [6224.543, 3676.4, 1845.4, 3676.4, 76.4],
            "sampleB1": [7787.676, 4023.0, 1005.2, 4023.0, 23.56],
            "sampleB2": [86457.5, 4684.4, 1309.0, 4684.4, 44.34]
        })
        df.index = ["id009", "id010", "id011", "id019", "id020"]
        var_df = pd.DataFrame({
            "ID": ["id009", "id010", "id011", "id019", "id020"],
            "compound_name": ["acCoA", "acCoA", "acCoA",
                              "unknown", "unknown"],
            "isotope": [0, 1, 2, 0, 1]})
        df.index = var_df["compound_name"].str.cat(
            var_df["isotope"].astype(str), sep="_m+")
        myobj = CompositeData("rule_tsv")
        myobj.frames_dict = {"utopik-info": df}
        config = {"isotopologues": "utopik-info", "mean_enrichment": None,
                  "abundances": None, "isotopologue_proportions": None}
        myobj.load_metabolite_to_isotopologue_df(config)
        myobj.fill_missing_data(config)
        result1 = myobj.frames_dict["abundances_computed"]
        result2 = myobj.frames_dict["mean_enrichment_computed"]
        result3 = myobj.frames_dict["isotopologue_props_computed"]

        self.assertListEqual(result1.loc["unknown", :].round(3).tolist(),
                             [2958.60, 3752.800, 4046.560, 4728.74])
        self.assertTrue(np.all(
            np.around(
                result2.loc["acCoA", :].to_numpy(), 6) == np.array(
                [0.320188, 0.313595, 0.235388, 0.039493])))
        self.assertAlmostEqual(result3.loc["unknown_m+0", "sampleB2"],
                               0.990623, places=6)
        self.assertAlmostEqual(result3.loc["acCoA_m+1", "sampleA2"],
                               0.312983, places=6)
        self.assertAlmostEqual(result3.loc["acCoA_m+1", "sampleA2"],
                               0.312983, places=6)
        self.assertAlmostEqual(result3.loc["acCoA_m+0", "sampleB1"],
                               0.607659, places=6)

    def test_true_key_value_available_frames(self):
        emulate_df = pd.DataFrame({"m0": [87, 85], "m1": [64, 37]})
        myobj = CompositeData("generic_xlsx")
        myobj.frames_dict = {"FractionsIsotopic": emulate_df,
                             "totalAbunds": emulate_df}
        config = {"isotopologues": None, "mean_enrichment": None,
                  "abundances": "totalAbunds",
                  "isotopologue_proportions": "FractionsIsotopic"}
        myobj.true_key_value_available_frames(config)
        result1 = myobj.available_frames
        result2 = myobj.reverse_available_frames
        self.assertListEqual(list(result1.values()), list(result2.keys()))
        self.assertListEqual(list(result2.values()), list(result1.keys()))
        self.assertListEqual(
            list(result1.values()),
            ["FractionsIsotopic", "totalAbunds"])
        self.assertListEqual(
            list(result2.values()),
            ["isotopologue_proportions", "abundances"])

    def test_stomp_fraction_values(self):

        df1 = pd.DataFrame({"X_m+0": [-0.013, 0.9], "X_m+1": [0.04, 1.025]}).T
        df2 = pd.DataFrame({"X": [0.1, 0.6], "W": [-0.43, 1.44]}).T
        myobj = CompositeData("VIBMEC_xlsx")
        myobj.frames_dict = {"IsotopicProps": {"cell": df1},
                             "FracContribs": {"cell": df2}}
        config = {
            "mean_enrichment": "FracContribs", "isotopologues": None,
            "isotopologue_proportions": "IsotopicProps", "abundances": None}
        emulate_arg = MimicArg()
        myobj.stomp_fraction_values(emulate_arg, config)
        self.assertTrue(np.all(
            myobj.frames_dict["IsotopicProps"]["cell"][0].to_numpy() == (
                np.array([0.00, 0.04])))
        )
        self.assertTrue(np.all(
            myobj.frames_dict["IsotopicProps"]["cell"][1].to_numpy() == (
                np.array([0.9, 1.0])))
        )
        self.assertTrue(np.all(
            myobj.frames_dict["FracContribs"]["cell"][0].to_numpy() == (
                np.array([0.1, 0.0])))
        )
        self.assertTrue(np.all(
            myobj.frames_dict["FracContribs"]["cell"][1].to_numpy() == (
                np.array([0.6, 1.0])))
        )

    def test_pull_internal_standard(self):
        df = pd.DataFrame({ "sample-a": [87, 64, 14],
                            "sample-b": [85, 37, 17]})
        df.index = ["X", "W_acid", "Z"]
        config = {
            "mean_enrichment": "FracContribs", "isotopologues": None,
            "isotopologue_proportions": "Prps", "abundances": "ToAbundances"}
        myobj = CompositeData("rule_tsv")
        myobj.frames_dict['ToAbundances'] = df
        emulate_arg = MimicArg()
        myobj.pull_internal_standard(config,  emulate_arg)
        result = myobj.internal_standard_abun_df
        self.assertListEqual(list(result.columns), ['sample', "W_acid"])
        self.assertTrue(np.all(np.array(result.shape) == np.array([2, 2])))
        self.assertTrue(np.all(np.array(
            result["W_acid"]) == np.array([64, 37])))
