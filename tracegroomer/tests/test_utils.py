from unittest import TestCase

import pandas as pd
import numpy as np

from tracegroomer import utils


class Test(TestCase):

    def test_isotopologues_meaning_df(self):
        my_list = ["acCoA_m+0", "acCoA_m+1", "acCoA_m+2",
                   "unknown_m+0", "unknown_m+1"]
        result = utils.isotopologues_meaning_df(
            isotopologues_full_list=my_list
        )
        self.assertListEqual(list(result.columns),
                             ["metabolite", "m+x", "isotopologue_name"])
        self.assertListEqual(
            list(result["metabolite"]),
            ["acCoA", "acCoA", "acCoA", "unknown", "unknown"])
        self.assertListEqual(
            list(result["isotopologue_name"]), my_list)
        self.assertListEqual(
            list(result["m+x"]), ["m+0", "m+1", "m+2", "m+0", "m+1"])

    def test_compute_abund_from_absolute_isotopol(self):
        df = pd.DataFrame({
            "sampleA1": [5769.87, 2879.3, 1956.0, 2879.3, 79.3],
            "sampleA2": [6224.543, 3676.4, 1845.4, 3676.4, 76.4],
            "sampleB1": [7787.676, 4023.0, 1005.2, 4023.0, 23.56],
            "sampleB2": [86457.5, 4684.4, 1309.0, 4684.4, 44.34]
        })
        df.index = ["acCoA_m+0", "acCoA_m+1", "acCoA_m+2",
                    "unknown_m+0", "unknown_m+1"]
        metabolites2isotopologues_df = pd.DataFrame({
            "isotopologue_name": ["acCoA_m+0", "acCoA_m+1", "acCoA_m+2",
                                  "unknown_m+0", "unknown_m+1"],
            "metabolite": ["acCoA", "acCoA", "acCoA", "unknown", "unknown"],
            "m+x": ["m+0", "m+1", "m+2", "m+0", "m+1"],
        })
        result = utils.compute_abund_from_absolute_isotopol(
            df, metabolites2isotopologues_df)
        self.assertAlmostEqual(result.loc["acCoA", "sampleA1"],
                               10605.17, places=4)
        self.assertAlmostEqual(result.loc["acCoA", "sampleB1"],
                               12815.876, places=4)
        self.assertAlmostEqual(result.loc["unknown", "sampleA2"],
                               3752.8, places=4)
        self.assertAlmostEqual(result.loc["unknown", "sampleB2"],
                               4728.74, places=4)

    def test_compute_isotopologues_proportions_from_absolute(self):
        df = pd.DataFrame({
            "sampleA1": [5769.87, 2879.3, 1956.0, 2879.3, 79.3],
            "sampleA2": [6224.543, 3676.4, 1845.4, 3676.4, 76.4],
            "sampleB1": [7787.676, 4023.0, 1005.2, 4023.0, 23.56],
            "sampleB2": [86457.5, 4684.4, 1309.0, 4684.4, 44.34]
        })
        df.index = ["acCoA_m+0", "acCoA_m+1", "acCoA_m+2",
                    "unknown_m+0", "unknown_m+1"]
        metabolites2isotopologues_df = pd.DataFrame({
            "isotopologue_name": ["acCoA_m+0", "acCoA_m+1", "acCoA_m+2",
                                  "unknown_m+0", "unknown_m+1"],
            "metabolite": ["acCoA", "acCoA", "acCoA", "unknown", "unknown"],
            "m+x": ["m+0", "m+1", "m+2", "m+0", "m+1"],
        })
        result = utils.compute_isotopologues_proportions_from_absolute(
            df, metabolites2isotopologues_df
        )

        self.assertAlmostEqual(result.loc['acCoA_m+0', 'sampleA2'],
                               0.529913, places=6)
        self.assertAlmostEqual(result.loc['acCoA_m+1', 'sampleB1'],
                               0.313908, places=6)
        self.assertAlmostEqual(result.loc['acCoA_m+2', 'sampleB2'],
                               0.014159, places=6)
        self.assertAlmostEqual(result.loc['unknown_m+0', 'sampleA1'],
                               0.973197, places=6)
        self.assertAlmostEqual(result.loc['unknown_m+1', 'sampleB2'],
                               0.009377, places=6)

    def test_compute_MEorFC_from_isotopologues_proportions(self):
        df = pd.DataFrame({
              "sampleA1": [0.54406, 0.2715,  0.18444, 0.9732,  0.0268],
              "sampleA2": [0.52991, 0.31298, 0.1571,  0.97964, 0.02036],
              "sampleB1": [0.60766, 0.31391, 0.07843, 0.99418, 0.00582],
              "sampleB2": [0.93517, 0.05067, 0.01416, 0.99062, 0.00938]})

        df.index = ["acCoA_m+0", "acCoA_m+1", "acCoA_m+2",
                    "unknown_m+0", "unknown_m+1"]

        metabolites2isotopologues_df = pd.DataFrame({
            "isotopologue_name": ["acCoA_m+0", "acCoA_m+1", "acCoA_m+2",
                                  "unknown_m+0", "unknown_m+1"],
            "metabolite": ["acCoA", "acCoA", "acCoA", "unknown", "unknown"],
            "m+x": ["m+0", "m+1", "m+2", "m+0", "m+1"],
        })
        result = utils.compute_MEorFC_from_isotopologues_proportions(
            df, metabolites2isotopologues_df
        )
        self.assertAlmostEqual(result.loc['acCoA', 'sampleA1'],
                               0.32019, places=5)
        self.assertAlmostEqual(result.loc['acCoA', 'sampleB1'],
                               0.235385, places=5)
        self.assertAlmostEqual(result.loc['unknown', 'sampleA2'],
                               0.02036, places=5)
        self.assertAlmostEqual(result.loc['unknown', 'sampleB2'],
                               0.00938, places=5)

    def test_compute_sums_isotopol_props(self):
        df = pd.DataFrame({"metabolite": ['a', 'a', 'b', 'b'],
                           "col1": [0.3, 0.7, 0.2, 0.8],
                           "col2": [0.7, 0.3, 0.4, 0.6]})
        result = utils.compute_sums_isotopol_props(df)
        self.assertEqual(result.loc['a', 'col1'], 1)
        self.assertEqual(result.loc['a', 'col2'], 1)
        self.assertEqual(result.loc['b', 'col1'], 1)
        self.assertEqual(result.loc['b', 'col2'], 1)

    def test_impute_custom_levels_to_df(self):
        melted_df = pd.DataFrame({
            "metabolite": ["AcCoA", "AcCoA", "AcCoA", "unknown", "unknown",
                           "gly", "gly", "gly", "AcCoA", "AcCoA", "AcCoA",
                           "unknown", "unknown", "gly", "gly", "gly"],
            "isotopologue_type": [0, 1, 2, 0, 1, 0, 1, 2, 0, 1, 2,
                                  0, 1, 0, 1, 2],
            "samples": ["s1", "s1", "s1", "s1", "s1", "s1", "s1", "s1", "s2",
                        "s2", "s2", "s2", "s2", "s2", "s2", "s2"],
            "value": [0.3, 0.6, 0.1, 0.4, 0.6, 0.2, 0.5, 0.3,
                      0.25, 0.62, 0.15, 0.5, 0.5, 0.23, 0.46, 0.31]
          })
        result = utils.impute_custom_levels_to_df(melted_df)

        self.assertEqual(result['metabolite'].dtype, 'category')
        self.assertFalse(np.all(
            result['metabolite'].cat.categories.to_numpy() == np.array(
                melted_df['metabolite'].unique()))
            )
        self.assertTrue(np.all(
            result['metabolite'].cat.categories.to_numpy() == np.array(
                ['unknown', 'gly', 'AcCoA'])
            )
        )

    def test_abund_divideby_internalStandard(self):
        df = pd.DataFrame({
            "sampleA1": [10605.17, 26438.33, 2958.6],
            "sampleA2": [11746.343, 21626.02, 3752.8],
            "sampleB1": [12815.876, 11652.1, 4046.56],
            "sampleB2": [92450.98, 18658.87, 4728.74]})
        df.index = ["acCoA", "unknown",  "W_acid"]
        internal_standard_df = pd.DataFrame({
            "samples": ["sampleA1", "sampleA2", "sampleB1", "sampleB2"],
            "W_acid": [2958.6, 3752.8, 4046.56, 4728.74]})
        internal_standard_df.index = internal_standard_df["samples"]
        frames_dict = {"MyAbunds": df}
        confdict = {"abundances": "MyAbunds"}

        result = utils.abund_divideby_internalStandard(
            frames_dict, confdict, internal_standard_df,
            use_internal_standard="W_acid")
        self.assertTrue(np.all(np.around(
            result['MyAbunds'].loc["acCoA", :].to_numpy(), 6) == np.array(
                [3.584523, 3.130021, 3.167104, 19.550870])))
        self.assertTrue(np.all(np.around(
            result['MyAbunds'].loc["unknown", :].to_numpy(), 6) == np.array(
                [8.936095, 5.762636, 2.879508,  3.945844])))
        self.assertTrue(np.all(np.around(
            result['MyAbunds'].loc["W_acid", :].to_numpy(), 6) == np.array(
                [1, 1, 1, 1])))

    def test_divide_by_amount_material(self):
        df = pd.DataFrame({
            "sampleA1": [5769.87, 2879.3, 1956.0, 2879.3, 79.3],
            "sampleA2": [6224.543, 3676.4, 1845.4, 3676.4, 76.4],
            "sampleB1": [7787.676, 4023.0, 1005.2, 4023.0, 23.56],
            "sampleB2": [86457.5, 4684.4, 1309.0, 4684.4, 44.34]
        })
        df.index = ["acCoA_m+0", "acCoA_m+1", "acCoA_m+2",
                    "unknown_m+0", "unknown_m+1"]

        micrograms_weight = pd.DataFrame({
            "0": [800, 1000, 1700, 2200]})  # amount_material
        micrograms_weight.index = [
            "sampleA1", "sampleA2", "sampleB1", "sampleB2"]
        frames_dict = {"MyIsotopes": df}
        confdict = {"isotopologues": "MyIsotopes"}

        result = utils.divide_by_amount_material(
            frames_dict, confdict, material_df=micrograms_weight,
            alternative_method=True, metric="isotopologues")

        witness = (df.loc["acCoA_m+1", "sampleA2"] / micrograms_weight.loc[
            "sampleA2", "0"]) * micrograms_weight["0"].mean()

        self.assertAlmostEqual(
            result['MyIsotopes'].loc["acCoA_m+1", "sampleA2"],
            witness, places=6)

        self.assertTrue(np.all(
            np.around(
                result['MyIsotopes'].loc['acCoA_m+0', :], 6
            ) == np.array(
                [10277.580938, 8869.973775, 6527.904882, 56000.880682]))
        )
        self.assertTrue(np.all(
            np.around(
                result['MyIsotopes'].loc['acCoA_m+2', :], 6
            ) == np.array(
                [3484.125000, 2629.695000,  842.594118,  847.875000]))
        )
        self.assertTrue(np.all(
            np.around(
                result['MyIsotopes'].loc['unknown_m+1', :], 6
            ) == np.array(
                  [141.253125, 108.870000,  19.748824,  28.720227]))
        )
