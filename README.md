# TraceGroomer

[![PyPI - Python Version](https://img.shields.io/pypi/v/TraceGroomer)](https://pypi.org/project/TraceGroomer/)
[![bioconda package](https://img.shields.io/conda/v/bioconda/TraceGroomer)](https://anaconda.org/bioconda/TraceGroomer)

TraceGroomer is a solution for formatting and normalising **Trace**r metabolomics given file(s), 
to produce the .csv files which are ready for [DIMet](https://github.com/cbib/DIMet) tool.

Currently, three **styles of format** of Tracer (or Isotope-labeled) metabolomics measurements files are accepted:

1. IsoCor results (.tsv measurments file).
2. Results provided by the VIB Metabolomics Expertise Center (El-Maven results are shaped by VIB MEC team into a multi-sheet .xlsx file).  
3. A 'generic' .xlsx measurements file.

For any type of these supported inputs, TraceGroomer generates an independent file for:
i) total metabolite abundances ii) Isotopologues iii) Isotopologues' proportions and iv) mean enrichment (a.k.a fractional contributions).

Automatic formatting is performed, as well as the normalization chosen by the user:
whether by the amount of material and/or by an internal standard.
Useful advanced options are offered (e.g. if the user has only Isotopologues' absolute values, TraceGroomer can generate all the other 
measurements files automatically).


_Note_ : this script does not correct for naturally occurring isotopologues. 
Your data must be already processed by another software that performs such correction.

--------------------------

## Requirements

TraceGroomer requires Python 3.10+.  Running in a virtual environment is highly recommended.

For installing locally, clone this repository, make sure you have activated 
your virtual environment with Python 3.10+
(`source MY_VIRTUAL_ENV/bin/activate`), with `poetry` installed.

Then install dependencies: locate yourself in `TraceGroomer` and run
```
poetry install
```

After this the tool is ready to use:
```
python -m tracegroomer --help
```


## How to use TraceGroomer

Its execution takes only few seconds. Below we explain how to proceed.


### Input files

**Compulsory files:**

- the measurements file (tsv, csv -tab delimited-, or xlsx)
  <details>
  <summary>
  Description
  </summary>
  
  The measurements file is given by a Metabolomics facility. 
  It is the result of the correction by software such as IsoCor, El-Maven, etc. 
  Some times the file can be further formatted in the Metabolomics facility, before
  being delivered to the end user, which is the case of VIB MEC delivered files.
     
  TraceGroomer also accepts a "generic" format, see more details in ['generic' type of data](#users-having-generic-data).
  <p>
	  
  </p>
  
   The user will find [here](#running-a-test-with-the-provided-examples) how get examples of IsoCor direct output, VIB MEC file, and generic file.
  We provide also more details in the sections [users having VIB results as input](#users-having-vib-results), [users having IsoCor results](#users-having-isocor-results) and  ['generic' type of data](#users-having-generic-data).


  </details>

- the **metadata** file, which describes the experimental setup. 

  <details>
  <summary>
  Description and example
  </summary>
   
   The metadata is a tab delimited .csv file provided by the user,
   which has to contain 6 columns named 
<code>name_to_plot</code>, <code>timepoint</code>, 
<code>timenum</code>, <code>condition</code>, 
<code>compartment</code>, <code>original_name</code>. 

   Here is the semantics of the columns:
   
   - <code>name_to_plot</code> is the string that will appear on the figures produced by DIMet
   - <code>condition</code> is the experimental condition
   - <code>timepoint</code> is the sampling time as it is defined in your experimental setup
     (it is an arbitary string that can contain non numerical characters)
   - <code>timenum</code> is the numerical encoding of the <code>timepoint</code>
   - <code>compartment</code> is the name of the cellular compartment for which the measuring
     has been done (e.g. "endo", "endocellular", "cyto", etc)
   - <code>original_name</code> contains the column names that are provided in the quantification files
   
   _Example_:
   
   | name_to_plot | condition | timepoint | timenum | compartment | original_name |
   |--------------|-----------|-----------|---------|-------------|---------------|
   | Cond1 T0     | cond1     | T0        | 0       | comp_name   | T0_cond_1     |
   | Cond1 T24    | cond1     | T24       | 24      | comp_name   | T24_cond_1    |
   | Cond2 T0     | cond2     | T0        | 0       | comp_name   | T0_cond_2     |
   | Cond3 T24    | cond2     | T24       | 24      | comp_name   | T24_cond_2    |

   The column `name_to_plot` is not used by `tracegroomer` 
      but it will be used by DIMet, so it is practical to set it from the start.

     _Note_: You can create this file with any spreadsheet program such as Excel or Google Sheets or LibreOfice Calc. At the moment of saving your file you specify that the delimiter must be a `tab` ("Tab delimiter" or similar option depending of your context), see https://support.microsoft.com/en-us/office/save-a-workbook-to-text-format-txt-or-csv-3e9a9d6c-70da-4255-aa28-fcacf1f081e6. 

  </details>


- the configuration file

  <details>
  <summary>
  Description and example
  </summary>
  This file contains basic needed information: the name of the metadata file,
  the names (but not the paths) of the output files, and the absolute path to the output folder.
  
  The comments (`#`) serve as guide. The user must fill after the colon of each field: 
  
   ```
   # coments start with #
   # -----------------------------
   
  # absolute path to output DIRECTORY : 
   groom_out_path :  ~/examples_TraceGrommer/data/example-isocor_data 
  
   metadata: metadata1   # file name, no extension. Must be in the output DIR 
   # names of the sheets in xlsx file that exist (null otherwise) 
   abundances : null  # total abundance
   mean_enrichment : null  # mean enrichment
   isotopologue_proportions : null  # isotopologue proportions
   isotopologues : isotopologuesCorrValues  # isotopologue absolute values  
   
   ```  
   
   When the fields `abundances`, `mean_enrichment`, and/or `isotopologue_proportions` are set `null`
  (as shown in the example above), 
  the respective output file will be automatically generated (when possible from existing quantifications).

   The user will find [here](#running-a-test-with-the-provided-examples) how to get examples.


   _Note_: There exist online editors for .yml files, such as https://yamlchecker.com/, just copy-paste and edit!

  
  </details>


**Facultative files**

   * the amount of material by sample (tab delimited csv file)
   * a file with metabolites to exclude (tab delimited csv file)

The facultative files are used through the command line, which is explained
in [Advanced options](#advanced-options)

You must organize your files as follows:
```
MyProject
├── data
│   ├── dataset1_data
│   │   ├── metadata_1.csv
│   │   └── TRACER_IsoCor_out_example.tsv
└── groom_files
    └── dataset1
              ├── amount_material_weightorcells.csv
              └── config-1-groom.yml

```
This structure is recommended to easily re-use the `data` folder for DIMet.

The generic command line is:

```
python3 -m tracegroomer --targetedMetabo_path $MEASUREMENTS \
    --type_of_file $MY_TYPE_OF_INPUT \
    $MY_BASIC_CONFIG
```
Where :
- `MEASUREMENTS` is the file that contains the measurements, in absolute path.
- `MY_TYPE_OF_INPUT` corresponds to one of: `IsoCor_out_tsv`, `VIBMEC_xlsx`, `generic_xlsx`

We recommend to run a test with the provided examples if this is the
first time you use TraceGroomer. Then re-use  the 
organization and the configurations, and modify the command line to be suitable to your data.

-------------------------------------

### Running a test with the provided examples

To perform a test using the examples we provide, please download
and uncompress our example from [Zenodo](https://zenodo.org/records/10660026). The structure of the folder is:

```
examples_TraceGroomer
├── data
│   ├── example-isocor_data
│   │   ├── metadata_1.csv
│   │   └── TRACER_IsoCor_out_example.tsv
│   ├── example-sheet_data
│   │   ├── metadata_3.csv
│   │   └── TRACER_generic_custom.xlsx
│   └── example-vib_data
│       ├── metadata_2.csv
│       └── TRACER_metabo_2.xlsx
└── groom_files
    ├── example-isocor
    │   ├── amount_material_weightorcells.csv
    │   └── config-1-groom.yml
    ├── example-sheet
    │   └── config-3-groom.yml
    └── example-vib
        ├── config-2-groom.yml
        ├── nbcells-or-amountOfMaterial.csv
        └── reject_list.csv

```

Pick the example most suited to your data:

1. IsoCor output (tsv file)
2. VIB MEC xlsx file
3. or a generic type of xlsx file


### Run the script

_Note_ : if the working folder is not the 'home' directory, modify accordingly the absolute paths in the .yml files and in the bash commands.
locate yourself in ``, then run:


**For IsoCor case**:

```
python3 -m tracegroomer \
   --targetedMetabo_path ~/examples_TraceGroomer/data/example-isocor_data/TRACER_IsoCor_out_example.tsv \
   --type_of_file IsoCor_out_tsv \
   ~/examples_TraceGroomer/groom_files/example-isocor/config-1-groom.yml
```

or, for VIB MEC case:

```
python3 -m tracegroomer \
  --targetedMetabo_path ~/examples_TraceGroomer/data/example-vib_data/TRACER_metabo_2.xlsx \
  --type_of_file VIBMEC_xlsx 
  ~/examples_TraceGroomer/groom_files/example-vib/config-2-groom.yml

```

or, for generic case:

```
python3 -m tracegroomer --targetedMetabo_path ~/examples_TraceGroomer/data/example-sheet_data/TRACER_generic_sheet.xlsx \
   --type_of_file generic_xlsx 
   ~/examples_TraceGroomer/groom_files/example-sheet/config-3-groom.yml
```


## The output

The output files are saved in the  folder that  
you specified in the config `.yml` file (`groom_out_path` field). The `data/[my_dataset]` location is recommended. 
A total of 4 output files are generated if the absolute isotopologues are provided, otherwise 3 files are generated.


In this way you simply copy the entire `data/` content to the folder structure that we want to run with  [DIMet](https://github.com/cbib/DIMet) !

The format of these output files is tab-delimited .csv.


--------------------

## Advanced options

We provide advanced options for this script, check the help:
```
python -m tracegroomer --help
```
they appear as 'optional arguments' in the help menu.


You can:

- normalize by the amount of material (number of cells, tissue weight), and/or
- normalize by an internal standard (present in your data) at choice

We have some indications that can slightly differ for [users having VIB results as input](#users-having-vib-results), [users having IsoCor results](#users-having-isocor-results) or users having ['generic' type of data](#users-having-generic-data). After consulting the one of your case, please visit [Advanced options](#advanced-options) section for the offered normalisations.


 
### Users having IsoCor results

Before explaining the advanced options for this kind of data, a short explanation about what TraceGroomer performs automatically as basic formatting:

A typical IsoCor results table is described in: https://isocor.readthedocs.io/en/latest/tutorials.html#output-files
 It consists of a .tsv file which has in columns the sample, metabolite, isotopologue and all quantifications, and the rows are in piled version (the samples are repeated vertically).
 
 Our script transforms specific columns of that file into tables. As the total metabolite abundance column is not present in the input data, the total abundance per metabolite is the automatic result of the sum, per metabolite, of Isotopologues' absolute values (see `AbundanceCorrected` below). So here the correspondances:

|column in the IsoCor file  | TraceGroomer output filename |
|--------|-------|
|corrected_area |IsotopologuesAbsolute |
|isotopologue_fraction | IsotopologuesProportions|
| mean_enrichment| MeanEnrichment|
|- | AbundanceCorrected |
 
We provide the example downloadable from Zenodo (see [here](#running-a-test-with-the-provided-examples))
        
Advanced options regarding to detection limit (LOD) and blanks will not have any effect on the IsoCor type of data: LOD is not provided in the data, and the same is true for blanks. 

All the other advanced options do have effect: those related to internal standard, amount of material, and isotopologues.
 
 
 
### Users having VIB results

As shown in the example, give the names of the sheets that are present in your excel file coherently in the .yml file. 
 
Our script performs, by default:
- the subtraction of the means of the blanks across all metabolites' abundance for each sample.
- seting to NaN the values of abundance that are under the limit of detection (LOD).
- excluding metabolites whose abundance values across all samples are under LOD (excluded then from all tables by default).
- stomping fractions values to be comprised between 0 and 1 (some negative and some superior to 1 values can occur after correction of naturally occurring isotopologues by certain software dedicated to such corrections)

You can modify all those options depending on your needs, they appear as 'optional arguments' in the help menu. 


### Users having generic data

We have created this option for those formats that are not the other two scenarios, so your data is expectd to be in the form of a .xlsx file with sheets similar as in the provided in `example-sheet`:

- this .xlsx file must **NOT** contain: formulas, symbols accompanying the numeric values, nor special characters.
- each sheet must correspond to one type of quantification, see _Notes_ below.
- the header (first row, with the molecules IDs), and the first column (with the samples) can contain non numeric values.
- the isotopologues names, in the header, must follow the convention `metaboliteID_labelX`: the substring `_label` is compulsory and is located between the metabolite name (or identifier) `metaboliteID` and the number of marked carbon atoms `X` (e.g. `glucose6phosphate_label0`).

**Important** : If you only have Isotopologue Absolute values, but not the other tables: put them as a single named sheet in your .xlsx file, and TraceGroomer automatically generate all the other types of tables for you ! 

_Notes_: 
- sheets corresponding to isotopologue Proportions (when available) and isotopologue Absolute values (compulsory if the proportions not available) must have isotopologues as columns and samples as rows.
- sheets corresponding to abundance and mean enrichment  (when available) must have metabolites as columns and samples as rows.
- the sheets corresponding to isotopologues measurements must be named with a name containing the string "isotopol". The names of the sheets must be unambiguous.


## Advanced options:

Available for all the types of supported input files:

- normalize by the amount of material (number of cells, tissue weight): setting the path to the file in `--amountMaterial_path` option. The file must be like [this csv file](groomexamples/toyp2/nbcells-or-amountOfMaterial.csv), and the first column must contain the same names as in metadata 'original\_name'.
- normalize by an internal standard (present in your data) at choice: using the advanced option `--use_internal_standard`.
- remove metabolites 
- print a preview of isotopologues values
  
