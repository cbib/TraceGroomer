# TraceGroomer

[![PyPI - Python Version](https://img.shields.io/pypi/v/tracegroomer)](https://pypi.org/project/tracegroomer/)
[![bioconda package](https://img.shields.io/conda/v/bioconda/tracegroomer)](https://anaconda.org/bioconda/tracegroomer)

TraceGroomer is a solution for formatting and normalising **Trace**r metabolomics given file(s), 
to produce the .csv files which are ready for [DIMet](https://github.com/cbib/DIMet) tool.

Currently, three **styles of format** of Tracer (or Isotope-labeled) metabolomics measurements files are accepted:

1. IsoCor results (.tsv measurments file).
2. Results provided by the VIB Metabolomics Expertise Center (El-Maven results are shaped by VIB MEC team into a multi-sheet .xlsx file).  
3. A 'generic' .xlsx measurements file, manually set by the user.

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

Install it via `pip`: `pip install tracegroomer`

Tracegroomer is also available as a [conda package](https://bioconda.github.io/conda-package_index.html)

Alternatively, if you are a developer, you can do a local install:
<details>
<summary>
Local install of TraceGroomer <sup><sub>(click to show/hide)</sub></sup>
</summary>
For a local install, clone this repository, make sure you have activated 
your virtual environment with Python 3.10+
(<code>source MY_VIRTUAL_ENV/bin/activate</code>), with <code>poetry</code> installed.

Then install dependencies: locate yourself in `TraceGroomer` and run
```
poetry install
```

After this, the tool is ready to use:
```
python -m tracegroomer --help
```
</details>
  

## Documentation

All the details about how to use TraceGroomer can be found on the dedicated [Wiki](https://github.com/cbib/TraceGroomer/wiki) page. 
This is where you will find the information of the supported formats, examples, and how to run TraceGroomer.


  
# Getting help

For any information or help running TraceGroomer, you can get in touch with: 

* [Johanna Galvis](mailto:deisy-johanna.galvis-rodriguez[AT]u-bordeaux.fr)
* [Macha Nikolski](mailto:macha.nikolski[AT]u-bordeaux.fr)
* [Benjamin Dartigues](mailto:benjamin.dartigues[AT]u-bordeaux.fr)

---

# LICENSE MIT

Copyright (c) 2024

    Johanna Galvis (1,2)    deisy-johanna.galvis-rodriguez@u-bordeaux.fr
    Benjamin Dartigues (2)	benjamin.dartigues@u-bordeaux.fr
    Slim Karkar (1,2)       slim.karkar@u-bordeaux.fr
    Helge Hecht (3,5)       helge.hecht@recetox.muni.cz
    Bjorn Gruening (4,5)    bjoern.gruening@gmail.com
    Macha Nikolski (1,2)    macha.nikolski@u-bordeaux.fr

    (1) CNRS, IBGC - University of Bordeaux,
    1, rue Camille Saint-Saens, Bordeaux, France

    (2) CBiB - University of Bordeaux,
    146, rue Leo Saignat, Bordeaux, France

    (3) RECETOX
    Faculty of Science, Masaryk University, Kotlářksá 2, 611 37 Brno, Czech Republic

    (4) University of Freiburg,
    Freiburg, Germany

    (5) Galaxy Europe

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
