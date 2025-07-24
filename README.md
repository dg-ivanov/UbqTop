# UbqTop 1.0

**Topologically-driven Top-down Mass Spectrometry Data Assesment for Ubiquitinated Proteins**

UbqTop is a centralized software package for automated processing of the TD-MS data for ubiquitinated proteins. This package is focused on the (1) graph-based approach for calculation of the theoretical fragmentation patterns for branched proteins
(2) fast experimantal fragment matching and (3) Bayesian-like scoring of the detected isoforms based on the experimental spectra.

UbqTop works in the vendor-independent manner and can process TD-MS data containig $a$, $b$, $c$, $x$, $y$, $z$ fragments, making it suitable
for analysis for TD-MS data using all common fragmentation techniques

Principle of operation for UbqTop is described in .... If you are using this software, please cite this paper

## Installation
UbqTop is completely written in Python 3 and requires several external libraries available in pip. Specifically, it uses pyOpenMS for calculation of the isotope patterns of fragment ions, and wxPython for GUI. Before running the scripts, all of the needed libraries need to be installed using pip.

UbqTop can be installed by cloning this repository. The whole code in splitted in two files: main.py and backend.py. backend.py contains all computational functions, while main.py contains the user interface.

_In order to start UbqTop, the user needs to load UbqTop GUI by running main.py script_

## Running the analysis

### 1. Spectra loading

UbqTop supports any TD-MS spectra presented in the x-y format. The example of spectral file can be found in the **EXAMPLES** folder. The loading dialog is tuned to load *.txt files, and supports both **comma-separated** and **tab-separated** files.


Users of Thermo instruments can export TD-MS data to clipboaed in the text format by right-click on the spectral window in QualBrowser and selecting "Copy to Clipboard (exact mass)" option. 
Users of the Bruker FT-ICR instruments can export file in the "Simple x-y text format" option in the Export Window in DataAnalysis.

**1. Load spectrum ...** button in the main window of UbqTop can be used to . 
If the spectrum is loaded correctly, it appears in the main spectrum window of UbqTop. In this case, status bar of UbqTop shows __Loading spectrum done__ message

### 2. Config file loading

In order to encode the branched proteins sequences, UbqTop uses the
