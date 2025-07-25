# UbqTop 1.0

**Topologically-driven Top-down Mass Spectrometry Data Assesment for Ubiquitinated Proteins**

UbqTop is a centralized software package for automated processing of the TD-MS data for ubiquitinated proteins. This package is focused on the (1) graph-based approach for calculation of the theoretical fragmentation patterns for branched proteins
(2) fast experimantal fragment matching and (3) Bayesian-like scoring of the detected isoforms based on the experimental spectra.

UbqTop works in the vendor-independent manner and can process TD-MS data containig $a$, $b$, $c$, $x$, $y$, $z$ fragments, making it suitable
for analysis for TD-MS data using all common fragmentation techniques

Principle of operation for UbqTop is described in the associated bioRxiv preprint: 

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

UbqTop uses internal FASTA-like format for encoding sequnces and topologies of branched proteins. The details on the format and principle of the encoding can be found the UbqTop paper. The examples of the config files are presented in **EXAMPLES** folder

The config file can be loaded into the UbqTop by clicking the **2. Load config from file ...** button or entered manually into the text window behind the button. The loaded text file appears in the same text window, and can be further modified. 
This content in the window is used for the calculation of fragment series and generation of the fragment series intervals

### 2. Isoform table generation

After the config file has been loaded the isoform table needs to be generated. In order to do this, the user needs to select fragment types that will be used in the analysis by multiple choice buttons, and then click **3. Generate Isoform Table** button. After that, the table with fragment series and fragment series intervals will be generated

### 3. Annotation of the MS2 spectrum

In order to initiate the annotation of the fragments in the MS2 spectrum according to the isoform table, click **3. Annotate spectrum**. The progress can be checked in the status bar of the software. Once the spectrum has been processed (full processing usually requires 1-10 s), the status bar shows message "Spectrum annotation processing done", and the series from the isoform table appears in the **Show fragment series** list. In order to manually examine the fragment coverage 
of particular series, it is possible to select the series of interest from the list and click **Show** button.

Once the specific series is selected, two diagrams appears behind the MS spectrum window: CP values diagram (on the left side) and series coverage diagram (on the right side). The details on the CP value meaning can be found in the associated paper. CP values diagram is a raster diagram where 
x-coordinate of the dot refers to the fragment number in series, and y-coordinate - to the charge state. The color of the raster dot is associated with the CP value of the particular charge state of the corresponding fragment. The dots in the CP values diagram are interactive: double-click on the dot shows isotopes of the particular fragment in the corresponding charge state in the MS2 spectrum window. The isotopes that were matched in the spectrum shows as a green grid lines, while the missing ones - as red grid lines. 

Parameters of the fragment matching can be adjusted in the parameters section. The __Intensity threshold__ is the minimum relative intensity that is used for peak picking. The __m/z accuracy__ is a maximum deviation (in Da) between experimental and theoretical isotope m/z values that can be used for matching fragments. The default values are adjusted for the Orbitrap 240k EThcD MS2 spectra.

Annotation of the sequence coverage is based on the CP1/CP2 criteria described in the paper. The CP1 and CP2 values can be adjusted in the parameters section within software. After changing the values, the annotation of spectra needs to be updated by clicking the **3. Annotate spectrum** button.

### 4. Calculation of scores

To open the bayesian scoring window the button **4. Calculate scores ...** needs to be clicked. Activation of this function will cause opening the second window that can be used for plotting of scores graphs.
In the property grid in the software it is possible to adjust the $\alpha$ and $\beta$ that are responsible for influece of the detected ($\alpha$) and missing ($\beta$) fragments. The pre-filled parameters were found to be the optimal for the characterization of the broad types of chains.

In order to generate scores, the **Calculate Scores** button needs to be clicked. The scores plots for each isoform will appear in the graph window, and the values of the scores after considering each fragments appears in the text window behind the **Calculate Scores** button. The captions on the scores plot can be edited by changing the parameters in the __Plot Parameters__ section of the parameters section in scores window. Scoring plots can be exported by clicking the **Save** icon on top of the save scores window.
