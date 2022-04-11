# Signal Processing and TA analysis

The content of this repository is the result of a project from the CAS in Applied Data Sceince which I completed at the Univeristy of Bern in 2019-2020. _This work stopped being maintained after the end of the project and so is out of date with the current TAS analysis_

The repository consists of a [conceptual design report](https://github.com/GarethJMoore/CAS_M1_GMoore/blob/master/CDR_Gmoore.pdf) and a [jupyter notebook](https://github.com/GarethJMoore/CAS_M1_GMoore/blob/master/CAS_M1_GMoore.ipynb), along with a data file.

### The conceptual design report (_'DCR_Gmoore.pdf'_) 
Describes the data flow for Transient Absorption Spectroscopic (TAS) data. This includes the processing, cleaning, presenting and analysis of the data in the form of a report.  

### The jupyter notebook notebook (_'CAS_M1_GMoore.ipynb'_) 
Presents and describes the code used to first process the signal, then produce the TA spectra. 

1. The singal processing esentially applies the work of [Anderson _et. al_](https://aip.scitation.org/doi/full/10.1063/1.2755391) to filter out noise from the measured TAS data.
2. Production of Spectra:
- Cleans the TAS data by averaging and background subtraction.
- Shows a the process of removing the 'Chirp' from the spectroscopic data.
- Presents the spectra in different forms.
3. Spectral analysis where spectral decomposition is discussed

