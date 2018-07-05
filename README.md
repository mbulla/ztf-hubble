# ztf-hubble

Tools for fitting ZTF Ia's and creating an Hubble diagram.

# Usage

All relevant functions are included in marshaltools.py (from Ulrich Feindt: https://github.com/ufeindt/marshaltools) and myfun.py

Salt2 fitting is performed in fitall.py. Fitting parameters are saved to the directory "runs" and the fits of all the individual supernovae are saved to the directory "fits".

Type Ia supernovae that are either (i) peculiar, (ii) with no spectroscopic redshift from the Marshall or (iii) with no photometry are not considered. Objects from the Marshall that are found to be duplicates (see duplicate function in myfun.py) are merged. 

Histograms of relevant parameters and Hubble diagrams are plotted in hubble.py. The user is asked to specify if cuts on (i) x1, (ii) earliest g/r data, (iii) salt2 reduced chi2 and (iv) number of data points are needed.
