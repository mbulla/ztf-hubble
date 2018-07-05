# ztf-hubble

Tools for fitting ZTF Ia's and creating an Hubble diagram.

# Usage

All relevant functions are included in marshaltools.py (from Ulrich Feindt: https://github.com/ufeindt/marshaltools) and myfun.py

Salt2 fitting is performed in fitall.py. Fitting parameters are saved to the directory "runs" and the fits of all the individual supernovae are saved to the directory "fits".

Type Ia supernovae in the Cosmology project that are either (i) peculiar, (ii) with no spectroscopic redshift or (iii) with no photometry are not considered in the fitting. Objects from the Marshall that are found to be duplicates (see duplicate function in myfun.py) are merged. 

Histograms of relevant parameters and Hubble diagrams are plotted in hdiag.py. The user is asked to specify if cuts on (i) x1, (ii) salt2 reduced chi2, (iii) maximum phase for g/r data and (iv) number of g/r data points are needed.
