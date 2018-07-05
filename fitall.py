import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from pylab import *
import sncosmo
import myfun as mf
import marshaltools as mt
import time
import datetime
from collections import OrderedDict as odict

##### OUTPUT FILE

save = 1

filename = "runs/fitall_"+datetime.datetime.now().strftime("%y-%m-%d-%H-%M-%S")

mt.filters("./")

alpha = 0.14
beta = 3.14

lcs = mt.ProgramList("Cosmology")
name = lcs.table['name']
ra = lcs.table['ra']
dec = lcs.table['dec']

nSN = len(name)

# Identify duplicates, dtheta = maximum separation (in arcsec)
dtheta = 2 
SN_dupl,theta_dupl = mf.duplicate(name,ra,dec,dtheta)

nobj = [1]*nSN
for i in range(0,len(SN_dupl)):
    nobj[np.where(name==SN_dupl[i][0])[0][0]] = 2
    nobj[np.where(name==SN_dupl[i][1])[0][0]] = -1


print("%d SNe"%nSN)
# Write header
if save == 1:
    out = open(filename+".txt","w")
    out.write("# %d Marshall SNe \n"%nSN)
    out.write("# LC fitter: salt \n")
    out.close()


###########
### FIT ###
###########

idupl = 0

for isn in range(0,nSN):

    # Are spectra available on the marshall?
    spec_flag = lcs.check_spec(name[isn]) 
    # Load light curve
    lc_dataSN = lcs.get_lightcurve(name[isn])
    # Table of data in sncosmo format (merging duplicates)
    lc_data, clas, idupl = mf.table_sncosmo_dupl(lc_dataSN,SN_dupl,nobj[isn],idupl)
    
    if len(lc_dataSN.table) > 0:

        # Removed duplicate (this will be merged below in table_sncosmo_dupl)
        if (nobj[isn] == -1):
            mf.print_and_save(name[isn],isn,nSN,lc_data,None,clas,spec_flag,save,-1,filename)

        else:
            
            # Retrieve redshift and ebv_mw
            z = lc_data.meta['z']
            ebv_mw = lc_data.meta['mwebv']
            # remove non-detection upper limits
            lc_data = lc_data[np.where(lc_data['flux']!=0)]
            # remove wrong early points in some objects
            lc_data = lc_data[np.where(lc_data['mjd']>58000)]

            # Peculiar Ia
            if (clas != "SN Ia" and clas != "SN Ia-norm"):
                mf.print_and_save(name[isn],isn,nSN,lc_data,None,clas,spec_flag,save,1,filename)

            # Missing redshift
            elif (z==None):
                mf.print_and_save(name[isn],isn,nSN,lc_data,None,clas,spec_flag,save,2,filename) 

            #Do fit
            else:

                dust = sncosmo.CCM89Dust()
                model=sncosmo.Model(source="salt2",effects=[dust],effect_names=['mw'],effect_frames=['obs']) 
                model.set(z=z,mwebv=ebv_mw,mwr_v=3.1)
               
                try:
                    fit = sncosmo.fit_lc(lc_data, model,['t0', 'x0', 'x1', 'c'],bounds={'x1':(-10,10)})
                    res = fit[0]
                    fitted_model = fit[1]
    
                    # Distance modulus + error
                    param = odict()
                    for n in res['vparam_names']:
                        param[n] = res['parameters'][res['param_names'].index(n)]
                    mB = mf._get_bandmag('bessellb', 'ab', **param)
                    mu = mB + alpha * fitted_model.parameters[3] - beta * fitted_model.parameters[4]
                    
                    grad = mf._get_bandmag_gradient('bessellb', 'ab', param, res['errors'], {})
                    J = np.array([[1,0,0,0],
                          [0,1,0,0],
                          [0,0,1,0],
                          [0,0,0,1],
                            grad])
                    cov = J.dot(res['covariance']).dot(J.T)
                    sig_mB = sqrt(cov[4,4])
                
                    c = np.delete( (np.delete(cov,(0,1),axis=0)) , (0,1), axis=1)
                    v = np.array([alpha,-beta,1.])
                    sig_mu = ( v.dot(c).dot(v.T) )**0.5
                    
                    ###### Plots and output files
                    mf.print_and_save(name[isn],isn,nSN,lc_data,fit,mu,sig_mu,save,0,filename)
                    
                except:
                   mf.print_and_save(name[isn],isn,nSN,lc_data,None,clas,spec_flag,save,3,filename)

    else:
        mf.print_and_save(name[isn],isn,nSN,lc_data,None,clas,spec_flag,save,4,filename)