import numpy as np
import sncosmo

#########################################
############ PEAK MAG FROM ULI ##########
#########################################

from collections import OrderedDict as odict

def _get_bandmag(band, magsys, t=0, rest_frame=True, **kwargs):
    """
    Returns mag at max for the model, band and magsys
    Arguments:
    model  -- sncosmo model, e.g. SALT2
    band   -- sncosmo band object or string, e.g. 'bessellb'
    magsys -- magnitude system, e.g. 'ab'
    Keyword arguments:
    t -- time relative to t0 (observer-frame), at which to evaluate
    rest_frame -- default: True, overrides the redshifts
    """
    model = sncosmo.Model(source='salt2')
    if rest_frame:
        kwargs['z'] = 0

    model.set(**kwargs)
    return model.bandmag(band,magsys,kwargs['t0'] + t)

def _get_bandmag_gradient(band, magsys, param, sig, fixed_param, 
                          t=0, rest_frame=True):
    """
    Return gradient of _get_bandmag as function of param
    param, sig must be dictionaries of means and uncertainties
    Best use odicts to make sure that the order of the components is correct
    """
    model = sncosmo.Model(source='salt2')
    out = []
    
    if rest_frame:
        if 'z' in param.keys():
            param['z'] = 0
        if 'z' in fixed_param.keys():
            fixed_param['z'] = 0

    model.set(**fixed_param)
    for key,val in param.items():
        model.set(**param)
        h = sig[key] / 100.
        
        model.set(**{key: val - h})
        m0 = model.bandmag(band, magsys, param['t0'] + t)

        model.set(**{key: val + h})
        m1 = model.bandmag(band, magsys, param['t0'] + t)
        
        out.append((m1 - m0) / (2. * h))

    return np.array(out)


################################################################
####### IDENTIFY SNe ON THE MARSHALL THAT ARE SAME SN ##########
################################################################

def duplicate(SN,ra,dec,dtheta):

    from itertools import combinations
    
    SN_comb = []
    ra_comb = []
    dec_comb = []

    for i in list(combinations(SN,2)):
        SN_comb.append(i)
    for i in list(combinations(ra,2)):
        ra_comb.append(i)
    for i in list(combinations(dec,2)):
        dec_comb.append(i)

    t = np.pi/180.
    
    theta = []

    for i in range(0,len(SN_comb)):
        costheta = np.sin(dec_comb[i][0]*t)*np.sin(dec_comb[i][1]*t)+np.cos(dec_comb[i][0]*t)*np.cos(dec_comb[i][1]*t)*np.cos(ra_comb[i][0]*t-ra_comb[i][1]*t)

        # Theta in arcsec
        theta.append( np.arccos(costheta) / t * 3600 )
    
    SN_dupl = np.array(SN_comb)[np.where(np.array(theta)<dtheta)]
    theta_dupl = np.array(theta)[np.where(np.array(theta)<dtheta)]

    return SN_dupl, theta_dupl


################################################################
##### Table of data in sncosmo format (merging duplicates) #####
################################################################

def table_sncosmo_dupl(data,SN_dupl,nobj,idupl):

    import marshaltools as mt
    from astropy.table import vstack

    ### No duplicates
    if (nobj==1):  

        # Light curves
        classification = data.classification
        z = data.redshift
        ebv_mw = data.mwebv

        lc_data = data.table_sncosmo
       
    ### Duplicates
    elif nobj==2:  

        lcs = mt.ProgramList("Cosmology")
        
        lc1 = lcs.get_lightcurve(SN_dupl[idupl][0])
        lc1_data = lc1.table_sncosmo

        lc2 = lcs.get_lightcurve(SN_dupl[idupl][1])
        lc2_data = lc2.table_sncosmo

        classification = lc1.classification
        z1 = lc1.redshift
        z2 = lc2.redshift
        if z1 == None:
            z = z2
        elif z2 == None:
            z = z1
        else:
            z = 0.5 * (float(z1) + float(z2))

        lc_data = vstack([lc1_data,lc2_data])
        lc_data.meta['z'] = z
        lc_data.meta['mwebv'] = lc1.mwebv   
        
        lc_data.sort('mjd')

        idupl+=1


    return lc_data, classification, idupl


################################################################
############### Print messages and saved parameters ############
################################################################

def print_and_save(sn,isn,nSN,data,fit,par1,par2,save,message,filename):

    if save == 1:
        out = open(filename+".txt","a+")

    # SN included in the analysis, good fit
    if message == 0:

        print "%s (%d out of %d)"%(sn,isn+1,nSN)

        if save == 1:

            out.write("%s "%sn)
            
            for i in range(0,len(fit[1].parameters)):
                out.write("%f "%fit[1].parameters[i])
            
            out.write("%f %f "%(par1,par2))
            
            lc_datag = data[np.where(data['band']=="p48g")]
            lc_datar = data[np.where(data['band']=="p48r")]
            
            if len(lc_datag)!=0 and len(lc_datar)!=0:
                out.write("%f "%(lc_datag['mjd'][0] - fit[1].parameters[1]))
                out.write("%f "%(lc_datar['mjd'][0] - fit[1].parameters[1]))
            elif len(lc_datag)==0:
                out.write("1e6 ")
                out.write("%f "%(lc_datar['mjd'][0] - fit[1].parameters[1]))
            elif len(lc_datar)==0:
                out.write("%f "%(lc_datag['mjd'][0] - fit[1].parameters[1]))
                out.write("1e6 ")
            else:
                out.write("1e6 ")
                out.write("1e6 ")  

            out.write("%d "%len(lc_datag) )
            out.write("%d "%len(lc_datar) )

            out.write("%f %d \n"%(fit[0].chisq,fit[0].ndof))

            fig = sncosmo.plot_lc(data, model=fit[1],errors=fit[0].errors,xfigsize=15,tighten_ylim=True)
            fig.savefig("fits/"+sn+".pdf",bbox_inches='tight')
                              
            
    elif message == -1:

        print "%s not considered, duplicate "%sn
        
        if save == 1:
            out.write("#%s not considered, duplicate \n"%sn) 


    # Not a normal Ia
    elif message == 1:

        print "%s not considered, not a SN Ia norm (%s)"%(sn,par1)
        
        if save == 1:
            out.write("#%s not considered, not a SN Ia norm (%s) \n"%(sn,par1))     

    # Missing redshift
    elif message == 2:

        if par2 == 1:
            print "%s not considered, missing redshift (spectra available) "%sn
            if save == 1:
                out.write("#%s not considered, missing redshift (spectra available) \n"%sn) 
        else:
            print "%s not considered, missing redshift "%sn
            if save == 1:
                out.write("#%s not considered, missing redshift \n"%sn)     

    # Bad fit
    elif message == 3:

        print "%s not considered, bad fit"%sn

        if save == 1:
            out.write("#%s not considered, bad fit \n"%sn)

    # No data points
    elif message == 4:

        print "%s not considered, no data points \n"%sn

        if save == 1:
            out.write("#%s not considered, no data points \n"%sn)


    if save == 1:
        out.close()
