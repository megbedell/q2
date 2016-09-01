import numpy as np
import q2
import pdb
import matplotlib.pyplot as plt

b_map_linear = {
# slope correction factor and error [dex/Gyr]
          'CI'  :   (0.030,0.003),
          'CH'  :   (0.030,0.003),
          'OI'  :   (0.016,0.002),
          'NaI' :   (0.015,0.003),
          'MgI' :   (0.010,0.0009),
          'MgII':   (0.010,0.0009),
          'AlI' :   (0.0147,0.0015),
          'SiI' :   (0.0051,0.0008),
          'SI'  :   (0.0046,0.0017),
          'CaI' :   (-0.0014,0.0008),
          'ScI' :   (0.0071,0.0018),
          'ScII':   (0.0071,0.0018),
          'TiI' :   (0.0061,0.0011),
          'TiII':   (0.0061,0.0011),
          'VI'  :   (-0.0034,0.0012),
          'CrI' :   (-0.0024,0.0005),
          'CrII':   (-0.0024,0.0005),
          'MnI' :   (0.0075,0.0019),
          'CoI' :   (0.012,0.003),
          'NiI' :   (0.0076,0.0017),
          'CuI' :   (0.020,0.003),
          'ZnI' :   (0.0122,0.0014),
          }

def getb_linear(species_id):
    try:
        b = b_map_linear[species_id]
    except:
        return (0.0, 0.0)
    return b


k_map_hyperbolic = {
          'CI'  :   16.3,
          'CH'  :   16.3,
          'OI'  :   48.0,
          'NaI' :    5.0,
          'MgI' :    8.0,
          'MgII':    8.0,
          'AlI' :    8.2,
          'SiI' :    5.9,
          'SI'  :    5.8,
          'CaI' :    5.0,
          'ScI' :    5.5,
          'ScII':    5.5,
          'TiI' :   15.0,
          'TiII':   15.0,
          'VI'  :    4.7,
          'CrI' :    6.2,
          'CrII':    6.2,
          'MnI' :    4.8,
          'CoI' :    5.2,
          'NiI' :    5.0,
          'CuI' :    7.1,
          'ZnI' :    6.3
          }

def getk_hyperbolic(species_id):
    try:
        k = k_map_hyperbolic[species_id]
    except:
        return 0.0
    return k


b_map_hyperbolic = {
          'CI'  :   22.8,
          'CH'  :   22.8,
          'OI'  :   45.8,
          'NaI' :   -6.3,
          'MgI' :   18.2,
          'MgII':   18.2,
          'AlI' :   15.5,
          'SiI' :   15.2,
          'SI'  :   15.1,
          'CaI' :   14.2,
          'ScI' :   10.0,
          'ScII':   10.0,
          'TiI' :   41.2,
          'TiII':   41.2,
          'VI'  :   12.6,
          'CrI' :   25.0,
          'CrII':   25.0,
          'MnI' :    9.8,
          'CoI' :    8.8,
          'NiI' :    8.6,
          'CuI' :   11.0,
          'ZnI' :   11.4
          }

def getb_hyperbolic(species_id):
    try:
        b = b_map_hyperbolic[species_id]
    except:
        return np.inf
    return b



def correct(Star, age, species_ids=None, method='linear', Ref=None, Ref_age=0.0, silent=True, errors=False):
    '''
    Apply GCE corrections and return adjusted abundances.
    method: 'linear' or 'hyperbolic', describes the GCE fit method used.
       linear: from Spina et al. 2016, arXiv:1606.04842
       hyperbolic: from Spina et al. 2016, arXiv:1606.04842
    '''
    if (Ref==None):
    	print "Ref not set; assuming it is the Sun (age = 4.6 Gyr)"
        Ref_age = 4.6
    if (Ref_age==0.0):
        # if Ref was set but no age given
        "Must set Ref_age keyword!"
        return None
    if species_ids == None:
        species_codes = sorted(set(Star.linelist['species']))
        species_ids = q2.abundances.getsp_ids(species_codes)
        if not silent:
            print '"species_ids" not provided'
            print 'Lines found for the following species: '+\
                  ','.join(species_ids)
            print ''
    Tc = []
    abund = []
    err = []
    for species_id in species_ids:
        # get abundances:
        Tc = np.append(Tc, q2.abundances.gettc(species_id))
        species_difab = np.array(getattr(Star, species_id)['difab'])
        species_difab = species_difab[species_difab != np.array(None)]  #remove Nones from where ref star was unavailable
        abund = np.append(abund, np.mean(species_difab))
        try:
            err = np.append(err, np.std(species_difab) / np.sqrt(len(species_difab) - 1) )
        except:
            err = np.append(err, 0.0)
        
    for t in set(Tc):
        # eliminate duplicate measurements of the same element
        ind = np.where(Tc == t)[0]
        if len(ind) == 2:  # won't take care of 3+ states of the same thing
            (abund[ind[0]], err[ind[0]]) = np.average(abund[ind], weights=err[ind], returned=True)
            abund = np.delete(abund, ind[1])
            err = np.delete(err, ind[1])
            Tc = np.delete(Tc, ind[1])
            species_ids = np.delete(species_ids, ind[1])
    
    for i,species_id in enumerate(species_ids):
        if not silent:
            print "*** Begin "+species_id+":"
        # fetch correction factors:
        if method is 'linear':
            (b, err_b) = getb_linear(species_id)
            if (b == 0.0):
                if not silent:
		            print "No GCE correction available."
                corr_factor = 0.0
                err_factor = 0.0
                continue
            corr_factor = (age - Ref_age)*b
            err_factor = (age - Ref_age)*err_b
        elif method is 'hyperbolic':
            k = getk_hyperbolic(species_id)
            b = getb_hyperbolic(species_id)
            if not np.isfinite(c):
    	        if not silent:
    		        print "No GCE correction available."
                corr_factor = 0.0
                err_factor = 0.0
    	        continue
            corr_factor = np.sqrt((age - k)**2/b**2 + 1.0) - np.sqrt((Ref_age - k)**2/b**2 + 1.0)
            err_factor = 0.0 # no errors available (yet)
        else:
            print "Correction method not recognized; no changes made."
            return
        
        # apply correction and error:
        abund[i] -= corr_factor
        err[i] = np.sqrt(err[i]**2 + err_factor**2)
        if not silent:
	        print "GCE correction of {0:6.3} +/- {1:6.3} dex made.".format(-corr_factor, np.abs(err_factor))
              
    
    return abund, err, Tc
    
        
        
        
