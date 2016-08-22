import numpy as np
import q2
import pdb
import matplotlib.pyplot as plt

b_map = {
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
          'ScI' :    5.6,
          'ScII':    5.6,
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

def getb(species_id):
    try:
        b = b_map[species_id]
    except:
        #logger.warning('no GCE correction available for: '+species_id)
        return 0.0
    return b


c_map = {
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
          'ScI' :   10.2,
          'ScII':   10.2,
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

def getc(species_id):
    try:
        c = c_map[species_id]
    except:
        #logger.warning('no GCE correction available for: '+species_id)
        return np.inf
    return c

def correct(Star, age, species_ids=None, Ref=None, Ref_age=0.0, silent=True, errors=False):
    if (hasattr(Ref, 'name') and Ref_age == 0.0):
    	print "Ref_age keyword must be set!"
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
        #species = q2.abundances.getsp(species_id)
        if not silent:
            print "*** Begin "+species_id+":"
        # get abundances:
        Tc = np.append(Tc, q2.abundances.gettc(species_id))
        species_difab = np.array(getattr(Star, species_id)['difab'])
        species_difab = species_difab[species_difab != np.array(None)]  #remove Nones from where ref star was unavailable
        abund = np.append(abund, np.mean(species_difab))
        err = np.append(err, np.std(species_difab) / np.sqrt(len(species_difab) - 1) )
        # make GCE correction:
        species_b = getb(species_id)
        species_c = getc(species_id)
        if (species_b == None or species_c == None):
	        if not silent:
		        print "No GCE correction available."
	        continue
        corr_factor = np.sqrt((age - species_b)**2/species_c**2 + 1.0) - np.sqrt((Ref_age - species_b)**2/species_c**2 + 1.0)
        abund[-1] -= corr_factor
        if not silent:
	        print "GCE correction of {0:6.3} dex made.".format(-corr_factor)

    for t in set(Tc):
         ind = np.where(Tc == t)[0]
         if len(ind) == 2:  #assumes there's only one other state
             (abund[ind[0]], err[ind[0]]) = np.average(abund[ind], weights=err[ind], returned=True)
             abund = np.delete(abund, ind[1])
             err = np.delete(err, ind[1])
             Tc = np.delete(Tc, ind[1])
              
    return abund, err, Tc
    
    
        
        
