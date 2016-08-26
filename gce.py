import numpy as np
import q2
import pdb
import matplotlib.pyplot as plt

b_map = {
# slope correction factor and error [10e-2 dex/Gyr]
          'CI'  :   (3.10,0.66),
          'CH'  :   (3.10,0.66),
          'OI'  :   (2.86,0.70),
          'NaI' :   (2.90,0.48),
          'MgI' :   (1.21,0.28),
          'MgII':   (1.21,0.28),
          'AlI' :   (1.77,0.41),
          'SiI' :   (0.95,0.27),
          'SI'  :   (1.90,0.65),
          'CaI' :   (-0.75,0.19),
          'ScI' :   (2.61,0.59),
          'ScII':   (2.61,0.59),
          'TiI' :   (0.95,0.22),
          'TiII':   (0.95,0.22),
          'VI'  :   (0.77,0.27),
          'CrI' :   (-0.38,0.16),
          'CrII':   (-0.38,0.16),
          'MnI' :   (1.84,0.32),
          'CoI' :   (2.41,0.40),
          'NiI' :   (1.89,0.30),
          'CuI' :   (3.21,0.52),
          'ZnI' :   (2.53,0.58),
          'SrI' :   (-4.84,0.83),
          'YII' :   (-3.72,0.69),
          'BaII':   (-5.93,0.89)
          }

def getb(species_id):
    try:
        b = b_map[species_id]
    except:
        #logger.warning('no GCE correction available for: '+species_id)
        return (0.0, 0.0)
    return b


def correct(Star, age, species_ids=None, Ref=None, Ref_age=0.0, silent=True, errors=False):
    if (Ref==None):
    	print "Ref not set; assuming it is the Sun (age = 4.6 Gyr)"
        Ref_age = 4.6
    if (Ref_age==0.0):
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
        err = np.append(err, np.std(species_difab) / np.sqrt(len(species_difab) - 1) )
    
    for t in set(Tc):
        # eliminate duplicate measurements of the same element
        ind = np.where(Tc == t)[0]
        if len(ind) == 2:  # won't take care of 3+ states of the same thing
            (abund[ind[0]], err[ind[0]]) = np.average(abund[ind], weights=err[ind], returned=True)
            abund = np.delete(abund, ind[1])
            err = np.delete(err, ind[1])
            Tc = np.delete(Tc, ind[1])
            species_ids = np.delete(species_ids[1])

    for i,species_id in enumerate(species_ids):
        # make corrections
        if not silent:
            print "*** Begin "+species_id+":"
        # make GCE correction:
        (b, b_err) = getb(species_id)
        if (b == 0.0):
	        if not silent:
		        print "No GCE correction available."
	        continue
        corr_factor = (age - Ref_age)*b*1e-2
        abund[i] -= corr_factor
        # adjust abundance error:
        err_factor = (age - Ref_age)*err_b*1e-2
        err[i] = np.sqrt(err[i]**2 + err_factor**2)
        if not silent:
	        print "GCE correction of {0:6.3} +/- {1:6.3} dex made.".format(-corr_factor, np.abs(err_factor))

              
    return abund, err, Tc
    
    
        
        
