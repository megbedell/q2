from . import moog
from .star import Star
import numpy as np
import datetime
import logging
from scipy import interpolate
import os
from .config import *
from .tools import read_csv
from collections import OrderedDict

logger = logging.getLogger(__name__)


def get_all(Data, output_file, species_ids=None, reference=None, grid='odfnew',
        errors=False, nlte=True):
    print('------------------------------------------------------')
    print('Initializing ...')
    start_time = datetime.datetime.now()
    print('- Date and time: '+start_time.strftime('%d-%b-%Y, %H:%M:%S'))
    print('- Model atmospheres: '+grid)
    print('- Star data: '+Data.star_data_fname)
    print('- Line list: '+Data.lines_fname)
    if reference:
        print('- Reference star: '+reference)
    print('------------------------------------------------------')
    if reference:
        ref = Star(reference)
        ref.get_data_from(Data)
        if hasattr(ref, 'feh_model'):           #####
            ref.feh = getattr(ref, 'feh_model') #####
        ref.get_model_atmosphere(grid)
    else:
        ref = None
    fout = open(output_file, 'w')
    header = 'id'
    if species_ids == None:
        species_codes = sorted(set(Data.lines['species']))
        species_ids = getsp_ids(species_codes)
        print('"species_ids" not provided')
        print('Lines found for the following species: '+\
              ','.join(species_ids))
        print('')
    for species_id in species_ids:
        header += ','+species_id+',e_'+species_id+',n_'+species_id
        if reference:
            header += ',['+species_id+'],e_['+species_id+\
                      '],n_['+species_id+']'
        if errors:
            header += ',err_'+species_id
    fout.write(header + '\n')
    for star_id in Data.star_data['id']:
        line = star_id
        print('')
        print('*'*len(star_id))
        print(star_id)
        print('*'*len(star_id))
        s = Star(star_id)
        try:
            s.get_data_from(Data)
            if hasattr(s, 'feh_model'):
                s.feh = getattr(s, 'feh_model')
            s.get_model_atmosphere(grid)
        except:
            print('No data available')
            logger.warning('Could not get all the necessary data')
            line += ','*(len(species_ids)*2)
            if reference:
                line += ','*(len(species_ids)*2)
            fout.write(line+'\n')
            continue

        print('Using [Fe/H] = {0:6.3f} for the model atmosphere'.format(s.feh))
        get_one(s, species_ids, ref, errors=errors, nlte=nlte)
        for species_id in species_ids:
            print('\n'+species_id+'\n'+'-'*len(species_id))
            if not hasattr(s, species_id):
                print('No data available')
                logger.warning('There are no '+species_id+' abundances '+\
                               'for this star')
                line += ',,'
                if reference:
                    line += ',,'
                continue
            mab = np.mean(getattr(s, species_id)['ab'])
            sab = np.std(getattr(s, species_id)['ab'])
            nab = len(getattr(s, species_id)['ab'])
            print("ABS = {0:6.3f} +/- {1:6.3f} , n = {2:.0f}".\
                  format(mab, sab, nab))
            line += ',{0:.3f},{1:.3f},{2:.0f}'.format(mab, sab, nab)
            if reference:
                da = getattr(s, species_id)['difab']
                da = np.array(da, dtype=np.float) #convert None to np.nan
                mda = np.ma.masked_array(da, np.isnan(da))
                mdifab = np.mean(mda)
                sdifab = np.std(mda)
                ndifab = mda.count()
                print("DIF = {0:6.3f} +/- {1:6.3f} , n = {2:.0f}".\
                      format(mdifab, sdifab, ndifab))
                line += ',{0:.3f},{1:.3f},{2:.0f}'.\
                        format(mdifab, sdifab, ndifab)
                if errors:
                    print("ERR = {0:5.3f} (DIF)".\
                          format(getattr(s, species_id)['err_difab']))
                    line += ',{0:.3f}'.\
                            format(getattr(s, species_id)['err_difab'])
            else:
                mdifab = 0
                if errors:
                    print("ERR = {0:5.3f} (ABS)".\
                          format(getattr(s, species_id)['err_ab']))
                    line += ',{0:.3f}'.\
                            format(getattr(s, species_id)['err_ab'])
            print('')
            llhd1 = 'Wavelength   ABS    RES '
            llhd2 = '----------  ----- ------'
            if reference:
                llhd1 += '   DIF    RES '
                llhd2 += '  -----  -----'
            print(llhd1+'\n'+llhd2)
            for wi, ab, difab in \
                zip(getattr(s, species_id)['ww'],
                    getattr(s, species_id)['ab'],
                    getattr(s, species_id)['difab']):
                if reference and difab != None:
                    print("{0:10.4f} {1:6.3f} {2:6.3f} {3:6.3f} {4:6.3f}".\
                          format(wi, ab, ab-mab, difab, difab-mdifab))
                else:
                    print("{0:10.4f} {1:6.3f} {2:6.3f}".\
                          format(wi, ab, ab-mab))

        fout.write(line+'\n')
    fout.close()
    print('')
    print('------------------------------------------------------')
    end_time = datetime.datetime.now()
    print('- Date and time: '+end_time.strftime('%d-%b-%Y, %H:%M:%S'))
    delta_t = (end_time - start_time).seconds
    hours, remainder = divmod(delta_t, 3600)
    minutes, seconds = divmod(remainder, 60)
    print('- Time elapsed: %sH %sM %sS' % (hours, minutes, seconds))
    print('Done!')
    print('------------------------------------------------------')
    print('')


def get_one(Star, species_ids=None, Ref=object, silent=True, errors=False, nlte=True):
    logger.info('Working on: '+Star.name)
    if species_ids == None:
        species_codes = sorted(set(Star.linelist['species']))
        species_ids = getsp_ids(species_codes)
        if not silent:
            print('"species_ids" not provided')
            print('Lines found for the following species: '+\
                  ','.join(species_ids))
            print('')
    for species_id in species_ids:
        species = getsp(species_id)
        if not silent:
            print("*** Begin "+species_id+":")
        if species == None:
            logger.warning('Not doing calculations for: '+species_id)
            continue
        logger.info('Working on: '+species_id)
        ab = moog.abfind(Star, species, species_id)
        if ab:
            setattr(Star, species_id, ab)
        else:            
            logger.warning('Did not calculate '+species_id+' abundances')
            continue            

        if (species_id == 'OI') & nlte:
            if not silent:
                print('777 nm oxygen abundances will be NLTE corrected')
            ao = []
            for wx in [7771.94, 7774.16, 7775.39]:
                k = np.where(abs(Star.OI['ww']-wx) < 0.05)
                if len(k[0]) == 1:
                    ao.append(np.mean(Star.OI['ab'][k]))
                else:
                    ao.append(0)
            aon = nlte_triplet(Star.teff, Star.logg, Star.feh, ao,
                               silent=silent)
            k= np.where(np.array(ao) > 0)
            getattr(Star, species_id)['ab'] = aon[k]

        getattr(Star, species_id)['ref'] = None
        if hasattr(Ref, 'name'):
            logger.info('Diferential analysis: '+Ref.name)
            if Star.name == Ref.name:
                logger.warning('Reference star object redefined!')
                Ref = Star
            if not hasattr(Ref, species_id):
                logger.info('Calculating reference star abundances: '+Ref.name)
                setattr(Ref, species_id, moog.abfind(Ref, species, species_id))
                

                if (species_id == 'OI') & nlte:
                    if not silent:
                        print('777 nm oxygen abundances will be NLTE '\
                              +'corrected (Reference)')
                    ao = []
                    for wx in [7771.94, 7774.16, 7775.39]:
                        k = np.where(abs(Ref.OI['ww']-wx) < 0.05)
                        if len(k[0]) == 1:
                            ao.append(np.mean(Ref.OI['ab'][k]))
                        else:
                            ao.append(0)
                    aon = nlte_triplet(Ref.teff, Ref.logg, Ref.feh, ao,
                                       silent=silent)
                    k= np.where(np.array(ao) > 0)
                    getattr(Ref, species_id)['ab'] = aon[k]
            else:
                logger.info('Reference star has '+species_id+\
                            ' abundances computed already: '+Ref.name)

            ws = getattr(Star, species_id)['ww']
            wr = getattr(Ref, species_id)['ww']
            ww = np.intersect1d(ws, wr)
            k  = [i for i, w in zip(range(len(ws)), ws) if w in ww]
            kr = [i for i, w in zip(range(len(wr)), wr) if w in ww]
            a = getattr(Star, species_id)['ab'][k] - \
                getattr(Ref, species_id)['ab'][kr]
            ax, ix = [], 0
            for wx in ws:
                if wx in ww:
                    ax.append(a[ix])
                    ix += 1
                else:
                    ax.append(None)
            getattr(Star, species_id)['difab'] = ax
            getattr(Star, species_id)['ref'] = Ref.name

        if not silent:
            ab = getattr(Star, species_id)['ab']
            difab = getattr(Star, species_id)['difab']

            aa = np.array(ab, dtype=np.float) #convert None to np.nan
            maa = np.ma.masked_array(aa, np.isnan(aa))

            da = np.array(difab, dtype=np.float) #convert None to np.nan
            mda = np.ma.masked_array(da, np.isnan(da))

            print("A({0})  = {1:6.3f} +/- {2:5.3f} (# of lines = {3})".\
                  format(species_id, np.mean(maa), np.std(maa), maa.count()))

            if hasattr(Ref, 'name'):
                print("[{0}/H] = {1:6.3f} +/- {2:5.3f} (# of lines = {3})".\
                      format(species_id, np.mean(mda), np.std(mda), mda.count()))

        if errors:
            error(Star, species_id, Ref=Ref, silent=silent)

        if not silent:
            print('---' + species_id + ' done')

    if not silent and len(species_ids) >= 1:
        print('All species completed')



def error(Star_in, species_id, Ref=object, silent=True):
    s = Star()
    s.__dict__ = Star_in.__dict__.copy()

    if not silent:
        print('-----------------------------')
        print('Error propagation for '+species_id+':')

    try:
        Ref.model_atmosphere_grid
        dab = getattr(Star_in, species_id)['difab']
        l2l_sct = np.std(dab)/np.sqrt(max([len(dab),2])-1)
        abx = 'difab'
    except:
        try:
            ab = getattr(Star_in, species_id)['ab']
            l2l_sct = np.std(ab)/np.sqrt(max([len(ab),2])-1)
            abx = 'ab'
        except:
            logger.error('Must calculate abundances before errors')
            return None

    if hasattr(s, 'err_teff'):
        if s.err_teff > 0:
            s.teff += s.err_teff
            s.get_model_atmosphere(s.model_atmosphere_grid)
            get_one(s, [species_id], Ref=Ref)
            ap = np.mean(getattr(s, species_id)[abx])
            s.teff -= 2*s.err_teff
            s.get_model_atmosphere(s.model_atmosphere_grid)
            get_one(s, [species_id], Ref=Ref)
            am = np.mean(getattr(s, species_id)[abx])
            a_teff = abs(ap-am)/2.
            s.teff += s.err_teff
        else:
            a_teff = 0.
    else:
        a_teff = 0.

    if hasattr(s, 'err_logg'):
        if s.err_logg > 0:
            s.logg += s.err_logg
            s.get_model_atmosphere(s.model_atmosphere_grid)
            get_one(s, [species_id], Ref=Ref)
            ap = np.mean(getattr(s, species_id)[abx])
            s.logg -= 2*s.err_logg
            s.get_model_atmosphere(s.model_atmosphere_grid)
            get_one(s, [species_id], Ref=Ref)
            am = np.mean(getattr(s, species_id)[abx])
            a_logg = abs(ap-am)/2.
            s.logg += s.err_logg
        else:
            a_logg = 0.
    else:
        a_logg = 0.

    if hasattr(s, 'err_feh'):
        if s.err_feh > 0:
            s.feh += s.err_feh
            s.get_model_atmosphere(s.model_atmosphere_grid)
            get_one(s, [species_id], Ref=Ref)
            ap = np.mean(getattr(s, species_id)[abx])
            s.feh -= 2*s.err_feh
            s.get_model_atmosphere(s.model_atmosphere_grid)
            get_one(s, [species_id], Ref=Ref)
            am = np.mean(getattr(s, species_id)[abx])
            a_feh = abs(ap-am)/2.
            s.feh += s.err_feh
        else:
            a_feh = 0.
    else:
        a_feh = 0.

    if hasattr(s, 'err_vt'):
        if s.err_vt > 0:
            s.vt += s.err_vt
            s.get_model_atmosphere(s.model_atmosphere_grid)
            get_one(s, [species_id], Ref=Ref)
            ap = np.mean(getattr(s, species_id)[abx])
            s.vt -= 2*s.err_vt
            s.get_model_atmosphere(s.model_atmosphere_grid)
            get_one(s, [species_id], Ref=Ref)
            am = np.mean(getattr(s, species_id)[abx])
            a_vt = abs(ap-am)/2.
            s.vt += s.err_vt
        else:
            a_vt = 0.
    else:
        a_vt = 0.

    a_tot = np.sqrt(a_teff**2+a_logg**2+a_feh**2+a_vt**2+l2l_sct**2)
    if not silent:
        print('Line to line scatter:  {0:.3f}'.format(l2l_sct))
        print('Error from Teff:       {0:.3f}'.format(a_teff))
        print('Error from logg:       {0:.3f}'.format(a_logg))
        print('Error from [Fe/H]:     {0:.3f}'.format(a_feh))
        print('Error from vt:         {0:.3f}'.format(a_vt))
        print('                      -------')
        print('Total abundance error: {0:.3f}'.format(a_tot))
        print('-----------------------------')

    try:
        Ref.model_atmosphere_grid
        getattr(Star_in, species_id)['err_difab'] = a_tot
    except:
        getattr(Star_in, species_id)['err_ab'] = a_tot

sp_map = {
          'LiI' :   3.0,
          'BeI' :   4.0,
          'BeII':   4.1,
          'BI'  :   5.0,
          'CI'  :   6.0,
          'CI2' :   6.1,
          'CH'  : 106.0,
          'CH2' : 106.1,
          'NI'  :   7.0,
          'OI'  :   8.0,
          'OI2' :   8.1,
          'FI'  :   9.0,
          'NaI' :  11.0,
          'MgI' :  12.0,
          'MgII':  12.1,
          'AlI' :  13.0,
          'SiI' :  14.0,
          'PI'  :  15.0,
          'SI'  :  16.0,
          'KI'  :  19.0,
          'CaI' :  20.0,
          'ScI' :  21.0,
          'ScII':  21.1,
          'TiI' :  22.0,
          'TiII':  22.1,
          'VI'  :  23.0,
          'CrI' :  24.0,
          'CrII':  24.1,
          'MnI' :  25.0,
          'FeI' :  26.0,
          'FeII':  26.1,
          'CoI' :  27.0,
          'NiI' :  28.0,
          'CuI' :  29.0,
          'ZnI' :  30.0,
          'RbI' :  37.0,
          'SrI' :  38.0,
          'SrII':  38.1,
          'YII' :  39.1,
          'ZrII':  40.1,
          'BaII':  56.1,
          'LaII':  57.1,
          'CeII':  58.1,
          'PrII':  59.1,
          'NdII':  60.1,
          'SmII':  62.1,
          'EuII':  63.1,
          'GdII':  64.1,
          'DyII':  66.1
          }

def getsp(species_id):
    try:
        species = sp_map[species_id]
    except:
        logger.warning('species id not recognized: '+species_id)
        return None
    return species

def getsp_ids(species_list):
    species_ids = []
    for species_code in species_list:
        try:
            species_id = [key for key in sp_map if sp_map[key] == species_code][0]
            species_ids.append(species_id)
        except:
            logger.warning('species_code '+str(species_code)+' not found')
    return species_ids

tc_map = {
          'LiI' :   1142.0,
          'BeI' :   1452.0,
          'BeII':   1452.0,
          'BI'  :    908.0,
          'CI'  :     40.0,
          'CH'  :     40.0,
          'CI2'  :    40.0,
          'CH2'  :    40.0,
          'NI'  :    123.0,
          'OI'  :    180.0,
          'OI2'  :   180.0,
          'FI'  :    734.0,
          'NaI' :    958.0,
          'MgI' :   1336.0,
          'MgII':   1336.0,
          'AlI' :   1653.0,
          'SiI' :   1310.0,
          'PI'  :   1229.0,
          'SI'  :    664.0,
          'KI'  :   1006.0,
          'CaI' :   1517.0,
          'ScI' :   1659.0,
          'ScII':   1659.0,
          'TiI' :   1582.0,
          'TiII':   1582.0,
          'VI'  :   1429.0,
          'CrI' :   1296.0,
          'CrII':   1296.0,
          'MnI' :   1158.0,
          'FeI' :   1334.0,
          'FeII':   1334.0,
          'CoI' :   1352.0,
          'NiI' :   1353.0,
          'CuI' :   1037.0,
          'ZnI' :    726.0,
          'RbI' :    800.0,
          'SrI' :   1464.0,
          'SrII':   1464.0,
          'YII' :   1659.0001,  # so it's not the same as Sc & Dy & Gd
          'ZrII':   1741.0,
          'BaII':   1455.0,
          'LaII':   1578.0,
          'CeII':   1478.0,
          'PrII':   1582.0001, # so it's not the same as Ti
          'NdII':   1602.0,
          'SmII':   1590.0,
          'EuII':   1356.0,
          'GdII':   1659.0002,  # so it's not the same as Sc & Dy & Y
          'DyII':   1659.0003  # so it's not the same as Sc & Y & Gd
          }

def gettc(species_id):
    # return 50% solar system gas condensation temperature from Table 8 of Lodders (2003)
    try:
        tc = tc_map[species_id]
    except:
        logger.warning('species id not recognized: '+species_id)
        return None
    return tc


def nlte_triplet(teff, logg, feh, ao, silent=True):
    if feh >= 0.4:
        feh = 0.4
    grid = read_csv(os.path.join(OTHER_PATH ,'nlte_triplet.csv'))

    t,g,f,dao0,dao1,dao2=[],[],[],[],[],[]
    for i in range(640):
        rg = range(i*7, i*7+7)
        x0 = interpolate.griddata(grid['ao'][rg], grid['dao0'][rg],\
                                  ao[0], method='cubic')
        x1 = interpolate.griddata(grid['ao'][rg], grid['dao1'][rg],\
                                  ao[1], method='cubic')
        x2 = interpolate.griddata(grid['ao'][rg], grid['dao2'][rg],\
                                  ao[2], method='cubic')
        x0, x1, x2 = float(x0), float(x1), float(x2)
        t.append(grid['teff'][rg[0]])
        g.append(grid['logg'][rg[0]])
        f.append(grid['feh'][rg[0]])
        dao0.append(x0)
        dao1.append(x1)
        dao2.append(x2)
    t = np.array(t)
    g = np.array(g)
    f = np.array(f)
    dao0 = np.array(dao0)
    dao1 = np.array(dao1)
    dao2 = np.array(dao2)

    tt,ff,dao00,dao11,dao22=[],[],[],[],[]
    for i in range(160):
        rg =range(i*4, i*4+4)
        x0 = interpolate.griddata(g[rg], dao0[rg], logg, method='cubic')
        x1 = interpolate.griddata(g[rg], dao1[rg], logg, method='cubic')
        x2 = interpolate.griddata(g[rg], dao2[rg], logg, method='cubic')
        x0, x1, x2 = float(x0), float(x1), float(x2)
        tt.append(t[rg[0]])
        ff.append(f[rg[0]])
        dao00.append(x0)
        dao11.append(x1)
        dao22.append(x2)
    tt = np.array(tt)
    ff = np.array(ff)
    dao00 = np.array(dao00)
    dao11 = np.array(dao11)
    dao22 = np.array(dao22)

    t,dao0,dao1,dao2=[],[],[],[]
    for i in range(16):
        rg =range(i*10, i*10+10)
        x0 = interpolate.griddata(ff[rg], dao00[rg], feh, method='cubic')
        x1 = interpolate.griddata(ff[rg], dao11[rg], feh, method='cubic')
        x2 = interpolate.griddata(ff[rg], dao22[rg], feh, method='cubic')
        x0, x1, x2 = float(x0), float(x1), float(x2)
        t.append(tt[rg[0]])
        dao0.append(x0)
        dao1.append(x1)
        dao2.append(x2)
    t = np.array(t)
    dao0 = np.array(dao0)
    dao1 = np.array(dao1)
    dao2 = np.array(dao2)

    x0 = interpolate.griddata(t, dao0, teff, method='cubic')
    x1 = interpolate.griddata(t, dao1, teff, method='cubic')
    x2 = interpolate.griddata(t, dao2, teff, method='cubic')
    x0, x1, x2 = float(x0), float(x1), float(x2)

    x0 = x0 - 0.0355
    x1 = x1 - 0.0180
    x2 = x2 - 0.0000

    if not silent:
        print('Wavelength (A) | A(O) LTE | Correction | A(O) NLTE')
        print("   7771.9      |  {0:6.3f}  |    {1:5.3f}   | {2:6.3f}".\
              format(ao[0], x0, ao[0]-x0))
        print("   7774.2      |  {0:6.3f}  |    {1:5.3f}   | {2:6.3f}".\
              format(ao[1], x1, ao[1]-x1))
        print("   7775.4      |  {0:6.3f}  |    {1:5.3f}   | {2:6.3f}".\
              format(ao[2], x2, ao[2]-x2))
    ax = [round(ao[0]-x0, 3),
          round(ao[1]-x1, 3),
          round(ao[2]-x2, 3)]

    aon = np.ma.masked_array(ax,np.isnan(ax))

    if not silent:
        print("A(O) LTE  = {0:6.3f} +/- {1:5.3f}".\
              format(np.mean(ao), np.std(ao)))
        print("A(O) NLTE = {0:6.3f} +/- {1:5.3f}".\
              format(np.mean(aon), np.std(aon)))

    return aon
