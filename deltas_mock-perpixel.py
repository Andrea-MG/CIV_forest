import numpy as np 
import matplotlib.pyplot as plt
import scipy as sp
from astropy.io import fits
from astropy.table import Table, join, vstack
import desispec.io
from desispec.coaddition import coadd_cameras,resample_spectra_lin_or_log
import fitsio
from random import seed
from random import random
import glob
import os

in_dir_spectra = '/global/cscratch1/sd/andreamg/v9/CIV-forest/quick-3.100/spectra-16/'# Where spectra-16 are
spectra_files = sorted(glob.glob(os.path.join(in_dir_spectra,'*/*/spectra*.fits')))

continuum_100 = []
z = []
fibermap_table = []
specobj_100_flux = []

wave_min = 3500
wave_max = 10000
wavelength = np.arange(wave_min,wave_max,0.19999384634318945)

for spectra in spectra_files[:2]:
#for spectra in spectra_files:
    base,spfile = os.path.split(spectra)
    
    # truth files
    truth_file = os.path.join(base,spfile.replace('spectra','truth'))
    hdu_truth_100 = fits.open(truth_file)
    cont_interp = sp.interpolate.interp1d(wavelength,hdu_truth_100[3].data['TRUE_CONT'])
    
    
    # zbest files
    zbest_file = os.path.join(base,spfile.replace('spectra','zbest'))
    hdu_zbest_100=fits.open(zbest_file)
    # We make two tables -  one for hdu[1] (z) and one for hdu[2] (fibermap)
    #z.append(hdu_zbest_100[1].data['Z'])
    z=hdu_zbest_100[1].data['Z']
    fibermap_table=hdu_zbest_100[2].data
 
    # coadd - spectra
    specobj_100 = desispec.io.read_spectra(spectra)
    specobj_100 = resample_spectra_lin_or_log(specobj_100,linear_step=0.8, wave_min =wave_min+1, wave_max =wave_max-1, fast = True)
    specobj_100 = coadd_cameras(specobj_100, cosmics_nsig=None)
    specobj_100_flux=specobj_100.flux['brz']
    specobj_100_wave=specobj_100.wave['brz']
    
    continuum_100=cont_interp(specobj_100.wave['brz'])
    
#fibermap_table_stack = np.hstack(fibermap_table)
    RA = fibermap_table['TARGET_RA']
    DEC = fibermap_table['TARGET_DEC']
    TARGETID = fibermap_table['TARGETID']
#specobj_100_flux = np.vstack(specobj_100_flux)
#continuum_100_stack = np.vstack(continuum_100)

#z=np.concatenate(z)
#print(z)

    z_min=np.min(z)
    z_max=np.max(z)
    zbins = 100
    sum_flux = np.zeros(zbins)
    lambda_min = 1420.0
    lambda_max = 1520.0
    lambda_min_obs = lambda_min*(1.0 + z)
    lambda_max_obs = lambda_max*(1.0 + z)
    wavelength_CIV_rf = [np.arange(lambda_min,lambda_max,1)]

    all_forests = []
    all_forests_rf = []
    all_fluxes_100 = []
    all_cont_100 = []
    all_trans_100 = []

    trans_100 = specobj_100_flux/continuum_100

    Flux_CIV_100 = []
    Flux_CIV_100 = np.ones(np.shape(specobj_100_flux))
    Cont_CIV_100 = np.ones(np.shape(continuum_100))

    for i in range(0,len(lambda_min_obs)):
        # Masks the CIV forest
        CIV_mask_100 = (specobj_100_wave < lambda_max_obs[i]) & (specobj_100_wave > lambda_min_obs[i])
        # Masks spectra
        Flux_CIV_100[i][CIV_mask_100] = specobj_100_flux[i][CIV_mask_100]
        # Masks continua
        Cont_CIV_100[i][CIV_mask_100] = continuum_100[i][CIV_mask_100]
        # Makes 1 whatever isn't in the CIV forest
        trans_100[i][~CIV_mask_100] = 1.

    mean_flux_transmission_100 = np.mean(trans_100,axis=0)
    mean_flux_transmission_100_matrix = np.ones(np.shape(continuum_100))

    for j in range(0,len(lambda_min_obs)):
        # Masks the CIV forest
        CIV_mask_100 = (specobj_100_wave < lambda_max_obs[j]) & (specobj_100_wave > lambda_min_obs[j])  
        # Masks mean_flux_transmission
        mean_flux_transmission_100_matrix[j][CIV_mask_100] = mean_flux_transmission_100[CIV_mask_100]

    Meanflux_CIV_100 = np.ones(np.shape(mean_flux_transmission_100))
    Meanflux_CIV_100[~CIV_mask_100] = 1.
    delta_100 = ((Flux_CIV_100) / (Cont_CIV_100 * mean_flux_transmission_100_matrix)) -1
    loglam_100 = np.log10(specobj_100_wave)

    #DELTA FILES
    seed(1)
    weights=np.ones(np.shape(loglam_100[CIV_mask_100]))
    thingid = [i for i in range(len(delta_100))]

    outdir = f"/global/cscratch1/sd/andreamg/v9/CIV-forest/quick-3.100/Picca/Delta_CIV/Delta/"
    delta_file_name = spfile.replace('spectra','delta')
    delta_file_100 = os.path.join(outdir,delta_file_name)

    results  = fitsio.FITS(delta_file_100, 'rw', clobber=True)
    #results  = fitsio.FITS(delta_file, 'rw', clobber=True)

    for i in range(len(delta_100)):

        header   = [{'name': 'RA',  'value': np.radians(RA[i]),  'comment': 'Right Ascension [rad]'},
                    {'name': 'DEC', 'value': np.radians(DEC[i]), 'comment': 'Declination [rad]'},
                    {'name': 'Z',   'value': z[i],                                 'comment': 'Redshift'},
                    {'name': 'TARGETID', 'value':TARGETID[i], 'comment': 'Target ID'},
                    {'name': 'THING_ID', 'value':thingid[i], 'comment': 'Thing ID (fake)'},
                    {'name': 'PLATE', 'value':thingid[i], 'comment': 'Plate (fake)'},
                    {'name': 'MJD', 'value':thingid[i], 'comment': 'MJD (fake)'},
                    {'name': 'FIBERID', 'value':thingid[i], 'comment': 'fiberid (fake)'},
                                    ]       
        CIV_mask_100 = (loglam_100 < np.log10(lambda_max_obs[i])) & (loglam_100 > np.log10(lambda_min_obs[i]))

        cols     = [loglam_100[CIV_mask_100], delta_100[i][CIV_mask_100], np.ones(np.shape(loglam_100[CIV_mask_100])), continuum_100[i][CIV_mask_100]]
        names    = ['LOGLAM', 'DELTA', 'WEIGHT', 'CONT']
        units    = ['LogAngstrom','', '','']
        comments = ['Log Wavelength', 'delta field', 'weight', 'continuum']


        results.write(cols,
                      names   = names,
                      header  = header,
                      comment = comments,
                      units   = units,
                      #extname  =  str(hdu_zbest_0[2].data["TARGETID"][i])
                      extname  =  str(TARGETID[i])
                     )

    results.close()
