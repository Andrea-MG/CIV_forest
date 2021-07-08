%pylab inline
import numpy as np 
import matplotlib.pyplot as plt
import scipy as sp
from astropy.io import fits
from astropy.table import Table, join
import desispec.io
from desispec.coaddition import coadd_cameras,resample_spectra_lin_or_log

hdu_truth_0=fits.open('/global/cscratch1/sd/andreamg/v9/CIV-forest/quick-2.0/truth-16-0.fits')
hdu_zbest_0=fits.open('/global/cscratch1/sd/andreamg/v9/CIV-forest/quick-2.0/zbest-16-0.fits')

hdu_truth_1=fits.open('/global/cscratch1/sd/andreamg/v9/CIV-forest/quick-2.1/truth-16-0.fits')
hdu_zbest_1=fits.open('/global/cscratch1/sd/andreamg/v9/CIV-forest/quick-2.1/zbest-16-0.fits')

hdu_truth_10=fits.open('/global/cscratch1/sd/andreamg/v9/CIV-forest/quick-2.2/truth-16-0.fits')
hdu_zbest_10=fits.open('/global/cscratch1/sd/andreamg/v9/CIV-forest/quick-2.2/zbest-16-0.fits')

hdu_truth_20=fits.open('/global/cscratch1/sd/andreamg/v9/CIV-forest/quick-2.3/truth-16-0.fits')
hdu_zbest_20=fits.open('/global/cscratch1/sd/andreamg/v9/CIV-forest/quick-2.3/zbest-16-0.fits')

hdu_truth_100=fits.open('/global/cscratch1/sd/andreamg/v9/CIV-forest/quick-2.4/truth-16-0.fits')
hdu_zbest_100=fits.open('/global/cscratch1/sd/andreamg/v9/CIV-forest/quick-2.4/zbest-16-0.fits')

z=hdu_zbest_0[1].data['Z']

wave_min = 3500
wave_max = 10000
wavelength = np.arange(wave_min,wave_max,0.19999384634318945)

specobj_100 = desispec.io.read_spectra('/global/cscratch1/sd/andreamg/v9/CIV-forest/quick-2.4/spectra-0.fits')
specobj_100 = resample_spectra_lin_or_log(specobj_100,linear_step=0.8, wave_min =wave_min+1, wave_max =wave_max-1, fast = True)
specobj_100 = coadd_cameras(specobj_100, cosmics_nsig=None)

specobj_20 = desispec.io.read_spectra('/global/cscratch1/sd/andreamg/v9/CIV-forest/quick-2.3/spectra-0.fits')
specobj_20 = resample_spectra_lin_or_log(specobj_20,linear_step=0.8, wave_min =wave_min+1, wave_max =wave_max-1, fast = True)
specobj_20 = coadd_cameras(specobj_20, cosmics_nsig=None)

specobj_10 = desispec.io.read_spectra('/global/cscratch1/sd/andreamg/v9/CIV-forest/quick-2.2/spectra-0.fits')
specobj_10 = resample_spectra_lin_or_log(specobj_10,linear_step=0.8, wave_min =wave_min+1, wave_max =wave_max-1, fast = True)
specobj_10 = coadd_cameras(specobj_10, cosmics_nsig=None)

specobj_1 = desispec.io.read_spectra('/global/cscratch1/sd/andreamg/v9/CIV-forest/quick-2.1/spectra-0.fits')
specobj_1 = resample_spectra_lin_or_log(specobj_1,linear_step=0.8, wave_min =wave_min+1, wave_max =wave_max-1, fast = True)
specobj_1 = coadd_cameras(specobj_1, cosmics_nsig=None)

flux_coadd_100=specobj_100.flux.values()
flux_coadd_20=specobj_20.flux.values()
flux_coadd_10=specobj_10.flux.values()
flux_coadd_1=specobj_1.flux.values()

cont_interp = sp.interpolate.interp1d(wavelength,hdu_truth_0[3].data['TRUE_CONT'])

continuum_100 = cont_interp(specobj_100.wave['brz'])
continuum_20 = cont_interp(specobj_20.wave['brz'])
continuum_10 = cont_interp(specobj_10.wave['brz'])
continuum_1 = cont_interp(specobj_1.wave['brz'])

"""
plt.figure(figsize=(20,10))
plt.plot(specobj_20.wave["brz"],specobj_100.flux['brz'][2]/continuum_100[2], linewidth=1, color="lightskyblue", alpha=0.7, label="x100")
plt.plot(specobj_20.wave["brz"],specobj_20.flux['brz'][2]/continuum_20[2], linewidth=1, color="cornflowerblue", alpha=0.7, label="x20")
plt.plot(specobj_20.wave["brz"],specobj_10.flux['brz'][2]/continuum_10[2], linewidth=1, color="royalblue", alpha=0.7, label="x10")
plt.plot(specobj_20.wave["brz"],specobj_1.flux['brz'][2]/continuum_1[2], linewidth=1, color="navy", alpha=0.7, label="x1")
plt.legend()
plt.xlabel("Wavelength [$\AA$]")
plt.ylabel("Transmitted Flux")
plt.savefig("CIV_transmission.png")

plt.figure(figsize=(20,10))
plt.plot(specobj_20.wave["brz"],specobj_100.flux['brz'][10]/continuum_100[10], linewidth=1, color="lightskyblue", alpha=0.7, label="x100")
plt.plot(specobj_20.wave["brz"],specobj_20.flux['brz'][10]/continuum_20[10], linewidth=1, color="cornflowerblue", alpha=0.7, label="x20")
plt.plot(specobj_20.wave["brz"],specobj_10.flux['brz'][10]/continuum_10[10], linewidth=1, color="royalblue", alpha=0.7, label="x10")
plt.plot(specobj_20.wave["brz"],specobj_1.flux['brz'][10]/continuum_1[10], linewidth=1, color="navy", alpha=0.7, label="x1")
plt.legend()
#plt.xlim(1000,2500)
plt.ylim(0.99,1.01)
"""

z_min=np.min(z)
z_max=np.max(z)
print("z_min = ", z_min)
print("z_max = ", z_max)
zbins = 100
sum_flux = np.zeros(zbins)

#CIV forest:

lambda_min = 1420.0
lambda_max = 1520.0
lambda_min_obs = lambda_min*(1.0 + z)
lambda_max_obs = lambda_max*(1.0 + z)
print("lambda min obs")
print(lambda_min_obs)
print("lambda max obs")
print(lambda_max_obs)

wavelength_CIV_rf = [np.arange(lambda_min,lambda_max,1)]
#np.shape(forests)
print(wavelength_CIV_rf)

all_forests = []
all_forests_rf = []
all_fluxes_100 = []
all_fluxes_20 = []
all_fluxes_10 = []
all_fluxes_1 = []
all_cont_100 = []
all_cont_20 = []
all_cont_10 = []
all_cont_1 = []
all_trans_100 = []
all_trans_20 = []
all_trans_10 = []
all_trans_1 = []

# TRANSMISSION

trans_100 = specobj_100.flux['brz']/continuum_100
trans_20 = specobj_20.flux['brz']/continuum_20
trans_10 = specobj_10.flux['brz']/continuum_10
trans_1 = specobj_1.flux['brz']/continuum_1

for i in range(0,len(lambda_min_obs)):
    CIV_mask_1 = (specobj_1.wave["brz"] < lambda_max_obs[i]) & (specobj_1.wave["brz"] > lambda_min_obs[i])
    CIV_mask_10 = (specobj_10.wave["brz"] < lambda_max_obs[i]) & (specobj_1.wave["brz"] > lambda_min_obs[i])
    CIV_mask_20 = (specobj_20.wave["brz"] < lambda_max_obs[i]) & (specobj_1.wave["brz"] > lambda_min_obs[i])
    CIV_mask_100= (specobj_100.wave["brz"] < lambda_max_obs[i]) & (specobj_1.wave["brz"] > lambda_min_obs[i])
  

    trans_100[i][~CIV_mask_100] = 1.
    trans_20[i][~CIV_mask_20] = 1.
    trans_10[i][~CIV_mask_10] = 1.
    trans_1[i][~CIV_mask_1] = 1.

# MEAN FLUX
    
mean_flux_transmission_100 = np.mean(trans_100,axis=0)
mean_flux_transmission_20 = np.mean(trans_20,axis=0)
mean_flux_transmission_10 = np.mean(trans_10,axis=0)
mean_flux_transmission_1 = np.mean(trans_1,axis=0)

"""
plt.figure(figsize=(20,10))
plt.plot(specobj_100.wave["brz"],mean_flux_transmission_100)
plt.xlabel("Wavelength [$\AA$]")
plt.ylabel("Mean Flux")
plt.savefig("meanflux.png")

plt.figure(figsize=(20,10))
plt.plot(specobj_100.wave["brz"],trans_100[30])
plt.plot(specobj_20.wave["brz"],trans_20[30])
plt.plot(specobj_10.wave["brz"],trans_10[30])
plt.plot(specobj_1.wave["brz"],trans_1[30])
plt.xlabel("Wavelength [$\AA$]")

plt.figure(figsize=(20,10))
plt.plot(specobj_100.wave["brz"],mean_flux_transmission_100, alpha=0.7, color="lightskyblue", label="x100")
plt.plot(specobj_20.wave["brz"],mean_flux_transmission_20, alpha=0.7, color= "cornflowerblue", label="x20")
plt.plot(specobj_10.wave["brz"],mean_flux_transmission_10, alpha=0.7, color="royalblue", label="x10")
plt.plot(specobj_1.wave["brz"],mean_flux_transmission_1, alpha=0.7, color="navy", label="x1")
plt.legend()
plt.xlabel("Wavelength [$\AA$]")
plt.savefig("meanflux_allCIV.png")
"""

# Deltas
delta_100 = ((specobj_100.flux['brz']) / (continuum_100 * (mean_flux_transmission_100))) - 1.0
delta_20  = ((specobj_20.flux['brz'])  / (continuum_20  * (mean_flux_transmission_20)))  - 1.0
delta_10  = ((specobj_10.flux['brz'])  / (continuum_10  * (mean_flux_transmission_10)))  - 1.0
delta_1   = ((specobj_1.flux['brz'])   / (continuum_1   * (mean_flux_transmission_1)))   - 1.0

plt.figure(figsize=(20,10))
plt.plot(specobj_100.wave["brz"],delta_100[30])
plt.plot(specobj_20.wave["brz"],delta_20[30])
plt.plot(specobj_10.wave["brz"],delta_10[30])
plt.plot(specobj_1.wave["brz"],delta_1[30])
