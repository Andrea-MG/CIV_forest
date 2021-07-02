# %pylab inline
import numpy as np 
import matplotlib.pyplot as plt
import scipy as sp
from astropy.io import fits
from astropy.table import Table, join
import desispec.io
from desispec.coaddition import coadd_cameras,resample_spectra_lin_or_log

# VARIABLES

mock_path = '/global/cscratch1/sd/andreamg/v9/CIV-forest/quick-2.0/'

# LOAD SPECTTRA

### hdu_spectra=fits.open(mock_path + 'spectra-0.fits')
hdu_truth=fits.open(mock_path + 'truth-16-0.fits')
hdu_zbest=fits.open(mock_path + 'zbest-16-0.fits')

z=hdu_zbest[1].data['Z']

wave_min = 3500
wave_max = 10000
wavelength = np.arange(wave_min,wave_max,0.19999384634318945)

# COADD SPECTRA

specobj = desispec.io.read_spectra(mock_path + 'spectra-0.fits')
specobj = resample_spectra_lin_or_log(specobj,linear_step=0.8, wave_min =wave_min+1, wave_max =wave_max-1, fast = True)
specobj = coadd_cameras(specobj, cosmics_nsig=None)

flux_coadd=specobj.flux.values()
cont_interp = sp.interpolate.interp1d(wavelength,hdu_truth[3].data['TRUE_CONT'])

continuum = cont_interp(specobj.wave['brz'])

z_min=np.min(z)
z_max=np.max(z)
print("z_min = ", z_min)
print("z_max = ", z_max)
zbins = 100
sum_flux = np.zeros(zbins)

# TRANSMISSION IN CIV FOREST

lambda_min = 1420.0
lambda_max = 1520.0
lambda_min_obs = lambda_min*(1.0 + z)
lambda_max_obs = lambda_max*(1.0 + z)
print("lambda min obs")
print(lambda_min_obs)
print("lambda max obs")
print(lambda_max_obs)

wavelength_CIV_rf = [np.arange(lambda_min,lambda_max,1)]

all_forests = []
all_forests_rf = []
all_fluxes = []
all_cont = []
all_trans = []


for i in range(0,len(lambda_min_obs)):
    CIV_mask = (specobj.wave["brz"] < lambda_max_obs[i]) & (specobj.wave["brz"] > lambda_min_obs[i])

    # interpolate continuum
    cont_interp = sp.interpolate.interp1d(wavelength,hdu_truth[3].data['TRUE_CONT'][i]) # Value cont_interp in specobj.wave["brz"]
    # arrays of continua of forests
    cont = cont_interp(specobj.wave["brz"])
 
    trans = specobj.flux['brz'][i]/cont
    trans[~CIV_mask] = 1.
    all_trans.append(trans)
    
# MEAN FLUX
    
mean_flux_transmission = np.mean(all_trans,axis=0)

# COMPUTATION OF DELTAS

delta = (specobj.flux['brz'])/(mean_flux_transmission)-1.
plt.figure(figsize=(20,10))
plt.plot(delta, color="red", label="x100", alpha=0.7)
plt.show()

# WRITE IN FILE

