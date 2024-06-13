#!/usr/bin/env python3

from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
import numpy as np

# Open the FITS file
hdu_list = fits.open('NGC7538_CII_subcube.fits')
data = hdu_list[0].data
header = hdu_list[0].header

# Create a WCS object for the FITS file
wcs = WCS(header)

# Check for NaN values in the data and handle them
intensity = np.nan_to_num(data)  # Convert NaNs to zero

# Convert intensity from K km/s to erg/s/cm2/sr
conversion_factor = 7e-6  # Conversion factor for [CII] 158 micron
intensity *= conversion_factor  # Apply conversion

# Extract the frequency axis information using the WCS object
# The third axis is typically the spectral axis
naxis3 = header['NAXIS3']
pix_coords = np.arange(naxis3)

# Get the frequency values for each pixel along the spectral axis
world_coords = wcs.all_pix2world([0]*naxis3, [0]*naxis3, pix_coords, 0)
frequency = world_coords[2] * u.Hz

# Print the frequency array for debugging
print("Frequency array (Hz):", frequency)

# Calculate the pixel area in steradians using absolute values to avoid negative areas
cdelt1 = np.abs(wcs.wcs.cdelt[0])  # Degrees per pixel along axis 1
cdelt2 = np.abs(wcs.wcs.cdelt[1])  # Degrees per pixel along axis 2
pixel_area_sr = (cdelt1 * u.deg).to(u.rad).value * (cdelt2 * u.deg).to(u.rad).value

# Integrate intensity over the frequency axis to get flux density
flux_density = np.nansum(intensity, axis=0)  # Summing over the frequency axis

# Calculate the total flux (flux density * pixel area in steradians)
total_flux = np.nansum(flux_density * pixel_area_sr)

# Debugging output for flux
print(f"Total Flux: {total_flux} erg/s/cm²")

# Convert total flux to luminosity
distance = 2.65 * 10**3 * u.pc  # Distance to NGC 7538 in parsecs
distance_cm = distance.to(u.cm).value  # Convert parsecs to cm

# Ensure total_flux is positive
if total_flux <= 0:
    print("Total flux is zero or negative. Please check the intensity values.")
    hdu_list.close()
    exit()
    
luminosity = total_flux * 4 * np.pi * distance_cm**2  # in erg/s

# Convert luminosity from erg/s to watts
luminosity_watts = luminosity * 1e-7

# Convert watts to solar luminosities (1 L_sun = 3.828 x 10^26 W)
solar_luminosity = luminosity_watts / (3.828 * 10**26)

# Print the result
print(f'[C II] Luminosity: {solar_luminosity:.2e} L_sun')

# Close the FITS file
hdu_list.close()
