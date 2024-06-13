# Intensity Conversion and [CII] (158 micron) Luminosity Calculation

This repository contains a Python script to convert the intensity data in a FITS file from K km/s to erg/s/cm²/sr at [C II] 158 micron, and then calculate the [C II] luminosity of the NGC 7538 region. The script performs unit conversions, handles the data, and calculates the luminosity in both Watts and solar luminosities.

## Requirements

- Python 3.x
- Astropy
- Numpy

You can install the necessary Python packages using:

```sh
pip install astropy numpy
```

## Script Overview

The script `calculate_luminosity.py` performs the following steps:

1. **Import Libraries**: Imports necessary libraries for file handling, data manipulation, and FITS file operations.
2. **Set Up Variables**: Defines the source name and file names for input and output FITS files.
3. **Open FITS File**: Reads the FITS file and extracts the data and header.
4. **Convert Intensity Units**: Converts the intensity data from K km/s to erg/s/cm²/sr using a given conversion factor.
5. **Extract Frequency Axis Information**: Uses the WCS information to extract the frequency axis information.
6. **Calculate Pixel Area**: Computes the pixel area in steradians.
7. **Integrate Intensity and Calculate Flux Density**: Integrates the intensity over the frequency axis to obtain flux density, excluding negative values.
8. **Calculate Total Flux**: Computes the total flux by summing the flux density values.
9. **Convert Total Flux to Luminosity**: Converts the total flux to luminosity in erg/s, then to watts, and finally to solar luminosities.
10. **Print the Result**: Outputs the calculated [C II] luminosity in solar luminosities.
11. **Close the FITS File**: Ensures the FITS file is properly closed after processing.

## Detailed Steps

### 1. Import Libraries

```python
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
import numpy as np
```

### 2. Set Up Variables

Define the source name and file names for the input and output FITS files.

### 3. Open FITS File

Reads the FITS file and extracts the data and header.

```python
# Open the FITS file
hdu_list = fits.open('NGC7538_CII_subcube.fits')
data = hdu_list[0].data
header = hdu_list[0].header
```

### 4. Convert Intensity Units

Convert intensity from K km/s to erg/s/cm²/sr using the conversion factor of 7.0e-6 for the [CII] 158 micrin line observations.

```python
# Convert intensity from K km/s to erg/s/cm²/sr
conversion_factor = 7e-6  # Conversion factor from Goicoechea et al. (2015)
intensity *= conversion_factor  # Apply conversion
```

### 5. Extract Frequency Axis Information

Use the WCS object to extract the frequency axis information.

```python
# Create a WCS object for the FITS file
wcs = WCS(header)

# Extract the frequency axis information using the WCS object
naxis3 = header['NAXIS3']
pix_coords = np.arange(naxis3)

# Get the frequency values for each pixel along the spectral axis
world_coords = wcs.all_pix2world([0]*naxis3, [0]*naxis3, pix_coords, 0)
frequency = world_coords[2] * u.Hz

# Print the frequency array for debugging
print("Frequency array (Hz):", frequency)
```

### 6. Calculate Pixel Area

Calculate the pixel area in steradians.

```python
# Calculate the pixel area in steradians using absolute values to avoid negative areas
cdelt1 = np.abs(wcs.wcs.cdelt[0])  # Degrees per pixel along axis 1
cdelt2 = np.abs(wcs.wcs.cdelt[1])  # Degrees per pixel along axis 2
pixel_area_sr = (cdelt1 * u.deg).to(u.rad).value * (cdelt2 * u.deg).to(u.rad).value
```

### 7. Integrate Intensity and Calculate Flux Density

Integrate the intensity over the frequency axis to get flux density.

```python
# Integrate intensity over the frequency axis to get flux density
flux_density = np.nansum(intensity, axis=0)  # Summing over the frequency axis

# Exclude negative values in flux_density
flux_density = np.where(flux_density > 0, flux_density, 0)
```

### 8. Calculate Total Flux

Calculate the total flux (flux density * pixel area in steradians).

```python
# Calculate the total flux (flux density * pixel area in steradians)
total_flux = np.nansum(flux_density * pixel_area_sr)

# Debugging output for flux
print(f"Total Flux: {total_flux} erg/s/cm²")
```

### 9. Convert Total Flux to Luminosity

Convert the total flux to luminosity in erg/s, then to watts, and finally to solar luminosities.

```python
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
```

### 10. Print the Result

Output the calculated [C II] luminosity in solar luminosities.

```python
# Print the result
print(f'[C II] Luminosity: {solar_luminosity:.2e} L_sun')
```

### 11. Close the FITS File

Ensure the FITS file is properly closed after processing.

```python
# Close the FITS file
hdu_list.close()
```

## Running the Script

Ensure that the FITS file `NGC7538_CII_subcube.fits` is in the same directory as the script. Run the script using:

```sh
python calculate_luminosity.py
```
## Result 

- Total Flux: 2.70e-08 erg/s/cm²
- [C II] Luminosity: 5.93e+03 L_sun

## License

This project is licensed under public licences. 
