
# 3D Plot of log(NH), log(C+), and Distance

This script generates an interactive 3D plot of logarithmic hyrodgen column density, logarithmic C+ 158 micron emission, and linear distance of each pixel to the ionizing star (IRS 6) in NGC 7538 values using data from FITS files. The final plot is saved as an HTML file for interactive viewing.

## Requirements

Ensure you have the following Python packages installed:

```sh
pip install numpy pandas plotly astropy reproject
```

## Usage

Place the following FITS files in the same directory as this script:

- `NGC7538_CII_moment0_cropped.fits`
- `NGC7538_Tdust_tau160_NH_220px.fits`

Run the script to generate the 3D plot and save it as an HTML file.

## Script Description

```python
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
import os

# Load the CII FITS file
cii_filename = "NGC7538_CII_moment0_cropped.fits"
cii_hdu = fits.open(cii_filename)
cii_data = cii_hdu[0].data
cii_header = cii_hdu[0].header
cii_wcs = WCS(cii_header)

# Load the NH FITS file
nh_filename = "NGC7538_Tdust_tau160_NH_220px.fits"
nh_hdu = fits.open(nh_filename)
nh_data = nh_hdu[0].data
nh_header = nh_hdu[0].header

# Extract the 4th slice from the 3D NH data
nh_data_2d = nh_data[3]  # Adjusted index as needed

# Reproject the 2D NH slice to match the CII WCS
nh_wcs = WCS(nh_header, naxis=2)
nh_reprojected, footprint = reproject_interp((nh_data_2d, nh_wcs), cii_wcs, shape_out=cii_data.shape)

# Calculate the distance of each pixel to the ionizing star
star_position = (130.7244, 82.69072) #taken from FITS, but use WCS instead. 
distance_kpc = 2.65  # Distance to NGC 7538 in kpc
pixel_scale = cii_header['CDELT2'] * 3600  # pixel scale in arcsec/pixel
arcsec_to_rad = np.pi / (180 * 3600)  # arcsec to rad conversion factor

# Convert pixel distances to physical distances in pc
y, x = np.indices(cii_data.shape)
distances = np.sqrt((x - star_position[0])**2 + (y - star_position[1])**2)
distances_rad = distances * pixel_scale * arcsec_to_rad
distances_pc = distances_rad * distance_kpc * 1000  # Convert to pc

# Flatten the data for plotting
NH_flat = nh_reprojected.flatten()
distance_flat = distances_pc.flatten()
C_plus_flat = cii_data.flatten()

# Mask invalid values
mask = ~np.isnan(NH_flat) & ~np.isnan(distance_flat) & ~np.isnan(C_plus_flat)
NH_flat = NH_flat[mask]
distance_flat = distance_flat[mask]
C_plus_flat = C_plus_flat[mask]

# Log-transform NH and C+
log_NH = np.log10(NH_flat)
log_C_plus = np.log10(C_plus_flat)

# Create an interactive 3D plot
fig = go.Figure(data=[go.Scatter3d(
    x=log_NH,
    y=log_C_plus,
    z=distance_flat,
    mode='markers',
    marker=dict(
        size=5,
        color=distance_flat,  # Set color to log(C+)
        colorscale='Viridis',
        opacity=0.8,
        colorbar=dict(title='Distance')
    )
)])

# Add labels
fig.update_layout(
    scene=dict(
        xaxis_title='log(NH)',
        yaxis_title='log(C+)',
        zaxis_title='Distance (pc)'
    ),
    title='3D Plot of log(NH), log(C+), and Distance',
)

# Save the plot as an HTML file
fig.write_html("3D_plot_logNH_logCplus_distance.html")
```

### Step-by-Step Explanation

1. **Import Libraries**:
   The script imports necessary libraries: `numpy`, `pandas`, `plotly`, `astropy.io`, and `reproject`.

2. **Load FITS Files**:
   The script opens the CII and NH FITS files using `astropy.io.fits`.

3. **Extract Data**:
   It extracts the 4th slice from the 3D NH data.

4. **Reproject Data**:
   The NH data is reprojected to match the WCS of the CII data using `reproject_interp`.

5. **Calculate Distance**:
   The script calculates the distance of each pixel to the ionizing star in parsecs (pc).

6. **Flatten Data**:
   The NH, distance, and C+ data are flattened for plotting.

7. **Mask Invalid Values**:
   Invalid (NaN) values are masked out.

8. **Log-transform Data**:
   The NH and C+ data are log-transformed.

9. **Create 3D Plot**:
   The script creates an interactive 3D scatter plot using `plotly`, with `log(NH)` on the x-axis, `log(C+)` on the y-axis, and distance on the z-axis.

10. **Save Plot**:
    The plot is saved as an HTML file for interactive viewing.
