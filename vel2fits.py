import numpy as np
from astropy.io import fits

# Read the ASCII data e.g. RT_vel_SILCC64_+z_CO10.dat)
data = np.loadtxt('RT_vel_FILENAME_DIRECTION_LINE.dat') #specify the input _vel_ file
ipix = 128 #specify the resolution
csize = 2 #size of cloud in pc
outputname = 'OUTPUT.fifs' #specify the name of the fits output
channels = 101 #number of channels (don't forget the velocity = 0 in your channel count!)

# Get the unique velocity values (assuming these are in the first column)
velocities = np.unique(data[:, 0])

# Prepare the 3D array (velocity, x, y) for the FITS file
# Create an empty 3D array with shape (200, 128, 128)
# 200 is the number of velocity channels (frames)
# 128x128 is the size of each image
image_shape = (channels, ipix, ipix)
image_data = np.zeros(image_shape)

# Fill the image_data array with the appropriate values from the ASCII file
for row in data:
    v, x, y, tr = row
    v_index = np.where(velocities == v)[0][0]  # Find the index of the velocity value
    x_index = int(x*float(ipix)/csize)  # Assuming x is in [0, 127] for the 128x128 grid
    y_index = int(y*float(ipix)/csize)  # Assuming y is in [0, 127] for the 128x128 grid
    image_data[v_index, x_index, y_index] = tr  # Fill in the temperature (Tr) value

# Create the FITS header
header = fits.Header()
header['SIMPLE'] = True
header['BITPIX'] = -32
header['NAXIS'] = 3
header['NAXIS1'] = ipix
header['NAXIS2'] = ipix
header['NAXIS3'] = channels
header['EXTEND'] = True
header['COMMENT'] = "FITS file created for velocity, x, y, Tr data"

# Create the Primary HDU (header data unit) from the 3D array
hdu = fits.PrimaryHDU(data=image_data, header=header)

# Write the FITS file
hdu.writeto(outputname, overwrite=True)

print("FITS file has been written")

