import rasterio
import numpy
import matplotlib.pyplot as plt
import normxcorr2 as norm
import copy
                
band_choice = 0 # 0-red, 1-green, 2-blue
valid_corr_value = 0.6 # sets the value for deciding if the kernel of defines symbol is valid

with rasterio.open('TM25_sk2.jpg', 'r') as ds:
    arr = ds.read()  # stores values like arr[band][row][column]
    rgb_c = 255 # RGB values from 0 to 255
    
    # normalize the values and store them to chosen band
    band = arr[band_choice]/rgb_c 

plt.imshow(band)
# used for later
#band_255 = band*255
#band_int8 = band_255.astype(numpy.uint8)            
#plt.imsave(f"valid_corr_points_on_image.jpg", band_int8) 


