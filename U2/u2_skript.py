import rasterio
import numpy
import matplotlib.pyplot as plt

band_choice = 2 # 0-red, 1-green, 2-blue

with rasterio.open('MMC_sk2.jpg', 'r') as ds:
    arr = ds.read()  # stores values like arr[band][row][column]
    rgb_c = 255 # RGB values from 0 to 255
    
    # normalize the values and store them to chosen band
    band = arr[band_choice]/rgb_c 

    
# symbol size
sym_width = 39
sym_height = 82

# reference symbol borders [symbol][left_hor, right_hor, up_vert, down_vert]
ref_sym_borders = [[1212, 1212 + sym_width, 4854, 4854 + sym_height], 
                   [1344, 1344 + sym_width, 5256, 5256 + sym_height],
                   [1433, 1433 + sym_width, 5327, 5327 + sym_height],
                   [5000, 5000 + sym_width, 5348, 5348 + sym_height],
                   [3768, 3768 + sym_width, 3905, 3905 + sym_height]
                   ]

# create a new averaged symbol
symbol_array = numpy.zeros([sym_height, sym_width])
for border in ref_sym_borders:
    r = 0
    for row in range(border[2]-1, border[3]-1):
        c = 0
        for col in range(border[0]-1, border[1]-1):
            symbol_array[r][c] += band[row][col]
            c += 1
        r += 1

symbol_array = symbol_array/len(ref_sym_borders) # divides -> average

# part for saving the reference image
symbol_array_255 = symbol_array*255
ref_symbol_int8 = symbol_array_255.astype(numpy.uint8)
plt.imsave(f"ref_symbol_average_{band_choice}.jpg", ref_symbol_int8)