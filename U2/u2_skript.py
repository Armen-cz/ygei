import rasterio
import numpy
import matplotlib.pyplot as plt
import normxcorr2 as norm
import copy

def draw_rect(array, row, col, width, height, border_size = 3):
    # draws a rectangle starting from left top corner 
    
    # draws the top line
    for c in range(col, col+width):
        for r in range(row, row+border_size):
            array[r][c] = 255
            
    # draws the bottom line
    for c in range(col, col+width):
        for r in range(row+height-border_size, row+height):
            array[r][c] = 255
            
    # draws left vertical line
    for c in range(col, col+border_size):
        for r in range(row+border_size, row+height-border_size):
            array[r][c] = 255
            
    # draws left vertical line
    for c in range(col+width-border_size, col+width):
        for r in range(row+border_size, row+height-border_size):
            array[r][c] = 255
            
            
def draw_rects(image, corr_array, border_size = 3):
    # input is an image to draw rects on, correlation array and border size of rectangle
    # output is a new drawn raster
    new_image = copy.deepcopy(image)
    
    # draws a rectangle if image[r][c] > valid correlation value
    for r, row in enumerate(corr_array):
        for c, value in enumerate(row):
            if value > valid_corr_value:
                draw_rect(new_image, r-int(sym_height/2), c-int(sym_width/2), sym_width, sym_height, border_size)
                
    return new_image
                
                
band_choice = 0 # 0-red, 1-green, 2-blue
valid_corr_value = 0.6 # sets the value for deciding if the kernel of defines symbol is valid

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
#plt.imsave(f"ref_symbol_average_{band_choice}.jpg", ref_symbol_int8)

# using normalized corralation function 
corr_image = norm.normxcorr2(symbol_array, band, "same")
#plt.imsave(f"corr_image.jpg", corr_image)

# creating a new image with only valid rectangles
valid_corr_rects = numpy.zeros([len(band), len(band[0])])
valid_corr_rects = draw_rects(valid_corr_rects, corr_image)
            
plt.imsave(f"valid_corr_rects.jpg", valid_corr_rects)

band_with_rects = draw_rects(band, corr_image)
band_255 = band_with_rects*255
band_int8 = band_255.astype(numpy.uint8)

            
plt.imsave(f"valid_corr_points_on_image.jpg", band_int8)


