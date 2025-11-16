import rasterio
import numpy
import matplotlib.pyplot as plt
import normxcorr2 as norm
import copy
import csv


def assign_value(array, r, c, value):
    # assins a value to array[r][c] if the indices are valid
    if 0 <= r < len(array) and 0 <= c < len(array[0]):
        array[r][c] = value


def draw_rect(array, row, col, width, height, border_size = 3):
    # draws a rectangle starting from left top corner 
    
    # draws the top line
    for c in range(col, col+width):
        for r in range(row, row+border_size):
            assign_value(array, r, c, 255)
            
    # draws the bottom line
    for c in range(col, col+width):
        for r in range(row+height-border_size, row+height):
            assign_value(array, r, c, 255)
            
    # draws left vertical line
    for c in range(col, col+border_size):
        for r in range(row+border_size, row+height-border_size):
            assign_value(array, r, c, 255)
            
    # draws left vertical line
    for c in range(col+width-border_size, col+width):
        for r in range(row+border_size, row+height-border_size):
            assign_value(array, r, c, 255)
            
            
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
          
                
def find_all_symbol_coor(corr_array):
    # finds all coordinates where the correlation value is higher than the set valid_corr_value
    coordinate_list = []
    for r, row in enumerate(corr_array):
        for c, value in enumerate(row):
            if value > valid_corr_value:
                coordinate_list.append([r, c])
    return coordinate_list
                
def filter_duplicates(coordinate_list):
    # filters the coordinate list to remove duplicates that are too close to each other
    # and returns new filtered list
    filtered_list = []
    
    for pair in coordinate_list:
        is_duplicate = False
        for filt_pair in filtered_list:
            if abs(pair[0]-filt_pair[0]) < sym_height/2 and abs(pair[1]-filt_pair[1]) < sym_width/2:
                is_duplicate = True
                break
        if not is_duplicate:
            filtered_list.append(pair)
    
    return filtered_list


def find_symbol_coor(coor_array):
    # returns a filtered list of found symbols
    coordinate_list = find_all_symbol_coor(coor_array)
    filtered_coordinate_list = filter_duplicates(coordinate_list)
    return filtered_coordinate_list


def coordinates_to_csv(coordinate_list, filename=None):
    # takes the coordinate list a transfers it to a csv file
    
    # checks if the list is None
    if coordinate_list is None:
        raise ValueError("coordinate_list is None")

    # checks if the filename was not set and assigns a default name if so
    if filename is None:    
        filename = 'found_coordinates.csv'
        
    # creates a header and writes the individual pairs by rows
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['x', 'y'])
        for pair in coordinate_list:
            x = int(pair[1])
            y = int(pair[0])
            
            writer.writerow([x, y])

    return filename
              
                
band_choice = 0 # 0-red, 1-green, 2-blue
valid_corr_value = 0.62 # sets the value for deciding if the kernel of defines symbol is valid

with rasterio.open('MMC_sk2.jpg', 'r') as ds:
    arr = ds.read()  # stores values like arr[band][row][column]
    rgb_c = 255 # RGB values from 0 to 255
    
    # normalize the values and store them to chosen band
    band = arr[band_choice]/rgb_c 

    
# symbol size
sym_width = 31 
sym_height = 78 

# reference symbol borders [symbol][left_hor, right_hor, up_vert, down_vert]
ref_sym_borders = [[1217, 1217 + sym_width, 4855, 4855 + sym_height], 
                   [1349, 1349 + sym_width, 5257, 5257 + sym_height],
                   [1438, 1438 + sym_width, 5328, 5328 + sym_height],
                   [5005, 5005 + sym_width, 5349, 5349 + sym_height],
                   [3773, 3773 + sym_width, 3906, 3906 + sym_height]
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
plt.imsave(f"ref_symbol_average_{band_choice}.jpg", ref_symbol_int8, cmap="gray")

# using normalized corralation function 
corr_image = norm.normxcorr2(symbol_array, band, "same")
#plt.imsave(f"corr_image.jpg", corr_image, cmap="gray")

# creating a new image with only valid rectangles
valid_corr_rects = numpy.zeros([len(band), len(band[0])])
valid_corr_rects = draw_rects(valid_corr_rects, corr_image)

# creating a csv of found symbol coordianates without repeats
found_symbol_coordinates = find_symbol_coor(corr_image)
output_file = f'found_coordinates_{band_choice}.csv'
coordinates_to_csv(found_symbol_coordinates, filename=output_file)
            
#plt.imsave(f"valid_corr_rects.jpg", valid_corr_rects)

band_with_rects = draw_rects(band, corr_image)
band_255 = band_with_rects*255
band_int8 = band_255.astype(numpy.uint8)

# draws the rectangles on the original raster            
#plt.imsave(f"found_symbols_on_image_{band_choice}.jpg", band_int8, cmap="pink")


