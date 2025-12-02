from graph_algorithms import *
from plot_algorithms import *
from file_to_graph_algorithms import *
import shapefile as shp  # Requires the pyshp package
import matplotlib.pyplot as plt
import random

"""Plots shortest paths from between two cities"""

plt.figure(figsize=(15,5))

plt.axis('equal')

file = 'silnice/s1_euklid.txt' # used for getting coordinates of road nodes
D, G, V, E = file_to_graph(file)

obce_file = "obce_souradnice.txt" # file with cities and its coordinates
obce, obce_jmena = obce_file_to_dict_list(obce_file) # transfers to dictionary (city name: (x, y)) and list of city names

#city1_name = "Svojšín" # start city 
#city2_name = "Kanina"  # end city

#city1_name = "Sokolov"
#city2_name = "Přelouč"

city1_name = "Úštěk"
city2_name = "Uhlířské Janovice"

start = city_to_node(city1_name, D, obce) # transfers city name to node
end = city_to_node(city2_name, D, obce)   # transfers city name to node

# plots additional polygons and roads
sf = shp.Reader("silnice/silnice_final_final.shp")
plot_shp_polygons(sf, 0.3, color="gray")

sf = shp.Reader("okresy/okresy.shp")
plot_shp_polygons(sf, 0.6)

sf = shp.Reader("kraje/kraje.shp")
plot_shp_polygons(sf, 1.3)

### Shortest path using euklid distance weight ###

file = 'silnice/s1_euklid.txt' # roads with euklid distances
D, G, V, E = file_to_graph(file)  # file to variables
pred, dist = shortest_path(G, start, end)          # finds shortest path
print(f"{round(dist/1000,3)} km: euklidovská vzdálenost")   # distance in kilometers

p = reconstPath(pred, start, end) # reconstructs path

C = my_dict2 = {y: x for x, y in D.items()}  # plots the path
plt.plot([C[x][0] for x in p], [C[y][1] for y in p], "r-", linewidth=1.5)

### Shortest path using euklid/max_speed weight ###

file = 'silnice/s2_rychlost.txt' # roads with euklid/max_speed distances
D, G, V, E = file_to_graph(file)  # file to variables
pred, dist = shortest_path(G, start, end)           # finds shortest path
print(f"{round(dist/1000,3)} hodin: rychlost + vzdálenost")  # distance in hours

p = reconstPath(pred, start, end) # reconstructs path

C = my_dict2 = {y: x for x, y in D.items()}  # plots the path
plt.plot([C[x][0] for x in p], [C[y][1] for y in p], "g-", linewidth=1.5)

### Shortest path using euklid/max_speed*curvature(klikatost) weight ###

file = 'silnice/s3_rychlost_klikatost.txt' # roads with euklid/max_speed*curvature distances
D, G, V, E = file_to_graph(file)  # file to variables
pred, dist = shortest_path(G, start, end)                       # finds shortest path
print(f"{round(dist/1000,3)} hodin: rychlost + klikatost + vzdálenost")  # distance in hours

p = reconstPath(pred, start, end) # reconstructs path

C = my_dict2 = {y: x for x, y in D.items()}  # plots the path
plt.plot([C[x][0] for x in p], [C[y][1] for y in p], "b-", linewidth=1.5)

plt.annotate(city1_name, obce[city1_name])   # plots the cities names
plt.annotate(city2_name, obce[city2_name])

plt.show() # shows the plot


"""Calculates minimal skeleton"""

# uses edges with weidht of distance/speed*curvature(klikatost)
file = 'silnice/s1_euklid.txt'
D, G, V, E = file_to_graph(file)

C = my_dict2 = {y: x for x, y in D.items()}

# plots additional polygons and roads
plt.figure(figsize=(15,5))
plt.axis('equal')

sf = shp.Reader("kraje/kraje.shp")
plot_shp_polygons(sf, 1.5)
    
sf = shp.Reader("silnice/silnice_final_final.shp")
plot_shp_polygons(sf, 0.3)

# sets what algorithm to use (can be used both)
boruvka = True
prime = False

if boruvka:# uses boruvka_kruskal algorithm for minimal skeleton
    skeleton_value, skeleton_tree = boruvka_kruskal(V, E) 
    plot_skeleton_tree(C, skeleton_tree)

if prime: # uses prime/jarnik algorithm for minimal skeleton
    skeleten_value_prime, skeleton_tree_prime = jarnik_prime(V, G)
    plot_skeleton_tree(C, skeleton_tree_prime, color="g")

plt.show() # plots minimal skeleton


"""calculates path from every point to every other point"""

# from every point to every point
file = 'silnice/s3_rychlost_klikatost.txt'
D, G, V, E = file_to_graph(file)
C = my_dict2 = {y: x for x, y in D.items()}
pred, _ = shortest_path(G)

# plot additional polygons
plt.figure(figsize=(15,5))
plt.axis('equal')

sf = shp.Reader("silnice/silnice_final_final.shp")
plot_shp_polygons(sf, 0.3)
    
sf = shp.Reader("okresy/okresy.shp")
plot_shp_polygons(sf, 0.6)
    
sf = shp.Reader("kraje/kraje.shp")
plot_shp_polygons(sf, 1.5)

# create random start and end points
for i in range(10): # creates 10 random paths
    city1_name = random.choice(obce_jmena) # chooses a random city name
    city2_name = random.choice(obce_jmena) # chooses a random city name
    
    start = city_to_node(city1_name, D, obce) # transfers city name to node
    end = city_to_node(city2_name, D, obce)   # transfers city name to node
    p = reconstPath(pred[start], start, end)  # reconstracts the path
    
    plt.plot([C[x][0] for x in p], [C[y][1] for y in p], "r-", linewidth=1.5)  # plots the path
    plt.annotate(city1_name, obce[city1_name])  # plots city names
    plt.annotate(city2_name, obce[city2_name])

plt.show()



