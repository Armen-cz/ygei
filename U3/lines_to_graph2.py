from collections import *
from queue import *
from numpy import *
from skript_ze_cviceni import *
import shapefile as shp  # Requires the pyshp package
import matplotlib.pyplot as plt
import time

def create_obce_dict(filename):
    obce = dict()
    with open(filename, "r", encoding="utf-8") as f:
        
        for line in f:
            parts = line.strip().split("\t")     # split by TAB
            name = parts[0]
            x = float(parts[1])
            y = float(parts[2])

            obce[name] = (x, y)
            
    return obce

def city_to_node(city_name, D, obce):
    try:
        city = obce[city_name]
    except KeyError:
        raise Exception(f"Město {city_name} není v databázi")
        
    shortest = 1000000
    shortest_key = 0
    
    for value, key in D.items():
        distance = math.sqrt(math.pow(city[0]-value[0],2) + math.pow(city[1]-value[1],2))
        if distance < shortest:
            shortest = distance
            shortest_key = key
    return shortest_key

def loadEdges(file_name):
    #Convert list of lines to the graph
    PS = []
    PE = []
    W = []
    with open(file_name) as f:
        for line in f:
            #Split
            x1, y1, x2, y2, w = line.split()
            
            #Add start, end points and weights to the list
            PS.append((float(x1), float(y1)))
            PE.append((float(x2), float(y2)))
            W.append(float(w))
    return PS, PE, W

def pointsToIDs(P):
    #Create a map: key = coordinates, value = id
    D = {}
    for i in range(len(P)):
        D[(P[i][0], P[i][1])] = i
        
    return D

def edgesToGraph(D, PS, PE, W):
    #Convert edges to undirected graph
    G = defaultdict(dict)

    for i in range(len(PS)):
        G[D[PS[i]]][D[PE[i]]] = W[i]
        G[D[PE[i]]][D[PS[i]]] = W[i]

    #print(G)
    return G


def graph_to_nodes(G):
    V = []
    for node, _ in G.items():
        V.append(node)
    return V


def graph_to_edges(G):
    E = set()
    n = len(G)
    for u in range(n):
        for v, wuv in G[u].items():
            E.add((u, v, wuv))
    return E


def file_to_graph(file):
    PS, PE, W = loadEdges(file)

    #Merge lists and remove unique points
    PSE = PS + PE
    PSE=unique(PSE,axis=0).tolist()
    PSE.insert(0, [1000000, 1000000])

    #Edges to graph
    D = pointsToIDs(PSE)
    G = edgesToGraph(D, PS, PE, W)

    V = graph_to_nodes(G)
    E = graph_to_edges(G)
    
    return D, G, V, E

#print(C)
plt.figure(figsize=(15,5))

plt.axis('equal')
#plot_graph(C, pred)

file = 'silnice/s1_bez_klipu.txt'
D, G, V, E = file_to_graph(file)

obce = create_obce_dict("obce_souradnice.txt")

city1_name = "Sokolov"
city2_name = "Přelouč"

start = city_to_node(city1_name, D, obce)
end = city_to_node(city2_name, D, obce)


sf = shp.Reader("silnice/silnice_final_final.shp")

for shape in sf.shapeRecords():
    x = [i[0] for i in shape.shape.points[:]]
    y = [i[1] for i in shape.shape.points[:]]
    plt.plot(x,y, "-", color="gray", linewidth=0.3)
    
sf = shp.Reader("okresy/okresy.shp")

for shape in sf.shapeRecords():
    x = [i[0] for i in shape.shape.points[:]]
    y = [i[1] for i in shape.shape.points[:]]
    plt.plot(x,y, "k-", linewidth=0.6)
    
sf = shp.Reader("kraje/kraje.shp")

for shape in sf.shapeRecords():
    x = [i[0] for i in shape.shape.points[:]]
    y = [i[1] for i in shape.shape.points[:]]
    plt.plot(x,y, "k-", linewidth=1.5)
    
    #print([C[x][0] for x in p])

#pred, dist = shortest_path(G, start, end)
pred, dist = bellman_ford(G, start, end)
print(f"{dist/1000} km: euklidovská vzdálenost")

p = reconstPath(pred, start, end)

C = my_dict2 = {y: x for x, y in D.items()}

plt.annotate(city1_name, obce[city1_name])
plt.annotate(city2_name, obce[city2_name])

plt.plot([C[x][0] for x in p], [C[y][1] for y in p], "r-", linewidth=1.5)

file = 'silnice/s2_bez_klipu.txt'
D, G, V, E = file_to_graph(file)

#pred, dist = shortest_path(G, start, end)
pred, dist = bellman_ford(G, start, end)
print(f"{dist/1000} hodin: rychlost + vzdálenost")

p = reconstPath(pred, start, end)


C = my_dict2 = {y: x for x, y in D.items()}

plt.plot([C[x][0] for x in p], [C[y][1] for y in p], "g-", linewidth=1.5)

file = 'silnice/s3_bez_klipu.txt'
D, G, V, E = file_to_graph(file)

#pred, dist = shortest_path(G, start, end)
pred, dist = bellman_ford(G, start, end)
print(f"{dist/1000} hodin: rychlost + klikatost + vzdálenost")

p = reconstPath(pred, start, end)


C = my_dict2 = {y: x for x, y in D.items()}

plt.plot([C[x][0] for x in p], [C[y][1] for y in p], "b-", linewidth=1.5)

plt.show()

# vykreslení minimální kostry kruskal/boruvka

file = 'silnice/s3_bez_klipu.txt'
D, G, V, E = file_to_graph(file)

C = my_dict2 = {y: x for x, y in D.items()}

plt.figure(figsize=(15,5))

plt.axis('equal')

for shape in sf.shapeRecords():
    x = [i[0] for i in shape.shape.points[:]]
    y = [i[1] for i in shape.shape.points[:]]
    plt.plot(x,y, "k-", linewidth=1.5)
    
sf = shp.Reader("silnice/silnice_final_final.shp")

for shape in sf.shapeRecords():
    x = [i[0] for i in shape.shape.points[:]]
    y = [i[1] for i in shape.shape.points[:]]
    plt.plot(x,y, "k-", linewidth=0.3)

starter_time = time.time()
skeleton_value, skeleton_tree = boruvka_kruskal(V, E)
print(skeleton_value)
plot_skeleton_tree(C, skeleton_tree)

skeleten_value_prime, skeleton_tree_prime = jarnik_prime(V, G)
print(skeleten_value_prime)
plot_skeleton_tree(C, skeleton_tree_prime, color="g")

plt.show()

"""calculates path from every point to every other point"""
# from every point to every point
pred, dist = shortest_path(G)

start = 100
end = 2982
p = reconstPath(pred[start], start, end)

start = 243
end = 7433
p2 = reconstPath(pred[start], start, end)

start = 2473
end = 11482
p3 = reconstPath(pred[start], start, end)

start = 5233
end = 9478
p4 = reconstPath(pred[start], start, end)
                
C = my_dict2 = {y: x for x, y in D.items()}
#print(C)
plt.figure(figsize=(15,5))

plt.axis('equal')
#plot_graph(C, pred)

sf = shp.Reader("silnice/silnice_final_final.shp")

for shape in sf.shapeRecords():
    x = [i[0] for i in shape.shape.points[:]]
    y = [i[1] for i in shape.shape.points[:]]
    plt.plot(x,y, "k-", linewidth=0.3)
    
sf = shp.Reader("okresy/okresy.shp")

for shape in sf.shapeRecords():
    x = [i[0] for i in shape.shape.points[:]]
    y = [i[1] for i in shape.shape.points[:]]
    plt.plot(x,y, "k-", linewidth=0.6)
    
    sf = shp.Reader("kraje/kraje.shp")

for shape in sf.shapeRecords():
    x = [i[0] for i in shape.shape.points[:]]
    y = [i[1] for i in shape.shape.points[:]]
    plt.plot(x,y, "k-", linewidth=1.5)
    
    #print([C[x][0] for x in p])
plt.plot([C[x][0] for x in p], [C[y][1] for y in p], "r-", linewidth=1.5)
plt.plot([C[x][0] for x in p2], [C[y][1] for y in p2], "r-", linewidth=1.5)
plt.plot([C[x][0] for x in p3], [C[y][1] for y in p3], "r-", linewidth=1.5)
plt.plot([C[x][0] for x in p4], [C[y][1] for y in p4], "r-", linewidth=1.5)

plt.show()
