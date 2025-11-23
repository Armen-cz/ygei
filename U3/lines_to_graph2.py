from collections import *
from queue import *
from numpy import *
from skript_ze_cviceni import *
import shapefile as shp  # Requires the pyshp package
import matplotlib.pyplot as plt
  
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


#Load edges
file = 'cesty_2015_euklid.csv'
#file = 'graph_disjkstra.txt'

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

pred = shortest_path(G, 48)

p = reconstPath(pred, 48, 12911)


C = my_dict2 = {y: x for x, y in D.items()}
#print(C)
plt.figure(figsize=(15,5))

plt.axis('equal')
#plot_graph(C, pred)

sf = shp.Reader("silnice/silnice_2015.shp")

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
plt.plot([-C[x][0] for x in p], [-C[y][1] for y in p], "r-", linewidth=1.5)

plt.show()

# vykreslení minimální kostry kruskal/boruvka

plt.figure(figsize=(15,5))

plt.axis('equal')

for shape in sf.shapeRecords():
    x = [i[0] for i in shape.shape.points[:]]
    y = [i[1] for i in shape.shape.points[:]]
    plt.plot(x,y, "k-", linewidth=1.5)

skeleton_value, skeleton_tree = boruvka_kruskal(V, E)

plot_skeleton_tree(C, skeleton_tree)

plt.show()
""" calculates path from every point to every other point
# from every point to every point
pred = shortest_path(G)

start = 100
end = 2982
p = reconstPath(pred[start], start, end)

start = 243
end = 7433
p2 = reconstPath(pred[start], start, end)

start = 2473
end = 12953
p3 = reconstPath(pred[start], start, end)

start = 5233
end = 9478
p4 = reconstPath(pred[start], start, end)
                
C = my_dict2 = {y: x for x, y in D.items()}
#print(C)
plt.figure(figsize=(15,5))

plt.axis('equal')
#plot_graph(C, pred)

sf = shp.Reader("silnice/silnice_2015.shp")

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
plt.plot([-C[x][0] for x in p], [-C[y][1] for y in p], "r-", linewidth=1.5)
plt.plot([-C[x][0] for x in p2], [-C[y][1] for y in p2], "r-", linewidth=1.5)
plt.plot([-C[x][0] for x in p3], [-C[y][1] for y in p3], "r-", linewidth=1.5)
plt.plot([-C[x][0] for x in p4], [-C[y][1] for y in p4], "r-", linewidth=1.5)

plt.show()"""
