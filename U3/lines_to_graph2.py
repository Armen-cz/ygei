from collections import *
from queue import *
from numpy import *
from skript_ze_cviceni import *
  
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

pred = dijkstra(G, 48, 12911)

p = reconstPath(pred, 48, 12911)


C = my_dict2 = {y: x for x, y in D.items()}
#print(C)
plt.figure(figsize=(15,5))

plt.axis('equal')
#plot_graph(C, pred)



import shapefile as shp  # Requires the pyshp package
import matplotlib.pyplot as plt

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

