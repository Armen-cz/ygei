from queue import *
from numpy import *
from collections import *
import math


def obce_file_to_dict_list(filename):
    # takes a filename and creates a cities dictionary with value = city, key = (X, Y) 
    # and a list of cities
    obce = dict()
    obce_names = list()
    with open(filename, "r", encoding="utf-8") as f: # opens and reads the file
        
        for line in f: # process each line
            parts = line.strip().split("\t")  # splits by blanc spaces
            name = parts[0]
            x = float(parts[1])
            y = float(parts[2])

            obce[name] = (x, y)
            obce_names.append(name)
            
    return obce, obce_names


def city_to_node(city_name, D, obce):
    # takes city name and transfers it to node number
    try: # tries if the city name exist
        city = obce[city_name]
    except KeyError:
        raise Exception(f"Město {city_name} není v databázi")
    
    shortest = 1000000
    shortest_key = 0
    
    for value, key in D.items(): # finds the closest node on the road to the city definition point
        distance = math.sqrt(math.pow(city[0]-value[0],2) + math.pow(city[1]-value[1],2)) # pythagoras
        if distance < shortest: # rewrites current shortest with the new shortest distance
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