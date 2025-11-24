# Graph

# výborná literatura - průvodcem labyrintem algoritmů od českých tvůrců

import graph_def
import matplotlib.pyplot as plt

G = graph_def.G
C = graph_def.C


# DU 
# vytvořit funkci, která cestu vizualizuje přes matplotlib
# bodová množina jako uzly
def plot_graph(C, pred):
    # plots the nodes
    text_dx = 3
    text_dy = 2
    for node, coor in C.items():
        plt.annotate(str(node), (-coor[0] - text_dx, -coor[1] - text_dy))
        
    plt.plot([-C[i][0] for i in range(1,len(C))], [-C[j][1] for j in range(1, len(C))], "ko")
    for node, coor in C.items():

        # checks if the predecessor exists
        if pred[node] != -1:
            node1_x, node1_y = -coor[0], -coor[1]
            node2_x, node2_y = -C[pred[node]][0], -C[pred[node]][1]
            dx, dy = node1_x - node2_x, node1_y - node2_y
            
            # draws an arrow between a node and its predecessor
            plt.arrow(node2_x, node2_y, dx, dy, width=0.5, length_includes_head = True, color="r", head_length = 7, head_width = 4)
    
    
def plot_skeleton_tree(C, T):
    for edge in T:
        node1_x, node1_y = -C[edge[0]][0], -C[edge[0]][1]
        node2_x, node2_y = -C[edge[1]][0], -C[edge[1]][1]
        plt.plot((node1_x, node2_x), (node1_y, node2_y), color="r")


# BFS - in(G,u), out(BF tree (predecesers))
def BFS(G, u):
    # define variables
    node_count = len(G)
    states = ["new"] * (node_count+1)
    pred = [-1] * (node_count+1)
    
    # add starting node and change status
    queue = [u]
    states[u] = "open"
    
    # while queue not empty
    while queue:
        # takes first node
        u = queue.pop(0)
        
        # browse adjacent nodes
        for v in G[u]:
            
            # node is new
            if states[v] == "new":
                
                # change status
                states[v] = "open"
                
                # update predecessor
                pred[v] = u
                
                # add v to queue
                queue.append(v)
                
        # after all marked children, states to closed
        states[u] = "closed"
        
    return pred


# BFS - in(G,u), out(BF tree (predecesers))
def DFSStack(G, u):
    # define variables
    node_count = len(G)
    states = ["new"] * (node_count+1)
    pred = [-1] * (node_count+1)
    
    # add starting node and change status
    stack = [u]
    
    
    # while queue not empty
    while stack:
        # takes last node
        u = stack.pop() 
        
        # change status
        states[u] = "open"
        
        # browse adjacent nodes
        for v in reversed(G[u]):
            
            # node is new
            if states[v] == "new":
                
                # update predecessor
                pred[v] = u
                
                # add v to queue
                stack.append(v)
                
        # after all marked children, states to closed
        states[u] = "closed"
        
    return pred
            

def reconstPath(pred, u, v):
    path = [v]
    
    # path shortening
    while v != u and v != -1:
        # update predecessor
        v = pred[v]
        
        # add to the list
        path.append(v)
        
    return path


def DFS(G, u):
    # define variables
    node_count = len(G)
    states = ["new"] * (node_count+1)
    pred = [-1] * (node_count+1)
    
    # browse all nodes
    for u, value in G.items():
        if states[u] == "new":
            DFSR(G, u, pred, states)
            
    return pred


def DFSR(G, u, pred, states):
    states[u] = "open"
    
    for v in G[u]:
        if states[v] == "new":
            pred[v] = u
            DFSR(G, v, pred, states)
            
    states[u] = "closed"


# # apply BFS
# p_bfs = BFS(G, 1)
# print(p_bfs)

# pred_path = reconstPath(p_bfs, 1, 9)
# print(pred_path)

# # apply DFS with recursion
# p = DFS(G, 1)
# print(p)

# pred_path = reconstPath(p, 1, 9)
# print(pred_path)

# # apply DFS with stack
# p_dfs = DFSStack(G, 1)
# print(p)

# pred_path = reconstPath(p_dfs, 1, 9)
# print(pred_path)

# plt.figure(figsize=(15,5))
# plt.subplot(121)
# plot_graph(C, p_bfs)
# plt.subplot(122)
# plot_graph(C, p_dfs)
# plt.suptitle("Left: BFS tree    Right: DFS tree")
# plt.show()


### 2. cvičení ###

# Priority queue: <prior: value>
# jarnik - modifikace dijkstra - dv = du
# alespoň 100 měst
# Dalnice: 130 km/h
# 1T-3T: 90 km/h
# sídla: 50 km/h
# t = s/v
# faktor klikatosti f = d(u,v)/||u-v||
# t' = f*(s/v)
# porovnat s mapy.cz, applemapy, google maps
# hranu exportovat jako u, v, w=||u-v|| (v arcgisu to jde nějak exporotovat)
# souřadnice na 3 desetinná místa

import graph_def3
import math
import queue

G = graph_def3.G
C = graph_def3.C

def dijkstra(G, start):
    # Dijkstra algorithmus
    
    n = len(G)
    pred = [-1]*(n+1) # Predecessors
    dist = [math.inf]*(n+1) # distances
    
    # create priority queue
    PQ = queue.PriorityQueue()
    
    # add start vertex
    PQ.put((0, start))
    
    dist[start] = 0
    
    # repeat until queue is empty
    while not PQ.empty():
        
        # get point with lowest du
        du, u = PQ.get()
        
        # browse all adjacent nodes
        for v, wuv in G[u].items():
            # relax
            if dist[v] > dist[u] + wuv:
                # update dv
                dist[v] = dist[u] + wuv
                
                # store predecessor
                pred[v] = u
                
                # add v to priority queue
                PQ.put((dist[v], v))
                
    return pred

def relaxace(G, start):
    # algoritm for finding shortest paths in graph with negative edges
    n = len(G)
    pred = [-1]*(n+1) # Predecessors
    dist = [math.inf]*(n+1) # distances
    states = ["N"]*(n-1) # states set to Not found

    states[start] = "O" # sets start node as Open
    dist[start] = 0 # sets distance of start node 0

    # create queue
    Q = queue.Queue()
    
    # add start vertex
    Q.put((start))

        # repeat until queue is empty
    while not Q.empty():
        
        # get point the is in the queue for longest
        u = Q.get()

        # browse all adjacent nodes
        for v, wuv in G[u].items():
            if dist[v] > dist[u] + wuv:
                dist[v] = dist[u] + wuv
                states[v] = "O"
                Q.put((v))
                pred[v] = u
        states[v] = "C" # state of original node is Closed
    return pred


def bellman_ford(G, start):
    # implementation of bellman_ford algorithm - for graph with negative edges
    n = len(G)
    pred = [-1]*(n+1) # Predecessors
    dist = [math.inf]*(n+1) # distances
    
    dist[start] = 0

    # relax edges n times
    for i in range(n+1):
        improved = False
        for u in range(n):
            for v, wuv in G[u].items():
                if dist[u] != math.inf and dist[v] > dist[u] + wuv:
                    dist[v] = dist[u] + wuv
                    pred[v] = u
                    improved = True
        # if no path is better, algorithm can stop sooner
        if not improved:
            break

        # if relaxation happends n-times, it has a negative cycle
        #if i == n:
        #    raise Exception("Graph contains negative cycle")

    return pred


def shortest_path(G, start = None):
    print("Searching graph")
    # if start is defined, calculates shortest paths from one point 
    # to all other points
    n = len(G)
    if start is not None:
        negative_edges = False
        for u in range(n):
            for v, wuv in G[u].items():
                if wuv < 0:
                    negative_edges = True
        if negative_edges is True:
            print("Found negative edges")
            print("Using Bellman-Ford algorithm")
            pred = bellman_ford(G, start)
        else:
            print("Found only positive edges")
            print("Using Dijkstra algorithm")
            pred = dijkstra(G, start)
        return pred
    else:
        # if no point is defined as start point
        # calculates shortest paths from every point to every point
        # usis dijkstra from each point and saves it to pred[point][pred]
        print("Calculates paths from every point to every point")
        pred = [None]*(n+1)
        for u in range(1, n):
            if u % 100 == 0:
                print(f"Calculated path for {u} points")
            pred[u] = dijkstra(G, u)
        return pred
            
# makes new set
def make_set(u, pred, rank):
    pred[u] = u
    rank[u] = 0


# find root
def find(pred, u, path_compression = None):
    
    if path_compression == "half": # half compression
        while pred[u] != u:
            pred[u] = pred[pred[u]] # moves to grandparent
            u = pred[u]             # predecessor is grandparent
        return u
    elif path_compression == "full": # full compression
        while pred[u] != u:
            u = pred[u]
        root = u
        while u != root:
            u_pred = pred[u]  #Store predecessor
            pred[u] = root    #Change predecessor to root
            u = u_pred        #Go to parent
        return u
    else: # does it without path compression
        while pred[u] != u:
            u = pred[u]
        return u



# weighted union
def union(pred, u, v, rank, path_compression=None):
    root_u = find(pred, u, path_compression=path_compression)
    root_v = find(pred, v, path_compression=path_compression)
    if root_u != root_v:
        if rank[root_u] > rank[root_v]:
            pred[root_v] = root_u
        elif rank[root_v] > rank[root_u]:
            pred[root_u] = root_v 
        else:
            pred[root_u] = root_v
            rank[root_v] = root_v + 1


def boruvka_kruskal(V, E, path_compression=None):
    if path_compression not in [None, "half", "full"]:
        print("Invalid path compression argument, path compression is not used")
        
    T=[] #Empty tree
    wt = 0 #Sum of weights of T
    pred = [-1] * (len(V) + 1) #List of roots
    rank = [math.inf] * (len(E) + 1) #Rank of the node
    for v in V: #Make set
        make_set(v, pred, rank) #Initilize p and r
    ES = sorted(E, key=lambda it:it[2]) #Sort edges by w
    for e in ES: #Process all edges
        u, v, w = e #Take an edge
        if (find(pred, u, path_compression=path_compression) != find(pred, v, path_compression=path_compression)): #roots u, v in different trees?
            union(pred, u, v, rank) #Create union
            T.append([u, v, w]) #Add edge to tree
            wt = wt + w #Compute weight of T
    return wt, T



# call dijkstra
pred = dijkstra(G, 1)

# path 
path = reconstPath(pred, 1, 9)

# print results
print(path)
    
# v graph_def2 jsou věci na minimání kostru

