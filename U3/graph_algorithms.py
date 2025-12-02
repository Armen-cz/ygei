# výborná literatura - průvodcem labyrintem algoritmů od českých tvůrců

import random
import math
import queue


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
    # reconstructs path from list of predecessors from start to end node
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
    for u, _ in G.items():
        if states[u] == "new":
            DFSR(G, u, pred, states) # uses recursive DFS
            
    return pred


def DFSR(G, u, pred, states):
    # recusive DFS, used in full DFS algorithm
    states[u] = "open"  # updates status to OPEN
    
    for v in G[u]:
        if states[v] == "new":
            pred[v] = u
            DFSR(G, v, pred, states) # calls itself
    
    states[u] = "closed"  # closes starting node


def dijkstra(G, start, end):
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
                
    return pred, dist[end]


def bellman_ford(G, start, end):
    # implementation of bellman_ford algorithm - for graph with negative edges
    n = len(G)
    pred = [-1]*(n+1) # predecessors
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
        if i == n:
            print("Found negative cycle")
            return None
        
    return pred, dist[end]


def shortest_path(G, start = None, end = 1):
    print("Searching graph")
    # if start is defined, calculates shortest paths from one point 
    # to all other points
    n = len(G)
    if start is not None: # if start node is defined
        negative_edges = False
        for u in range(n): # finds negative edges
            for v, wuv in G[u].items():
                if wuv < 0:
                    negative_edges = True
        if negative_edges is True: # if negative edges were found, uses bellman-ford
            print("Found negative edges")
            print("Using Bellman-Ford algorithm")
            pred, dist = bellman_ford(G, start, end)
        else: # if graph does not have negative edges, uses dijkstra
            print("Found only positive edges")
            print("Using Dijkstra algorithm")
            pred, dist = dijkstra(G, start, end)
        return pred, dist
    else:
        # if no point is defined as start point
        # calculates shortest paths from every point to every point
        # usis dijkstra from each point and saves it to pred[point][pred]
        print(f"Calculates paths from every point to every point [{n}: point total]")
        pred = [None]*(n+1)
        for u in range(1, n):
            if u % 100 == 0:
                print(f"Calculated path for {u} points")
            pred[u], dist = dijkstra(G, u, end)
        return pred, dist
            

def make_set(u, pred, rank):   # makes new set
    pred[u] = u
    rank[u] = 0


def find(pred, u, path_compression = None):  # find root
    
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
    # finds roots of start and end nodes
    root_u = find(pred, u, path_compression=path_compression)
    root_v = find(pred, v, path_compression=path_compression)
    if root_u != root_v: # if the roots are different, they are joined together
        if rank[root_u] > rank[root_v]:
            pred[root_v] = root_u
        elif rank[root_v] > rank[root_u]:
            pred[root_u] = root_v 
        else: # if they are the same size, v is joined to u
            pred[root_u] = root_v
            rank[root_v] = root_v + 1


def boruvka_kruskal(V, E, path_compression=None):
    if path_compression not in [None, "half", "full"]:
        print("Invalid path compression argument, path compression is not used")
        
    T=[] # empty tree
    wt = 0 # sum of weights of T
    pred = [-1] * (len(V) + 1) # list of predecessors
    rank = [math.inf] * (len(E) + 1) # ranks of the node
    for v in V: # makes set
        make_set(v, pred, rank) # initilizes p and r
    ES = sorted(E, key=lambda it:it[2]) # sorts edges by weight
    for e in ES: # processes all edges
        u, v, w = e # takes an edge
        if (find(pred, u, path_compression=path_compression) != find(pred, v, path_compression=path_compression)): # roots u, v in different trees
            union(pred, u, v, rank) # unites nodes
            T.append([u, v, w]) # adds edge to tree
            wt += w # adds weight to tree
    return wt, T


def jarnik_prime(V, G): 
    T = [] # empty tree 
    wt = 0 # sum of weights of T 
    u = random.choice(V) # chooses random node 
    tree_nodes = set([u]) # set of found nodes
    n = len(G) 
    while len(tree_nodes) != n: # while there are nodes to add
        if len(tree_nodes) % 1000 == 0: # for ensurance the code is working
            print(f"{len(tree_nodes)} bodů ve stromu") 
        min_wuv = math.inf 
        min_v = math.inf 
        for node in tree_nodes: 
            for v, wuv in G[node].items(): 
                if v in tree_nodes: # if node is already added to the tree -> continues
                    continue 
                if wuv < min_wuv: # selects node and edge with the least weight 
                    min_wuv = wuv 
                    min_v = v 
                    u = node 
        
        if min_wuv == math.inf: # if no other node was found, this breaks the cycle
            break

        T.append([u, min_v, min_wuv])  # adds node (edge) to the tree
        tree_nodes.add(min_v)  
        wt += min_wuv # adds weight to tree

    return wt, T