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
        plt.annotate(str(node), (coor[0] + text_dx, coor[1] + text_dy))
        
    plt.plot([C[i][0] for i in range(1,len(C)+1)], [C[j][1] for j in range(1, len(C)+1)], "ko")
    for node, coor in C.items():

        # checks if the predecessor exists
        if pred[node] != -1:
            node1_x, node1_y = coor[0], coor[1]
            node2_x, node2_y = C[pred[node]][0], C[pred[node]][1]
            dx, dy = node1_x - node2_x, node1_y - node2_y
            
            # draws an arrow between a node and its predecessor
            plt.arrow(node2_x, node2_y, dx, dy, width=0.5, length_includes_head = True, color="r", head_length = 7, head_width = 4)
    
    
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


# apply BFS
p_bfs = BFS(G, 1)
print(p_bfs)

pred_path = reconstPath(p_bfs, 1, 9)
print(pred_path)

# apply DFS with recursion
p = DFS(G, 1)
print(p)

pred_path = reconstPath(p, 1, 9)
print(pred_path)

# apply DFS with stack
p_dfs = DFSStack(G, 1)
print(p)

pred_path = reconstPath(p_dfs, 1, 9)
print(pred_path)

plt.figure(figsize=(15,5))
plt.subplot(121)
plot_graph(C, p_bfs)
plt.subplot(122)
plot_graph(C, p_dfs)
plt.suptitle("Left: BFS tree    Right: DFS tree")
plt.show()


