import matplotlib.pyplot as plt
import shapefile as shp  # Requires the pyshp package

def plot_graph(C, pred):    # plots the nodes into graph, used for DFS and BFS tree
    # plots the node numbers
    text_dx = 3 
    text_dy = 2
    for node, coor in C.items():
        plt.annotate(str(node), (coor[0] + text_dx, coor[1] + text_dy))
    # plots the nodes
    plt.plot([C[i][0] for i in range(1,len(C))], [C[j][1] for j in range(1, len(C))], "ko")
    
    for node, coor in C.items(): # draws arrow with direction from predecessor
        # checks if the predecessor exists
        if pred[node] != -1:
            node1_x, node1_y = coor[0], coor[1]
            node2_x, node2_y = C[pred[node]][0], C[pred[node]][1]
            dx, dy = node1_x - node2_x, node1_y - node2_y
            
            # draws an arrow between a node and its predecessor
            plt.arrow(node2_x, node2_y, dx, dy, width=0.5, length_includes_head = True, color="r", head_length = 7, head_width = 4)
    

def plot_skeleton_tree(C, T, color="r"): # plots a skeleton tree
    # input: coordinates of all nodes, tree, color (voluntary)
    for edge in T:
        node1_x, node1_y = C[edge[0]][0], C[edge[0]][1]
        node2_x, node2_y = C[edge[1]][0], C[edge[1]][1]
        plt.plot((node1_x, node2_x), (node1_y, node2_y), color=color)
        
        
def plot_shp_polygons(sf, line_width, color="k"): # plots shapefiles 
    # plots each shape (polygon, line) separately
    for shapeRec in sf.shapeRecords():
        shp = shapeRec.shape # takes one shape
        # each "part" is a separate polygon ring
        parts = shp.parts.tolist() + [len(shp.points)]
        for i in range(len(parts) - 1):
            start_i, end_i = parts[i], parts[i+1]
            x = [p[0] for p in shp.points[start_i:end_i]]
            y = [p[1] for p in shp.points[start_i:end_i]]
            plt.plot(x, y, "-", linewidth=line_width, color=color)
            
    
