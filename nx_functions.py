'''
Functions in this library:
    1. shortest_path
    2. find_paths
    3. simple_min_cut
    4. augmented_min_cut
'''

class BreakException(Exception):
    pass

def path_not_found():
    raise BreakException()

def shortest_path(g, u, v):
    
    '''
    It is a "BFS" method and return the first path with minimum distance from u to v 
    
    '''
    if u == v:
        return 'same node'  
    
    visited = [] #list of visited nodes
    queue = [[u]] #queue structure
    
    visited.append(u) #add u to visited list
    
    nodes = [] #list of nodes --> path
    
    #calculate all paths with n distance
    #when i have v in nodes return nodes
    while queue:
        s = queue.pop(0)

        for neighbour in g[s[-1]]:
            if neighbour not in visited:
                nodes = list(s)
                nodes.append(neighbour)
                queue.append(nodes)
                visited.append(neighbour)
                
                if neighbour == v:
                    return nodes
                
    return path_not_found()

def find_paths(g, u, v):
    
    '''
    1. Find the shortest path
    2. Remove a random edge of the s.p. from the graph
    3. Repeat 1-2 until there are no paths
    '''
    
    graph = g.copy()
    
    paths = []
    
    while True:
        try:
            path = shortest_path(graph, u, v)
        except BreakException:
            break

        paths.append(path)
        index = random.choice(range(1, len(path)))
        graph.remove_edge(path[index-1], path[index])    
                
    return paths  

def simple_min_cut(g, u, v, starter):
    
    '''
    1. Chose one path from u to v --> starter
    2. Remove all edges of the path
    3. Reapet 1 and 2 until there's no path
    4. return the number of path found
    '''

    #print(f'random_path: {starter}')
    index=1
    
    my_g = g.copy()
    for j in range(1, len(starter)):
        my_g.remove_edge(starter[j-1], starter[j])
        #print(f'{index} removed: {starter[j-1]}, {starter[j]}')

    while True:
        try:
            nodes_list = shortest_path(my_g, u, v)
            #print('----')
        except BreakException:
            #print(f'Another path not found! MinCut: {index}')
            break

        index+=1
        for j in range(1, len(nodes_list)):
            my_g.remove_edge(nodes_list[j-1], nodes_list[j])
            #print(f'{index} removed: {nodes_list[j-1]}, {nodes_list[j]}')

    return index

def augmented_min_cut(g, u, v):
    
    '''
    1) if u and v are the same node -> return 'same node'
    2) if deg of u or v is 1 or 0 -> return 1 or 0
    3) paths = call find_paths to calculate different paths
    4) if paths have just one path -> return 1
    5) use simple_min_cut for each path in paths
    6) return the max of the min cuts
    '''

    if u == v:
        return 'Same node'
    
    deg_u = g.degree(u)
    deg_v = g.degree(v)
    
    if deg_u == 1 or deg_v == 1:
        return 1
    
    if deg_u == 0 or deg_v == 0:
        return 0

    paths = find_paths(g, u, v)
    print('Starter Paths: ', paths)
    lenght = len(paths)
    if lenght == 1:
        return 1
    results = []
    
    for i in range(lenght):
        index = simple_min_cut(g, u, v, paths[i]) 
        results.append(index)
    
    print('All min-cut: ', results)
    return max(results)

