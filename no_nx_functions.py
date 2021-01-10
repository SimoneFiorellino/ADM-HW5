'''
Functions in this library:
    1. explore
    2. shortest_path
    3. get_reversed_path
    4. get_articles
    4. categories_distance
    5. transform
    6. pagerank
'''
import pandas as pd
from matplotlib.pyplot import figure
import matplotlib.pyplot as plt
from random import randint
from collections import Counter 
import collections
import statistics
import time
from itertools import chain 
import random
import operator
import copy

#################################### general ########################################

def explore(node,limit, out_links):
   
    ''' This function takes as an input one node and a click limit.
    Limit can be set to None if we don't want any limit.
    The function returns a dictionary where keys are distances and values are all the nodes 
    which are at distance d (key) from the input node'''
    
    # initializing variables
    distance_tree = {}
    explored = set()
    clicks = 1 
    
    # checking if node has at list one outlink
    try:
        out_nodes = set(out_links[node])
    except:
        print('No out-link from', node)
        return distance_tree
    
    '''
    1. We iterate until the new set of out-nodes to explore is non-empty 
    2. For each iteration we only consider unexplored nodes
    3. If the node is a key in out_links, we add the set of nodes associated to it to the new set of nodes
    4. We update all variables and increase distance by one'''

    while out_nodes != set(): 
        
        if limit != None:
            if clicks > limit+1:        
                return distance_tree
        
        distance_tree[clicks] = out_nodes
        new_out_nodes = set() 
        
        for node in out_nodes: 
            try:
                new_out_nodes = new_out_nodes.union(set(out_links[node])) 
            except:
                continue
        
        explored = explored.union(out_nodes) 
        out_nodes = new_out_nodes.difference(explored)
        clicks += 1 
    
    return distance_tree

def shortest_path (node1,node2, out_links):
    
    ''' This function takes in input two nodes and returns the length of the shortest path.
    If node 1 does not have any outlink or if no path is found between node 1 and node 2, None is returned. '''
    
    # exeptions
    if node1 == node2:
        return 0, []
    
    # set of outnodes of node 1
    try:
        out_nodes = set(out_links[node1]) 
    except:
        print('No out-link from', node1)
        return 'Inf', []
    
    # initializing explored nodes set and minimum distance
    explored = set() 
    clicks = 1
    parents =  {node : {1:[node1]} for node in out_nodes}
    
    ''' ITERATIONS: we iterate until we find node 2. In the case in which we have not found node 2 
    but we don't have any new node to explore, we assume the two nodes are disconnected.
    
    1. We check the new set of out-nodes to explore is non-empty 
    2. We iterate through unexplored nodes
    3. If the node is a key in out_links, we add the set of nodes associated to it to the 
    new set of nodes in order to proceed with exploration
    4. In order to retrieve the actual set of nodes met to get from node1 to node2,
    for each explored node we save its parent (or parents) at each iteration and retrace the path
    backwards from node2 to node 1.'''
   
    while node2 not in parents.keys(): 
        clicks += 1 
        
        if out_nodes != set(): 
            new_out_nodes = set() 

            for node in out_nodes:  
                try:
                    new_out_nodes = new_out_nodes.union(set(out_links[node]))
                    for new_node in out_links[node]:
                        if new_node in parents.keys():
                            if clicks in parents[new_node]:
                                parents[new_node][clicks].append(node)
                            else:
                                parents[new_node].update({clicks:[node]})
                        else:
                            parents[new_node] = {clicks:[node]}
                except:
                    continue

            explored = explored.union(out_nodes)
            out_nodes = new_out_nodes.difference(explored) 
    
        else: 
            print('No path found')
            return 'Inf', []
    
    path = [node2]
    while path[-1] != node1:
        path = get_reversed_path(parents,clicks,node2)
    
    return clicks, [node2]+path

def get_reversed_path(parents, dist, node):
    
    '''we iterate backwards from node2 to node1 through the parents dictionary until we find node1 
    at distance "clicks" from node 2. In the case of branched path we explore them randomly untile we 
    find the one we are interested in.'''
    
    path = []
    for i in reversed(range(1,dist+1)):
        parent_list = parents[node][i]
        if len(parent_list) > 1:
            index = random.randint(0,len(parent_list)-1)
            parent = parent_list[index]
        else:
            parent = parent_list[0]
        path.append(parent)
        node = parent
    
    return path

####################################### RQ2 ##########################################

def get_articles(node,n_clicks,out_links):
    '''This function works as the explore function defined above, but with a different output format.'''
    tree = explore(node,n_clicks,out_links)
    articles = []
    for d in tree.keys():
        articles+=tree[d]
    return articles

####################################### RQ3 ##########################################

def visit_all(C, category_nodes, in_links, out_links):
    
    ''' PRIMARY STEPS :
    1. We retrieve the set of nodes in C
    2. We find the center and check that it exists (meaning that at list one node has a non zero in-degree)
    3. We check that at most one node has outdegree 0'''
    
    p = copy.copy(category_nodes[C])
    max_indeg = 0
    out_deg_p = {}
    center = ''
    
    for node in p:
        if node in in_links.keys():    
            if len(in_links[node]) > max_indeg:
                center = node
                max_indeg = len(in_links[node])
        try:
            out_deg_p[node] = len(out_links[node])
        except:
            out_deg_p[node] = 0
    
    zero_outdeg = list(out_deg_p.values()).count(0)
    
    p = set_of_pages(p)
    print('Nodes to visit', p)
    
    if center == '' or zero_outdeg > 1:
        print('More than one node has outdegree zero')
        return 'Not possible'
    elif center in p:
        p.remove(center)
    
    print('Center is node', center, 'with indegree', max_indeg)
    
    ''' SECONDARY STEPS :
    Now we apply the "get_next" function until we either find all nodes or we are stuck for two consecutive
    iterations at the same node'''
    
    start = center
    exploration = [center]
    to_visit = set(copy.copy(p))
    
    while len(to_visit) >= 1:
        tree = explore(start, None, out_links)
        next_node, path, visited = get_next(start,to_visit,tree,out_deg_p)
        if next_node == start or next_node == None or path == []:
            return 'Not possible'
        exploration += path + [next_node]
        to_visit =  to_visit - visited
        start = next_node
    
    return exploration, p+[center]

def set_of_pages(l):
    l_index = randint(0, len(l) - 1)
    r_index = randint(l_index + 1, len(l))
    return l[l_index:r_index]

def get_next(start, p, tree, out_deg_p):
        
    ''' We want to sort all nodes left to visit by their distance from the starting node, in decreasing order.
    The nodes are sorted in this way to be more time efficient, since we avoid to look two times for the same 
    shortest path.
    1. We create a set of nodes that we want to find and generate the exploration tree of our start node.
    2. We iterate all the distances in the tree and if we find any of the nodes of interest we save them 
    in a dictionary.
    3. We get the list of nodes nodes to be visited sorted by furthest to the nearest. '''
    
    p_decreasing = sort_nodes(tree, p)
    if p_decreasing ==[]:
        print('Nodes', p, 'cannot be reached from node', start)
        return None, [], set()
    
    print('Nodes to visit', len(p_decreasing))
        
    ''' Now we want to get all shortest paths from start node to any other node left to visit. 
    1. We have a set of nodes to visit to be iterated in a specific order. We always pick the one which is 
    further away from the start node (element 0).
    2. We save the shortest path to this node in a dictionary, remove the node from the list of nodes to visit 
    and save the nodes that we have met along the path if they are in the set of nodes to visit.
    3. We remove all the other nodes of interest met along the way from the set of nodes to visit. 
    
    notice: we alredy know that these nodes will not have the greatest number of crossed nodes '''
    
    rel_path, crossed = get_paths(p_decreasing, p, start, out_deg_p)
    
    ''' Among all nodes that are left to visit we pick as "next" according to this features:
    1. We check the node that maximises the number of crossed nodes. However if the outdegree of this node is 
    zero we check if we have any othe node to visit, if yes, we discard this node.
    2. If two nodes cross the same number of nodes, we pick the one with higher outdegree.
    3. If two nodes cross the same number of nodes and have equal outdegree, we pick the closer one'''
    
    best_node = p_decreasing[0]
    max_crossed = len(crossed[best_node])
    
    for node in crossed.keys():
        if len(crossed[node]) > max_crossed:
            if out_deg_p[node] == 0:
                if p - set(rel_path[node]+[best_node]) == set():
                    best_node = node
                    max_crossed = len(crossed[node])
            else:
                best_node = node
                max_crossed = len(crossed[node])

        elif len(crossed[node]) == max_crossed:          
            if out_deg_p[node] > out_deg_p[best_node] : 
                    best_node = node
            elif len(crossed[node]) > len(crossed[best_node]):
                best_node = node

    if rel_path[best_node] == []:
        print('Cannot reach', best_node)
    else:
        print('Path to', best_node, ':', rel_path[best_node])
    
    return best_node, rel_path[best_node], set(rel_path[best_node]+[best_node])

def sort_nodes(tree, p):
    
    distances = {}
    to_find = set(copy.copy(p))
    for d in tree.keys():
        visited = set(tree[d]).intersection(to_find)
        distances[d] = visited                 
    
    p_decreasing = collections.OrderedDict(sorted(distances.items()))
    p_decreasing = list(chain.from_iterable(list(p_decreasing.values())))
        
    return p_decreasing

def get_paths(p_decreasing, p, start, out_links):
    
    rel_path = {node:[] for node in p}
    crossed = {node:[] for node in p}
    to_visit = copy.copy(p_decreasing)
    
    while len(to_visit) > 1:
        node = to_visit[0]
        clicks,path = shortest_path(start, node, out_links)
        rel_path[node] = path
        crossed_list = set(path).intersection(set(to_visit)) 
        crossed[node] = crossed_list 
        to_visit.remove(node)
        for node_a in crossed_list:
            if node_a in to_visit:
                to_visit.remove(node_a)
                
    return rel_path, crossed

####################################### RQ5 ##########################################

def categories_distance(c1, category_nodes, out_links):
    
    c1_nodes = [int(i) for i in category_nodes[c1]] # int casting needed for comparison
    distances_dict = {k: [] for k in category_nodes.keys()}
    t = []
    
    '''
    1. We consider the exploration tree of each node in c1 with the exploration function defined above. 
    2. For each category c2 we save a list of distances between any node in c2 (node2) and node1.
    3. For a new node in c1 (node1), we will keep updating the list of distances for each category c2.
       If any node in c2 is disconnected from node1 (not in node1 tree) we consider it to be at a distance
       from equal to the max distance of node1 tree + 1.
    4. Finally, we will take the median of distances for each category and sort the dictionary.
    '''
    
    #step 1
    for node1 in c1_nodes:
        if node1 in out_links.keys():
            distance_tree = explore(node1,None, out_links)
            t.append(max(list(distance_tree.keys())) + 1)
    #step 2
            for cat in category_nodes.keys():  
                dist = []
                disconnected = []
                for node1_dist in distance_tree.keys():
                    c2_nodes = [int(i) for i in category_nodes[cat]]
                    nodes_in_tree = len(set(c2_nodes).intersection(set(distance_tree[node1_dist])))
                    dist += [node1_dist]*nodes_in_tree # number of nodes at distance d x distance d
    # step 3
                disconnected = [max(distance_tree.keys())+1]*(len(c2_nodes)-len(dist)) 
                distances_dict[cat] += dist + disconnected
    # step 4
    for key in distances_dict.keys():
        distances_dict[key] = statistics.median(distances_dict[key])
    
    distances_dict = {k: v for k, v in sorted(distances_dict.items(), key=lambda item: item[1])}
    threshold = statistics.median(t)
    d_filt= {k: v for k,v in distances_dict.items() if distances_dict[k] < threshold}
    
    return d_filt 

####################################### RQ6 ##########################################


def transform(links, category_nodes, nodes_category):
    
    ''' We transform our graph into one which has categories as nodes. 
    In order to do this we consider the inlinks and outlinks of each node separately.
    1. For each category we concatenate the nodes associated to any of the nodes of that category.
    2. We weight links by counting how many nodes in a category are linked to another node (or other nodes)
    in another category.'''
    
    links_cat = {k: [] for k in category_nodes.keys()}
    
    for node in nodes_category.keys():
        if int(node) in links.keys():
            for link in links[int(node)]:
                if link in nodes_category.keys():
                    if links_cat[nodes_category[node]] != nodes_category[link]:    # to avoid selflinks
                        links_cat[nodes_category[node]].append(nodes_category[link])

    for cat in links_cat.keys():
        categories = {}
        to_transform = set(links_cat[cat])
        for node in to_transform:
            categories[node] = links_cat[cat].count(node)
        links_cat[cat] = categories
                                     
    return links_cat

def pagerank(n_iter, in_links_cat, out_links_cat):
    
    ''' INITIALIZATION
    This page rank algorithm initialize the ranking vector assuming that all pages have the same importance.
    1. we store a list of category names with fixed index that we use for the ranking vector
    2. we initalize all values as 1 / number of nodes
    '''
    categories_names = list(set(in_links_cat.keys()).union(set(out_links_cat.keys())))
    init = [1/len(categories_names)]*len(categories_names)
    
    
    ''' ITERATIONS
    1. At each iteration we initalize a new vector to store the ranking score of each node.
    2. For each node we initialize a 0 score and save the nodes (In) that have an outlink to that node (node1).
    3. For each node (node2) that has an outlink to the analyzed node (node1), we define:
        a) Probability : probability that node1 is reached through node2
        b) Rank : rank score of node2 in the previous iteration
    4. We multiply probability of node 2 by its ranking score
    5. We sum this score for all nodes in In and assign this final score to the position of node1 in the rank vector.
    6. We update init vector.
     '''
    for i in range(n_iter):
        new_iter = [0]*len(categories_names)
        for node1 in categories_names:
            score = 0
            pos1 = categories_names.index(node1)
            if node1 in in_links_cat.keys():  
                In = in_links_cat[node1]
                for node2 in In:
                    pos2 = categories_names.index(node2)
                    prob2 = out_links_cat[node2][node1]/sum(out_links_cat[node2].values())
                    rank2 = init[pos2]
                    score += rank2*prob2
            new_iter[pos1] = score
        init = [float(i)/sum(new_iter) for i in new_iter]
        
    return dict(zip(categories_names,init))