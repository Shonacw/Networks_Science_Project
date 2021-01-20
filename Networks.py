#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 10:53:24 2020

@author: ShonaCW
"""

import networkx as nx
import random
import matplotlib.pyplot as plt 
import logbin230119 as logbin
import numpy as np


def Add(G, Og, Number_Edges):
    """
    Function which adds integer 'Number_Edges' to the node 'Og'.
    """
    G.add_node(Og)
    edges_added = 0 #keep track of how many edges added
    Num = len(G.edges) #only consider previously existing nodes/edges 
    while True:
        #randomly choose an existing edge
        #randomly choose one of the nodes it is connected to
        node = list(G.edges)[random.randint(0, Num-1)][random.randint(0,1)] 
        #check that this edge does not already exist
        Binary = (node, Og) in G.edges
        if Binary == False:
            #if this edge does not already exist, create the edge
            G.add_edge(node, Og)
            edges_added +=1 #keep track of number of edges added

        if edges_added == Number_Edges:
            #break out loop once required number of edges have been added
            break



def Function(m, Num_nodes):
    """
    Function for the BA model. Inputs are 'm': the number of edges we wish to 
    add to each new node, and 'Num_nodes': the total number of nodes we would 
    like to add to the system.
    Returns the nx graph 'G'.
    """
    G = nx.Graph()
    
    #First, initialise. 
    for i in range(m+1): #if only adding one edge need two nodes so if N=1 this will add node 0 and node 1
        G.add_node(i)
    for i in range(m):
        G.add_edge(i, i+1) #making a line
        
    #now add remaining nodes + edges between them
    for i in range(m+1, Num_nodes):
        Add(G, i, m) #add N edges to node i 

    return G

def p_inft_theor(m, k):
    """
    Function to calculate the statistic probability distribution expected from
    theory. 
    """
    A = 2 * m * (m + 1)
    result = A / (k * (k + 1) * (k + 2))
    
    return result
    
def segment_list(List, n = int): 
    """
    Function to split up a list 'List' into sublists of length 'n'
    """
    segmented = [List[x:x+n] for x in range(0, len(List), n)]
    return segmented
        
#%%

G = nx.Graph() #.add_node(i) for i in range(4)]
for i in range(4):
   G.add_node(i)
    
G.add_edge(0, 1)
G.add_edge(1, 2)
G.add_edge(2, 3)
print('original', G.edges)

for i in range(4, 20):
    Add(i , 3) #must have at least the same number of nodes ALREADY DEFINED as how many edges to add
#%%
x=[0,1,2,3,4,5,6,7,8,9]
res = segment_list(x, 2)
print(res)
#%%    
M = [1, 2, 4] # 8, 16, 32, 64]
C = ['r', 'y', 'pink', 'g', 'b', 'c', 'k', 'magenta']
#%%
#degree_list = [] #will contain list of lists
#Define a number of functions to speed up processing
Append = np.append
Degree_hist = nx.degree_histogram
Linspace = np.linspace
Logbin = logbin.logbin
Transpose = np.transpose
Mean = np.mean
Degree = []
Degree_save = []

repeats = 10 
for m in M:
    Degree = []
    for i in range(repeats):
        G = Function(m, 1000) #evetually do 100,000 and also repeats 
        data = list(zip(*G.degree))[1] #extracts the degrees 
        Degree = Append(Degree, data) #append data to 
    print('DONE')
    split_data = segment_list(Degree, int((len(Degree) / int(repeats)))) #split list into 10 sublists so can take average
    Degree_save = Append(Degree_save, Mean(split_data, axis = 0)) #save average of the lists. this should be same length as M

#save data to file
#data is one long list that will require splitting where necessary
np.savetxt('degrees.csv', Degree_save, delimiter=',')

#%%
"""GRAPH: DEGREES AGAINST K"""
#load required data
Degrees = np.loadtxt('degrees.csv', delimiter=',')
        
degree_list = list(segment_list(Degrees, 1000)) #split data into each m 

plt.figure(2)
plt.xlabel('k')
plt.ylabel(r'$p_{\infty}$(k)')

for i in [0,1,2]:
    #plot data
    x, y = Logbin(list(degree_list[i]), 1.3)
    plt.plot(x, y, 'x', color = C[i])
    
    #plot theoretical expectation
    x_theor = np.linspace(min(x), max(x))
    y_theor = [p_inft_theor(i+1, i) for i in x_theor]
    plt.plot(x_theor, y_theor, label= f'm={M[i]}', color = C[i])

plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()
#%% 

#old code (works!) for just first two m sizes    
G_1 = Function(1, 1000)
G_2 = Function(2, 1000)

data_1 = list(zip(*G_1.degree))[1] #this is the data we want to save 
data_2 = list(zip(*G_2.degree))[1]

#one = nx.degree_histogram(G_1)
#two = nx.degree_histogram(G_2)

x, y = logbin.logbin(data_1, 1.3)
x_2, y_2 = logbin.logbin(data_2, 1.3)


x_theor_1 = np.linspace(min(x), max(x))
x_theor_2 = np.linspace(min(x), max(x))
y_theor_1 = [p_inft_theor(1, i) for i in x_theor_1]
y_theor_2 = [p_inft_theor(2, j) for j in x_theor_2]


plt.figure(3)
plt.plot(x, y, 'x', color = 'r')
plt.plot(x_2, y_2, 'x', color = 'b')
plt.plot(x_theor_1, y_theor_1, color = 'r', label='m=1')
plt.plot(x_theor_2, y_theor_2, color = 'b', label='m=2')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('k')
plt.ylabel(r'$p_{\infty}$(k)')
plt.legend()
