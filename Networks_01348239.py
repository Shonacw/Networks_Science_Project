#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 10:46:43 2020

@author: ShonaCW
"""

import random
import matplotlib.pyplot as plt 
import logbin230119 as logbin
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import chi2_contingency
from scipy.stats import ks_2samp

#CODE IS SPLIT INTO THREE SECTIONS:
#1: BA Model with Random Attachment
#2: BA Model with Preferential Attachment
#3: BA Model with Random Walk
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""FUNCTIONS FOR SECTION 1"""
def Function(m, Num_nodes, LK = False):
    """
    Function for the BA model. Inputs are 'm': the number of edges we wish to 
    add to each new node, and 'Num_nodes': the total number of nodes we would 
    like to add to the system.
    Returns the graph 'G'.
    
    Also has option to measure the maximum degree of the Graph. This is measured
    at N(t)= Num_nodes, i.e. the max degree of the graph at the point when all 
    nodes have been added.
    """
    G = [[]] #each sublist will contain the nodes the element is connected to
    G_connections = [] #to keep track of which nodes have been connected, for prob
    
    #to increase efficiency
    Append_Gconn = G_connections.append
    Append_node = G.append

    
    #First, initialise.                                                          
    for i in range(2*m): #so that 2m+1 nodes are addd                      \
        #add node
        numnodes = len(G) #store number of nodes so far
        Append_node([])
        for j in range(numnodes):
            G[i+1].append(j)
            G[j].append(i+1)
            #add straight line edges
            Append_Gconn(i+1)
            Append_Gconn(j)

       
    #now add remaining nodes + edges between them
    for i in range(2*m + 1, Num_nodes):                                        
        Add(G, G_connections, i, m) #add m edges to node i
        
        #un-hashtag the following if need to check the code is working correctly
        #print('G:', len(G))
        #print('G_CONN:', len(G_connections))
    
    Largest_K = len(max(G, key=len))
    if LK==True:
        return Largest_K
    
    return G

def Function_K(m, Num_nodes):
    """
    Function for the BA model. Inputs are 'm': the number of edges we wish to 
    add to each new node, and 'Num_nodes': the total number of nodes we would 
    like to add to the system.
    This function measures the maximum degree of the graph at various points 
    during the simulation run, namely when N(t) = a square number.
    """
    G = [[]] #each sublist will contain the nodes the element is connected to
    G_connections = [] #to keep track of which nodes have been connected, for prob
    Largest_K = [] #store largest k values (k_1) per iteration
    
    #to increase efficiency
    Append_Gconn = G_connections.append
    Append_node = G.append
    
    #Points at which the largest k will be measured (exponentially spaced)
    numbers = np.linspace(1,1000,1000).tolist()
    points_to_measure = [i ** 2 for i in numbers]

    points_to_measure_updated = [x for x in points_to_measure if x>= (2*m+1)]
    
    #First, initialise.                                
    for i in range(2*m): #so that 2m+1 nodes are added
        #add node
        numnodes = len(G) #store number of nodes so far
        Append_node([])
        for j in range(numnodes):
            G[i+1].append(j)
            G[j].append(i+1)
            #add straight line edges
            Append_Gconn(i+1)
            Append_Gconn(j)
    
       
    #now add remaining nodes + edges between them
    for i in range(2*m + 1, Num_nodes):                
        Add(G, G_connections, i, m) #add m edges to node i
        if i+1 in points_to_measure_updated:   #change back to i
            Largest_K.append(len(max(G, key=len)))
    
    return Largest_K

def Add(G, G_connections, Og, m):
    """
    Function which adds integer 'Number_Edges' to the node 'Og'.
    """
    Random = random.randint
    Append_Gconn = G_connections.append
    Append_node = G.append
    
    new_connections = [] #of max length m. nodes to which 'Og' is being connected
    Num = len(G_connections) #only consider previously existing nodes/edges
    Append_node([]) #add one new node
    while True:
        #randomly choose a node to connect edge with
        node = G_connections[Random(0, Num-1)]
        #check that this edge does not already exist
        if node not in new_connections:
            #if this edge does not already exist, create the edge
            G[Og].append(node) 
            G[node].append(Og) 
            new_connections.append(node)
            Append_Gconn(Og)
            Append_Gconn(node)

        if len(new_connections) == m:
            #break out loop once required number of edges have been added
            break

def Degrees_G(Graph):
    """
    Function to find the degree list of Graph data
    """
    degrees_list = [len(i) for i in Graph]
    return degrees_list


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
    Function to split up a list 'List' into sublists of length 'n'.
    """
    segmented = [List[x:x+n] for x in range(0, len(List), n)]
    return segmented

def Largest_Degree(N, m):
    """
    Function which returns the theoretical value of k_1, for N iterations of 
    a system size m. 
    """
    result = (-1 + (1 + 4*N*m*(m + 1))**0.5)/ 2
    return result


#define some general functions that are used repeatedly
Mean = np.mean
Logbin=logbin.logbin

#define the m values I will be testing
M = [ 2, 4, 8, 16, 32, 64]

#define the colours which i will use for plotting different m, for consistency
C = ['r', 'y', 'pink', 'g', 'b', 'c', 'k', 'magenta']
#%%
"""PRELIMINARY TESTS"""
#for checking that the system is working correctly
#un-hashtag the "print(len(G))" commands from 'Function(m, Num_nodes)' above.
X = Function(2, 20)
print(X)

#%%
"""COLLECTION DATA SECTION 1 but without taking mean"""
M = [2, 4, 8, 16, 32, 64]
C = ['r', 'y', 'pink', 'g', 'b', 'c', 'k', 'magenta']

Degree_save_1 = []

repeats = 50 
N = 100000
for m in M:
    Degree = []
    for i in range(repeats):
        G = Function(m, N)
        data_1 = Degrees_G(G) #extracts the degrees 
        Degree.append(data_1) #append data to 'Degree' list
    print('DONE')
    
    deg = [val for sublist in Degree for val in sublist]
    #save data
    np.savetxt('degrees_1_nomean{}.csv'.format(m), deg, delimiter=',')
    Degree_save_1.append(deg) 

#%%
"""GRAPH: DEGREES AGAINST K"""
Degrees = Degree_save_1

data_y = []
theory = []

plt.figure()
plt.xlabel('k')
plt.ylabel('$p_{N}(k)$')

for i in range(len(M)):
    #plot data
    x, y = Logbin(list(Degrees[i]), 1.3)
    plt.plot(x, y, 'x', color = C[i])
    data_y.append(y)
    
    #plot theoretical expectation
    y_theor = [p_inft_theor(M[i], j) for j in x]
    theory.append(y_theor)
    plt.plot(x, y_theor, label= f'm={M[i]}', color = C[i])

plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()

"""DATA VS. THEORY FOR SECTION 1"""
plt.figure()
plt.xlabel(r'$p_{\infty}$(k)')
plt.ylabel(r'$p_{N}$(k)')
for i in range(len(M)):

    plt.plot(theory[i], data_y[i], 'x', color = C[i], label = f'm={M[i]}')
    
    x = np.linspace(min(theory[i]), max(theory[i]), len(theory[i]))
    y = x
    plt.plot(x, y, C= 'k')     #this is only fitting to the last set of data shoudl really fit to all
    
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()

"""CHI SQUARED GOODNESS OF FIT CODE"""
dataa = []
chi_vals = []
p_vals = []
for i in range(len(M)):
    Contingency = []
    Contingency.append(data_y[i])
    Contingency.append(theory[i])
    dataa.append(chi2_contingency(Contingency))
    
    chi_vals.append(dataa[i][0])
    p_vals.append(dataa[i][1])

print('Chi_vals:', chi_vals)
print('P_vals:', p_vals)

"""KS STATS GOODNESS OF FIT CODE"""
dataa = []
KS_vals = []
p_vals_KS = []
for i in range(len(M)):
    dataa.append(ks_2samp(data_y[i], theory[i]))
    
    KS_vals.append(dataa[i][0])
    p_vals_KS.append(dataa[i][1])
print('KS_vals:', KS_vals)
print('P_vals_ks:', p_vals_KS)

#%%
"""LARGEST AVERAGE DEGREE RUN: K1 MEASURED JUST AT END OF ITERATION"""

m = 2 
repeats = 100

N_l = np.linspace(100, 100000, 21) #testing 21 different N's
N_int = [int(i) for i in N_l] #N's must be integers not floats

#define required lists
midi_ks_all = []
ks_all = []
ks_error = []

for i in N_int:
    midi_ks_all = []
    for j in range(repeats):
        K1 = Function(m, int(i), LK = True)
        midi_ks_all.append(K1) #K1 here is a value, measured at end of each run
    ks_all.append(np.mean(midi_ks_all)) #store the mean values of k1 for each N
    ks_error.append(np.std(midi_ks_all)) #store the std of k1's for each N
    print('DONE')

#save data
np.savetxt('degrees_1_all.csv', ks_all, delimiter=',')
np.savetxt('errors_degrees_1_all.csv', ks_error, delimiter=',')
#%% 
"""LARGEST AVERAGE DEGREE RUN: K1 MEASURED MULTIPLE TIMES DURING ITERATION"""

m = 2
runs = 30
N = 1000000 #KEEP AS A MILLION

#define list to store the lists of max degrees during each run
ks_all_ = []

for j in range(runs):
    K1 = Function_K(m, N)
    ks_all_.append(K1) #K1 here is a list, measured when N(t) is a sqr numb.
    print('DONE')

#save data
np.savetxt('degrees_1_all_MIA_2.csv', ks_all_, delimiter=',')
#%%
"""PLOTTING: K1 MEASURED JUST AT END OF ITERATION, variable N"""
#import data
data = np.loadtxt('degrees_1_all.csv', delimiter=',')
errors_IMP = np.loadtxt('errors_degrees_1_all.csv', delimiter=',')
repeats = 10

#Standard Error: standard deviation / sqrt(number of repeats)
Sqrt_Rep = np.sqrt(repeats)
errors = [i / Sqrt_Rep for i in errors_IMP]

#as k1 scales with sqrt(N), calculate sqrt(N)'s
SQRT_N = [(i)**0.5 for i in N_int]


"""Calculating the error on the best-fit gradient for K1 vs. SQRT(N)"""
data_upper = [data[i] + errors[i] for i in range(len(data))]
P_1 = np.polyfit(SQRT_N, data_upper, 1)
Poly_Info_1 = np.poly1d(P_1)
gradient_1 = round(Poly_Info_1[1], 4)

data_lower = [(data[i] - errors[i]) for i in range(len(data))]
P_2 = np.polyfit(SQRT_N, data_lower, 1)
Poly_Info_2 = np.poly1d(P_2)
gradient_2 = round(Poly_Info_2[1], 4)

error_in_grad = (gradient_1 - gradient_2)/2
print("Error in gradient", error_in_grad)


"""PLOTTING K1 AGAINST SQRT(N)"""
plt.figure()
P, cov = np.polyfit(SQRT_N, data, 1, cov=True)
Poly_Info = np.poly1d(P)
gradient = round(Poly_Info[1], 4)
#plot data
plt.errorbar(SQRT_N, data, yerr=errors, fmt='o')
#plot linear fit
plt.plot(SQRT_N, Poly_Info(SQRT_N), label=f'Linear Fit Slope={gradient} $\pm$ {round(error_in_grad,4)}')
#also plot theory
y = [np.sqrt(m*(m + 1))* i - 0.5 for i in SQRT_N]
plt.plot(SQRT_N, y, '--', label='Theoretical Slope = 2.4494')
plt.xlabel(r'$\sqrt{N}$')
plt.ylabel(r'$ \langle k_1 \rangle$')
plt.legend()
plt.show()
#find Pearsons fit
print("Pearsons test for SQRTN vs ks_all:", pearsonr(SQRT_N, data))


"""PLOTTING ERRORS AGAINST SQRT(N)"""
plt.figure()
#ks_error = [0.5*i for i in errors]
P_std = np.polyfit(SQRT_N, errors_IMP, 1)
Poly_Info_std = np.poly1d(P_std)
gradient_std = round(Poly_Info_std[1], 4)
print(gradient_std)
#plot data
plt.plot(SQRT_N, errors_IMP, 'o')
#plot linear fit
plt.plot(SQRT_N, Poly_Info_std(SQRT_N), label=f'best fit, gradient={gradient_std}')
plt.xlabel(r'$\sqrt{N}$')
plt.ylabel(r'$\sigma$')
plt.title('Errors in K_1 against SQRT(N)')
plt.legend()
plt.show()
print("Pearsons test for SQRTN vs std devs:", pearsonr(SQRT_N, errors_IMP))

#%%
"""PLOTTING: K1 MEASURED MULTIPLE TIMES DURING ITERATION"""
m=2
repeats = 10

#import data
midi_ks = np.loadtxt('degrees_1_all_MIA.csv', delimiter=',') #10 runs (better)
#midi_ks = np.loadtxt('degrees_1_all_MIA_2.csv', delimiter=',') #30 runs

#re-define the N(t) values at which the Graph's max degree is being measured
numbers = np.linspace(1,1000,1000).tolist()
points_to_measure = [i ** 2 for i in numbers]
Ns = [x for x in points_to_measure if x>= (2*m+1)]
SQRT_N = [j**0.5 for j in Ns]

#calculate means standard deviations for plotting
mean_ks_imported = Mean(midi_ks, axis = 0) #save average of the lists
standard_devs = list(np.std(midi_ks, axis = 0) )

#select fewer data points to plot to improve visualisation
SQRT_N_fewer = SQRT_N[100::100]
mean_ks = mean_ks_imported[100::100]
standard_devs_fewer = standard_devs[100::100]

#Standard Error: standard deviation / sqrt(number of repeats)
errors_proper = [i / np.sqrt(repeats) for i in standard_devs_fewer]

"""PLOT OF ALL K1s against N"""
plt.figure(4)
for i in range(len(midi_ks)):
    plt.plot(Ns, midi_ks[i], 'x-', markersize=2, linewidth=2)
plt.xlabel('N')
plt.ylabel(r'$k_{1}$')
plt.show()

"""PLOT OF ALL K1s against SQRT(N)"""
plt.figure(5)
for i in range(len(midi_ks)):
    plt.plot(SQRT_N, midi_ks[i],  'x-', markersize=2, linewidth=2)
plt.xlabel(r'$\sqrt{N}$')
plt.ylabel(r'$k_{1}$')
plt.show()

"""PLOTTING STANDARD DEVIATIONS AGAINST SQRT(N)"""
plt.figure(7)
P_std = np.polyfit(SQRT_N_fewer, standard_devs_fewer, 1)
Poly_Info_std = np.poly1d(P_std)
gradient_1 = round(Poly_Info_std[1], 4) #shoudl be 1 
#plot data
plt.plot(SQRT_N_fewer, standard_devs_fewer, 'o')
#plot linear fit
plt.plot(SQRT_N_fewer, Poly_Info_std(SQRT_N_fewer), label=f'Slope={gradient_1}')
plt.xlabel(r'$\sqrt{N}$')
plt.ylabel(r'$\sigma$')
plt.legend()
plt.show()


"""PLOT TO FIND ERROR IN K1"""
plt.figure()
#errors is a list of the standard devs
y = [mean_ks_imported[i]/ standard_devs[i] for i in range(len(mean_ks_imported))]
plt.plot(SQRT_N, y, 'x')
tends = np.mean(y[-100:]) #take average value of last three points 
print(tends)
plt.plot([0, max(SQRT_N)], [tends, tends], '--', label=f'y={round(tends,4)}')
plt.xlabel(r'$\sqrt{N}$')
plt.ylabel(r'$\langle k_1 \rangle$/ $\sigma$')
plt.legend()
plt.show()

"""Calculating the error on the best-fit gradient for K1 vs. SQRT(N)"""
data_upper = [mean_ks[i] + errors_proper[i] for i in range(len(mean_ks))]
P_1 = np.polyfit(SQRT_N_fewer, data_upper, 1)
Poly_Info_1 = np.poly1d(P_1)
gradient_1 = round(Poly_Info_1[1], 4)

data_lower = [(mean_ks[i] - errors_proper[i]) for i in range(len(mean_ks))]
P_2 = np.polyfit(SQRT_N_fewer, data_lower, 1)
Poly_Info_2 = np.poly1d(P_2)
gradient_2 = round(Poly_Info_2[1], 4)

error_in_grad = (gradient_1 - gradient_2)/2
print("Error in gradient", error_in_grad)

"""AVERAGE K1 AGAINST SQRT(N)""" #including the standard deviation on each k
plt.figure()
P, cov = np.polyfit(SQRT_N_fewer, mean_ks, 1, cov=True)
Poly_Info = np.poly1d(P)
gradient_2 = round(Poly_Info[1], 4)
#plot data
plt.errorbar(SQRT_N_fewer, mean_ks, yerr=errors_proper, fmt='o')
#plot linear fit
plt.plot(SQRT_N, Poly_Info(SQRT_N), label=f'Linear Fit Slope={gradient_2} $\pm$ {round(error_in_grad, 4)}')
#plot theory
sqrt_N_fewer = SQRT_N_fewer
sqrt_N_fewer.append(1000)
sqrt_N_fewer.append(0)
y = [np.sqrt(m*(m + 1))* i - 0.5 for i in SQRT_N_fewer]
plt.plot(SQRT_N_fewer, y, '--', label='Theoretical Slope = 2.4494')

plt.xlabel(r'$\sqrt{N}$')
plt.ylabel(r'$\langle k_1 \rangle$')
plt.legend()
plt.show()
print("Pearsons test for SQRTN vs ks_all:", pearsonr(SQRT_N, mean_ks_imported))
#%%
"""COLLECTING DATA FOR DATA COLLAPSE"""
m = 2
N = [100, 1000, 10000, 100000, 1000000]
repeats = 100 

#define lists
largest_k = [] #will store largest k for one run-through 
Degree_save_N = [] #list of mean degrees

Deg_N_Append = Degree_save_N.append

for j in N:
    Degree_N = [] #for gathering lists of degrees for each N
    for i in range(repeats):
        G = Function(m, j)
        data_N = Degrees_G(G) #extract the degrees 
        
        #save lists 
        Degree_N.append(data_N)
    
    #save average of the lists
    print('DONE')
    deg = [val for sublist in Degree_N for val in sublist]
    np.savetxt('degrees_N{}.csv'.format(j), deg, delimiter=',') #save data
    Degree_save_N.append(deg)


#%%
"""LOAD COLLECTED DATA"""
N = [1000, 10000, 100000, 1000000]

Degree_save_N_ = [] 
for i in N:
    Degree_save_N_.append(np.loadtxt('degrees_N{}.csv'.format(i), delimiter=','))
#%%
"""PLOTTING THE DEGREE DATA FOR DIFFERENT N"""
#Logbin data and plot
N = [100, 1000, 10000, 100000, 1000000]
for i in range(len(N)):
    x, y = Logbin(list(Degree_save_N[i]), 1.3)
    plt.plot(x,y, 'x-', label=f'N={N[i]}')
y_theor = [p_inft_theor(2, ll) for ll in x]
plt.plot(x, y_theor)
plt.xlabel('k')
plt.ylabel(r'$p_{N}$(k)')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()

"""PLOTTING THE DATA COLLAPSE"""
plt.figure(8)
plt.xlabel(r'k/$\sqrt{N}$')
plt.ylabel(r'$p_{N}$(k) / $p_{\infty}$(k)')

N = [100, 1000, 10000, 100000, 1000000]

for i in range(len(N)):
    xx =[]
    yy = []
    #logbin data
    x, y = Logbin(list(Degree_save_N[i]), 1.3)

    #calculate theoretical expectation from logbinned data
    y_theor = [p_inft_theor(2, ll) for ll in x] #as m=2, x or x_theor? 
    SQRTN = np.sqrt(N[i])
    for j in range(len(x)):
        #calculate manipulated data for data collapse
        yy.append(y[j] / y_theor[j])
        xx.append(x[j] / SQRTN) #k/ k1
    
    plt.plot(xx, yy, '.', markersize=8, label=f'N={N[i]}')
    
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""SECTION TWOO"""
def Function_2(m, Num_nodes):
    """
    Function for the BA model. Inputs are 'm': the number of edges we wish to 
    add to each new node, and 'Num_nodes': the total number of nodes we would 
    like to add to the system.
    Returns the graph 'G'.
    """
    G = [[]] #each sublist will contain the nodes the element is connected to

    
    #to increase efficiency
    Append_node = G.append
    
    #First, initialise. 
    for i in range(2*m): #so that 2m+1 nodes are added
        #add node
        numnodes = len(G) #store number of nodes so far
        Append_node([])
        for j in range(numnodes):
            G[i+1].append(j)
            G[j].append(i+1)
            
    #now add remaining nodes + edges between them
    for i in range(2*m +1, Num_nodes): 
        Add_2(G, i, m) #add m edges to node i
    
    return G

def Function_K_2(m, Num_nodes):
    """
    Function for the BA model. Inputs are 'm': the number of edges we wish to 
    add to each new node, and 'Num_nodes': the total number of nodes we would 
    like to add to the system.
    This function measures the maximum degree of the graph at various points 
    during the simulation run, namely when N(t) = a square number.
    """
    G = [[]] #each sublist will contain the nodes the element is connected to
    Largest_K = [] #store largest k values (k_1) per iteration
    
    #to increase efficiency
    Append_node = G.append
    
    #Points at which the largest k will be measured (exponentially spaced)
    numbers = np.linspace(1,1000,1000).tolist()
    points_to_measure = [i ** 2 for i in numbers]

    points_to_measure_updated = [x for x in points_to_measure if x>= (2*m+1)]
    
    #First, initialise.                                
    for i in range(2*m): #so that 2m+1 nodes are added
        #add node
        numnodes = len(G) #store number of nodes so far
        Append_node([])
        for j in range(numnodes):
            G[i+1].append(j)
            G[j].append(i+1)
    
       
    #now add remaining nodes + edges between them
    for i in range(2*m + 1, Num_nodes):
        Add_2(G, i, m) #add m edges to node i
        if i+1 in points_to_measure_updated: 
            Largest_K.append(len(max(G, key=len)))
    
    return Largest_K

def Add_2(G, Og, m):
    """
    Function which adds integer 'Number_Edges' to the node 'Og'.
    """
    Random = random.randint
    Append_node = G.append
    
    Append_node([]) #add one new node
    new_connections = [] #of max length m/nodes to which 'Og' is being connected
    while True:
        #randomly choose a node to connect edge with
        node = Random(0, len(G)-1) 
        #check that this edge does not already exist
        if node not in new_connections:
            #if this edge does not already exist, create the edge
            G[Og].append(node)
            G[node].append(Og)
            new_connections.append(node)

        if len(new_connections) == m:
            #break out loop once required number of edges have been added
            break
        
def p_inft_theor_2(m, k):
    """
    Function to calculate the statistic probability distribution expected from
    theory. 
    """

    result = (1/ m) * ((m / (m+1))**(k-m))
    
    return result

def Largest_Degree_2(N, m):
    """
    Function which returns the theoretical value of k_1, for N iterations of 
    a system size m. 
    """
    result = m + (-1/(np.log(m/(m+1))))*np.log(N)
    return result

#%%
"""COLLECTION DATA SECTION 2"""
M = [2, 4, 8, 16, 32, 64]
C = ['r', 'y', 'pink', 'g', 'b', 'c', 'k', 'magenta']

Degree_save_2 = []

repeats = 50 
N = 100000
for m in M:
    Degree = []
    for i in range(repeats):
        G = Function_2(m, N) 
        data_2 = Degrees_G(G) #extracts the degrees 
        Degree.append(data_2) #append data to 'Degree' list
    print('DONE')
    
    deg = [val for sublist in Degree for val in sublist]
    np.savetxt('degrees_2{}.csv'.format(m), deg, delimiter=',') #save data
    Degree_save_2.append(deg) 
#%%
"""COLLECTION DATA FOR SECTION 3 COMPARISON"""
m = 2
N = 100000
G = Function_2(m, N)
data_2 = Degrees_G(G)

#%%
"""LOAD DEGREE DATA"""
M = [1, 2, 4, 8, 16, 32, 64]
C = ['r', 'y', 'pink', 'g', 'b', 'c', 'k', 'magenta']
Degree_save_2_ = []

for i in M:
    Degree_save_2_.append(np.loadtxt('degrees_2{}.csv'.format(i), delimiter=','))
    print('DONE')
#%%
"""p(k) vs k"""
data_y_2 = []
theory_2 = []

plt.figure()
plt.xlabel('k')
plt.ylabel(r'$p_{N}(k)$')

for i in range(len(M)):
    #plot data
    x, y = Logbin(list(Degree_save_2[i]))
    plt.plot(x, y, 'x', markersize=4, color = C[i])
    
    #plot theoretical expectation
    x_theor = np.linspace(min(x), max(x), len(x))
    y_theor = [p_inft_theor_2(M[i], j) for j in x]
    data_y_2.append(y)
    theory_2.append(y_theor)
    plt.plot(x, y_theor, label= f'm={M[i]}', color = C[i])

plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()

#%%
"""DATA VS THEORY P VALUES FOR SECTION 2"""
plt.figure()
plt.xlabel(r'$p_{\infty}$(k)')
plt.ylabel(r'$p_{N}$(k)')
for i in range(len(M)-1):
    #plot data
    print(i)
    plt.plot(data_y_2[i], theory_2[i], 'x', label = f'm={M[i]}', color = C[i])
    #big_theory = [p_inft_theor_2(M[i], j) for j in Degrees_2[i]]
    #plt.plot(big_theory, Degrees_2[i], 'x', label= f'm={M[i]}', color = C[i])
    x = np.linspace(min(theory_2[i]), max(theory_2[i]), len(theory_2[i]))
    y = x
    plt.plot(x, y, C= 'k')
    
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()


"""CHI SQUARED GOODNESS OF FIT CODE PART 2"""
dataa_2 = []
chi_vals_2 = []
p_vals_2 = []
for i in range(len(M)-1):
    Contingency_2 = []
    Contingency_2.append(data_y_2[i])
    Contingency_2.append(theory_2[i])
    dataa_2.append(chi2_contingency(Contingency_2))
    
    chi_vals_2.append(dataa_2[i][0])
    p_vals_2.append(dataa_2[i][1])

print('Chi_vals:', chi_vals_2)
print('P_vals:', p_vals_2)      

"""KS STATS GOODNESS OF FIT CODE"""
dataa = []
KS_vals = []
p_vals_KS = []
for i in range(len(M)):
    dataa.append(ks_2samp(data_y_2[i], theory_2[i]))
    
    KS_vals.append(dataa[i][0])
    p_vals_KS.append(dataa[i][1])
print('KS_vals:', KS_vals)
print('P_vals_ks:', p_vals_KS)

#%%
"""COLLECTING DATA: LARGEST K SECTION 2"""

m = 2
N_l = np.linspace(100, 100000, 21)
repeats = 30

midi_ks_all_2 = [] #for gathering lists
ks_all = []
ks_error = []
ks_nonmean = []

for i in N_l:
    midi_ks_all_2 = []
    for j in range(repeats):
        K1 = Function_K_2(m, int(i), LK = True)
        midi_ks_all_2.append(K1)

    ks_all.append(np.mean(midi_ks_all_2))
    ks_error.append(np.std(midi_ks_all_2))
    print('DONE')

np.savetxt('sec2_ degrees_1_all.csv', ks_all, delimiter=',')
np.savetxt('sec_2errors_degrees_1_all.csv', ks_error, delimiter=',')

#%%
"""LARGEST AVERAGE DEGREE RUN: K1 MEASURED MULTIPLE TIMES DURING ITERATION""" 

m = 2
runs = 30
N = 1000000 #KEEP AS A MILLION

#define list to store the lists of max degrees during each run
ks_all_MIA_2 = []

for j in range(runs):
    K1 = Function_K_2(m, N)
    ks_all_MIA_2.append(K1) #K1 here is a list, measured when N(t) is a sqr numb.
    print('DONE')

#save data
np.savetxt('sec2_ degrees_1_all_MIA.csv', ks_all_MIA_2, delimiter=',')
#%%
"""PLOT OF ALL K1s against N SECTION 2"""
midi_ks_2 =  np.loadtxt('sec2_ degrees_1_all.csv', delimiter=',')
errors_2 = np.loadtxt('sec_2errors_degrees_1_all.csv', delimiter=',')

errors_proper = [i / np.sqrt(repeats) for i in errors_2]

N_l = np.linspace(100, 100000, 21) #testing 21 different N's 
N_int = [int(i) for i in N_l]

plt.figure(10)
plt.plot(N_l, midi_ks_2,  'x-', markersize=2, linewidth=2)
plt.xlabel('N')
plt.ylabel(r'$k_{1}$')
plt.title(f'N={N}, m={m}, Repeats={repeats}')
plt.show()


"""AVERAGE K1 AGAINST log(N) SECTION 2"""#including the std on each k
plt.figure(11) 
x = [np.log(n) for n in N_l]
y = [Largest_Degree_2(i, m) for i in N_l]

P_2, cov_2 = np.polyfit(x, midi_ks_2, 1, cov=True)
Poly_Info_2 = np.poly1d(P_2)
plt.plot(x, Poly_Info_2(x), label=f'Linear Fit Slope={gradient_2}')
gradient_2 = round(Poly_Info_2[1], 4)
plt.errorbar(x, midi_ks_2, yerr=errors_proper, fmt='o', zorder=1)

P_2, cov_2 = np.polyfit(x, y, 1, cov=True)
Poly_Info_2 = np.poly1d(P_2)
gradient_22 = round(Poly_Info_2[1], 4)
plt.plot(x, y, '--', label=f'Theoretical Slope={gradient_22}') #plotting theory
plt.xlabel(r'log(N)')
plt.ylabel(r'$k_{1}$')
plt.legend()
plt.show()
#%%
"""PLOTTING: K1 MEASURED MULTIPLE TIMES DURING ITERATION"""
m = 2
#import data
midi_ks_2 = np.loadtxt('sec2_ degrees_1_all_MIA.csv', delimiter=',') #10 runs
#midi_ks = np.loadtxt('degrees_1_all_MIA_2.csv', delimiter=',') #30 runs

#re-define the N(t) values at which the Graph's max degree is being measured
numbers = np.linspace(1,1000,1000).tolist()
points_to_measure = [i ** 2 for i in numbers]
Ns = [x for x in points_to_measure if x>= (2*m+1)]
LOG_N = [np.log(j) for j in Ns]

#calculate means standard deviations for plotting
mean_ks_imported = Mean(midi_ks_2, axis = 0) #save average of the lists
standard_devs = list(np.std(midi_ks_2, axis = 0) ) 

#select fewer data points to plot to improve visualisation
Ns_fewer = Ns[100::100]
LOG_N_fewer_2 = LOG_N[100::100]
mean_ks_2 = mean_ks_imported[100::100]
standard_devs_fewer_2 = standard_devs[100::100]

repeats = 30
#Standard Error: standard deviation / sqrt(number of repeats)
errors_proper_2 = [i / np.sqrt(repeats) for i in standard_devs_fewer_2]

"""PLOT OF ALL K1s against N"""
plt.figure(4)
for i in range(len(midi_ks_2)):
    plt.plot(Ns, midi_ks_2[i], 'x-', markersize=2, linewidth=2)
plt.xlabel('N')
plt.ylabel(r'$k_{1}$')
plt.show()

"""PLOT OF ALL K1s against SQRT(N)"""
plt.figure(5)
for i in range(len(midi_ks_2)):
    plt.plot(LOG_N, midi_ks_2[i],  'x-', markersize=2, linewidth=2)
plt.xlabel('log(N)')
plt.ylabel(r'$k_{1}$')
plt.show()

"""Calculating the error on the best-fit gradient for K1 vs. SQRT(N)"""
data_upper = [mean_ks_2[i] + errors_proper_2[i] for i in range(len(mean_ks_2))]
P_1 = np.polyfit(LOG_N_fewer_2, data_upper, 1)
Poly_Info_1 = np.poly1d(P_1)
gradient_1 = round(Poly_Info_1[1], 4)

data_lower = [(mean_ks_2[i] - errors_proper_2[i]) for i in range(len(mean_ks_2))]
P_2 = np.polyfit(LOG_N_fewer_2, data_lower, 1)
Poly_Info_2 = np.poly1d(P_2)
gradient_2 = round(Poly_Info_2[1], 4)

error_in_grad = (gradient_2 - gradient_1)/2
print("Error in gradient", error_in_grad)

"""AVERAGE K1 AGAINST SQRT(N)""" #including the standard deviation on each k
plt.figure()
P, cov = np.polyfit(LOG_N_fewer_2, mean_ks_2, 1, cov=True)
Poly_Info = np.poly1d(P)
gradient_2 = round(Poly_Info[1], 4)
#plot data
plt.errorbar(LOG_N_fewer_2, mean_ks_2, yerr=errors_proper_2, fmt='o')
#plot linear fit
plt.plot(LOG_N_fewer_2, Poly_Info(LOG_N_fewer_2), label=f'Linear Fit Slope={gradient_2} $\pm$ {round(error_in_grad, 4)}')
#plot theory

y = [Largest_Degree_2(i, 2) for i in Ns_fewer]

P, cov = np.polyfit(LOG_N_fewer_2, y, 1, cov=True)
Poly_Info = np.poly1d(P)
gradient_2 = round(Poly_Info[1], 4)
#y_log = [np.log(i) for i in y]
plt.plot(LOG_N_fewer_2, y, '--', label=f'Theoretical Slope = {gradient_2}')

plt.xlabel('log(N)')
plt.ylabel(r'$\langle k_1 \rangle$')
plt.legend()
plt.show()

#%%
"""COLLECTING DATA FOR DATA COLLAPSE"""
#collect data
m = 2
N = [100, 1000, 10000, 100000, 1000000]
repeats = 100 

largest_k = [] #will store largest k for one run-through 
Degree_save_N_2 = [] #list of mean degrees

Deg_N_Append = Degree_save_N_2.append

for j in N:
    Degree_N = [] #for gathering lists of degrees for each N
    for i in range(repeats):
        G = Function_2(m, j) #evetually do 100,000 and also repeats  
        data_N = Degrees_G(G) #extracts the degrees 
        
        #save lists 
        Degree_N.append(data_N)

    print('DONE')
    deg = [val for sublist in Degree_N for val in sublist]
    np.savetxt('degrees_N_2{}.csv'.format(j), deg, delimiter=',') #save data
    Degree_save_N_2.append(deg)
#%%
"""PLOTTING THE DEGREE DATA FOR DIFFERENT N"""

N = [100, 1000, 10000, 100000, 1000000]

#load required data
Degree_save_N_2 = []
for i in range(len(N)):
    Degree_save_N_2.append(np.loadtxt(f'degrees_N_2{i}.csv', delimiter=','))
for i in range(1, len(Degree_save_N_2)):
    x, y = Logbin(list(Degree_save_N_2[i]))
    plt.plot(x,y, 'x-', label=f'N={N[i]}')
plt.xlabel('k')
plt.ylabel(r'$p_{N}$(k)')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()

"""PLOTTING THE DATA COLLAPSE SECTION 2"""
plt.figure(13)
plt.xlabel(r'k/log(N)')
plt.ylabel(r'$p_{N}$(k) / $p_{\infty}$(k)')

xx = []
yy = []
for i in range(len(N)):
    xx =[]
    yy = []
    #logbin data
    #print('full data before;', len(Degree_save_N[i]))
    x, y = Logbin(list(Degree_save_N_2[i]), 1)
    #print('data after logbin, x:', len(x))
    #plt.plot(x,y, 'x')
    #calculate theoretical expectation from logbinned data
    #x_theor = np.linspace(min(x), max(x), len(x))
    y_theor = [p_inft_theor_2(2, ll) for ll in x] #as m=2, x or x_theor? 
    for j in range(len(x)):
        #calculate manipulated data for data collapse
        #yy.append(y[j] / y_theor[j]) #Pn / Pinft
        yy.append(y[j] / y_theor[j])
        xx.append(x[j] / np.log(N[i])) #k/ k1
    
    plt.plot(xx, yy, '.', markersize=8, color = C[i], label=f'N={N[i]}')

plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""SECTION 3"""
    
def Function_NOK_3(m, Num_nodes, q):
    """
    Function for the BA model. Inputs are 'm': the number of edges we wish to 
    add to each new node, and 'Num_nodes': the total number of nodes we would 
    like to add to the system.
    Returns the graph 'G'.
    """
    G = [[]] #each sublist will contain the nodes the element is connected to
    #G_connections = [] #to keep track of which nodes have been connected, for prob
    #Largest_k = [] #store largest k values (k_1) per iteration
    
    #to increase efficiency
    #Append_Gconn = G_connections.append
    Append_node = G.append
    #Append_K = Largest_k.append
    
    #First, initialise. 
    for i in range(2*m): #so that 2m+1 nodes are added
        #add node
        numnodes = len(G) #store number of nodes so far
        Append_node([])
        for j in range(numnodes):
            G[i+1].append(j)
            G[j].append(i+1)
    
    
    walk_list = []
    #now add remaining nodes + edges between them
    for i in range(2*m + 1, Num_nodes): 
        WALK = Add_3(G, i, m, q) #add m edges to node i
        walk_list.append(WALK)
        

    return G, walk_list
    
def Add_3(G, Og, m, q):
    """
    Function which adds integer 'Number_Edges' to the node 'Og'.
    """
    Random = random.randint
    RandUnif = random.uniform
    RandChoice = random.choice
    
    maxx = len(G)-1
    G.append([]) #add one new node
    new_connections = [] #of max length m. nodes to which 'Og' is being connected
    walk_length_list = []
    while True:
        #randomly choose a node to connect edge with
        node_start = Random(0, maxx) 
        
        #decide whether to start random walk
        prob_accept = (1 - q)
        #choose a random number between 0 and 1
        random_number = RandUnif(0,1)
        d = 0 #counter for how many random-walk-steps are taken 
        while random_number < prob_accept: #loop if random number < p_acc 
            #randomly choose of the nodes node_start is connected to to move to
            node_next = RandChoice(G[node_start])
            d += 1 
            
            node_start = node_next
            random_number = RandUnif(0,1) #choose new random number
            
            #will keep looping until finds a node_start with correct prob

        node = node_start
        
        if node not in new_connections:
            #if this edge does not already exist, create the edge
            G[Og].append(node)
            G[node].append(Og)
            new_connections.append(node)
            #only save random-walk list if node is accepted
            walk_length_list.append(d)
        
        if len(new_connections) == m:
            #break out loop once required number of edges have been added
            break
            
    return walk_length_list

    
def l_theor_3(q):
    """
    Function to calculate the statistic probability distribution expected from
    theory. 
    """
    result = q / (1-q)
    return result

#%%    
"""COLLECTION DATA SECTION 3 for WALKIES"""
#M = [2, 4, 16] # 8, 16] #, 32, 64]
C = ['r', 'y', 'pink', 'g', 'b', 'c', 'k', 'magenta']
#m = 2
Degree_save_3 = []
Walks_save_3 = []

m = 2
repeats = 5
k = 100000
Q = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

#for q in Q:
for q in Q:
    Degree = []
    Walk_list = []
    for i in range(repeats):
        G, walks = Function_NOK_3(m, k, q)  
        data_3 = Degrees_G(G) #extracts the degrees 
        Degree.append(data_3) #append data to 'Degree' list
        Walk_list.append(walks)
        
    print('DONE')    
    deg = [val for sublist in Degree for val in sublist]
    np.savetxt('degrees_3{}.csv'.format(m), deg, delimiter=',') #save data
    Degree_save_3.append(deg)
    
    deg_w = [val for sublist in Walk_list for val in sublist]
    np.savetxt('RANDwalk_3{}.csv'.format(m), deg, delimiter=',')
    Walks_save_3.append(deg_w)
    
#%%
"""DATA COLECTION FOR P(K) VS K GRAPH"""
#M = [2, 4, 16] # 8, 16] #, 32, 64]
C = ['r', 'y', 'pink', 'k', 'magenta', 'g', 'b', 'c']

Degree_save_3 = []
Walks_save_3 = []

#repeats = 5
m = 8
k = 100000
Q = [1, 0.5, 0.1] #so q = 0 (RA), 0.5, 0.9 (PA)

#for q in Q:
for q in Q:
    Degree = []
    Walk_list = []
    G, walks = Function_NOK_3(m, k, q) 
    data_3 = Degrees_G(G) #extracts the degrees 
    Degree.append(data_3)  #append data to 'Degree' list
    Walk_list.append(walks)
    print('DONE')
   
    deg = [val for sublist in Degree for val in sublist]
    np.savetxt('degrees_3{}.csv'.format(m), deg, delimiter=',') #save data
    Degree_save_3.append(deg)
    
    deg_w = [val for sublist in Walk_list for val in sublist]
    np.savetxt('RANDwalk_3{}.csv'.format(m), deg, delimiter=',')
    Walks_save_3.append(deg_w)
#%%
"""p(k) vs k"""
#load data
#FOR2M = Degree_save_3
#FOR2M = np.loadtxt('degrees_3.csv', delimiter=',')

Q = [1, 0.5, 0.1]
C = ['r', 'b', 'g']
data_y_3 = []
theory_2 = []

plt.figure()
plt.xlabel('k')
plt.ylabel('p(k)')
x, y = Logbin(list(Degree_save_2[2]), 1.3)
#x = [i-1 for i in x]
plt.plot(x, y, label='PRA', color='orange')
x, y = Logbin(list(Degree_save_1[2]), 1.3)
#x = [i-1 for i in x]
plt.plot(x, y, label='PPA', color='c')

for i in range(len(Q)):
    #plot data
    x, y = Logbin(list(FOR2M[i]), 1.3)
    q = 1 - Q[i]
    plt.plot(x, y, 'x', color = C[-i], label= f'q={q}')
    y_theor = [p_inft_theor_2(M[i], j) for j in x]
    data_y_3.append(y)
    theory_2.append(y_theor) #for chi-squared comparison
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()

#%%
"""KS STATS GOODNESS OF FIT CODE: Prefferential Attachment"""
dataa = []
KS_vals = []
p_vals_KS = []
for i in range(len(M)):
    dataa.append(ks_2samp(data_y_2[i], theory_2[i]))
    KS_vals.append(dataa[i][0])
    p_vals_KS.append(dataa[i][1])
print('KS_vals:', KS_vals)
print('P_vals_ks:', p_vals_KS)

"""CHI SQUARED GOODNESS OF FIT CODE PART 3:  Random Attachment"""
dataa_3 = []
chi_vals_3 = []
p_vals_3 = []
for i in range(len(M)):
    Contingency_3 = []
    Contingency_3.append(data_y_3[i])
    Contingency_3.append(theory_2[i])
    dataa_3.append(chi2_contingency(Contingency_3))
    
    chi_vals_3.append(dataa_3[i][0])
    p_vals_3.append(dataa_3[i][1])

print('Chi_vals:', chi_vals_3)
print('P_vals:', p_vals_3)   
  

x_theor = np.linspace(min(x), max(x), len(x))
y_theor = [p_inft_theor_2(M[i], j) for j in x]

#%%
"""PLOTTING AVERAGE RANDOM WALK VS THEORY FOR VARIOUS q"""
#load data
#walkies = Walks_save_3

Q = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]  #this is p_acc
real_q = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0]

errors_walks = [np.std(walkies[i])/ np.sqrt(len(walkies[i])) for i in range(len(walkies))] 

plt.figure()
QQ= np.linspace(0, 0.9, 100)
theory = [l_theor_3(i) for i in QQ]
plt.plot(QQ, theory, '-', color='g', label='Theory')
for i in range(len(Q)-1):
    plt.errorbar(real_q[i], np.mean(walkies[i]), yerr=errors_walks[i], fmt='o', color='b')
    #plt.plot(Q[i], l_theor_3(Q[i]), 'o', color='o')
plt.plot(real_q[-1], np.mean(walkies[-1]), 'o', color='b', label='Data')
plt.xlabel('q')
plt.ylabel(r'$\langle d \rangle$')
#theory.insert(0,0)
plt.yscale('log')
plt.legend()