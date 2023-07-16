# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 21:21:27 2020

@author: argdi
"""

import networkx as nx
import random
import numpy as np
import sys

class SIR():
    def __init__(self, alpha, beta, tI, k, k0, S, I, X, pWatts, nei, maxtime):
        self.maxtime = maxtime
        #initialize the parameter values
        self.tI = tI
        self.beta = beta
        self.k = k
        self.k0 = k0
        self.N = S + I + X       
        self.alpha = alpha
        
        #define the underlying structure, the theoretical network model that will be used for the simulation
        
        #self.graph = nx.watts_strogatz_graph(self.N, nei, pWatts) 
        #p = nei/(self.N - 1.) #probability for link creations. It will be used for random networks
        #self.graph = nx.fast_gnp_random_graph(self.N, p)
        #self.graph = nx.gaussian_random_partition_graph(N, 100, 100, 10.*p, 1.5*p)
        #sizes = list(np.zeros(100, dtype=int) + 100)
        #self.graph = nx.random_partition_graph(sizes, 2*p, 2.5*p)
        #self.graph = nx.barabasi_albert_graph(N, int(nei/2))
        self.graph = nx.dual_barabasi_albert_graph(self.N, 1, 2, 0.5)
        self.neighbours = [] #neighbors of each node
        
        nodes = np.array(nx.nodes(self.graph)) #nodes
        
        avgDegree = 0
        for i in nodes:
            self.neighbours.append(list(self.graph.neighbors(i)))
           
        
        self.tau = alpha/nei 
        #nodes of the compartments of the model
        self.S_node = nodes #initially all nodes are S
        self.I_node = np.array([], int)
        self.R_node = np.array([], int)
        self.X_node = np.array([], int)
        self.tIarray = np.zeros(self.N) + 2*self.maxtime
        
        #the initial X nodes (confirmed cases)
        for i in range(X):
            xTemp = random.choice(self.S_node)
            xIndTpl = np.where(self.S_node==xTemp)
            if len(xIndTpl) != 1:
               sys.exit("exit")
            xInd = xIndTpl[0][0]
            self.S_node = np.delete(self.S_node, xInd)
            self.X_node = np.append(self.X_node, xTemp)
            
        #the initial infected nodes
        for i in range(I):
            self.sick = random.choice(self.S_node)
            sickTmp = np.where(self.S_node==self.sick)
            if len(sickTmp) != 1:
               sys.exit("exit")
            sickInd = sickTmp[0][0]
            
            self.S_node = np.delete(self.S_node, sickInd) #infected node
            self.I_node = np.append(self.I_node, self.sick)       
            self.tIarray[self.sick] = self.tI
         
        
       
        self.Infected = np.array([])
        self.Infected = np.append(self.Infected, I)
        self.t = 1 #the initial process is at t = 1
        
        #time series of the compartments
        self.susceptible = np.array([])
        self.recovered = np.array([])
        self.susceptible = np.append(self.susceptible, len(self.S_node))
        self.recovered = np.append(self.recovered, 0)
        self.newInf = np.array([len(self.I_node)])
        self.newS = np.zeros(maxtime)
        self.X = np.array([len(self.X_node)])
     
    
    def run(self):
        #runs the stochastic Monte Carlo simulation of SIRX
        self.t = 2 #the MC simulation starts at t = 2       
        e = 0.2
        while self.t<=self.maxtime:
            #each timestep we initially check all the infection events 
            newI = np.array([], int)
            temp = 0
            #loop through each S - I link
            for node in self.I_node:
                for j in self.neighbours[node]: 
                    #if a random number 0-1 is smaller than transmissibility 
                    #we move the node it from S to a temporary list newI                    
                    if (random.random() < self.tau): 
                        if j in self.S_node and j not in newI:                            
                            newI = np.append(newI, j) #it is infected now
                            sickTpl = np.where(self.S_node==j)
                            if len(sickTpl) != 1:
                                sys.exit("error")
                            sickInd = sickTpl[0][0]
                            self.S_node = np.delete(self.S_node, sickInd) 
                            
                            temp += 1

            #the order of the next events is random
            randomArray = [1, 2, 3, 4]
            random.shuffle(randomArray)
                       
           
            for randNum in randomArray:
                if abs(randNum - 1) < e:
                    #recovery
                    for node in self.I_node:
                        if random.random() < self.beta : 
                            
                            recover_tpl = np.where(self.I_node==node)
                            if len(recover_tpl) != 1:
                                sys.exit("exit")
                            recover_ind = recover_tpl[0][0]
                            self.I_node = np.delete(self.I_node, recover_ind)
                            self.R_node = np.append(self.R_node, node)  
    
                #I to X events
                if abs(randNum - 2) < e:     
                    for node in self.I_node:
                       
                        if random.random() < self.k0:                        
                            Ik0tpl = np.where(self.I_node==node)
                            if len(Ik0tpl) != 1:
                                sys.exit("exit")
                            Ik0index = Ik0tpl[0][0]
                            self.I_node = np.delete(self.I_node, Ik0index)
                            self.X_node = np.append(self.X_node, node)
                
                if abs(randNum - 3) < e:
                    for node in self.I_node:
                        #if len(self.I_node)>0:
                        if random.random() < self.k:                        
                            Iktpl = np.where(self.I_node==node)
                            if len(Iktpl) != 1:
                                sys.exit("exit")
                            Ikindex = Iktpl[0][0]
                            self.I_node = np.delete(self.I_node, Ikindex)
                            self.X_node = np.append(self.X_node, node)
                #S to R events
                if abs(randNum - 4) < e:
                  for node in self.S_node:
                    if random.random() < self.k0:
                        Sk0tpl = np.where(self.S_node==node)
                        if len(Sk0tpl) != 1:
                            sys.exit("exit")
                        Sk0index = Sk0tpl[0][0]
                        self.S_node = np.delete(self.S_node, Sk0index)
                        self.R_node = np.append(self.R_node, node)       
            
            #renew the lists of time series
            self.I_node = np.append(self.I_node, newI)
            self.newInf = np.append(self.newInf, len(newI))  
            self.Infected = np.append(self.Infected, len(self.I_node))
            self.susceptible = np.append(self.susceptible, len(self.S_node))
            self.recovered = np.append(self.recovered, len(self.R_node))
            self.X = np.append(self.X, len(self.X_node))
            
            self.t += 1
            
            random.shuffle(self.I_node)

            
        return [self.Infected, self.susceptible, self.recovered, self.newInf, self.X]
    

