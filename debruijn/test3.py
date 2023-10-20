import argparse
import os
import sys
import networkx as nx
from pathlib import Path
from networkx import DiGraph, all_simple_paths, lowest_common_ancestor, has_path, random_layout, draw, spring_layout
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
from typing import Iterator, Dict, List
#matplotlib.use("Agg")


def read_fastq(path):
    fastq_dict = {}
    num_seq = 0
    liste = []
    with open(path, "r") as fastq_file:
        seq_id = ""
        for line in fastq_file:
            line = line.strip()
            if line.startswith("@"):
                seq_id = line[1:].split()[0]
                fastq_dict[seq_id] = {"seq_ordre": "", "seq_val": ""}
            elif line.startswith("+"):
                continue
            elif "seq_ordre" in fastq_dict[seq_id] :  
                fastq_dict[seq_id]["seq_val"] += line.strip()
            else : 
                fastq_dict[seq_id]["seq_ordre"] += line.strip()
        for id, output in fastq_dict.items():
            num_seq += 1
            #print(f'ID de la sequence : {id}')
            #print(f'Séquence ({num_seq}) :  {output["seq_ordre"][:30]}')
            #print(f'{output["seq_val"][:30]}')
            #print(fastq_dict[seq_id]["seq_val"]) affiche toute les sequences 
            liste.append(output["seq_val"][:30])
    liste = str(liste)
    return liste
    
            


def divide(texte) : 
    liste_kmer = []
    reponse2 = input("entrez le nombre correspondant a la longueur des kmers : ")
    n = int(reponse2)
    for i in range(len(texte) - n +1): 
        kmer = texte[i:i+n]
        if kmer.startswith("'") : 
            continue
        elif kmer.endswith("'") : 
            continue
        elif kmer.startswith(","):
            continue
        elif kmer.endswith(","):
            continue
        elif kmer.startswith("["):
            continue
        elif kmer.endswith("]"):
            continue
        elif kmer.startswith(" "):
            continue
        elif kmer.endswith(" "):
            continue
        else : 
            liste_kmer.append(kmer)
    for j in liste_kmer : 
        print(j)
    #print(liste_kmer)
    return liste_kmer
    
    
def kmer_dict(path,N) : 
    fastq_dict = {}
    num_seq = 0
    liste = []
    dico_kmer = {}

    with open(path, "r") as fastq_file:
        seq_id = ""
        for line in fastq_file:
            line = line.strip()
            if line.startswith("@"):
                seq_id = line[1:].split()[0]
                fastq_dict[seq_id] = {"seq_ordre": "", "seq_val": ""}
            elif line.startswith("+"):
                continue
            elif "seq_ordre" in fastq_dict[seq_id] :  
                fastq_dict[seq_id]["seq_val"] += line.strip()
            else : 
                fastq_dict[seq_id]["seq_ordre"] += line.strip()
        for id, output in fastq_dict.items():
            num_seq += 1
            #print(f'ID de la sequence : {id}')
            #print(f'Séquence ({num_seq}) :  {output["seq_ordre"][:30]}')
            #print(f'{output["seq_val"][:30]}')
            #print(fastq_dict[seq_id]["seq_val"]) affiche toute les sequences 
            liste.append(output["seq_val"][:30])
    liste = str(liste)
    for i in range(len(liste)-1) : 
        mot = liste[i:i+N] #ici N va determiner la taille du kmer
        if mot not in A :
            continue
        else : 
            if mot in dico_kmer : 
                dico_kmer[mot] += 1 
            else : 
                dico_kmer[mot] = 1
    return dico_kmer
          
        
    
def build_graph(Dico) :  
 
    liste_edge = []
    Debruijn = nx.DiGraph()
    for kmer, count in Dico.items() : 
        debut, fin = kmer[:-1], kmer[1:]
        Debruijn.add_edge(debut,fin, weight = count)
        liste_edge.append(Debruijn)

    return Debruijn
    
def remove_paths(Digraph, path_list, delete_entry_node, delete_sink_node):

    for i in path_list : 
        if delete_entry_node and not delete_sink_node : 
            Digraph.remove_node_from(i[:-1])
        if delete_entry_node and delete_sink_node : 
            DiGraph.remove_node_from(i)
        if delete_sink_node and not delete_entry_node : 
            Digraph.remove_node_from(i[1:])
        if not delete_entry_node and not delete_sink_node : 
            Digraph.remove_node_from(i[1:-1])
    return Digraph

def select_best_paths(Digraph, path_list, path_length, weight_avg_list, delete_entry_node = False, delete_sink_node = False):
    
    max_weight_avg = max(weight_avg_list)
    max_length = max(path_length)
    to_remove = []
    
    if statistics.stdev(weight_avg_list) > 0 :
    
        for i in range(len(path_list)) :
            if weight_avg_list[i] < max_weight_avg :
                to_remove.append(path_list[i])
            else : 
                best_path_W = path_list[i]
            
    if statistics.stdev(path_length) > 0 :
        
        for i in range(len(path_list)):
            if path_length[i] < max_length : 
                to_remove.append(path_list[i])
        else : 
            best_path_L = path_list[i]
    
    remove = path_list[to_remove]
    Digraph = remove_paths(Digraph, remove, delete_entry_node, delete_sink_node)
    return Digraph
        
   
def path_average_weight(graph,path):

    return statistics.mean([d["weight"] for (u,v,d) in graph.subgraph(path).edges(data = True)])


def solve_bubble(graph,ancestor_node,descendant_node): 
    
    path_list = list(nx.all_simple_paths(graph, source = ancestor_node, target = descendant_node))
    path_length = []
    weight_avg_list = []
    
    for i in path_list : 
        path_length.append(len(i))
        weight_avg_list.append(path_average_weight(graph,i))
    graph = select_best_paths(Digraph, path_list, path_length, weight_avg_list)
    return graph 
    
    
def simplify_bubble(graph) : 

    bubble = False  
    
    for i in graph.nodes() : 
        predecessor = [p for p in graph.predecessors(i)]
        if len(list(predecessor)) > 1 : 
            for combi in combinations(predecessor,2) : 
                ancestor = nx.lowest_common_ancestor(graph,combi[0],combi[1])
                if ancestor is not None : 
                    bubble = True
                    break
        if bubble == True : 
            break
    if bubble : 
        graph = simplify_bubble(solve_bubble(graph,ancestor,i))
    return graph 
    
def solve_entry_tips(graph, starting_nodes):

    list_path = []
    weight_path = []
    len_path = []
    for node in starting_nodes : 
        for descendant in nx.descendants(graph, node) : 
            predecessor = list(graph.predecessors(descendant))
            if len(predecessor) > 1:
                for path in nx.all_simple_paths(kmer_graph, node, descendant) : 
                    list_path.append(path)
                    len_path.append(len(path))
                    weight_path.append(path_average_weight(graph, path))
    graph = select_best_paths(graph, list_path, len_path, weight_path, delete_entry_node = True, delete_sink_node = False)
    return graph

def solve_out_tips(graph,ending_nodes):
    
    path_list = []
    path_length = []
    weight_avg_list = []
    for node in graph.nodes() : 
        list_successor = list(graph.succesors(node))
            if len(list_successor) > 1 :
                for end in ending_nodes : 
                    for path in nx.all_simple_paths(graph,node;end)
                        path_list.append(path)
                    
                    
        if len(path_list) != 0 : 
            for path in path_list : 
                path_length.append(len(path))
                weight_avg_list.append(path_average_weight(graph,path))
            graph = select_best_paths(graph, path_list, path_length, weight_avg_list, delete_entry_node = False, delete_sink_node = True)
        
        return graph
    
    
    
def get_starting_nodes(graph):
    
    starting_nodes = []
    
    for node in graph.nodes() : 
        if not any (graph.predecessors(node)): #????
            starting_nodes.append(node)
            
    return starting_nodes
    
def get_sink_nodes(graph) : 
    
    ending_nodes = []
    
    for node in graph.nodes() : 
        if not any (graph.succesors(node)) : 
            ending_nodes.append(node)
            
    return ending_nodes
    
    
def get_contigs(graph, starting_nodes, ending_nodes) : 
    contigs_list = []
    for input_node in starting_nodes : 
        for output_node in ending_nodes : 
            for path in nx.all_simple_paths(graph, source = input_node, target = output_node) : 
                contig = path[0]
                for comp in path[1:] :
                    contig += comp[1:]
                contigs_list.append((contig,len(contig)))
    return contigs_list

def save_contigs(contigs_list, output_file) : 

    with open(output_file, "w") as file : 
        count = 0 
        for i in contigs_list : 
            file.write(f"contig number {count} len = {i[1]}\n") 
            file.write(f"{i[0]}")
            count += 1
    file.close()
    
    
#A CONTINUER A PARTIR DES AUTRES GIT    
    
    
    
    



 


#main programme

#1
path = "eva71_plus_perfect.fq"
read_fastq(path)

#2
texte = read_fastq(path)
divide(texte)

#3
reponse3 = input("entrez le nombre correspondant a la longueur des kmers : ")
N = int(reponse3)
A = divide(texte)
print(kmer_dict(path,N)) 

#remarque : demande 3 fois la valeur de N au lieu de 2 ://

#4
Dico = kmer_dict(path,N)
graph = build_graph(Dico)
pos = nx.spring_layout(graph)
nx.draw(graph,pos,with_labels = True, node_size = 1000, font_size = 10, font_color = 'white')
edge_labels = {edge : graph[edge[0]][edge[1]]['weight'] for edge in graph.edges}
nx.draw_networkx_edge_labels(graph,pos,edge_labels=edge_labels)
plt.show()
print(build_graph(Dico))

#5


