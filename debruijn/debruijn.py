#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
from itertools import combinations
import os
import sys
from pathlib import Path
import networkx as nx
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
#matplotlib.use("Agg") affichage du graphe ne marche pas avec cette commande initiee

__author__ = "HOLO Amon-Donovan"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["HOLO Amon-Donovan"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "HOLO Amon-Donovan"
__email__ = "donovanholo@gmail.com"
__status__ = "Developpement"

def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=Path,
                        default=Path(os.curdir + os.sep + "contigs.fasta"),
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=Path,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(Path) :
    """Extract reads from fastq files.

    :param fastq_file: (Path) Path to the fastq file.
    :return: A generator object that iterate the read sequences. 
    """
    fastq_dict = {}
    num_seq = 0
    liste = []
    with open(Path, "r") as fastq_file:
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
            print(f'ID de la sequence : {id}')
            print(f'Séquence ({num_seq}) :  {output["seq_ordre"][:30]}')
            print(f'{output["seq_val"][:30]}')
            print(fastq_dict[seq_id]["seq_val"]) #affiche toute les sequences 
            liste.append(output["seq_val"])
    liste = str(liste)
    return liste



def cut_kmer(read, kmer_size) : #read = texte
    """Cut read into kmers of size kmer_size.
    
    :param read: (str) Sequence of a read.
    :return: A generator object that provides the kmers (str) of size kmer_size.
    """
    liste_kmer = []
    reponse2 = input("entrez le nombre correspondant a la longueur des kmers : ") #kmer_size
    n = int(reponse2) #kmer_size
    for i in range(len(read) - n +1): 
        kmer = read[i:i+n]
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

    return liste_kmer
    


def build_kmer_dict(Path, kmer_size) :
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """

    fastq_dict = {}
    num_seq = 0
    liste = []
    dico_kmer = {}

    with open(Path, "r") as fastq_file:
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
            liste.append(output["seq_val"][:30])
    liste = str(liste)
    for i in range(len(liste)-1) : 
        mot = liste[i:i+kmer_size] #ici N va determiner la taille du kmer
        if mot not in list_verif :
            continue
        else : 
            if mot in dico_kmer : 
                dico_kmer[mot] += 1 
            else : 
                dico_kmer[mot] = 1
    return dico_kmer


def build_graph(kmer_dict):
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    liste_edge = []
    graph = nx.DiGraph()
    for kmer, count in kmer_dict.items() : 
        debut, fin = kmer[:-1], kmer[1:]
        graph.add_edge(debut,fin, weight = count)
        liste_edge.append(graph)

    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node) :
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    for i in path_list : 
        if delete_entry_node and not delete_sink_node : 
            graph.remove_node_from(i[:-1])
        if delete_entry_node and delete_sink_node : 
            graph.remove_node_from(i)
        if delete_sink_node and not delete_entry_node : 
            graph.remove_node_from(i[1:])
        if not delete_entry_node and not delete_sink_node : 
            graph.remove_node_from(i[1:-1])
    return graph


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False) :
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length_list: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
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
    graph = remove_paths(graph, remove, delete_entry_node, delete_sink_node)
    return graph


def path_average_weight(graph, path) :
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])


def solve_bubble(graph, ancestor_node, descendant_node):
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph 
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    path_list = list(nx.all_simple_paths(graph, source = ancestor_node, target = descendant_node))
    path_length = []
    weight_avg_list = []
    
    for i in path_list : 
        path_length.append(len(i))
        weight_avg_list.append(path_average_weight(graph,i))
    graph = select_best_path(graph, path_list, path_length, weight_avg_list)

    return graph 


def simplify_bubbles(graph) -> DiGraph:
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
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
        graph = simplify_bubbles(solve_bubble(graph,ancestor,i))
    return graph


def solve_entry_tips(graph, starting_nodes) :
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :param starting_nodes: (list) A list of starting nodes
    :return: (nx.DiGraph) A directed graph object
    """
    list_path = []
    weight_path = []
    len_path = []
    for node in starting_nodes : 
        for descendant in nx.descendants(graph, node) : 
            predecessor = list(graph.predecessors(descendant))
            if len(predecessor) > 1:
                for path in nx.all_simple_paths(graph, node, descendant) : 
                    list_path.append(path)
                    len_path.append(len(path))
                    weight_path.append(path_average_weight(graph, path))
    graph = select_best_path(graph, list_path, len_path, weight_path, delete_entry_node = True, delete_sink_node = False)
    return graph


def solve_out_tips(graph, ending_nodes) :
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :param ending_nodes: (list) A list of ending nodes
    :return: (nx.DiGraph) A directed graph object
    """
    path_list = []
    path_length = []
    weight_avg_list = []
    for node in graph.nodes() : 
        list_successor = list(graph.succesors(node))
        if len(path_list) != 0 : 
            for path in path_list : 
                path_length.append(len(path))
                weight_avg_list.append(path_average_weight(graph,path))
            if len(list_successor) > 1 :
                for end in ending_nodes : 
                    for path in nx.all_simple_paths(graph,node,end) : 
                        path_list.append(path)
                    
                    
        
            graph = select_best_path(graph, path_list, path_length, weight_avg_list, delete_entry_node = False, delete_sink_node = True)
        
        return graph


def get_starting_nodes(graph):
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    starting_nodes = []
    
    for node in graph.nodes() : 
        if not any (graph.successors(node)): 
            starting_nodes.append(node)
            
    return starting_nodes


def get_sink_nodes(graph):
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    ending_nodes = []
    
    for node in graph.nodes() : 
        if not any (graph.succesors(node)) : 
            ending_nodes.append(node)
            
    return ending_nodes


def get_contigs(graph, starting_nodes, ending_nodes):
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object 
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
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
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, "w") as file : 
        count = 0 
        for i in contigs_list : 
            file.write(f"contig number {count} len = {i[1]}\n") 
            file.write(f"{i[0]}")
            count += 1
    file.close()


def draw_graph(graph: DiGraph, graphimg_file: Path) -> None: # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (Path) Path to the output file
    """                                   
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    print(elarge)
    # Draw the graph with networkx
    pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file.resolve())


#==============================================================
# Main program
#==============================================================
def main() :
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    Path = "eva71_plus_perfect.fq"

    read = read_fastq(Path)

    selection = input("entrez le nombre correspondant a la taille des kmers : ")
    kmer_size = int(selection)

    list_verif = cut_kmer(read,kmer_size)

    kmer_dict = build_kmer_dict(Path, kmer_size)
    graph = build_graph(kmer_dict)
    pos = nx.spring_layout(graph)
    nx.draw(graph,pos,with_labels = True, node_size = 1000, font_size = 10, font_color = 'white')
    edge_labels = {edge : graph[edge[0]][edge[1]]['weight'] for edge in graph.edges}
    nx.draw_networkx_edge_labels(graph,pos,edge_labels=edge_labels)
    plt.show()

    starting_nodes = get_starting_nodes(graph)
    graph = solve_entry_tips(graph, starting_nodes)

    ending_nodes = get_sink_nodes(graph)
    graph = solve_out_tips(graph, ending_nodes)


    graph = simplify_bubbles(graph)
    contigs_list = get_contigs(graph, starting_nodes, ending_nodes)

    save_contigs(contigs_list, args.output_file)





    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)


if __name__ == '__main__': # pragma: no cover
    main()
