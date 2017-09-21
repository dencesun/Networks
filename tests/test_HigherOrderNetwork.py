# -*- coding: utf-8 -*-
from CoreAlogorithm import HigerOrderNetwork as hon
from time import time
import networkx as nx
from scipy.io import savemat

def paper_figure1():
    """

    graph data g come from Jure Leskovec(2016) Higher-order organization of complex networks Fig 1
    :return: None

    """
    g = nx.DiGraph()
    g.add_edge(1, 2)
    g.add_edge(1, 3)
    g.add_edge(1, 4)
    g.add_edge(1, 5)
    g.add_edge(1, 6)
    g.add_edge(2, 1)
    g.add_edge(2, 3)
    g.add_edge(2, 4)
    g.add_edge(2, 5)
    g.add_edge(6, 2)
    g.add_edge(6, 7)
    g.add_edge(7, 2)
    g.add_edge(7, 6)
    g.add_edge(8, 2)
    g.add_edge(8, 6)
    g.add_edge(8, 9)
    g.add_edge(8, 10)
    g.add_edge(9, 6)
    g.add_edge(9, 7)
    g.add_edge(9, 8)
    g.add_edge(9, 10)

    W = hon.motif_m7(g)
    cluster, condv, condc, order = hon.spectral_partitioning(W)
    print("paper figure1's result\n")
    print('condc: ', condc)
    print('cluser: ', cluster)


def paper_figure2():
    """

    graph data g come from Jure Leskovec(2016) Higher-order organization of complex networks Fig 2
    :return: None

    """

    data = '..\data\C-elegans-frontal.txt'
    DG = hon.create_network(data)
    A = hon.create_adjacency(DG)
    W = hon.motif_bifan(DG)
    W = hon.largest_connect_component(W)
    cluster, condv, condc, order = hon.spectral_partitioning(W)

    print("largest_connect_component's shape: ", W.shape)
    print("C-elegans's result")
    print('condc: ', condc)

    # save as matlab file format
    savemat('..\data\c-elegans.mat', {'W': W})
    print('complete save motif adjacency matrix in data')


def main():
    paper_figure1()
    paper_figure2()

if __name__ == '__main__':
    start = time()
    main()
    end = time()
    print('program runtime = ', str(end - start) + 's')
