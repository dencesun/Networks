# -*- coding: utf-8 -*-
from CoreAlogorithm import *
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

    W = HigherOrderNetwork.motif_m7(g)
    cluster, condv, condc, order = HigherOrderNetwork.spectral_partitioning(W)
    print("\n\npaper figure1's result")
    print('condc: ', condc)
    print('cluster\n', cluster)


def paper_figure2():
    """

    graph data g come from Jure Leskovec(2016) Higher-order organization of complex networks Fig 2
    :return: None

    """
    # path for mac
    # data = '/Users/dencesun/Desktop/Networks/data/C-elegans-frontal.txt'
    # path for linux
    data = '/home/sun/Desktop/Networks-master/data/C-elegans-frontal.txt'
    # data = '../data/C-elegans-frontal.txt'
    DG = HigherOrderNetwork.create_network(data)
    A = HigherOrderNetwork.create_adjacency(DG)
    W = HigherOrderNetwork.motif_bifan(DG)
    W = HigherOrderNetwork.largest_connect_component(W)
    cluster, condv, condc, order = HigherOrderNetwork.spectral_partitioning(W)

    print("\n\npaper_figure2's result")
    print("largest_connect_component's shape: ", W.shape)
    print("C-elegans's result")
    print('condc: ', condc)
    print('cluster\n', cluster)
    # save as matlab file format
    # savemat('/Users/dencesun/Desktop/Networks/data/c-elegans.mat', {'W': W})
    # save path for linux
    savemat('/home/sun/Desktop/Networks-master/data/c-elegans.mat', {'W': W})
    print('complete save motif adjacency matrix in data')


def higher_order_network():
    paper_figure1()
    paper_figure2()


