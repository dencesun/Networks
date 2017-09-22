from CoreAlogorithm import Modularity as ml
import networkx as nx
from time import time


def test_modularity():
    g = nx.karate_club_graph()
    q, s = ml.modularity(g)
    print('modularity:  ', q)

