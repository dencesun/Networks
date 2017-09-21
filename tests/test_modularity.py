from CoreAlogorithm import Modularity as ml
import networkx as nx
from time import time

def main():
    g = nx.karate_club_graph()
    q, s = ml.modularity(g)
    print('modularity:  ', q)


if __name__ == '__main__':
    start = time()
    main()
    end = time()
    print('program runtime: ', str(end - start) + 's')