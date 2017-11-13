from CoreAlogorithm import Modularity as ml
import networkx as nx


def test_modularity():
    g = nx.karate_club_graph()
    q, s = ml.modularity(g)
    print('modularity:  ', q)
    ans = list()
    for i in range(len(s)):
        if s[i] == 1:
            ans.append(i+1)

    print(ans)
    print('community size: ', len(ans))

