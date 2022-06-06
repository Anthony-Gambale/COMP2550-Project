import networkx as nx
import numpy as np

def max_cut(graph):
    bestCutSize = -1
    bestSets = []
    bestBitStrings = []
    for x in range(1, 2 << (graph.number_of_nodes()-1), 2): #Step by two so that last digit is always the same (to handle symmetry) .ljust(graph.number_of_nodes(), "1")
        setOfNodes = binaryToSet(x)
        cutSize = nx.cut_size(graph, setOfNodes)
        if(cutSize>bestCutSize):
            bestSets = [setOfNodes]
            bestCutSize = cutSize
            bestBitStrings = [str(bin(notBin(x, graph.number_of_nodes())))[:1:-1].ljust(graph.number_of_nodes(), "0")]
        elif(cutSize==bestCutSize):
            bestSets.append(setOfNodes)
            bestBitStrings.append(str(bin(notBin(x, graph.number_of_nodes())))[:1:-1].ljust(graph.number_of_nodes(), "0"))
    return bestCutSize, bestSets, bestBitStrings

def notBin (x, length):
    return (2 << (length-1)) - 1 - x

def binaryToSet(binInt):
        rtn = set()
        for x in range(int(np.log2(binInt))+1):
            if ((binInt >> x) & 1) == 1:
                rtn.add(x)
        return rtn
