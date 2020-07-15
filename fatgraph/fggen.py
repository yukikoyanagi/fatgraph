import itertools
from fatgraph import Fatgraph, FatgraphB

def generateall(*args, l=0):
    '''
    Generate a list of Fatgraphs with a given number of vertices and
    valences and l marked points. e.x. generateall(2,3,4) returns a 
    list of all Fatgraphs with three vertices, one 2-valent, 
    one 3-valent and one 4-valent with no marked points.
    Note the generated list includes disconnected graphs. Use
    Fatgraph.isconnected() to check for connectedness.
    '''
    vertices = []
    start = 1
    for i in args:
        end = start + i
        vertices.append(tuple(range(start, end)))
        start = end
    halfedges = [i for j in vertices for i in j]
    fatgraphs = []
    for exclude in itertools.combinations(halfedges, l):
        newedges = [i for i in halfedges if i not in exclude]
        alledges = [pairs for pairs in all_pairs(newedges)]
        fatgraphs.extend([Fatgraph(vertices, edges)
                          for edges in alledges])
    return fatgraphs

def generateallB(*args, l=0):
    '''
    Generate a list of valid FatgraphBs with a given number of 
    vertices and valences.
    The generated graphs are always connected.
    '''
    fatgraphBs =[]
    fatgraphs = generateall(*args, l=l)
    for graph in fatgraphs:
        for ie in itertools.product(
                *[paralleledges(vertex)
                  for vertex in graph.vertices]):
            int_edges = [p for ps in ie for p in ps]
            fg = FatgraphB.from_fatgraph(graph, int_edges)
            if fg.isvalid():
                fatgraphBs.append(fg)
    return fatgraphBs

def paralleledges(vertex):
    ids = [(i, -(i+1)) for i in range(int(len(vertex)/2))]
    for j in range(int(len(vertex)/2)):
        yield [(vertex[id[0]+j], vertex[id[1]+j]) for id in ids]
        
            

def all_pairs(lst):
    if not lst:
        yield []
        return
    if len(lst) % 2 == 1:
        # Handle odd length list
        for i in range(len(lst)):
            for result in all_pairs(lst[:i] + lst[i+1:]):
                yield result
    else:
        a = lst[0]
        for i in range(1,len(lst)):
            pair = (a,lst[i])
            for rest in all_pairs(lst[1:i]+lst[i+1:]):
                yield [pair] + rest
