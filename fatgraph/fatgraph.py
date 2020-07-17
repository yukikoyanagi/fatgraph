import re
from itertools import accumulate
from operator import mul
from permutation import Permutation

import numpy as np

DONOR_COL = 12
ACCPT_COL = 13
FLAGS_COL = 21


class Fatgraph(object):
    """
    A fatgraph object representing a fatgraph. A fatgraph is stored
    as a pair of permutations, one representing vertices and the
    other edges. For example v=[(1,2,3,4)] and e=[(1,2),(3,4)] is
    a figure-of-eight fargraph on a single 4-valent vertex.
    Constructor takes two arguments, each an iterable of iterables,
    for example vertices=((1,2,3),(4,5,6)), edges=((1,2),(3,4)(5,6))
    """

    def __init__(self, vertices, edges):
        # Check the arguments define a valid fatgraph
        self.checkvalidity(vertices, edges)
        self.vertices = vertices
        self.edges = edges
        self._vs = Permutation.from_cycles(*vertices)
        self._es = Permutation.from_cycles(*edges)
        vs = [i for j in vertices for i in j]
        es = [i for j in edges for i in j]
        self.unpaired = set(vs) - set(es)

    def __repr__(self):
        return ('{0.__module__}.{0.__name__}(\n'
                'vertices={1}\n'
                'edges={2})').format(type(self),
                                     self.vertices,
                                     self.edges)

    def __eq__(self, other):
        return all([self._vs == other._vs,
                    self._es == other._es,
                    self.unpaired == other.unpaired])

    def __ne__(self, other):
        return not self == other

    @property
    def boundaries(self):
        bs = self._vs * self._es
        cycles = bs.to_cycles()
        # to_cycles omits one-vertex boundaries as identity, so
        # need to add them.
        fixed = set(i for v in self._vs.to_cycles() for i in v) - \
            set(j for c in cycles for j in c)
        for f in fixed:
            cycles.append((f,))
        return cycles

    @property
    def genus(self):
        if len(self._vs.to_cycles()) == 0:
            return 0
        else:
            return int((2 - len(self._vs.to_cycles())
                        + len(self._es.to_cycles())
                        - len(self.boundaries)) / 2)

    @classmethod
    def from_hbonds(cls, hbfile, bbtype='alpha'):
        """
        Create fatgraph object from PDB Hbond file. The resulting
        fatgraph is simplified; that is, non-bonded half-edges are
        not recorded. Note this does not change genus.
        TODO: return non-simplified fatgraph
        """
        bonds = []
        with open(hbfile, 'r') as fh:
            for k, line in enumerate(fh):
                cols = line.split()
                # we only consider intra-chain bonds (indicated by __),
                # and the last two characters must be either U (unique)
                # or S (strong)
                if not re.search("__[US][US]$", cols[FLAGS_COL]):
                    continue
                don = int(cols[DONOR_COL])
                acc = int(cols[ACCPT_COL])
                bonds.append((don, acc))

        try:
            dons, accs = zip(*bonds)
        except ValueError as e:
            if len(bonds) == 0:
                return cls([()], [()])
            else:
                raise e

        if bbtype == 'alpha':
            # We assume no backbone twist along the entire protein
            dons = sorted(dons)
            accs = sorted(accs, reverse=True)
            dons_dic = {v: k+1 for k, v in enumerate(dons)}
            dons_max = dons_dic[dons[-1]]
            accs_dic = {v: k+1+dons_max for k, v in enumerate(accs)}
        elif bbtype == 'beta':
            # Assume backbone twist for every c-alpha linkage
            upper, lower = [], []
            for d in dons:
                if d % 2 == 1:
                    upper.append(d)
                else:
                    lower.append(d)
            for a in accs:
                if a % 2 == 1:
                    upper.append(a + 1)
                else:
                    lower.append(a + 1)
            upper = sorted(upper)
            lower = sorted(lower, reverse=True)
            upper_dic = {v: k+1 for k, v in enumerate(upper)}
            upper_max = upper_dic[upper[-1]]
            lower_dic = {v: k+1+upper_max for k, v in enumerate(lower)}
            dons_dic, accs_dic = {}, {}
            for d in dons:
                if d % 2 == 1:
                    dons_dic[d] = upper_dic[d]
                else:
                    dons_dic[d] = lower_dic[d]
            for a in accs:
                if a % 2 ==1:
                    accs_dic[a] = upper_dic[a + 1]
                else:
                    accs_dic[a] = lower_dic[a + 1]
            

        edges = []
        for bond in bonds:
            edges.append((dons_dic[bond[0]], accs_dic[bond[1]]))

        vertices = [list(range(1, len(dons_dic)+len(accs_dic)+1))]

        return cls(vertices, edges)

    @classmethod
    def checkvalidity(self, vertices, edges):
        # Check supplied vertices and edges do not contain duplicates
        vs = [i for j in vertices for i in j]
        if len(vs) != len(set(vs)):
            raise ValueError('Vertices are not unique')
        es = [i for j in vertices for i in j]
        if len(es) != len(set(es)):
            raise ValueError('Edges are not unique')
        # Check all half-edges are present in vertices
        if not set(es) <= set(vs):
            raise ValueError('Edges do not connect to vertices')

    def isconnected(self):
        '''
        Return True if this graph is connected. False otherwise
        '''
        remaining = self.vertices.copy()
        visited = {}
        visit = [remaining[0], ]
        remaining.remove(remaining[0])
        while remaining:
            # get all vertices connected to the first vertex in
            # remaining list
            tovisit = []
            for v in visit:
                for h in v:
                    edge = [e for e in self.edges if h in e]
                    if not edge:
                        # h was an unpaired half-edge.
                        continue
                    edge = edge[0]
                    if edge[0] == h:
                        other = edge[1]
                    else:
                        other = edge[0]
                    tv = [r for r in remaining if other in r]
                    if len(tv) > 1:
                        raise ValueError('Half-edges are not unique')
                    tovisit.extend(tv)
            if not tovisit:
                # We've gone through all vertices in visit list,
                # and there is no more to visit even though remaining
                # list is not empty.
                return False
            visit = list(set(tovisit)) # de-duplicate
            for v in visit:
                remaining.remove(v)

        # Nothing remaining. It is connected.
        return True


class FatgraphB(Fatgraph):
    '''
    Fatgraph class for beta vertex. This class, in addition to the
    Fatgraph class, contains information about vertex interior; which
    edge is conneted to which *inside* the vertices. 
    Note the beta vertices we consider has two unpaired half-edges.
    We assume the half-edge 1 and one connected to it to be the two
    unpaired half-edges.
    '''

    def __init__(self, vertices, edges, interiors):
        self._is = Permutation.from_cycles(*interiors)
        super().__init__(vertices, edges)

    def __repr__(self):
        return '\n'.join([super().__repr__(),
                          'interiors={}'.format(self.interiors)])

    def __eq__(self, other):
        return all([super().__eq__(other),
                    self.interiors == other.interiors])

    def __ne__(self, other):
        return not self == other

    @property
    def interiors(self):
        return self._is.to_cycles()

    @classmethod
    def from_fatgraph(cls, fatgraph, interiors):
        assert isinstance(fatgraph, Fatgraph)
        return cls(fatgraph.vertices,
                   fatgraph.edges,
                   interiors)

    @classmethod
    def from_pmat(cls, mat):
        "Make fatgraph from pairing matrix"
        '''
        [[0,0,0],     vertices: [(1,2,3,4,5,6),]
         [1,0,0], --> edges: [(2,3), (5,6)]
         [0,1,0]]     i_edges: [(1,6), (2,5), (3,4)]
        '''

        def findsheetat(k, mat):
            "Find sheet that starts with k'th strand"
            sheet = [k]
            row = mat[k]
            if np.count_nonzero(row) > 1:
                raise ValueError('Sheet must start at an edge strand.')
            elif np.count_nonzero(row) < 1:
                return sheet
            prevk = k
            thisk = np.flatnonzero(row)[0]
            thisrow = mat[thisk]
            while np.count_nonzero(thisrow) != 1:
                if thisk in sheet:
                    raise ValueError('Pairing matrix has barrel(s).')
                sheet.append(thisk)
                i, j = np.flatnonzero(thisrow)
                if i == prevk:
                    nextk = j
                else:
                    nextk = i
                prevk, thisk = thisk, nextk
                thisrow = mat[thisk]
            else:
                if thisk in sheet:
                    raise ValueError('Pairing matrix has barrel(s).')
                sheet.append(thisk)
            return sheet

        def findsheets(p_mat):
            "Sheets from (symmetric) pairting matrix"
            '''
            findsheets([[0,0,1,0],
                        [0,0,1,0],
                        [1,1,0,0],
                        [0,0,0,0]])
            --> {(0,2,1),(3,)}
            '''
            wmat = np.copy(p_mat)
            edges = [i for i in range(len(wmat))
                     if np.count_nonzero(wmat[i])<2]
            if not edges:
                raise ValueError('Pairing matrix has no edge strand.')
            sheets = set()
            while edges:
                e = edges.pop(0)
                sheet = tuple(findsheetat(e, wmat))
                if sheet[0] > sheet[-1]:
                    sheet = sheet[::-1]
                sheets.add(sheet)
                if sheet[-1] in edges:
                    edges.remove(sheet[-1])
            return sheets

        def makevertices(sheets, p_mat):
            "Vertices from sheets and pairing matrix."
            '''
            [[0,1,2]], [[0, 0, 0],
                        [1, 0, 1], --> [[(0,1),(11,10),(21,20)]]
                        [0, 0, 0]]
            '''
            vertices = []
            for sheet in sheets:
                oseq = []
                for i, j in zip(sheet, sheet[1:]):
                    if i > j:
                        i, j = j, i
                    if p_mat[i,j]:
                        oseq.append(1)
                    else:
                        oseq.append(-1)
                    oseq = list(accumulate(oseq, mul))
                vertex = []
                i = sheet[0]
                vertex.append((i*10, i*10+1))
                for i, s in enumerate(sheet[1:]):
                    if oseq[i]>0:
                        vertex.append((s*10, s*10+1))
                    else:
                        vertex.append((s*10+1, s*10))
                first = min(vertex)
                if first[0] > first[1]:
                    vertex = [(s[1],s[0]) for s in vertex]
                vertices.append(vertex)
            return vertices

        symmat = mat + mat.T
        if np.where(np.sum(symmat, axis=1)>2)[0].size:
            raise ValueError('Pairing matrix has bifurcation(s).')
        sheets = findsheets(symmat)
        if len(mat) != len([i for sheet in sheets for i in sheet]):
            raise ValueError('Pairing matrix has barrel(s).')
        vertices = makevertices(sheets, mat)
        vdict = {}
        for vertex in vertices:
            #We want anti-clockwise ordering of half-edges on each vertex
            l = [v[0] for v in vertex] + [v[1] for v in vertex][::-1]
            r = len(vdict)
            for i in range(r, r+len(l)):
                vdict[l[i - r]] = i + 1

        #Now we find edges
        halfedges = sorted([i for vertex in vertices
                            for v in vertex
                            for i in v])
        edges = [(i, j) for i, j in zip(halfedges, halfedges[1:])]
        #edges = [(halfedges[i], halfedges[i+1])
        #         for i in range(1, len(halfedges)-2, 2)]

        #Translate vertices and edges according to vdict
        v = []
        p = 1
        for vertex in vertices:
            v.append(tuple(range(p, 2*len(vertex)+p)))
            p += 2*len(vertex)

        iv = []
        iedges = [e for vertex in vertices for e in vertex]
        for e in iedges:
            iv.append((vdict[e[0]], vdict[e[1]]))

        e = []
        for d in edges:
            edge = (vdict[d[0]], vdict[d[1]])
            if edge[0] > edge[1]:
                edge = edge[::-1]
            if edge not in iv:
                e.append(edge)

        return cls(v, e, iv)

    def isvalid(self):
        '''
        Return true if the graph has no path smaller than the maximal
        (i.e. all edges are traversed), when following the exterior
        and interior edges.
        '''
        h_edges = [i for j in self.vertices for i in j]
        edges = self.interiors + self.edges
        marked = set(h_edges) - set([i for e in self.edges for i in e])
        if len(marked)//2 == 0: # no marked points
            last = 1
        else:
            last = list(marked)[0]
        while edges:
            e = [l for l in edges if last in l]
            if e:
                edge = e[0]
                edges.remove(edge)
                if edge[0] == last:
                    last = edge[1]
                else:
                    last = edge[0]
            else:
                # No match for last, we have a loop or a path
                return False
        # We emptied the edges list. So we have a maximal loop.
        return True
