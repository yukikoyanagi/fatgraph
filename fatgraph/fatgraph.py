import re
from permutation import Permutation

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
    for example verticex=((1,2,3),(4,5,6)), edges=((1,2),(3,4)(5,6))
    """

    def __init__(self, vertices, edges):
        self._vs = Permutation.from_cycles(*vertices)
        self._es = Permutation.from_cycles(*edges)

    def __repr__(self):
        return ('{0.__module__}.{0.__name__}(\n'
                'vertices={1}\n'
                'edges={2})').format(type(self),
                                     self.vertices,
                                     self.edges)

    def __eq__(self, other):
        return all([self.vertices == other.vertices,
                    self.edges == other.edges])

    def __ne__(self, other):
        return not self == other

    @property
    def vertices(self):
        return self._vs.to_cycles()

    @property
    def edges(self):
        return self._es.to_cycles()

    @property
    def boundaries(self):
        bs = self._vs * self._es
        cycles = bs.to_cycles()
        # to_cycles omits one-vertex boundaries as identity, so
        # need to add them.
        fixed = set(i for v in self.vertices for i in v) - \
            set(j for c in cycles for j in c)
        for f in fixed:
            cycles.append((f,))
        return cycles

    @property
    def genus(self):
        return int((2 - len(self.vertices)
                    + len(self.edges)
                    - len(self.boundaries)) / 2)

    @classmethod
    def from_hbonds(cls, hbfile):
        """
        Create fatgraph object from PDB Hbond file. The resulting
        fatgraph is simplified; that is, non-bonded half-edges are
        not recorded. Note this does not change genus.
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
        dons = sorted(dons)
        dons_dic = {v: k+1 for k, v in enumerate(dons)}
        dons_max = dons_dic[dons[-1]]
        accs = sorted(accs, reverse=True)
        accs_dic = {v: k+1+dons_max for k, v in enumerate(accs)}

        edges = []
        for bond in bonds:
            edges.append((dons_dic[bond[0]], accs_dic[bond[1]]))

        vertices = [list(range(1, accs_dic[accs[-1]]+1))]

        return cls(vertices, edges)
