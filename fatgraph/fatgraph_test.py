import pytest
from permutation import Permutation
import numpy as np
from .fatgraph import Fatgraph
from .fatgraph import FatgraphB

class TestFatgraphB(object):
    def test_from_fatgraph(self):
        fg = Fatgraph([(1,2,3,4),], [(1,2), (3,4)])
        i = [(1,3), (2,4)]
        fgb = FatgraphB.from_fatgraph(fg, i)
    
    def test_isvalid(self):
        fg = FatgraphB([(1,2,3,4),],
                       [(1,2), (3,4)],
                       [(1,4), (2,3)])
        assert fg.isvalid()
        fg._is = Permutation.from_cycles((1,2), (3,4))
        assert not fg.isvalid()
        fgb = FatgraphB([tuple(range(1,9)),],
                        [(1,8), (2,3), (4,7), (5,6)],
                        [(2,7), (3,6), (4,5), (1,8)])
        assert not fgb.isvalid()
        fgb = FatgraphB([(1,2,3,4), (5,6,7,8)],
                        [(1,8), (2,7), (3,5), (4,6)],
                        [(1,2), (3,4), (5,6), (7,8)])
        assert not fgb.isvalid()
        fgb._is = Permutation.from_cycles((1,4), (2,3), (5,6), (7,8))
        assert fgb.isvalid()
        fg = FatgraphB([(1,2,3,4,5,6),],
                       [(3,6), (4,5)],
                       [(1,4), (2,3), (5,6)])
        assert fg.isvalid()
        fg = FatgraphB([(1,2,3,4,5,6),],
                       [(3,6), (4,5)],
                       [(1,2), (6,3), (5,4)])
        assert not fg.isvalid()
        fg = FatgraphB([(1,2,3,4,5,6),],
                       [(1,4), (3,5)],
                       [(1,6), (2,5), (3,4)])
        assert fg.isvalid()

    def test_from_pmat(self):
        mat = np.array([[0,1,0,0],
                        [0,0,1,1],
                        [0,0,0,0],
                        [0,0,1,0]])
        with pytest.raises(ValueError):
            fgb = FatgraphB.from_pmat(mat)
        mat = np.array([[0,1,1],
                        [0,0,1],
                        [0,0,0]])
        with pytest.raises(ValueError):
            fgb = FatgraphB.from_pmat(mat)
        mat = np.array([[0,1,0,0,0],
                        [0,0,0,0,0],
                        [0,0,0,0,1],
                        [0,0,1,0,0],
                        [0,0,0,1,0]])
        with pytest.raises(ValueError):
            fgb = FatgraphB.from_pmat(mat)
        mat = np.array([[0,0,0],
                        [1,0,0],
                        [0,1,0]])
        fg = FatgraphB([(1,2,3,4,5,6),],
                      [(2,3), (5,6)],
                      [(1,6), (2,5), (3,4)])
        assert FatgraphB.from_pmat(mat) == fg
        mat = np.array([[0,0,0],
                        [1,0,0],
                        [1,0,0]])
        fg = FatgraphB([(1,2,3,4,5,6),],
                      [(1,4), (5,6)],
                      [(1,6), (2,5), (3,4)])
        assert FatgraphB.from_pmat(mat) == fg
        mat = np.array([[0,0,0],
                        [1,0,0],
                        [0,0,0]])
        fg = FatgraphB([(1,2,3,4),(5,6)],
                      [(3,4), (2,5)],
                      [(1,4), (2,3), (5,6)])
        assert FatgraphB.from_pmat(mat) == fg

class TestFatgraph(object):
    def test_Fatgraph(self):
        g1 = Fatgraph([(1,2,3),], [(1,2),])
        assert g1.unpaired == {3} # set
        g2 = Fatgraph([(2,3,1),], [(2,1),])
        assert g1 == g2
        g2 = Fatgraph([(1,2),], [(2,1),])
        assert g1 != g2
        g1 = Fatgraph([(1,2,3),(4,5,6)], [(1,2), (3,4), (5,6)])
        assert len(g1.boundaries) == 3
        g1 = Fatgraph([(1,2,3),(4,5,6)], [(1,4), (2,6), (3,5)])
        assert len(g1.boundaries) == 3
        g2 = Fatgraph([(1,2,3),(4,5,6)], [(1,4), (2,5), (3,6)])
        assert len(g2.boundaries) == 1
        assert g1.genus == 0
        assert g2.genus == 1
        g3 = Fatgraph([], [])
        assert g3.genus == 0
        g4 = Fatgraph([(1,2,3,4),(5,6)], [(3,4), (2,5)])
        assert len(g4.boundaries) == 2
        assert g4.genus == 0

        with pytest.raises(ValueError):
            assert Fatgraph([(1,2,3), (3,4)], [(1,2), (3,4)])
            assert Fatgraph([(1,2,3,4),], [(1,2), (2,3)])
            assert Fatgraph([(1,2,3),], [(1,2), (3,4)])

    def test_isconnected(self):
        fg = Fatgraph([(1,2), (3,4)], [(1,3), (2,4)])
        assert fg.isconnected()
        fg = Fatgraph([(1,2), (3,4)], [(1,2), (3,4)])
        assert not fg.isconnected()
        fg = Fatgraph([(1,2), (3,4)], [(1,3), ])
        assert fg.isconnected()
        fg = Fatgraph([(1,2), (3,4)], [(1,2), ])
        assert not fg.isconnected()
        fg = Fatgraph([(1,2,3), (4,5,6), (7,8,9)],
                      [(1,3), (2,5), (6,7), (8,9)])
        assert fg.isconnected()
        fg = Fatgraph([(1,2,3), (4,5,6), (7,8,9)],
                      [(1,3), (2,5), (6,4), (8,9)])
        assert not fg.isconnected()
