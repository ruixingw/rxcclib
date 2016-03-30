#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import rx.molecules as rxmol
import rx.chemfiles as rxccf
import unittest
from io import StringIO

class TestConnect(unittest.TestCase):
    def test_xyzfile(self):
        benzene=rxmol.molecule("benzene")
        benfile=rxccf.File("ben")
        benfile.readfchk()
        xyzfile=benfile.xyzfile()
        xyzfile=StringIO(xyzfile)
        benzene.readfromxyz(xyzfile)
        for atom in benzene:
            print(atom.name)


if __name__=='__main__':
    unittest.main()
