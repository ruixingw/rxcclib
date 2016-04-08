#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import rx.molecules as rxmol
import rx.chemfiles as rxccf
import unittest
from io import StringIO

class TestConnect(unittest.TestCase):
    def test_xyzfile(self):
        benzene=rxmol.Molecule("benzene")
        benfile=rxccf.File("samples/bencom")
        benfile.readfchk()
        xyzfile=benfile.xyzfile
        xyzfile=StringIO(xyzfile)
        benzene.readfromxyz(xyzfile)
        self.assertEqual(benzene[1].name,'C1')
        self.assertEqual(benzene[12].name,'H12')


if __name__=='__main__':
    unittest.main()
