#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import rx.molecules as rxmol
import rx.chemfiles as rxccfile
import unittest
from io import StringIO

class TestConnect(unittest.TestCase):
    def test_xyzfile(self):
        benzene=rxmol.Molecule("benzene")
        benfile=rxccfile.File("samples/bencom")
        benfile.readfchk()
        xyzfile=benfile.xyzfile
        xyzfile=StringIO(xyzfile)
        benzene.readfromxyz(xyzfile)
        self.assertEqual(benzene[1].name,'C1')
        self.assertEqual(benzene[12].name,'H12')
    def test_mmfile(self):
        mmfile=rxccfile.File('samples/mmfile')
        mmfile.com.read()
        xyz=StringIO(mmfile.com.xyz)
        benmol=rxmol.Molecule('benzene')
        benmol.readfromxyz(xyz)
        benmol.readtypefromlist(mmfile.com.atomtypelist)
        benmol.readchargefromlist(mmfile.com.atomchargelist)
        for atom in benmol:
            print(atom,atom.atomtype,atom.atomcharge)
        benmol.readconnectivity(mmfile.com.connectivity)
        print("Bond: ",len(benmol.bondlist.keys()))
        print("Angle: ",len(benmol.anglelist.keys()))
        print("Dihd: ",len(benmol.dihdlist.keys()))


if __name__=='__main__':
    unittest.main()
