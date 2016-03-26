from rx.molecules import *
import unittest


class TestMole(unittest.TestCase):
    def setUp(self):
        self.mole=molecule("H2O2")
        self.mole.addatom('O',[0.0,0.0,0.0],unit='angstrom')
        self.mole.addatom('H',[0.0,1.0,0.0],unit='angstrom')
        self.mole.addatom('O',[0.0,0.0,1.0])
        self.mole.addatom('H',[1.0,0.0,0.0],unit='angstrom')
    def test_atom(self):
        self.assertEqual(self.mole[1].xyz,[0.0,0.0,0.0])
        self.assertEqual(self.mole[2].xyz,[0.0,1.0,0.0])
        self.assertEqual(self.mole[3].xyz,[0.0,0.0,0.5291772086])
        self.assertEqual(self.mole[4].xyz,[1.0,0.0,0.0])
        self.assertEqual(str(self.mole[1].name),'O1')
        self.assertEqual(self.mole.atomlist[1].elementid,8)
        self.assertEqual(str(self.mole.atomlist[1].mymolecule.name),'H2O2')
        ite=iter(self.mole)
        self.assertEqual(next(ite).name,'O1')
        self.assertEqual(next(ite).name,'H2')
        self.assertEqual(next(ite).name,'O3')
        self.assertEqual(next(ite).name,'H4')



if __name__=='__main__':
    unittest.main()
