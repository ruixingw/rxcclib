import rx.molecules as rxmol
import unittest


class TestMole(unittest.TestCase):
    def setUp(self):
        self.mole=rxmol.molecule("H2O2")
        self.assertTrue(isinstance(self.mole,rxmol.molecule))

        self.mole.addatom('O',np.array([0.0,0.0,0.0]),unit='angstrom')
        self.mole.addatom('H',np.array([0.0,1.0,0.0]),unit='angstrom')
        self.mole.addatom('O',np.array([0.0,0.0,1.0]))
        self.mole.addatom('H',np.array([1.0,0.0,0.0]),unit='angstrom')
        self.assertTrue(isinstance(self.mole[1],rxmol.Atom))

        self.mole.addbond(1,2)
        self.mole.addbond(1,3)
        self.mole.addbond(3,4)

        self.mole.addangle(2,1,3)
        self.mole.addangle(4,3,1)
        self.assertTrue(isinstance(self.mole.angle(2,1,3),rxmol.Angle))

        self.mole.adddihd(2,1,3,4)
        self.assertTrue(isinstance(self.mole.dihd(4,3,1,2),rxmol.Dihd))


        # @property natoms
        self.assertEqual(self.mole.natoms,4)

        self.assertEqual(self.mole.atom(1).elementid,8)
        # def __iter__
        ite=iter(self.mole)
        # def __next__
        self.assertEqual(next(ite).name,'O1')
        self.assertEqual(next(ite).name,'H2')
        self.assertEqual(next(ite).name,'O3')
        self.assertEqual(next(ite).name,'H4')
        # def __getitem__
        self.assertEqual(str(self.mole[1].name),'O1')

    def test_mole(self):
        # @property xyz
        np.testing.assert_array_equal(self.mole[1].xyz,np.array([0.0,0.0,0.0]))
        np.testing.assert_array_equal(self.mole[2].xyz,np.array([0.0,1.0,0.0]))
        np.testing.assert_array_equal(self.mole[3].xyz,np.array([0.0,0.0,0.5291772086]))
        np.testing.assert_array_equal(self.mole[4].xyz,np.array([1.0,0.0,0.0]))
        # @property atomnum
        self.assertEqual(self.mole[1].atomnum,1)
        # @property elementid
        self.assertEqual(self.mole[1].elementid,8)
        # @property atomsym
        self.assertEqual(self.mole[1].atomsym,'O')
        # @property name
        self.assertEqual(self.mole[1].name,'O1')
        # @property mymolecule
        self.assertEqual(str(self.mole.atom(1).mymolecule.name),'H2O2')
    # test_bond(self):
        # Bondlength
        self.assertEqual(self.mole.bond(1,3).length,0.5291772086)
    # test_angle(self):
        self.assertEqual(self.mole.angle(2,1,3).anglevalue,90.0)
    # test dihd
        self.assertEqual(self.mole.dihd(4,3,1,2).dihdvalue,90.0)
    # test_bondfunc(self):
        self.mole[1].atomtype='oh'
        self.mole[2].atomtype='ho'
        self.mole[3].atomtype='oh'
        self.mole[4].atomtype='ho'
        self.mole.addbondfunc(self.mole.bond(1,2))
        self.mole.addbondfunc(self.mole.bond(1,3))
        self.mole.addbondfunc(self.mole.bond(3,4))
        self.mole.addanglefunc(self.mole.angle(2,1,3))
        self.mole.addanglefunc(self.mole.angle(4,3,1))
        self.mole.adddihdfunc(self.mole.dihd(2,1,3,4))

    def test_readfile(self):
        benz=rxmol.molecule("Benzene")
        benz.readfromxyz("ben.xyz")
        for atom in benz:
            print(atom.name)
if __name__=='__main__':
    import numpy as np
    unittest.main()
