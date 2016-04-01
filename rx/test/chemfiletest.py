import rx.chemfiles as rxccfile
import rx.molecules as rxmol
import unittest,os


class TestFile(unittest.TestCase):
    def test_File(self):
        self.file=rxccfile.File('benfreq')
        self.assertIsInstance(self.file,rxccfile.File)
        self.assertIsInstance(self.file.com,rxccfile.gauCOM)
        self.assertIsInstance(self.file.log,rxccfile.gauLOG)
        self.assertIsInstance(self.file.fchk,rxccfile.gauFCHK)
        self.assertIsInstance(self.file.ac,rxccfile.amberAC)
    def test_fchk(self):
        self.file=rxccfile.File('benfreq')
        self.assertEqual(self.file.readfchk(),None)

        self.assertEqual(self.file.natoms,12)
        self.assertEqual(self.file.multiplicity,1)
        self.assertEqual(self.file.totalcharge,0)
        self.assertIsInstance(self.file.xyzfile,str)
        tmp=self.file.xyzfile[::-1]
        tmp=float(tmp.split()[0][::-1])
        self.assertEqual(tmp,-0.00242749278)
        hess=self.file.find33Hessian(3,5)
        self.assertEqual(hess[0][0],-2.62467513e-2)
        self.assertEqual(hess[1][1],3.38596612e-2)
        self.assertEqual(hess[2][2],7.20824333e-3)
    def test_com(self):
        self.file=rxccfile.File('bencom')
        self.file.rung09()
        self.file.isover()
        self.file.runformchk()
        self.file.com.read()
    def test_ac(self):
        self.file=rxccfile.File('benac')
        self.file.readac()
        self.assertEqual(self.file.atomtypelist[0],None)
        self.assertEqual(self.file.atomtypelist[1],'ca')
        self.assertEqual(self.file.atomtypelist[12],'ha')
        self.assertEqual(self.file.atomchargelist[0],None)
        self.assertEqual(self.file.atomchargelist[1],-0.117738)
        self.assertEqual(self.file.atomchargelist[12],0.117738)
    def test_log(self):
        self.file=rxccfile.File('benresp')
        self.file.runantecham()
    def tearDown(self):
        os.system('rm A* q* Q* p* gaussian*')

#        self.assertEqual()


if __name__=='__main__':
    import numpy as np
    unittest.main()
