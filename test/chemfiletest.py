import rx.chemfiles as rxccfile
import rx.molecules as rxmol
import unittest,os
from io import StringIO

rxccfile.gauCOM.g09rt='myg09boon'
rxccfile.gauCOM.g09a2rt='myg09a2boon'
os.system('rm A* q* Q* p* esout *gaussian* samples/bencom.fchk samples/bencom.chk samples/bencom.log')
class TestFile(unittest.TestCase):
    def test_comfchk(self):
        self.file=rxccfile.File('samples/bencom')
        self.assertIsInstance(self.file,rxccfile.File)
        self.assertIsInstance(self.file.com,rxccfile.gauCOM)
        self.assertIsInstance(self.file.log,rxccfile.gauLOG)
        self.assertIsInstance(self.file.fchk,rxccfile.gauFCHK)
        self.assertIsInstance(self.file.ac,rxccfile.amberAC)
        self.file.rung09()
        self.file.isover()
        self.file.runformchk()
        self.file.com.read()
        self.assertEqual(self.file.readfchk(),None)
        self.assertEqual(self.file.natoms,12)
        self.assertEqual(self.file.multiplicity,1)
        self.assertEqual(self.file.totalcharge,0)
        self.assertIsInstance(self.file.xyzfile,str)
        tmp=self.file.xyzfile[::-1]
        tmp=float(tmp.split()[0][::-1])
        self.assertEqual(tmp,-0.00242749278)
        hess=self.file.find33Hessian(3,5)
        self.assertEqual(hess[0][0],-2.62909045e-2)
        self.assertEqual(hess[1][1],3.38743754e-2)
        self.assertEqual(hess[2][2],7.19580040e-3)
    def test_logac(self):
        self.file=rxccfile.File('samples/benresp')
        self.file.runantecham()
        self.file.readac()
        self.assertEqual(self.file.atomtypelist[0],None)
        self.assertEqual(self.file.atomtypelist[1],'ca')
        self.assertEqual(self.file.atomtypelist[12],'ha')
        self.assertEqual(self.file.atomchargelist[0],None)
        self.assertEqual(self.file.atomchargelist[1],-0.117738)
        self.assertEqual(self.file.atomchargelist[12],0.117738)
    # def test_MMcom(self):
    #     mmfile=rxccfile.File('samples/mmfile')
    #     mmfile.com.read()
    #     xyz=StringIO(mmfile.com.xyz)




#        self.assertEqual()


if __name__=='__main__':
    import numpy as np
    unittest.main()
