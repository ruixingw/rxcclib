#!/usr/bin/env python3
import rxcclib.chemfiles as rxccfile
import rxcclib.molecules as rxmol
import unittest,os,logging
from io import StringIO

rxccfile.gauCOM.g09rt='myg09boon'
rxccfile.gauCOM.g09a2rt='myg09a2boon'
os.system('rm A* q* Q* p* esout *gaussian* samples/bencom.fchk samples/bencom.chk samples/bencom.log')
rxccfile.gauCOM.rung09=lambda self,x='myg09boon': rxccfile.gauCOM.rung09(self,x)
class TestFile(unittest.TestCase):
    def test_comfchk(self):
        file=rxccfile.File('samples/bencom')
        self.assertIsInstance(file,rxccfile.File)
        self.assertIsInstance(file.com,rxccfile.gauCOM)
        self.assertIsInstance(file.log,rxccfile.gauLOG)
        self.assertIsInstance(file.fchk,rxccfile.gauFCHK)
        self.assertIsInstance(file.ac,rxccfile.amberAC)
        file.com.rung09()
        file.com.isover()
        file.runformchk()
        self.assertEqual(file.com.read(),True)
        self.assertEqual(file.fchk.read(),True)
        self.assertEqual(file.natoms,12)
        self.assertEqual(file.multiplicity,1)
        self.assertEqual(file.totalcharge,0)
        self.assertIsInstance(file.xyzfile,str)
        tmp=file.xyzfile[::-1]
        tmp=float(tmp.split()[0][::-1])
        hess=file.fchk.find33Hessian(3,5)
        self.assertEqual(hess[0][0],-2.62909045e-2)
        self.assertEqual(hess[1][1],3.38743754e-2)
        self.assertEqual(hess[2][2],7.19580040e-3)
    def test_logac(self):
        file=rxccfile.File('samples/benresp')
        file.log.runantecham()
        file.ac.read()
        self.assertEqual(file.atomtypelist[0],None)
        self.assertEqual(file.atomtypelist[1],'ca')
        self.assertEqual(file.atomtypelist[12],'ha')
        self.assertEqual(file.atomchargelist[0],None)
        self.assertEqual(file.atomchargelist[1],-0.117738)
        self.assertEqual(file.atomchargelist[12],0.117738)
    # def test_MMcom(self):
    #     mmfile=rxccfile.File('samples/mmfile')
    #     mmfile.com.read()
    #     xyz=StringIO(mmfile.com.xyz)




#        self.assertEqual()


if __name__=='__main__':
    import numpy as np
    logging.basicConfig(level=logging.INFO)
    logging.info('hahaha')
    unittest.main()
