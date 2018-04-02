import os
import time
import logging
import numpy as np
import rxcclib.utils as utils
import rxcclib.geometry.molecules as rxmol
import rxcclib.utils.periodictable as periodictable


class PDBFile(object):
    def __init__(self,filename):
        self.abspath = os.path.abspath(filename)

    def read(self):
        logging.info('Read PDB:' + self.abspath)
        self.atomlist = []
        self.coordslist = []
        self.atomtypelist = []
        self.connectivity = ''
        # PDB parser
        with open(self.abspath, 'r') as f:
            for line in f:
                if line.find('CRYST1') == 0:
                    # a,b,c,alpha,beta,gamma
                    self.crystalparm = list(map(float,[line[6:15],line[15:24],line[24:33],line[33:40],line[40:47],line[47:54]]))

                    self.sGroup = line[55:66].strip()
                    self.z = int(line[65:70])
                if line.find('ATOM')==0 or line.find('HETATM')==0:
                    self.atomlist.append(line[76:78].strip())
                    self.coordslist.append(float(line[30:38]))
                    self.coordslist.append(float(line[38:46]))
                    self.coordslist.append(float(line[46:54]))
                    self.atomtypelist.append(line[12:16].strip())
                if line.find('CONECT')==0:
                    self.connectivity += line[6:]

        return

