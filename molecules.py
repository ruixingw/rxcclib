# Molecules, atoms definition
from __future__ import print_function
import numpy as np
import cclib.parser.utils as cclibutils
import inspect, logging
from io import StringIO


class MoleDefError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class dihdforceconst(object):
    def __init__(self, value, dihd):
        self.value = value
        self.dihd = dihd
        self.repr = dihd.repr

    def __repr__(self):
        return repr(self.value)

    def __str__(self):
        return str(self.value)

    def __call__(self, value):
        self.value = value

    @property
    def forceconst(self):
        return self.value

    @forceconst.setter
    def forceconst(self, value):
        self.value = value


class Molecule(object):
    def __init__(self, moleculename):
        if not isinstance(moleculename, str):
            logging.error('Error: Molecule name must be str')
            raise MoleDefError('Error: Molecule name must be str')

        self.name = moleculename
        self.__atomlist = [None]
        self.__bondlist = {}
        self.__anglelist = {}
        self.__dihdlist = {}
        self.__improperlist = {}

    @property
    def natoms(self):
        return len(self.__atomlist) - 1

    @property
    def bondlist(self):
        return self.__bondlist

    @property
    def anglelist(self):
        return self.__anglelist

    @property
    def dihdlist(self):
        return self.__dihdlist

    @property
    def improperlist(self):
        return self.__improperlist
    # Add atom and internal coordinates
    def addatom(self, idorsym, coords, unit='angstrom'):
        self.__atomlist.append(Atom(self, idorsym, coords, unit))

    def addbond(self, atomnum1, atomnum2):
        if atomnum1 > atomnum2:
            atomnum1, atomnum2 = atomnum2, atomnum1
        self.__bondlist.update({str(atomnum1) + '-' + str(atomnum2): Bond(
            self, atomnum1, atomnum2)})

    def addangle(self, atomnum1, atomnum2, atomnum3):
        if atomnum1 > atomnum3:
            atomnum1, atomnum3 = atomnum3, atomnum1
        self.__anglelist.update(
            {str(atomnum1) + '-' + str(atomnum2) + '-' + str(atomnum3): Angle(
                self, atomnum1, atomnum2, atomnum3)})

    def adddihd(self, atomnum1, atomnum2, atomnum3, atomnum4):
        if atomnum2 > atomnum3:
            atomnum2, atomnum3 = atomnum3, atomnum2
            atomnum1, atomnum4 = atomnum4, atomnum1
        elif atomnum2 == atomnum3:
            if atomnum1 > atomnum4:
                atomnum1, atomnum4 = atomnum4, atomnum1
        self.__dihdlist.update({
            str(atomnum1) + '-' + str(atomnum2) + '-' + str(atomnum3) + '-' +
            str(atomnum4): Dihd(self, atomnum1, atomnum2, atomnum3, atomnum4)
        })

    def addimproper(self, atomnum1, atomnum2, atomnum3, atomnum4):
        self.__improperlist.update(
            {str(atomnum1) + '-' + str(atomnum2) + '-' + str(atomnum3) + '-' +
             str(atomnum4): Improper(self, atomnum1, atomnum2, atomnum3,
                                     atomnum4)})

    # Get atom and internal coordiantes
    def atom(self, atomnum):
        return self.__atomlist[atomnum]

    def bond(self, atomnum1, atomnum2):
        if atomnum1 > atomnum2:
            atomnum1, atomnum2 = atomnum2, atomnum1
        return self.__bondlist[str(atomnum1) + '-' + str(atomnum2)]

    def angle(self, atomnum1, atomnum2, atomnum3):
        if atomnum1 > atomnum3:
            atomnum1, atomnum3 = atomnum3, atomnum1
        return self.__anglelist[str(atomnum1) + '-' + str(atomnum2) + '-' +
                                str(atomnum3)]

    def dihd(self, atomnum1, atomnum2, atomnum3, atomnum4):
        if atomnum2 > atomnum3:
            atomnum2, atomnum3 = atomnum3, atomnum2
            atomnum1, atomnum4 = atomnum4, atomnum1
        elif atomnum2 == atomnum3:
            if atomnum1 > atomnum4:
                atomnum1, atomnum4 = atomnum4, atomnum1
        return self.__dihdlist[str(atomnum1) + '-' + str(atomnum2) + '-' + str(
            atomnum3) + '-' + str(atomnum4)]

    def improper(self, atomnum1, atomnum2, atomnum3, atomnum4):
        return self.__improperlist[str(atomnum1) + '-' + str(atomnum2) + '-' +
                                   str(atomnum3) + '-' + str(atomnum4)]

    # getitem and iteration
    def __getitem__(self, key):
        if isinstance(key, int):
            return self.__atomlist[key]
        if isinstance(key, slice):
            start = key.start
            stop = key.stop
            if start is None:
                start = 1
            L = []
            for x in range(start, stop):
                L.append(self.__atomlist[x])
            return L

    def __iter__(self):
        for i in range(self.natoms):
            yield self[i + 1]

    # Read structure from coords
    def readfromxyz(self, string):
        f = StringIO(string)
        for line in f:
            tmp = line.split()
            atom = tmp[0]
            if atom.isdigit():
                atom = int(atom)
            coords = np.array([tmp[1], tmp[2], tmp[3]], dtype=float)
            self.addatom(atom, coords, unit='angstrom')
    # Read connectivity
    def readconnectivity(self, conntystring):
        f = StringIO(conntystring)
        for line in f:
            try:
                tmp = line.split()
            except:
                logging.debug('End of connectivity, return.')
                return
            ite = iter(tmp)
            item0 = next(ite)
            if not item0.isdigit():
                logging.debug('End of connectivity, return.')
                return
            a = int(item0)
            try:
                b = int(next(ite))
            except StopIteration:
                continue
            self.addbond(a, b)
            while True:
                try:
                    b = next(ite)
                    b = next(ite)
                    b = int(b)
                    self.addbond(a, b)
                except StopIteration:
                    break
        self.geninlcoords()

    def geninlcoords(self):
        angles = []
        dihds = []
        for atom1 in self:  # atom1: atom obj
            for atom2 in atom1.neighbor:  # atom2: obj
                if atom2.atomnum == atom1.atomnum:
                    continue
                for atom3 in atom2.neighbor:  #atom3: obj
                    if atom3.atomnum == atom2.atomnum or atom3.atomnum == atom1.atomnum:
                        continue
                    a = atom1.atomnum
                    b = atom2.atomnum
                    c = atom3.atomnum
                    if a > c:
                        a, c = c, a
                    angles.append(str(a) + '-' + str(b) + '-' + str(c))
                    for atom4 in atom3.neighbor:
                        if atom4.atomnum == atom3.atomnum or atom4.atomnum == atom2.atomnum or atom4.atomnum == atom1.atomnum:
                            continue
                        a = atom1.atomnum
                        b = atom2.atomnum
                        c = atom3.atomnum
                        d = atom4.atomnum
                        if b > c:
                            b, c = c, b
                            a, d = d, a
                        elif b == c:
                            if a > d:
                                a, d = d, a
                        dihds.append(str(a) + '-' + str(b) + '-' + str(c) + '-'
                                     + str(d))
        angles = list(set(angles))
        dihds = list(set(dihds))
        for item in angles:
            tmp = [int(x) for x in item.split('-')]
            self.addangle(*tmp)
        for item in dihds:
            tmp = [int(x) for x in item.split('-')]
            self.adddihd(*tmp)

    def readtypefromlist(self, L):
        if len(L) != len(self.__atomlist):
            logging.error(
                "Error when reading atomtype from list: length is not consistent with natoms")
            raise MoleDefError(
                "Error when reading atomtype from list: length is not consistent with natoms")
        ite = iter(L)
        next(ite)
        for atom in self:
            atom.atomtype = next(ite)

    def readchargefromlist(self, L):
        if len(L) != len(self.__atomlist):
            logging.error(
                "Error when reading atomcharge from list: length is not consistent with natoms")
            raise MoleDefError(
                "Error when reading atomcharge from list: length is not consistent with natoms")
        ite = iter(L)
        next(ite)
        for atom in self:
            atom.atomcharge = next(ite)


class Atom(object):
    '''
    Class "Atom" for atom object.
    Never add an atom directly. Always use molecule.addatom method instead, otherwise the molecule.__atomlist is not correct.

    >>> mole=molecule("CO")  # Define a molecule
    >>> mole.addatom(6,np.array([0.0,0.0,0.0]))
    >>> mole.addatom('O',np.array([0.0,0.0,2.0]))
    >>> mole.name
    'CO'
    >>> mole.atom(1).name
    'C1'
    >>> mole.atom(2).name
    'O2'
    >>> mole.atom(1).mymolecule.name
    'CO'
    >>> mole.atom(1).coords
    array([ 0.,  0.,  0.])
    >>> setattr(mole.atom(1),'atomtype','c2')
    >>> mole.atom(1).atomtype
    'c2'
    '''

    periotable = cclibutils.PeriodicTable()

    def __init__(self, mole, idorsym,
                 coords,
                 unit='angstrom'):  # molecule object,int,[float,float,float]
        # Assertion
        callername = inspect.stack()[1][3]
        assert callername == 'addatom', "Atom must be added via Molecule.addatom method"
        assert isinstance(
            mole,
            Molecule), "First argument must be a molecule object!. Use Molecule.addatom method to avoid this problem."
        assert unit != 'bohr' or unit != 'angstrom', "Coordinate unit must be bohr or angstrom"

        self.__mymolecule = mole
        if isinstance(idorsym, int):
            self.elementid = idorsym
            try:
                self.atomsym = Atom.periotable.element[self.elementid]  #str
            except KeyError:
                logging.critical(
                    "Error when adding atom: Idtosym not defined for atomic no:"
                    + str(self.elementid))
                raise MoleDefError(
                    "Error when adding atom: Idtosym not defined for atomic no:"
                    + str(self.elementid))
        elif isinstance(idorsym, str):
            self.atomsym = idorsym
            try:
                self.elementid = Atom.periotable.number[self.atomsym]
            except KeyError:
                logging.critical(
                    "Error when adding atom: Idtosym not defined for atomic symbol:"
                    + self.atomsym)
                raise MoleDefError(
                    "Error when adding atom: Idtosym not defined for atomic symbol:"
                    + self.atomsym)
        else:
            logging.critical(
                "Error when adding atom: Expected atomic NO(int) or symbol(str) for input, received a"
                + str(type(idorsym)))
            raise MoleDefError(
                "Error when adding atom: Expected atomic NO(int) or symbol(str) for input, received a"
                + str(type(idorsym)))

        if unit == 'bohr':
            self.coords = cclibutils.convertor(coords, "bohr", "Angstrom")
        elif unit == 'angstrom':
            self.coords = coords

        self.atomnum = mole.natoms + 1
        self.atomtype = self.name
        self.__neighbor = []
        self.atomcharge = None

    @property
    def neighbor(self):
        return self.__neighbor

    def addneighbor(self, atomobj):
        if isinstance(atomobj, Atom):
            self.__neighbor.append(atomobj)
            self.__neighbor = list(set(self.__neighbor))
        else:
            logging.error(
                "Error when adding neighbor: atomobj must be an integer.")
            raise MoleDefError(
                "Error when adding neighbor: atomobj must be an integer.")

    def delneighbor(self, atomobj):
        if isinstance(atomobj, Atom):
            self.__neighbor.remove(atomobj)
        else:
            logging.error(
                "Error when deleting neighbor: atomobj must be an integer.")
            raise MoleDefError(
                "Error when deleting neighbor: atomobj must be an integer.")

    @property
    def name(self):
        return str(self.atomsym) + str(self.atomnum)

    @property
    def mymolecule(self):
        return self.__mymolecule

    def __str__(self):
        return "Atom object for atom " + self.name


class Bond(object):
    def __init__(self, mole, a, b):  # self, atomid a, atomid b
        if a > b:
            a, b = b, a
        self.__a = mole[a]
        self.__b = mole[b]
        self.vec = self.__a.coords - self.__b.coords
        self.repr = self.__a.name + ' ' + self.__b.name
        self.__a.addneighbor(mole[b])
        self.forceconst = 0.0
        self.__b.addneighbor(mole[a])

    def __getitem__(self, value):
        if value == 1:
            return self.__a
        elif value == 2:
            return self.__b
        else:
            raise MoleDefError("Index for bond object must be 1 or 2.")

    @property
    def length(self):
        return np.linalg.norm(self.vec)

    def __str__(self):
        return "Bond object of bond " + self.__a.name + '-' + self.__b.name


class Angle(object):
    def __init__(self, mole, a, b, c):
        if a > c:
            a, c = c, a
        self.__a = mole[a]
        self.__b = mole[b]
        self.__c = mole[c]
        self.__ab = mole[a].coords - mole[b].coords
        self.__bc = mole[b].coords - mole[c].coords
        self.repr = self.__a.name + ' ' + self.__b.name + ' ' + self.__c.name
        self.forceconst = 0.0

    def __getitem__(self, value):
        if value == 1:
            return self.__a
        elif value == 2:
            return self.__b
        elif value == 3:
            return self.__c
        else:
            raise MoleDefError("Index for angle object must be 1, 2 or 3.")

    @property
    def anglevalue(self):
        v1u = self.__ab / np.linalg.norm(self.__ab)
        v2u = self.__bc / np.linalg.norm(self.__bc)
        angle = 180.0 - np.arccos(np.dot(v1u, v2u)) * 180.0 / np.pi
        if np.isnan(angle):
            if (v1u == v2u).all():
                return 0.0
            else:
                return 180.0
        return angle

    def __str__(self):
        return "Angle object of angle " + self.__a.name + '-' + self.__b.name + '-' + self.__c.name


class Dihd(object):
    def __init__(self, mole, a, b, c, d):
        if b > c:
            b, c = c, b
            a, d = d, a
        elif b == c:
            if a > d:
                a, d = d, a
        self.__a = mole[a]
        self.__b = mole[b]
        self.__c = mole[c]
        self.__d = mole[d]
        self.repr = self.__a.name + ' ' + self.__b.name + ' ' + self.__c.name + ' ' + self.__d.name
        self.forceconst = [dihdforceconst(0.0, self),
                           dihdforceconst(0.0, self),
                           dihdforceconst(0.0, self),
                           dihdforceconst(0.0, self)]
        self.phase = [0, 0, 0, 0]
        self.npaths = 1.0

    def __getitem__(self, value):
        if value == 1:
            return self.__a
        elif value == 2:
            return self.__b
        elif value == 3:
            return self.__c
        elif value == 4:
            return self.__d
        else:
            raise MoleDefError(
                "Index for dihedral object must be 1, 2 ,3 or 4.")

    @property
    def dihdvalue(self):
        v1 = self.__a.coords - self.__b.coords
        v2 = self.__b.coords - self.__c.coords
        v3 = self.__c.coords - self.__d.coords
        v1u = v1 / np.linalg.norm(v1)
        v2u = v1 / np.linalg.norm(v2)
        v3u = v1 / np.linalg.norm(v3)
        n1 = np.cross(v1u, v2u)
        n2 = np.cross(v2u, v3u)
        dihd = np.arccos(np.dot(n1, n2)) * 180.0 / np.pi
        if np.isnan(dihd):
            if (n1 == n2).all():
                return 0.0
            else:
                return 180.0
        return dihd

    def __str__(self):
        return "Dihedral object of dihedral " + self.__a.name + '-' + self.__b.name + '-' + self.__c.name + '-' + self.__d.name


class Improper(object):
    def __init__(self, mole, a, b, c, d):
        self.__a = mole[a]
        self.__b = mole[b]
        self.__c = mole[c]
        self.__d = mole[d]
        self.forceconst = 0.0
        self.phase = 180.0
        self.npaths = 2.0
        self.repr = self.__a.name + ' ' + self.__b.name + ' ' + self.__c.name + ' ' + self.__d.name

    def __getitem__(self, value):
        if value == 1:
            return self.__a
        elif value == 2:
            return self.__b
        elif value == 3:
            return self.__c
        elif value == 4:
            return self.__d
        else:
            raise MoleDefError(
                "Index for improper object must be 1, 2 ,3 or 4.")

    @property
    def impropervalue(self):
        v1 = self.__a.coords - self.__b.coords
        v2 = self.__b.coords - self.__c.coords
        v3 = self.__c.coords - self.__d.coords
        v1u = v1 / np.linalg.norm(v1)
        v2u = v1 / np.linalg.norm(v2)
        v3u = v1 / np.linalg.norm(v3)
        n1 = np.cross(v1u, v2u)
        n2 = np.cross(v2u, v3u)
        improper = np.arccos(np.dot(n1, n2)) * 180.0 / np.pi
        if np.isnan(improper):
            if (n1 == n2).all():
                return 0.0
            else:
                return 180.0
        return improper

    def __str__(self):
        return "Improper object of improper " + self.__a.name + '-' + self.__b.name + '-' + self.__c.name + '-' + self.__d.name


if __name__ == '__main__':
    import doctest
    doctest.testmod()
