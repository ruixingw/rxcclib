# Molecules, atoms definition
from __future__ import print_function
import math
import numpy as np

class MoleDefError(Exception):
    def __init__(self,value):
        self.value=value
    def __str__(self):
        return repr(self.value)


class molecule(object):
    def __init__(self,moleculename):
        if not isinstance(moleculename,str):
            raise MoleDefError('Error: Molecule name must be str')

        self.name=moleculename
        self.__atomlist=[0]
        self.bonds={}
        self.angles={}
        self.dihedrals={}
        self.bondfunc={}
        self.anglefunc={}
        self.dihdfunc={}
        self.pointer=0
    @property
    def natoms(self):
        return len(self.__atomlist)-1
    @natoms.setter
    def natoms(self,value):  # Should only be called by other method
        assert isinstance(value,int), "Added an non-int to natoms"
        self.__natoms=value

    def addatom(self,idorsym,xyz,unit='bohr'):
        self.__atomlist.append(Atom(self,idorsym,xyz,unit))
    @property
    def atomlist(self):
        return self.__atomlist

    def addbond(self,atomid1,atomid2):
        bonds(self,atomid1,atomid2)
    def addangle(self,atomid1,atomid2,atomid3):
        angles(self,atomid1,atomid2,atomid3)
    def adddihd(self,atomid1,atomid2,atomid3,atomid4):
        dihedrals(self,atomid1,atomid2,atomid3,atomid4)
    def addbondfunc(self,value):
        bondfunc(self,value)
    def addanglefunc(self,value):
        anglefunc(self,value)
    def adddihdfunc(self,value):
        dihdfunc(self,value)
    def __getitem__(self,key):
        if isinstance(key,int):
            return self.__atomlist[key]
        if isinstance(key,slice):
            start=key.start
            stop=key.stop
            if start is None:
                start=1
            L=[]
            for x in range(start,stop):
                L.append(self.__atomlist[x])
            return L
    def __iter__(self):
        self.pointer=0
        return self
    def __next__(self):
        self.pointer+=1
        if self.pointer <=self.natoms:
            return self[self.pointer]
        else:
            raise StopIteration()


class Atom(object):
    '''
    Class "Atom" for atom object.
    Never add an atom directly. Always use molecule.addatom method instead, otherwise the molecule.atomlist is not correct.

    >>> mole=molecule("CO")  # Define a molecule
    >>> mole.addatom(6,[0.0,0.0,0.0])
    >>> mole.addatom('O',[0.0,0.0,2.0])
    >>> mole.name
    'CO'
    >>> mole.atom[1]
    C1
    >>> mole.atom[2]
    O2
    '''
    bohr=0.5291772086 # bohr to angstrom
    __idtosym={1:'H',5:'B',6:'C',7:'N',8:'O',9:'F',13:'Al',14:'Si',15:'P',16:'S',17:'Cl',26:'Fe',28:'Ni',29:'Cu',30:'Zn'}
    __symtoid={v:k for k,v in __idtosym.items()}

    def __init__(self,mole,idorsym,xyz,unit='bohr'):  # molecule object,int,[float,float,float]
        assert isinstance(mole,molecule),"First argument must be a molecule object!. Use molecule.addatom method to avoid this problem."
        assert unit!='bohr' or unit!='angstrom', "Coordinate unit must be bohr or angstrom"
        self.__molecule=mole
        if isinstance(idorsym,int):
            self.__elementid=idorsym
            try:
                self.__atomsym=Atom.__idtosym[self.__elementid] #str
            except KeyError:
                print("Error when adding atom: Idtosym not defined for atomic no:",self.__elementid)
                quit()
        elif isinstance(idorsym,str):
            self.__atomsym=idorsym
            try:
                self.__elementid=Atom.__symtoid[self.__atomsym]
            except KeyError:
                print("Error when adding atom: Idtosym not defined for atomic symbol:",self.__atomsym)
                quit()
        else:
                print("Error when adding atom: Expected atomic NO(int) or symbol(str) for input, received a",type(idorsym))
                quit()
        if unit=='bohr':
            self.__x=xyz[0]*Atom.bohr
            self.__y=xyz[1]*Atom.bohr
            self.__z=xyz[2]*Atom.bohr
        elif unit=='angstrom':
            self.__x=xyz[0]
            self.__y=xyz[1]
            self.__z=xyz[2]

        self.__atomnum=mole.natoms+1
        return
    @property
    def xyz(self):
        return [self.__x,self.__y,self.__z]
    @property
    def name(self):
        return str(self.atomsym)+str(self.atomnum)
    @property
    def atomnum(self):
        return self.__atomnum
    @property
    def elementid(self):
        return self.__elementid
    @property
    def atomsym(self):
        return self.__atomsym
    @property
    def mymolecule(self):
        return self.__molecule
    # @mymolecule.setter
    # def mymolecule(self,value):
    #     if not isinstance(value,molecule):
    #         print("Error when changing molecule of,",self.atomsym,self.atomnum,"from",self.molecule,"to",value,": Expected a molecule object, received a",type(value))
    #         quit()
    #     self.__molecule=value


class bonds(object):
    def __init__(self,molecule,a,b): # self, atomid a, atomid b
        if a>b:
            a,b=b,a
        self.a=molecule.atoms[a]
        self.b=molecule.atoms[b]
        molecule.bonds.update({str(a)+','+str(b):self})
        self.vec=vector(self.a.x-self.b.x,self.a.y-self.b.y,self.a.z-self.b.z)
    def length(self):
        return self.vec.length()
class angles(object):
    def __init__(self,molecule,a,b,c):
        if a>c:
            a,c=c,a
        self.a=molecule.atoms[a]
        self.b=molecule.atoms[b]
        self.c=molecule.atoms[c]
        self.ab=vector(self.a.x-self.b.x,self.a.y-self.b.y,self.a.z-self.b.z)
        self.bc=vector(self.b.x-self.c.x,self.b.y-self.c.y,self.b.z-self.c.z)
        molecule.angles.update({str(a)+','+str(b)+','+str(c):self})

    def angle(self):
        innerproduct=self.ab.x*self.bc.x+self.ab.y*self.bc.y+self.ab.z*self.bc.z
        product=self.ab.length()*self.bc.length()
        cos=innerproduct/product
        angle=180-math.acos(cos)*180/math.pi
        return angle
class dihedrals(object):
    def __init__(self,molecule,a,b,c,d):
        if b>c:
            b,c=c,b
            a,d=d,a
        elif b==c:
            if a>d:
                a,d=d,a
        self.a=molecule.atoms[a]
        self.b=molecule.atoms[b]
        self.c=molecule.atoms[c]
        self.d=molecule.atoms[d]

        molecule.dihedrals.update({str(a)+','+str(b)+','+str(c)+','+str(d):self})


class bondfunc(object):

    def __init__(self,molecule,bondobj):
        a=bondobj.a.atomtype
        b=bondobj.b.atomtype
        if a>b:
            a,b=b,a
        bondobj.func=self
        self.link=a+' '+b
        molecule.bondfunc.update({self.link:self})
class anglefunc(object):

    def __init__(self,molecule,angleobj):
        a=angleobj.a.atomtype
        b=angleobj.b.atomtype
        c=angleobj.c.atomtype
        if a>c:
            a,c=c,a
        angleobj.func=self
        self.link=a+' '+b+' '+c
        molecule.anglefunc.update({self.link:self})
class dihdfunc(object):

    def __init__(self,molecule,dihdobj):
        a=dihdobj.a.atomtype
        b=dihdobj.b.atomtype
        c=dihdobj.c.atomtype
        d=dihdobj.d.atomtype
        if b>c:
            a,d=d,a
            b,c=c,b
        elif b==c:
            if a>d:
                a,d=d,a
        dihdobj.func=self
        self.link=a+' '+b+' '+c+' '+d
        molecule.dihdfunc.update({self.link:self})



if __name__=='__main__':
    import doctest
    doctest.testmod()
