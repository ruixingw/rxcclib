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
        self.__bondlist={}
        self.__anglelist={}
        self.__dihdlist={}
        self.bondfunc={}
        self.anglefunc={}
        self.dihdfunc={}
        self.__pointer=0
    @property
    def natoms(self):
        return len(self.__atomlist)-1
    def addatom(self,idorsym,xyz,unit='bohr'):
        self.__atomlist.append(Atom(self,idorsym,xyz,unit))
    def addbond(self,atomnum1,atomnum2):
        self.__bondlist.update({str(atomnum1)+','+str(atomnum2):Bond(self,atomnum1,atomnum2)})
    def addangle(self,atomnum1,atomnum2,atomnum3):
        self.__anglelist.update({str(atomnum1)+','+str(atomnum2)+','+str(atomnum3):Angle(self,atomnum1,atomnum2,atomnum3)})
    def adddihd(self,atomnum1,atomnum2,atomnum3,atomnum4):
        self.__dihdlist.update({str(atomnum1)+','+str(atomnum2)+','+str(atomnum3)+','+str(atomnum4):Dihd(self,atomnum1,atomnum2,atomnum3,atomnum4)})

    def atom(self,atomnum):
        return self.__atomlist[atomnum]
    def bond(self,atomnum1,atomnum2):
        if atomnum1>atomnum2:
            atomnum1,atomnum2=atomnum2,atomnum1
        return self.__bondlist[str(atomnum1)+','+str(atomnum2)]
    def angle(self,atomnum1,atomnum2,atomnum3):
        if atomnum1>atomnum3:
            atomnum1,atomnum3=atomnum3,atomnum1
        return self.__anglelist[str(atomnum1)+','+str(atomnum2)+','+str(atomnum3)]
    def dihd(self,atomnum1,atomnum2,atomnum3,atomnum4):
        if atomnum2>atomnum3:
            atomnum2,atomnum3=atomnum3,atomnum2
            atomnum1,atomnum4=atomnum4,atomnum1
        elif atomnum2==atomnum3:
            if atomnum1>atomnum4:
                atomnum1,atomnum4=atomnum4,atomnum1
        return self.__dihdlist[str(atomnum1)+','+str(atomnum2)+','+str(atomnum3)+','+str(atomnum4)]


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
        self.__pointer=0
        return self
    def __next__(self):
        self.__pointer+=1
        if self.__pointer <=self.natoms:
            return self[self.__pointer]
        else:
            raise StopIteration("Exceeded all atoms")


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
    >>> mole.atom(1).xyz
    array([ 0.,  0.,  0.])
    >>> setattr(mole.atom(1),'atomtype','c2')
    >>> mole.atom(1).atomtype
    'c2'
    '''
    bohr=0.5291772086 # bohr to angstrom
    __idtosym={1:'H',5:'B',6:'C',7:'N',8:'O',9:'F',13:'Al',14:'Si',15:'P',16:'S',17:'Cl',26:'Fe',28:'Ni',29:'Cu',30:'Zn'}
    __symtoid={v:k for k,v in __idtosym.items()}

    def __init__(self,mole,idorsym,xyz,unit='bohr'):  # molecule object,int,[float,float,float]
        assert isinstance(mole,molecule),"First argument must be a molecule object!. Use molecule.addatom method to avoid this problem."
        assert unit!='bohr' or unit!='angstrom', "Coordinate unit must be bohr or angstrom"
        self.__mymolecule=mole
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
            self.__xyz=xyz*Atom.bohr
        elif unit=='angstrom':
            self.__xyz=xyz

        self.__atomnum=mole.natoms+1
        return
    @property
    def xyz(self):
        return self.__xyz
    @property
    def name(self):
        return str(self.__atomsym)+str(self.__atomnum)
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
        return self.__mymolecule
    # @mymolecule.setter
    # def mymolecule(self,value):
    #     if not isinstance(value,molecule):
    #         print("Error when changing molecule of,",self.atomsym,self.atomnum,"from",self.molecule,"to",value,": Expected a molecule object, received a",type(value))
    #         quit()
    #     self.__molecule=value


class Bond(object):
    def __init__(self,mole,a,b): # self, atomid a, atomid b
        if a>b:
            a,b=b,a
        self.__a=mole[a]
        self.__b=mole[b]
        self.__vec=self.__a.xyz-self.__b.xyz

    @property
    def length(self):
        return np.linalg.norm(self.__vec)
class Angle(object):
    def __init__(self,mole,a,b,c):
        if a>c:
            a,c=c,a
        self.__a=mole[a]
        self.__b=mole[b]
        self.__c=mole[c]
        self.__ab=mole[a].xyz-mole[b].xyz
        self.__bc=mole[b].xyz-mole[c].xyz
    @property
    def anglevalue(self):
        v1u=self.__ab/np.linalg.norm(self.__ab)
        v2u=self.__bc/np.linalg.norm(self.__bc)
        angle=np.arccos(np.dot(v1u,v2u))*180.0/np.pi
        if np.isnan(angle):
            if(v1u==v2u).all():
                return 0.0
            else:
                return 180.0
        return angle
class Dihd(object):
    def __init__(self,mole,a,b,c,d):
        if b>c:
            b,c=c,b
            a,d=d,a
        elif b==c:
            if a>d:
                a,d=d,a
        self.__a=mole[a]
        self.__b=mole[b]
        self.__c=mole[c]
        self.__d=mole[d]
    @property
    def dihdvalue(self):
        v1=self.__a.xyz-self.__b.xyz
        v2=self.__b.xyz-self.__c.xyz
        v3=self.__c.xyz-self.__d.xyz
        v1u=v1/np.linalg.norm(v1)
        v2u=v1/np.linalg.norm(v2)
        v3u=v1/np.linalg.norm(v3)
        n1=np.cross(v1u,v2u)
        n2=np.cross(v2u,v3u)
        dihd=np.arccos(np.dot(n1,n2))*180.0/np.pi
        if np.isnan(dihd):
            if (n1==n2).all():
                return 0.0
            else:
                return 180.0
        return dihd

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
