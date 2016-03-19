from __future__ import print_function
import math

class geometry(object):

    def __init__(self):
        self.natoms=0
        self.atoms=[0]
        self.atomsnumber=0
        self.bonds={}
        self.angles={}
        self.dihedrals={}
        self.bondfunc={}
        self.anglefunc={}
        self.dihdfunc={}
    def addatom(self,elementid,xyz):
        atoms(self,elementid,xyz)
        self.natoms+=1
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

class atoms(object):

    bohr=0.5291772086 # bohr to angstrom
    idtoname={1:'H',5:'B',6:'C',7:'N',8:'O',9:'F',13:'Al',14:'Si',15:'P',16:'S',17:'Cl',26:'Fe',28:'Ni',29:'Cu',30:'Zn'}

    def __init__(self,geom,elementid,xyz):  # int,[float,float,float]
        self.elementid=elementid   # int
        self.elementname=atoms.idtoname[elementid] #str
        self.atomtype=''
        self.charge=''  # char
        self.x=xyz[0]*atoms.bohr
        self.y=xyz[1]*atoms.bohr
        self.z=xyz[2]*atoms.bohr
        geom.atomsnumber+=1
        self.atomid=geom.atomsnumber
        geom.atoms.append(self)
    def xyz(self):
        return [self.x,self.y,self.z]

class vector(object):
    def __init__(self,x,y,z):
        self.x=x
        self.y=y
        self.z=z
    def length(self):
        return (self.x**2+self.y**2+self.z**2)**(0.5)

class bonds(object):
    def __init__(self,geom,a,b): # self, atomid a, atomid b
        if a>b:
            a,b=b,a
        self.a=geom.atoms[a]
        self.b=geom.atoms[b]
        geom.bonds.update({str(a)+','+str(b):self})
        self.vec=vector(self.a.x-self.b.x,self.a.y-self.b.y,self.a.z-self.b.z)
    def length(self):
        return self.vec.length()
class angles(object):
    def __init__(self,geom,a,b,c):
        if a>c:
            a,c=c,a
        self.a=geom.atoms[a]
        self.b=geom.atoms[b]
        self.c=geom.atoms[c]
        self.ab=vector(self.a.x-self.b.x,self.a.y-self.b.y,self.a.z-self.b.z)
        self.bc=vector(self.b.x-self.c.x,self.b.y-self.c.y,self.b.z-self.c.z)
        geom.angles.update({str(a)+','+str(b)+','+str(c):self})

    def angle(self):
        innerproduct=self.ab.x*self.bc.x+self.ab.y*self.bc.y+self.ab.z*self.bc.z
        product=self.ab.length()*self.bc.length()
        cos=innerproduct/product
        angle=180-math.acos(cos)*180/math.pi
        return angle
class dihedrals(object):
    def __init__(self,geom,a,b,c,d):
        if b>c:
            b,c=c,b
            a,d=d,a
        elif b==c:
            if a>d:
                a,d=d,a
        self.a=geom.atoms[a]
        self.b=geom.atoms[b]
        self.c=geom.atoms[c]
        self.d=geom.atoms[d]

        geom.dihedrals.update({str(a)+','+str(b)+','+str(c)+','+str(d):self})


class bondfunc(object):

    def __init__(self,geom,bondobj):
        a=bondobj.a.atomtype
        b=bondobj.b.atomtype
        if a>b:
            a,b=b,a
        bondobj.func=self
        self.link=a+' '+b
        geom.bondfunc.update({self.link:self})
class anglefunc(object):

    def __init__(self,geom,angleobj):
        a=angleobj.a.atomtype
        b=angleobj.b.atomtype
        c=angleobj.c.atomtype
        if a>c:
            a,c=c,a
        angleobj.func=self
        self.link=a+' '+b+' '+c
        geom.anglefunc.update({self.link:self})
class dihdfunc(object):

    def __init__(self,geom,dihdobj):
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
        geom.dihdfunc.update({self.link:self})
