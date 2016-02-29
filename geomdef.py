from __future__ import print_function
import math

def addatom(elementid,xyz):
    at.append(atoms(elementid,xyz))
    return len(at)-1
def addbond(atomid1,atomid2):
    bonds(atomid1,atomid2)
def addangle(atomid1,atomid2,atomid3):
    angles(atomid1,atomid2,atomid3)
def adddihd(atomid1,atomid2,atomid3,atomid4):
    dihedrals(atomid1,atomid2,atomid3,atomid4)
def swap(a,b):
    return b,a

class atoms(object):
    number=0
    bohr=0.5291772086 # bohr to angstrom
    idtoname={'1':'H','5':'B','6':'C','7':'N','8':'O','9':'F','13':'Al','14':'Si','15':'P','16':'S','17':'Cl','26':'Fe','28':'Ni','29':'Cu','30':'Zn'}

    def __init__(self,elementid,xyz):  # int,[float,float,float]
        self.elementid=elementid   # int
        self.elementname=atoms.idtoname[str(elementid)] #str
        self.atomtype=''
        self.charge=''  # char
        self.x=xyz[0]*atoms.bohr
        self.y=xyz[1]*atoms.bohr
        self.z=xyz[2]*atoms.bohr
        atoms.number+=1
        self.atomid=atoms.number
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
    list={}
    number=0
    def __init__(self,a,b): # self, atomid a, atomid b
        if a>b:
            a,b=swap(a,b)
        self.a=at[a]
        self.b=at[b]
        self.bf=''
        bonds.number+=1
        self.number=bonds.number
        bonds.list.update({str(a)+','+str(b):self})
        bonds.list.update({str(b)+','+str(a):self})
        self.vec=vector(self.a.x-self.b.x,self.a.y-self.b.y,self.a.z-self.b.z)
    def length(self):
        return self.vec.length()
class angles(object):
    list={}
    number=0
    def __init__(self,a,b,c):
        if a>c:
            a,c=swap(a,c)
        self.a=at[a]
        self.b=at[b]
        self.c=at[c]
        self.af=''
        self.ab=vector(self.a.x-self.b.x,self.a.y-self.b.y,self.a.z-self.b.z)
        self.bc=vector(self.b.x-self.c.x,self.b.y-self.c.y,self.b.z-self.c.z)
        angles.number+=1
        self.number=angles.number
        angles.list.update({str(a)+','+str(b)+','+str(c):self})
        angles.list.update({str(c)+','+str(b)+','+str(a):self})

    def angle(self):
        innerproduct=self.ab.x*self.bc.x+self.ab.y*self.bc.y+self.ab.z*self.bc.z
        product=self.ab.length()*self.bc.length()
        cos=innerproduct/product
        angle=180-math.acos(cos)*180/math.pi
        return angle
class dihedrals(object):
    list={}
    number=0
    def __init__(self,a,b,c,d):
        if b>c:
            b,c=swap(b,c)
            a,d=swap(a,d)
        elif b==c:
            if a>d:
                a,d=swap(a,d)
        self.a=at[a]
        self.b=at[b]
        self.c=at[c]
        self.d=at[d]
        self.df=''
        dihedrals.number+=1
        dihedrals.list.update({str(a)+','+str(b)+','+str(c)+','+str(d):self})
        dihedrals.list.update({str(d)+','+str(c)+','+str(b)+','+str(a):self})

class bondfunc(object):
    list={}
    def __init__(self,bondobj):
        a=bondobj.a.atomtype
        b=bondobj.b.atomtype
        if a>b:
            a,b=swap(a,b)
        bondobj.func=self
        self.link=a+' '+b
        bondfunc.list.update({self.link:self})
class anglefunc(object):
    list={}
    def __init__(self,angleobj):
        a=angleobj.a.atomtype
        b=angleobj.b.atomtype
        c=angleobj.c.atomtype
        if a>c:
            a,c=swap(a,c)
        angleobj.func=self
        self.link=a+' '+b+' '+c
        anglefunc.list.update({self.link:self})
class dihdfunc(object):
    list={}
    def __init__(self,dihdobj):
        a=dihdobj.a.atomtype
        b=dihdobj.b.atomtype
        c=dihdobj.c.atomtype
        d=dihdobj.d.atomtype
        if b>c:
            a,d=swap(a,d)
            b,c=swap(b,c)
        elif b==c:
            if a>d:
                a,d=swap(a,d)
        dihdobj.func=self
        self.link=a+' '+b+' '+c+' '+d
        dihdfunc.list.update({self.link:self})




at=['0']
xyz=['']
