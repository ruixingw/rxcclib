# From fchk read Charge, Multiplicity, Coordinates
from __future__ import print_function
import os,time


class File(object):

    def __init__(self,name):
        self.__name=name
        self.__comname=name+'.com'
        self.__logname=name+'.log'
        self.__chkname=name+'.chk'
        self.__fchkname=name+'.fchk'
        self.__acname=name+'.ac'
        self.__com=gauCOM(self)
        self.__log=gauLOG(self)
        self.__fchk=gauFCHK(self)
        self.__ac=amberAC(self)
    def getname(self):
        return self.__name
    def comname(self):
        return self.__comname
    def logname(self):
        return self.__logname
    def chkname(self):
        return self.__chkname
    def fchkname(self):
        return self.__fchkname
    def acname(self):
        return self.__acname
    def com(self):
        return self.__com
    def log(self):
        return self.__log
    def fchk(self):
        return self.__fchk
    def ac(self):
        return self.__ac
    def assignatomtogeom(self,geom):
        self.fchk().assignatom(geom)
    def assignchargetogeom(self,geom):
        state=self.ac().assigncharge(geom)
        return state
    def assigntypetogeom(self,geom):
        state=self.ac().assigntype(geom)
        return state
    def readfchk(self):
        if not self.fchk().read():
            print("Error in reading fchk:",self.fchkname())
            quit()
    def readac(self):
        state=self.ac().readAC()
        return state
    def getnatoms(self):
        return self.fchk().getnatoms()
    def readHessianFromFchk(self):
        state=self.fchk().readhessian()
        if not state:
            print("Error in reading hessian from:",self.fchkname())
        return state
    def find33Hessian(self,i,j):
        return self.fchk().find33Hessian(i,j)
    def rung09(self):
        state=self.com().rung09()
        return state
    def rung09a2(self):
        state=self.com().rung09a2()
        return state
    def isover(self):
        state=self.com().isover()
        return state
    def gethessian(self):
        return self.fchk().gethessian()
    def formchk(self):
        string='formchk '+self.chkname()+' '+self.fchkname()
        print('  ',string)
        iferror=os.popen(string)
        if iferror.read().find('Error')>=0:
            print('   Error in formatting',self.chkname())
            return False
        return True
    def getcharge(self):
        return self.fchk().getcharge()
    def getspin(self):
        return self.fchk().getspin()
    def antecham(self):
        print(self)
        print(self.ac())
        self.ac().antecha()


class gauFCHK(File):
    def __init__(self,father):
        self.__father=father
        self.__filename=self.__father.getname()+'.fchk'
        self.__xyzlist=[]
        self.__atomlist=['0']
        self.__hessian=False
        self.__charge=''
        self.__spin=''
        self.__natoms=0
    def read(self):
        print('Read fchk:',self.__father.fchkname())
        with open(self.__father.fchkname(),'r') as f:
            while True:
                string=f.readline()
                if string.find('Charge')>=0:
                    self.__charge=string.split('I')[1].strip(' \n')
                if string.find('Multiplicity')>=0:
                    self.__spin=string.split('I')[1].strip(' \n')
                if string.find('Atomic numbers')>=0:
                    self.__natoms=int(string.split('=')[1].strip(' \n'))
                    while True:
                        string=f.readline()
                        if string.find('Nuclear charges')>=0:
                            break
                        for x in string.split():
                            self.__atomlist.append(x.strip(' '))

                if string.find('Current cartesian')>=0:
                    while True:
                        string=f.readline()
                        if string.find('Force Field')>=0:
                            break
                        for x in string.split():
                            self.__xyzlist.append(x.strip(' '))
                    break
        self.__xyzlist=[float(x) for x in self.__xyzlist]
        self.__atomlist=[int(x) for x in self.__atomlist]
        return True
    def assignatom(self,geom):
        xyz=['']
        for i in range(0,len(self.__xyzlist),3):
            xyz.append(self.__xyzlist[i:i+3])
        for i in range(1,self.__natoms+1):
            geom.addatom(self.__atomlist[i],xyz[i])
    def getnatoms(self):
        return self.__natoms
    def xyzlist(self):
        return self.__xyzlist
    def readhessian(self):
        if self.__hessian:
            print("Hessian is already read from:",self.name)
            return False
        else:
            hess=[]
            with open(self.__father.fchkname(),'r') as f:
                ifread=False
                for line in f.readlines():
                    if line.find('Dipole')>=0:
                        ifread=False
                    if ifread:
                        hess.extend([float(x) for x in line.split()])
                    if line.find('Cartesian Force Constants')>=0:
                        ifread=True
            self.__hessian=hess
            return True

    def gethessian(self):
        if self.__hessian==False:
            print('Hessian has not been read from:',self.__father.fchkname())
            return False
        else:
            return self.__hessian
    def getcharge(self):
        return self.__charge
    def getspin(self):
        return self.__spin
    def findHessianElement(self,i,j): #i, j: coordinate number
        if i<j:
            i,j=j,i
        num=i*(i-1)/2+j
        num=int(num)
        return self.__hessian[num-1]
    def find33Hessian(self,i,j):    #i, j: atom number
        if i<j:
            i,j=j,i
        tthess=[]
        i1=3*(i-1)+1
        i2=3*(i-1)+2
        i3=3*i
        j1=3*(j-1)+1
        j2=3*(j-1)+2
        j3=3*j
        tthess.append([self.findHessianElement(i1,j1),self.findHessianElement(i1,j2),self.findHessianElement(i1,j3)])
        tthess.append([self.findHessianElement(i2,j1),self.findHessianElement(i2,j2),self.findHessianElement(i2,j3)])
        tthess.append([self.findHessianElement(i3,j1),self.findHessianElement(i3,j2),self.findHessianElement(i3,j3)])
        return tthess

class amberAC(File):
    antecommand='antechamber -c resp'
    def __init__(self,father):
        self.__father=father
        self.__atomtype=['']
        self.__charge=['']
    def readAC(self):
        with open(self.__father.acname(),'r') as f:
            string=f.readline()
            string=f.readline()
            for string in f.readlines():
                if string.find('BOND')>=0:
                    break
                ac=string.split()
                self.__atomtype.append(ac[len(ac)-1])
                self.__charge.append(ac[len(ac)-2])
    def assigntype(self,geom):   #!!!!!!!!!!!
        for i in range(1,geom.natoms+1):
            geom.atoms[i].atomtype=self.__atomtype[i]
    def assigncharge(self,geom): #!!!!!!!!!!!!!
        for i in range(1,geom.natoms+1):
            geom.atoms[i].charge=self.__charge[i]
    def antecha(self):
        print('Runing antechamber: \n')
        command=amberAC.antecommand+' -i '+self.__father.logname()+' -fi gout -o '+self.__father.acname()+' -fo ac'
        print(command)
        os.system(command)

class gauCOM(File):
    g09rt='g09boon'
    g09a2rt='g09a2boon'
    def __init__(self,father):
        self.__father=father

    def rung09(self):
        ifchk=1 # if no chk, add.
        with open(self.__father.comname(),'r') as f:
            for line in f.readlines():
                if line.find('%chk')>=0:
                    line='%chk='+self.__father.chkname()+'\n'
                    ifchk=0
            if ifchk==1:
                f.seek(0)
                content=f.read()
        if ifchk==1:
            with open(self.__father.comname(),'w') as f:
                f.write('%chk='+self.__father.chkname()+'\n')
                f.write(content)
        print('Run g09 : '+gauCOM.g09rt+' '+self.__father.comname())
        os.system(gauCOM.g09rt+' '+self.__father.comname())

    def rung09a2(self):
        ifchk=1 # if no chk, add.
        with open(self.__father.comname(),'r') as f:
            for line in f.readlines():
                if line.find('%chk')>=0:
                    line='%chk='+self.__father.chkname()+'\n'
                    ifchk=0
            if ifchk==1:
                f.seek(0)
                content=f.read()
        if ifchk==1:
            with open(self.__father.comname(),'w') as f:
                f.write('%chk='+self.__father.chkname()+'\n')
                f.write(content)
        print('Run g09a2 : '+gauCOM.g09a2rt+' '+self.__father.comname())
        os.system(gauCOM.g09a2rt+' '+self.__father.comname())

    def isover(self):
        print('Waiting g09 for',self.__father.comname(),'...')
        while True:
            output=''
            if not os.path.isfile(self.__father.logname()):
                continue
            with open(self.__father.logname(),'r') as f:
                for x in f.readlines()[:-6:-1]:
                    output+=x
            if output.find('Normal termination')>=0:
                print('    ..normal termination')
                return True
            if output.find('Error termination')>=0:
                print('Error termination in '+self.__father.comname())
                return False
            time.sleep(2)
class gauLOG(File):
    def __init__(self,father):
        self.__father=father
        self.__freq=0
    def getnatoms(self):
        with open(self.__father.getlogname(),'r') as f:
            for x in f.readlines():
                if x.find('NAtoms')>=0:
                    self.natoms=int(x.split()[1])
                    break

    def xyzlast(self): #Uncomplete
        with open(self.__father.getlogname(),'r') as f:
            for x in list(reversed(f.readlines())):
                if x.find('orientation')>=0:
                    self.orn=x
                    break
        return self.orn
    def getfreq(self):
        freq=[]
        with open(self.__father.logname(),'r') as f:
            for line in f.readlines():
                if line.find('Frequencies -- ')>=0:
                    freq.extend(line.split()[2:])
            if freq!=[]:
                freq=list(map(float,freq))
                # if (len(freq)>3*self.natoms-5):
                #     freq=freq[:len(freq)/2]
            else:
                print("No Freq found")
                return ['No Freq found']
        self.__freq=freq
        return self.__freq
