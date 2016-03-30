# From fchk read Charge, Multiplicity, Coordinates
from __future__ import print_function
import os,time
import numpy as np

class chemfileParseError(Exception):
    def __init(self,value):
        self.value=value
    def __repr__(self):
        return repr(self.value)

class File(object):

    def __init__(self,name):
        self.__name=name
        self.comname=name+'.com'
        self.logname=name+'.log'
        self.chkname=name+'.chk'
        self.fchkname=name+'.fchk'
        self.acname=name+'.ac'
        self.__com=gauCOM(self)
        self.__log=gauLOG(self)
        self.__fchk=gauFCHK(self)
        self.__ac=amberAC(self)
        self.default='fchk'
        self.__natoms=None
        self.__mlpty=None
        self.__totalcharge=None
    # molecule name
    @property
    def name(self):
        return self.__name

    # subfile objects
    @property
    def com(self):
        return self.__com
    @property
    def log(self):
        return self.__log
    @property
    def fchk(self):
        return self.__fchk
    @property
    def ac(self):
        return self.__ac

    # molecule properties; setter should only be called by subfile methods
    @property
    def natoms(self):
        return self.__natoms
    @natoms.setter
    def natoms(self,value):
        if self.__natoms==None:
            self.__natoms=value
        else:
            if self.__natoms!=value:
                raise chemfileParseError("Error: natoms is already read, and not consistent with new value.")
    @property
    def multiplicity(self):
        return self.__mlpty
    @multiplicity.setter
    def multiplicity(self,value):
        if self.__mlpty==None:
            self.__mlpty=value
        else:
            if self.__mlpty!=value:
                raise chemfileParseError("Error: multiplicity is already read, and not consistent with new value")
    @property
    def totalcharge(self):
        return self.__totalcharge
    @totalcharge.setter
    def totalcharge(self,value):
        if self.__totalcharge==None:
            self.__totalcharge=value
        else:
            if self.__totalcharge!=value:
                raise chemfileParseError("Error: totalcharge is already read, and not consistent with new value")
    @property
    def xyzfile(self):
        souc='fchk'
        if self.default!='fchk':
            souc=self.default
        if souc=='fchk':
            return self.__fchk.xyz


    # fetch data as function
    def find33Hessian(self,i,j):
        return self.fchk.return33Hessian(i,j)


    # Parse
    def readfchk(self):
        if not self.fchk.read():
            raise chemfileParseError("Error in reading fchk:"+self.fchkname)

    def readac(self):
        if not self.ac.read():
            raise chemfileParseError("Error in reading ac:"+self.acname)

    # File operation
    def rung09(self):
        state=self.com.rung09()

    def rung09a2(self):
        state=self.com.rung09a2()

    def isover(self):
        state=self.com.isover()

    def formchk(self):
        string='formchk '+self.chkname+' '+self.fchkname
        print('  ',string)
        iferror=os.popen(string)
        if iferror.read().find('Error')>=0:
            print('   Error in formatting',self.chkname)
            iferror.close()
            return False
        iferror.close()
        return True
    def antecham(self):
        print(self.ac())
        self.ac().antecha()






class gauFCHK(object):
    def __init__(self,father):
        self.__father=father
        self.__filename=self.__father.name+'.fchk'
        self.__coordslist=[]
        self.__atomlist=[None]
        self.__readstate=None
        self.__totalcharge=None
        self.__mlpty=None
        self.__natoms=None
        self.__hessian=[]
        self.__xyz=''
    def read(self):
        if self.__readstate==True:
            print("Warning in fchk.read(): already read")
            return True
        print('Read fchk:',self.__father.fchkname)
        # FCHK parser
        with open(self.__father.fchkname,'r') as f:
            string=next(f)
            for string in f:
                if string.find('Charge')==0:
                    self.__totalcharge=int(string.split('I')[1])
                    self.__father.totalcharge=self.__totalcharge
                if string.find('Multiplicity')==0:
                    self.__mlpty=int(string.split('I')[1])
                    self.__father.multiplicity=self.__mlpty
                if string.find('Atomic numbers')==0:
                    self.__natoms=int(string.split('=')[1])
                    self.__father.natoms=self.__natoms
                    string=next(f)
                    while string.find('Nuclear charges')<0:
                        self.__atomlist.extend([int(x) for x in string.split()])
                        string=next(f)

                if string.find('Current cartesian coordinates')==0:
                    string=next(f)
                    while string.find('Force Field')<0:
                        self.__coordslist.extend([float(x) for x in string.split()])
                        string=next(f)
                # Read Hessian
                if string.find('Cartesian Force Constants')==0:
                    string=next(f)
                    while string.find('Dipole')<0:
                        self.__hessian.extend([float(x) for x in string.split()])
                        string=next(f)
                #Stop
        self.__coordslist=np.array(self.__coordslist)
        self.__atomlist=np.array(self.__atomlist)
        self.__hessian=np.array(self.__hessian)
        for i in range(0,len(self.__atomlist)-1):
            tmp=str(self.__atomlist[i+1])+'   '+str(self.__coordslist[3*i])+'   '+str(self.__coordslist[3*i+1])+'   '+str(self.__coordslist[3*i+2])+'\n'
            self.__xyz+=tmp

        return True
    @property
    def xyz(self):
        return self.__xyz
    @property
    def coordslist(self):
        return self.__coordslist
    @property
    def hessian(self):
        if self.__hessian==False:
            raise chemfileParseError('Hessian has not been read from:'+self.__father.fchkname)
        else:
            return self.__hessian
    @property
    def hessian(self):
        return self.__hessian
    def findHessianElement(self,i,j): #i, j: coordinate number
        if i<j:
            i,j=j,i
        num=i*(i-1)/2+j
        num=int(num)
        return self.__hessian[num-1]
    def return33Hessian(self,i,j):    #i, j: atom number
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
        tthess=np.array(tthess)
        return tthess

class amberAC(object):

    def __init__(self,father):
        self.__father=father
        self.__atomtype=[None]
        self.__charge=[None]

    def read(self):
        with open(self.__father.acname,'r') as f:
            string=f.readline()
            string=f.readline()
            for string in f.readlines():
                if string.find('BOND')>=0:
                    break
                ac=string.split()
                self.__atomtype.append(ac[len(ac)-1])
                self.__charge.append(ac[len(ac)-2])
    @property
    def atomtypelist(self):
        return self.__atomtype
    @property
    def atomchargelist(self):
        return self.__charge

class gauCOM(object):
    g09rt='g09boon'
    g09a2rt='g09a2boon'
    def __init__(self,father):
        self.__father=father
        self.__xyzfile=''
        self.__atomlist=[None]
        self.__coordslist=[]
        self.__connectivity=''
    # Parse
    def read(self):
        with open(self.__father.comname,'r') as f:
            line=next(f)
            counter=0
            command=''
            ifconnect=False
            for line in f:
                if line=='\n':
                    counter+=1
                    continue
                if line.find('#')>=0:
                    while True:
                        command+=line
                        line=next(f)
                        if line=='\n':
                            counter+=1
                            break
                    if command.find('connectivity')>=0:
                        ifconnect=True
                def molespecs(line):
                    self.__xyzfile+=line
                    tmp=line.split()[0]
                    if tmp.find('-')>=0:
                        self.__atomlist.append(tmp.split('-')[0])
                    else:
                        self.__atomlist.append(tmp)
                    self.__coordslist.extend(line.split()[1:4])

                if counter==2:
                    self.__father.multiplicity=int(line.split()[0])
                    self.__father.totalcharge=int(line.split()[1])
                    print("mlpty=",self.__father.multiplicity)
                    print('charge=',self.__father.totalcharge)
                    while counter==2:
                        line=next(f)
                        if line=='\n':
                            counter+=1
                            break
                        molespecs(line)
                def connect(line):
                    self.__connectivity+=line
                if counter==3:
                    if ifconnect:
                        line=next(f)
                        while counter==3:
                            if line=='\n':
                                counter+=1
                                break
                            connect(line)
                            line=next(f)

        self.__coordslist=np.array(self.__coordslist)
        print(self.__atomlist)
        print(self.__coordslist)
        print(self.__connectivity)



    # File operation
    def rung09(self):
        ifchk=1 # if no chk, add.
        with open(self.__father.comname,'r') as f:
            for line in f.readlines():
                if line.find('%chk')>=0:
                    line='%chk='+self.__father.chkname+'\n'
                    ifchk=0
            if ifchk==1:
                f.seek(0)
                content=f.read()
        if ifchk==1:
            with open(self.__father.comname,'w') as f:
                f.write('%chk='+self.__father.chkname+'\n')
                f.write(content)
        print('Run g09 : '+gauCOM.g09rt+' '+self.__father.comname)
        os.system(gauCOM.g09rt+' '+self.__father.comname)
    def rung09a2(self):
        ifchk=1 # if no chk, add.
        with open(self.__father.comname,'r') as f:
            for line in f.readlines():
                if line.find('%chk')>=0:
                    line='%chk='+self.__father.chkname+'\n'
                    ifchk=0
            if ifchk==1:
                f.seek(0)
                content=f.read()
        if ifchk==1:
            with open(self.__father.comname,'w') as f:
                f.write('%chk='+self.__father.chkname+'\n')
                f.write(content)
        print('Run g09a2 : '+gauCOM.g09a2rt+' '+self.__father.comname)
        os.system(gauCOM.g09a2rt+' '+self.__father.comname)


    def isover(self):
        print('Checking g09 termination for',self.__father.comname,'...')
        while True:
            output=''
            if not os.path.isfile(self.__father.logname):
                continue
            with open(self.__father.logname,'r') as f:
                for x in f.readlines()[:-6:-1]:
                    output+=x
            if output.find('Normal termination')>=0:
                print('    ..normal termination')
                return True
            if output.find('Error termination')>=0:
                print('Error termination in '+self.__father.comname)
                return False
            time.sleep(2)
class gauLOG(object):
    antecommand='antechamber -c resp'
    def __init__(self,father):
        self.__father=father
        self.__freq=0
    def getnatoms(self):
        with open(self.__father.getlogname,'r') as f:
            for x in f.readlines():
                if x.find('NAtoms')>=0:
                    self.natoms=int(x.split()[1])
                    break

    def coordslast(self): #Uncomplete
        with open(self.__father.getlogname,'r') as f:
            for x in list(reversed(f.readlines())):
                if x.find('orientation')>=0:
                    self.orn=x
                    break
        return self.orn
    def getfreq(self):
        freq=[]
        with open(self.__father.logname,'r') as f:
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
    def runantechm(self):
        print('Runing antechamber: \n')
        command=gauLOG.antecommand+' -i '+self.__father.logname+' -fi gout -o '+self.__father.acname+' -fo ac'
        print(command)
        os.system(command)
