# From fchk read Charge, Multiplicity, Coordinates
from __future__ import print_function
import os,time
import numpy as np
from io import StringIO
import rx.molecules as rxmol

class rxccError(Exception):
    def __init(self,value):
        self.value=value
    def __repr__(self):
        return repr(self.value)

class File(object):

    def __init__(self,name):
        pwd=os.path.abspath('.')
        self.__name=os.path.join(pwd,name)
        self.comname=self.__name+'.com'
        self.logname=self.__name+'.log'
        self.chkname=self.__name+'.chk'
        self.fchkname=self.__name+'.fchk'
        self.acname=self.__name+'.ac'
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
                raise rxccError("Error: natoms is already read, and not consistent with new value: Now is "+str(self.__natoms)+", New value is "+str(value))
    @property
    def multiplicity(self):
        return self.__mlpty
    @multiplicity.setter
    def multiplicity(self,value):
        if self.__mlpty==None:
            self.__mlpty=value
        else:
            if self.__mlpty!=value:
                raise rxccError("Error: multiplicity is already read, and not consistent with new value: Now is "+str(self.__mlpty)+", New value is "+str(value))
    @property
    def totalcharge(self):
        return self.__totalcharge
    @totalcharge.setter
    def totalcharge(self,value):
        if self.__totalcharge==None:
            self.__totalcharge=value
        else:
            if self.__totalcharge!=value:
                raise rxccError("Error: totalcharge is already read, and not consistent with new value: Now is "+str(self.__totalcharge)+", New value is "+str(value))
    @property
    def xyzfile(self):
        souc='fchk'
        if self.default!='fchk':
            souc=self.default
        if souc=='fchk':
            return self.__fchk.xyz
        elif souc=='com':
            return self.__com.xyz
    @property
    def atomtypelist(self):
        return self.__ac.atomtypelist
    @property
    def atomchargelist(self):
        return self.__ac.atomchargelist


    # fetch data as function
    def find33Hessian(self,i,j):
        return self.fchk.return33Hessian(i,j)


    # Parse
    def readfchk(self):
        if not self.fchk.read():
            raise rxccError("Error in reading fchk:"+self.fchkname)

    def readac(self):
        if not self.ac.read():
            raise rxccError("Error in reading ac:"+self.acname)

    # File operation
    def rung09(self):
        state=self.com.rung09()

    def rung09a2(self):
        state=self.com.rung09a2()

    def isover(self):
        state=self.com.isover()

    def runformchk(self):
        string='formchk '+self.chkname+' '+self.fchkname
        print('  ',string)
        iferror=os.popen(string)
        if iferror.read().find('Error')>=0:
            raise rxccError('   Error in formatting'+self.chkname)
            iferror.close()
            return False
        iferror.close()
        return True
    def runantecham(self):
        self.log.runantecham()





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
            raise rxccError('Hessian has not been read from:'+self.__father.fchkname)
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
                self.__charge.append(float(ac[len(ac)-2]))
        return True
    @property
    def atomtypelist(self):
        return self.__atomtype
    @property
    def atomchargelist(self):
        return self.__charge

class gauCOM(object):
    g09rt='g09'
    g09a2rt='g09'
    def __init__(self,father):
        self.__father=father
        self.__xyzfile=''
        self.__atomlist=[None]
        self.__atomtypelist=[None]
        self.__atomchargelist=[None]
        self.__coordslist=[]
        self.__connectivity=''
        self.__dihdfunc=[]
        self.__anglefunc=[]
        self.__bondfunc=[]
        self.additionfunc=[]
        self.nozomuvdw=[]
        self.__xyz=''
        self.commandline=''
    @property
    def connectivity(self):
        return self.__connectivity
    @property
    def atomlist(self):
        return self.__atomlist
    @property
    def atomtypelist(self):
        return self.__atomtypelist
    @property
    def atomchargelist(self):
        return self.__atomchargelist
    @property
    def xyz(self):
        return self.__xyz
    @property
    def nozomudihdfunc(self):
        return self.__dihdfunc
    @property
    def nozomuanglefunc(self):
        return self.__anglefunc
    @property
    def nozomubondfunc(self):
        return self.__bondfunc

    # Parse
    def read(self):
        with open(self.__father.comname,'r') as f:
            counter=0
            ifconnect=False
            ifamber=False
            line=''
            for line in f:
                if line=='\n':
                    counter+=1
                    continue
                if counter==0:
                    while True:
                        self.commandline+=line
                        line=next(f)
                        if line=='\n':
                            counter+=1
                            break
                    if self.commandline.find('connectivity')>=0:
                        ifconnect=True
                    if self.commandline.find('amber')>=0:
                        ifamber=True
                def molespecs(line):
                    self.__xyzfile+=line
                    tmp=line.split()[0]
                    if tmp.find('-')>=0:
                        self.__atomlist.append(tmp.split('-')[0])
                        if tmp.count('-')==2:
                            tmp=tmp.split('-')
                            self.__atomtypelist.append(tmp[1])
                            self.__atomchargelist.append(tmp[2])
                        elif tmp.count('-')==3:
                            tmp=tmp.split('-')
                            self.__atomtypelist.append(tmp[1])
                            self.__atomchargelist.append(-float(tmp[3]))
                    else:
                        self.__atomlist.append(tmp)
                    self.__coordslist.extend(line.split()[1:4])

                if counter==2:
                    self.__father.multiplicity=int(line.split()[1])
                    self.__father.totalcharge=int(line.split()[0])
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

                if counter==4:
                    if ifamber:
                        line=next(f)
                        while counter==4:
                            if line=='\n':
                                counter+=1
                                break
                            thisline=mmfunction(line)
                            if thisline.type=='dihd':
                                self.__dihdfunc.append(thisline)
                            elif thisline.type=='angle':
                                self.__anglefunc.append(thisline)
                            elif thisline.type=='bond':
                                self.__bondfunc.append(thisline)
                            elif thisline.type=='else':
                                self.additionfunc.append(thisline)
                            elif thisline.type=='vdw':
                                self.nozomuvdw.append(thisline)
                            line=next(f)




        self.__coordslist=np.array(self.__coordslist)
        for i in range(0,len(self.__atomlist)-1):
            tmp=str(self.__atomlist[i+1])+'   '+str(self.__coordslist[3*i])+'   '+str(self.__coordslist[3*i+1])+'   '+str(self.__coordslist[3*i+2])+'\n'
            self.__xyz+=tmp


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
        with open(self.__father.logname,'r') as f:
            for x in f.readlines():
                if x.find('NAtoms')>=0:
                    self.natoms=int(x.split()[1])
                    return self.natoms
                    break

    def coordslast(self): #Uncomplete
        with open(self.__father.logname,'r') as f:
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
    def runantecham(self):
        print('Runing antechamber: \n')
        command=gauLOG.antecommand+' -i '+self.__father.logname+' -fi gout -o '+self.__father.acname+' -fo ac'
        print(command)
        os.system(command)
class mmfunction(object):
    magicnum='XXXXXX'
    def __init__(self,line):
        tmp=line.split()
        fun=[]
        self.type=None
        self.repr=None
        for item in tmp:
            fun.append(item)
        def newfloat(value):
            if value=='XXXXXX':
                return mmfunction.magicnum
            else:
                return float(value)
        if fun[0]=='AmbTrs':
            self.type='dihd'
            self.a=fun[1]
            self.b=fun[2]
            self.c=fun[3]
            self.d=fun[4]
            self.npaths=float(fun[13])
            for i,paras in enumerate(fun[9:13]):
                if newfloat(paras)!=0.000:
                    i+=1
                    self.value=newfloat(paras)
                    break
            self.periodicity=i
            self.phase=int(fun[4+self.periodicity])
            self.repr=self.a+' '+self.b+' '+self.c+' '+self.d
        elif fun[0]=='HrmBnd1':
            self.type='angle'
            self.a=fun[1]
            self.b=fun[2]
            self.c=fun[3]
            self.value=newfloat(fun[4])
            self.eqvalue=newfloat(fun[5])
            self.repr=self.a+' '+self.b+' '+self.c
        elif fun[0]=='HrmStr1':
            self.type='bond'
            self.a=fun[1]
            self.b=fun[2]
            self.value=newfloat(fun[3])
            self.eqvalue=newfloat(fun[4])
            self.repr=self.a+' '+self.b
        # elif fun[0]=='ImpTrs':
        #     self.type='improper'
        #     self.a=fun[1]
        #     self.b=fun[2]
        #     self.c=fun[3]
        #     self.d=fun[4]
        #     self.value=newfloat(fun[5])
        #     self.eqvalue=newfloat(fun[6])
        #     self.repr=self.type+self.a+self.b+self.c+self.d
        elif fun[0]=='VDW':
            self.type='vdw'
            self.content=line
        else:
            self.type='else'
            self.content=line
