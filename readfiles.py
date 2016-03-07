# From fchk read Charge, Multiplicity, Coordinates
import os,time
class readfchk(object):

    def __init__(self,file):
        self.xyzlist=[]
        self.atomlist=['0']
        print('Read fchk:',file)
        with open(file,'r') as f:
            while True:
                string=f.readline()
                if string.find('Charge')>=0:
                    self.charge=string.split('I')[1].strip(' \n')
                if string.find('Multiplicity')>=0:
                    self.spin=string.split('I')[1].strip(' \n')
                if string.find('Atomic numbers')>=0:
                    self.natoms=int(string.split('=')[1].strip(' \n'))
                    while True:
                        string=f.readline()
                        if string.find('Nuclear charges')>=0:
                            break
                        for x in string.split():
                            self.atomlist.append(x.strip(' '))

                if string.find('Current cartesian')>=0:
                    while True:
                        string=f.readline()
                        if string.find('Force Field')>=0:
                            break
                        for x in string.split():
                            self.xyzlist.append(x.strip(' '))
                    break
        self.xyzlist=[float(x) for x in self.xyzlist]
        self.atomlist=[int(x) for x in self.atomlist]
    def assignatom(self,geom):
        xyz=['']
        for i in range(0,len(self.xyzlist),3):
            xyz.append(self.xyzlist[i:i+3])
        for i in range(1,self.natoms+1):
            geom.addatom(self.atomlist[i],xyz[i])
class readac(object):
    def __init__(self,file):
        self.atomtype=['']
        self.charge=['']
        with open(file,'r') as f:
            string=f.readline()
            string=f.readline()
            for string in f.readlines():
                if string.find('BOND')>=0:
                    break
                ac=string.split()
                self.atomtype.append(ac[len(ac)-1])
                self.charge.append(ac[len(ac)-2])
    def assigntype(self,geom):
        for i in range(1,geom.natoms+1):
            geom.atoms[i].atomtype=self.atomtype[i]
    def assigncharge(self,geom):
        for i in range(1,geom.natoms+1):
            geom.atoms[i].charge=self.charge[i]

class readlog(object):
    def __init__(self,file):
        self.filename=file
        with open(file,'r') as f:
            for x in f.readlines():
                if x.find('NAtoms')>=0:
                    self.natoms=int(x.split()[1])
                    break


    def xyzlast(self):
        with open(self.filename,'r') as f:
            for x in list(reversed(f.readlines())):
                if x.find('orientation')>=0:
                    self.orn=x
                    break
        return self.orn
    def freq(self):
        freq=[]
        with open(self.filename,'r') as f:
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
        return freq


class gauFile(object):
    g09rt='g09boon'
    g09a2rt='g09a2boon'
    antechamber='antechamber -c resp'
    def __init__(self,name):
        self.name=name
        self.com=name+'.com'
        self.log=name+'.log'
        self.chk=name+'.chk'
        self.fchk=name+'.fchk'
        self.ac=name+'.ac'
    def rung09(self):
        ifchk=1
        with open(self.com,'r') as f:
            for line in f.readlines():
                if line.find('%chk')>=0:
                    ifchk=0
            if ifchk==1:
                f.seek(0)
                content=f.read()
        if ifchk==1:
            with open(self.com,'w') as f:
                f.write('%chk='+self.chk+'\n')
                f.write(content)
        print('Run g09 : '+gauFile.g09rt+' '+self.com)
        os.system(gauFile.g09rt+' '+self.com)
    def rung09a2(self):
        print('Run g09a2 : '+gauFile.g09a2rt+' '+self.com)
        os.system(gauFile.g09a2rt+' '+self.com)
    def isover(self):
        print('Waiting g09 for',self.com,'...')
        while True:
            output=''
            if not os.path.isfile(self.log):
                continue
            with open(self.log,'r') as f:
                for x in f.readlines()[:-6:-1]:
                    output+=x
            if output.find('Normal termination')>=0:
                print('    ..normal termination')
                return True
            if output.find('Error termination')>=0:
                print('Error termination in '+self.com)
                return False
            time.sleep(2)
    def getresp(self):
        print('Runing antechamber: \n')
        os.system(gauFile.antechamber+' -i '+self.log+' -fi gout -o '+self.ac+' -fo ac')
