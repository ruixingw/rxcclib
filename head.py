# From fchk read Charge, Multiplicity, Coordinates

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
