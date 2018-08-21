#!/usr/bin/env python3
import pymatgen as mg
import shutil
import os

class moveSlab():
    def __init__(self, filename):
        t1 = filename.find('POSCAR')
        t2 = filename.find('poscar')
        t3 = filename.find('vasp')
        if t1 >= 0 or t2 >= 0 or t3 >= 0:
            if filename != 'POSCAR':
                shutil.copy(filename, 'POSCAR')
                self.structure = mg.Structure.from_file('POSCAR')
                os.remove('POSCAR')
            else:
                self.structure = mg.Structure.from_file(filename)

        self.lattice = self.structure.lattice

    def group(self,func=1):
        if func==1:
            self.up = sorted(
                filter(lambda x:x.frac_coords[-1]>0.5, self.structure),
                key=lambda x: x.frac_coords[-1])
            self.down = sorted(
                filter(lambda x:x.frac_coords[-1]<0.5, self.structure),
                key=lambda x: x.frac_coords[-1])
        else:
            func(self)

        self.all=sorted(self.up+self.down, key=lambda x: x.frac_coords[-1])
        assert len(self.all) == len(self.structure)


    def center(self):
        diff = abs(
            1 - self.up[-1].frac_coords[-1]) - self.down[0].frac_coords[-1]
        if abs(diff) < 1e-10:
            return
        self.move(self.all, diff / 2)

    def move(self, group, dist):
        for atom in group:
            atom._fcoords[-1] += dist
            if atom.frac_coords[-1] > 1:
                atom._fcoords[-1] -= 1
            elif atom.frac_coords[-1]<-1:
                atom._fcoords[-1] += 1
            atom._coords = atom._lattice.get_cartesian_coords(atom._fcoords)
        self.center()

    def goal(self, ingoal,outgoal=12):
        while abs(self.indist[1] - ingoal) >= 1e-5 or abs(self.outdist[1]-outgoal)>= 1e-5:
            outdiff = self.outdist[1]-outgoal
            self.changelattice(-outdiff/5)
            indiff = self.indist[1] - ingoal
            frac_indiff=indiff/self.lattice.c
            self.move(self.up,-frac_indiff/2)

    @property
    def indist(self):
        cdist = self.up[0].coords[-1] - self.down[-1].coords[-1]
        fdist = self.up[0].frac_coords[-1] - self.down[-1].frac_coords[-1]
        return fdist, cdist

    @property
    def outdist(self):
        fdist = 1 - self.up[-1].frac_coords[-1] + self.down[0].frac_coords[-1]
        cdist = fdist*self.lattice.c
        return fdist, cdist

    def tofile(self, filename):
        with open(filename, 'w') as f:
            f.write(self.structure.to('POSCAR'))


    def changelattice(self, dist):
        matrix = self.structure.lattice.matrix.copy()
        matrix[-1][-1] += dist
        newlattice = mg.core.lattice.Lattice(matrix)
        self.structure._lattice = newlattice
        self.lattice = newlattice
        for atom in self.structure:
            atom._lattice = newlattice
            atom._fcoords = atom._lattice.get_fractional_coords(atom.coords)
        self.center()

    def fix(self):
        shutil.move('POSCAR','bak.vasp')
        newpos=''
        coords=[]
        with open('bak.vasp') as f:
            for line in f:
                if line.find('direct')>=0:
                    newpos+='S\n'
                    newpos+='Direct\n'
                    break
                newpos+=line
            for line in f:
                coords.append(list(map(float,line.split()[:3])))
        for item in coords:
            if item[2] < self.all[120].frac_coords[-1]+0.001 or item[2]>self.all[-40].frac_coords[-1]-0.01:
                item.extend(['F','F','F'])
            else:
                item.extend(['T','T','T'])
        for item in coords:
            newpos+='   '.join(list(map(str,item))) + '\n'
        with open('POSCAR','w') as f:
            f.write(newpos)

    def batch(self,goals,basename=''):
        for item in goals:
            if os.path.exists(str(item)):
                shutil.rmtree(str(item))
            os.mkdir(str(item))
            os.chdir(str(item))
            self.goal(item)
            self.tofile('POSCAR')
            shutil.copy('../INCAR','.')
            shutil.copy('../KPOINTS','.')
            shutil.copy('../POTCAR','.')
            shutil.copy('../sub.slurm_dt2','.')
            os.system('sed -i "" "s/NAME/'+str(item)+basename+'/g" sub.slurm_dt2')
            os.chdir('..')



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('poscar')
    parser.add_argument('--start',type=float,default=2.0)
    parser.add_argument('--stop',type=float,default=6.0)
    parser.add_argument('--interval',type=float,default=0.5)
    parser.add_argument('--name',default='')

    args=parser.parse_args()
    goals=[]
    num=args.start
    while num<=args.stop:
        goals.append(num)
        num+=args.interval

    a = moveSlab(args.poscar)
    a.group()
    a.batch(goals,args.name)




