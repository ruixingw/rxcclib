import pymatgen as mg
a=mg.Structure.from_file('POSCAR')

print(mg.io.vasp.inputs.Kpoints.automatic_density(a,1000))
while 1:
    num=input('Input density:')
    print(mg.io.vasp.inputs.Kpoints.automatic_density(a,int(num)))
