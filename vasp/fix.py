import shutil
shutil.move('POSCAR','bakPOSCAR')
newpos=''
coords = []
with open('bakPOSCAR') as f:
    for line in f:
        if line.find('Direct')>=0:
            newpos+='S\n'
            newpos+='Direct\n'
            break
        newpos+=line
    for line in f:
        coords.append(list(map(float,line.split())))

for item in coords:
    if item[2] < 0.2:
        item.append('F')
        item.append('F')
        item.append('F')
    else:
        item.append('T')
        item.append('T')
        item.append('T')

for item in coords:
    newpos+='  '.join(list(map(str,item))) + '\n'

with open('POSCAR','w') as f:
    f.write(newpos)


