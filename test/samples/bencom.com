%chk=bencom.chk
%mem=12gb
%nproc=6
#p b3lyp/genecp nosymm geom=connectivity
opt freq iop(7/33=1)

opt-title

0 1
C     -0.196066032299      0.083954610270      0.000060359368
C      1.202825076336      0.083884471006      0.000506018868
C      1.902331654047      1.295325035437     -0.000047136416
C      1.202946707718      2.506835148041     -0.001045596240
C     -0.195944409385      2.506905316939     -0.001491091778
C     -0.895450997150      1.295464738220     -0.000938299116
H     -0.739751345816     -0.857625934164      0.000489987063
H      1.746416421211     -0.857750015636      0.001281451623
H      2.989606822804      1.295270625436      0.000299297106
H      1.746632277886      3.448415533193     -0.001475037473
H     -0.739536124154      3.448539598789     -0.002266379146
H     -1.982726171199      1.295519560979     -0.001284573853

1 2 1.5 6 1.5 7 1.0
2 3 1.5 8 1.0
3 4 1.5 9 1.0
4 5 1.5 10 1.0
5 6 1.5 11 1.0
6 12 1.0
7
8
9
10
11
12

C H 0
6-31+g*
****


