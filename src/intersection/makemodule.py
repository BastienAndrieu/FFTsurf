import os

src = '/home/bandrieu/FFTsurf/src/'

f = open( src + 'modules/mod_intersection.f90', 'w' )
f.write('module mod_intersection\n')
f.write('implicit none\n')
f.write('contains\n')

for filename in os.listdir( src + 'intersection/' ):
    if filename.endswith(".f90"): 
        print '\tmod_intersection.f90', '<---', filename
        f.write('\n\n\n')
        f.writelines( [l for l in open(src + 'intersection/' + filename).readlines()] )

f.write('end module mod_intersection\n')

f.close()
