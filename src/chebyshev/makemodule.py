import os

VERBOSE = False


nom_module = 'chebyshev'
src = '/d/bandrieu/GitHub/FFTsurf/src/'
#src = '/home/bastien/GitHub/FFTsurf/src/'

f = open( src + 'modules/mod_' + nom_module + '.f90', 'w' )
f.write('module mod_' + nom_module + '\n')

f.writelines( [l for l in open(src + nom_module + '/header').readlines()] )


f.write('\n\ncontains\n')

for filename in os.listdir( src + nom_module + '/' ):
    if filename.endswith(".f90"): 
        if VERBOSE:
            print '\tmod_' + nom_module + '.f90', '<---', filename
        f.write('\n\n\n')
        f.writelines( [l for l in open(src + nom_module + '/' + filename).readlines()] )

f.write('end module mod_' + nom_module + '\n')

f.close()
