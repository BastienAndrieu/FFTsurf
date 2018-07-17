VERBOSE = True

progname = 'polylines'
f = open( '../' + progname + '.f90', 'w' )
f.write('program ' + progname + '\n')

f.writelines( [l for l in open('main.f90').readlines()] )


#f.write('\n\ncontains\n')

listfiles = ['intersect_2Dpolylines.f90']


for filename in listfiles:
    if VERBOSE:
        print '\t * <---', filename
    f.write('\n\n\n\n\n')
    f.writelines( [l for l in open('../dev_intersection/' + filename).readlines()] )

f.write('end program ' + progname + '\n')

f.close()
