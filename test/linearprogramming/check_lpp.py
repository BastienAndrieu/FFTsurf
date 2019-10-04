import numpy as np
import matplotlib.pyplot as plt

f = open( 'result.dat', 'r' )
fname = f.readline().rstrip('\n').strip()
stat = int( f.readline() )
if stat == 0:
    x = np.array( [ float(a) for a in f.readline().split() ] )
f.close()


f = open( fname, 'r' )


f.readline() #'dimension'
dim = int( f.readline() )
f.readline() #'objectif'
c = np.array( [ float(a) for a in f.readline().split() ] )
f.readline() #'contraintes'
nc = int( f.readline() )
A = np.zeros( (nc,dim+1) )
for i in range(nc):
    A[i,:] = [ float(a) for a in f.readline().split() ]
f.close()


# vertices
#v = np.linalg.lstsq( A[:,0:dim], -A[:,dim] )
#print v
v = np.empty((0,2))
for i in range(nc-1):
    for j in range(i+1,nc):
        det = A[i,0]*A[j,1] - A[i,1]*A[j,0]
        if abs(det) > 1.e-15:
            vij = np.array( [(-A[i,2]*A[j,1] + A[j,2]*A[i,1])/det , 
                       (-A[i,0]*A[j,2] + A[i,2]*A[j,0])/det ] )
            inside = True
            for k in range(1,nc):
                if k != i and k != j:
                    if A[k,0]*vij[0] + A[k,1]*vij[1] + A[k,2] < 0.0:
                        inside = False
                        break
            if inside:
                #v.append( vij )
                v = np.vstack( (v, vij) )

print v

# figure
fig, ax = plt.subplots()
ax.plot( v[:,0], v[:,1], 'k.' )
#ax.quiver( 0.0, 0.0, c[0], c[1] )
if stat == 0:
    ax.plot( x[0], x[1], 'r*' )

ax.quiver(x[0], x[1], c[0], c[1], color='r')

ax.set_aspect('equal', adjustable='box')
plt.show()
