import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys

scl = 7.0
ratio = 16.0/9.0


args = sys.argv
if len(args) < 2:
    nsurf = 2
else:
    nsurf = int(args[1])


nrow = np.sqrt(float(nsurf)/ratio)
nrow = max(1, int(round(nrow)))
#ncol = (nsurf + 1)/nrow
ncol = int( np.ceil( float(nsurf)/float(nrow) ) )
m = nrow*ncol


#print 'nsurf = ', nsurf, ', nrow = ', nrow, ', ncol = ', ncol

fig, axes = plt.subplots(nrows=nrow, ncols=ncol)

for isurf in range(m):
    j = isurf%ncol
    i = (isurf-j)/ncol
    if nrow > 1:
        ax = axes[i,j]
    else:
        ax = axes[j]
    if isurf < nsurf:
        tree = np.loadtxt( 'treessi_' + str(isurf+1) + '.dat' )
        if len(tree.shape) < 2:
            rect = patches.Rectangle( ( tree[0], tree[2] ), # lower-left corner
                               tree[1] - tree[0],    # u-extent
                               tree[3] - tree[2],    # v-extent
                               fill = False,
                               edgecolor = 'k',
                               linewidth = 1 )
            ax.add_patch( rect )
        else:
            for region in tree:
                rect = patches.Rectangle( ( region[0], region[2] ), # lower-left corner
                                          region[1] - region[0],    # u-extent
                                          region[3] - region[2],    # v-extent
                                          fill = False,
                                          edgecolor = 'k',
                                          linewidth = 1 )
                ax.add_patch( rect )
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim( -1, 1 )
        ax.set_ylim( -1, 1 )
        ax.set_title( 'Surface #' + str(isurf+1),
                               size = 12, weight = 'bold' )
        ax.set_xticks( np.linspace(-1.0,1.0,5) )
        ax.set_yticks( np.linspace(-1.0,1.0,5) )
    else:
        fig.delaxes(ax)
plt.show()

