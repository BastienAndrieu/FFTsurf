import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys

fig, ax = plt.subplots()

tree = np.loadtxt('tree_surf37.dat')

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
#ax.set_xticks( np.linspace(-1.0,1.0,5) )
#ax.set_yticks( np.linspace(-1.0,1.0,5) )
plt.show()
