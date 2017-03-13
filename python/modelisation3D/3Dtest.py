import mpl_toolkits.mplot3d as a3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#import scipy as sp
# Vertex data
verts= [    
      (-1, -1, -1), (-1, -1, 1), (-1, 1, 1), (-1, 1, -1),
      (1, -1, -1), (1, -1, 1), (1, 1, 1), (1, 1, -1)
        ]

# Face data
faces = np.array([ 
   [0, 1, 2, 3], [4, 5, 6, 7], [0, 3, 7, 4], [1, 2, 6, 5],     
   [0, 1, 5, 4], [2, 3, 7, 6]
   ])

ax = a3.Axes3D(plt.figure())
ax.dist=30
ax.azim=-140
#ax.elev=20['enter image description here'][1]
ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.set_zlim([-1,1])

for i in np.arange(len(faces)):
    square=[verts[faces[i,0]], verts[faces[i,1]], verts[faces[i, 2]], verts[faces[i, 3]]]
    face = a3.art3d.Poly3DCollection([square])
    #face.set_color(colors.rgb2hex(sp.rand(3)))
    face.set_edgecolor('k')
    face.set_alpha(0.5)
    ax.add_collection3d(face)

plt.show()