"""
Visualize wavefield
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import animation
import matplotlib.patches as patches
from sys import argv
from pylab import *

nx = 208
nz = 208 #308
nt = 1000
nnt =1000

Bcwavefield = np.fromfile("Boundaryfield.dat", dtype=np.float64).reshape(nt, 100)
Bcwavefield_back = np.fromfile("Boundaryfield_back_mute.dat", dtype=np.float64).reshape(nt, 100)
force_test = np.fromfile("force_test.dat", dtype=np.float64).reshape(1000, 100)

wavefield = np.fromfile("Fullwavefield.dat", dtype=np.float64).reshape(nt, nx, nz)
wavefield_back = np.fromfile("Fullwavefield_back_mute.dat", dtype=np.float64).reshape(nt, nx, nz)
#wavefield_old = np.fromfile("Fullwavefield_back_old.dat", dtype=np.float64).reshape(nt, nx, nz)
#wavefield_cur = np.fromfile("Fullwavefield_back_cur.dat", dtype=np.float64).reshape(nt, nx, nz)


#te = np.fromfile("force_0.799.dat", dtype=np.float64).reshape(nx, nz)
#test = np.fromfile("force_test.dat", dtype=np.float64).reshape(nt, 100)

t = np.linspace(0.2, 1.0, num=nnt)
tt = np.linspace(0., 1.0, num=nt)

rect = patches.Rectangle((50,50),100,100,linewidth=1,edgecolor='r',facecolor='none')


fig, ax = plt.subplots(1,1)
ax.add_patch(rect)
img = ax.imshow(np.zeros((nx, nz), dtype=np.float32), vmin=-10, vmax=10)

def animate(i):
    img.set_data(wavefield_back[i,:,:])

anim = animation.FuncAnimation(fig, animate, frames=nt, interval=2)

fig.colorbar(img)
plt.ylim([208, 50])
ymin, ymax = plt.ylim()
#plt.show()

def forceAspect(ax,aspect):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)


plt.figure(2)
#plt.subplot(1,2,1)
a = imshow(Bcwavefield, extent=[0,100,nt,0],aspect='auto')
#a = imshow(Bcwavefield, extent=[0,200,1000,0],aspect='auto')
colorbar(a)
plt.xlabel("Distance (m)")
plt.ylabel("Time (s)")

rect1 = patches.Rectangle((50,50),100,100,linewidth=1,edgecolor='r',facecolor='none')


plt.figure(4)
b = imshow(Bcwavefield_back, extent=[0,100,nt,0],aspect='auto')
colorbar(b)
plt.xlabel("Distance (m)")
plt.ylabel("Time (s)")
plt.title("Backward Boundary Wavefield Muted the Top Wave")


plt.figure(5)
c = imshow(force_test, extent=[0,100,2000,0],aspect='auto')
colorbar(c)
plt.xlabel("Distance (m)")
plt.ylabel("Time (s)")
plt.title("Force Test")


migrated_image = np.zeros((nx, nz), dtype=np.float32)
for i in range(nx):
    for j in range(nz):
        migrated_image[i, j] = np.sum(wavefield_back[:, i, j]*wavefield[::-1, i, j])

rect1 = patches.Rectangle((51,51),100,100,linewidth=1,edgecolor='r',facecolor='none')
fig1, ax1 = plt.subplots()
ax1.add_patch(rect1)
ax1.axhline(114, color='red', lw=1)
img1 = ax1.imshow(migrated_image, cmap='coolwarm_r', vmin=-150, vmax=150)
plt.title("Cross_correlation Image")
#plt.imshow(migrated_image,cmap = 'coolwarm_r', vmin = -150, vmax = 150, interpolation='bilinear')


#plt.figure(3)
#b = imshow(test, aspect='auto')
#colorbar(b)
'''
plt.figure(4)
c = imshow(te, aspect='auto')
colorbar(c)
'''
plt.show()

#plt.subplot(1,2,2)

