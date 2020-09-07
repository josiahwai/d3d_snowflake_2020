
#-------------------------------------
# PLOT: NSTX-U divertor configurations

import h5py
import numpy as np 
import math
import matplotlib.pyplot as plt
import scipy.io
from   matplotlib.patches import Rectangle
from   matplotlib.collections import PatchCollection

#---------------------------------------
# Define global parameters for the plots

plt.rc('text', usetex = True)

#-------------------------------------------------------
# Create figure object and define the size of the figure

fig = plt.figure(figsize = (3.10, 3.10))

ax1 = fig.add_axes([0.010, 0.500, 0.490, 0.495])
ax2 = fig.add_axes([0.500, 0.500, 0.490, 0.495])
ax3 = fig.add_axes([0.500, 0.015, 0.490, 0.485])
ax4 = fig.add_axes([0.010, 0.015, 0.490, 0.485])

ax1.set_xticks([])
ax1.set_yticks([])

ax2.set_xticks([])
ax2.set_yticks([])
 
ax3.set_xticks([])
ax3.set_yticks([])

ax4.set_xticks([])
ax4.set_yticks([])

ax1.axis([0.1, 1.1, -2.1, -1.1])
ax2.axis([0.1, 1.1, -2.1, -1.1])
ax3.axis([0.1, 1.1, -2.1, -1.1])
ax4.axis([0.1, 1.1, -2.1, -1.1])

ax1.spines['top']   .set_linewidth(0.5)
ax1.spines['bottom'].set_linewidth(0.5)
ax1.spines['left']  .set_linewidth(0.5)
ax1.spines['right'] .set_linewidth(0.5)

ax2.spines['top']   .set_linewidth(0.5)
ax2.spines['bottom'].set_linewidth(0.5)
ax2.spines['left']  .set_linewidth(0.5)
ax2.spines['right'] .set_linewidth(0.5)

ax3.spines['top']   .set_linewidth(0.5)
ax3.spines['bottom'].set_linewidth(0.5)
ax3.spines['left']  .set_linewidth(0.5)
ax3.spines['right'] .set_linewidth(0.5)

ax4.spines['top']   .set_linewidth(0.5)
ax4.spines['bottom'].set_linewidth(0.5)
ax4.spines['left']  .set_linewidth(0.5)
ax4.spines['right'] .set_linewidth(0.5)

#----------------------------------------
# Define limiter data
#----------------------------------------

rlim = [ 0.315000,  0.315000,  0.415000,  0.415000,  0.435000,  0.571200,  0.617000,  1.043300,  1.043300,  1.319200,  1.335800, 
         1.485100,  1.489000,  1.563800,  1.570000,  1.573700,  1.575000,  1.573700,  1.570000,  1.563800,  1.489000,  1.485100, 
         1.335800,  1.319200,  1.043300,  1.043300,  0.617000,  0.571200,  0.435000,  0.415000,  0.415000,  0.315000,  0.315000]

zlim = [ 0.000000,  1.050000,  1.270000,  1.578000,  1.623400,  1.623400,  1.628000,  1.460287,  1.430000,  1.039700,  0.997600, 
         0.545000,  0.490000,  0.114100,  0.076360,  0.038261,  0.000000, -0.038261, -0.076360, -0.114100, -0.490000, -0.545000, 
        -0.997600, -1.039700, -1.430000, -1.460287, -1.628000, -1.623400, -1.623400, -1.578000, -1.270000, -1.050000,  0.000000]

ax1.plot(rlim, zlim, 'gray', linewidth = 2.0)
ax2.plot(rlim, zlim, 'gray', linewidth = 2.0)
ax3.plot(rlim, zlim, 'gray', linewidth = 2.0)
ax4.plot(rlim, zlim, 'gray', linewidth = 2.0)

#----------------------------------------
# Load the equilibria and grid parameters
#----------------------------------------

rg129 = scipy.io.loadmat('rg129.mat')['rg129']
zg129 = scipy.io.loadmat('zg129.mat')['zg129']

[RGG129, ZGG129] = np.meshgrid(rg129,zg129)

psizr_standard  = scipy.io.loadmat('psizr_standard.mat') ['psizr_standard']
psizr_exactSnow = scipy.io.loadmat('psizr_exactSnow.mat')['psizr_exactSnow']
psizr_snowMinus = scipy.io.loadmat('psizr_snowMinus.mat')['psizr_snowMinus']
psizr_snowPlus  = scipy.io.loadmat('psizr_snowPlus.mat') ['psizr_snowPlus']

psibry_standard  = scipy.io.loadmat('psibry_standard.mat') ['psibry_standard'] [0]
psibry_exactSnow = scipy.io.loadmat('psibry_exactSnow.mat')['psibry_exactSnow'][0]
psibry_snowMinus = scipy.io.loadmat('psibry_snowMinus.mat')['psibry_snowMinus'][0]
psibry_snowPlus  = scipy.io.loadmat('psibry_snowPlus.mat') ['psibry_snowPlus'] [0]

#--------------------------------------------
# Plot the primary separatrices
#--------------------------------------------

sepP_standard  = ax1.contour(RGG129, ZGG129, psizr_standard,  levels = [ psibry_standard  ], colors = 'k') 
sepP_exactSnow = ax2.contour(RGG129, ZGG129, psizr_exactSnow, levels = [ psibry_exactSnow ], colors = 'k')
sepP_snowPlus  = ax4.contour(RGG129, ZGG129, psizr_snowPlus,  levels = [ psibry_snowPlus  ], colors = 'k')
sepP_snowMinus = ax3.contour(RGG129, ZGG129, psizr_snowMinus, levels = [ psibry_snowMinus ], colors = 'k')

plt.setp(sepP_standard. collections, linestyle = 'solid', linewidths = 1.0)
plt.setp(sepP_exactSnow.collections, linestyle = 'solid', linewidths = 1.0)
plt.setp(sepP_snowPlus. collections, linestyle = 'solid', linewidths = 1.0)
plt.setp(sepP_snowMinus.collections, linestyle = 'solid', linewidths = 1.0)

#--------------------------------------------
# Plot open flux surfaces
#--------------------------------------------

fluxOpenStandard = psibry_standard  - np.arange(0.001, 0.100, 0.010)
fluxOpenExact    = psibry_exactSnow - np.arange(0.001, 0.100, 0.010)
fluxOpenPlus     = psibry_snowPlus  - np.arange(0.001, 0.100, 0.010)
fluxOpenMinus    = psibry_snowMinus - np.arange(0.001, 0.100, 0.010)

fluxOpenStandard = fluxOpenStandard[::-1]
fluxOpenExact    = fluxOpenExact[::-1]
fluxOpenPlus     = fluxOpenPlus[::-1]
fluxOpenMinus    = fluxOpenMinus[::-1]

fOpenStandard = ax1.contour(RGG129, ZGG129, psizr_standard,  levels = fluxOpenStandard)
fOpenExact    = ax2.contour(RGG129, ZGG129, psizr_exactSnow, levels = fluxOpenExact)
fOpenPlus     = ax4.contour(RGG129, ZGG129, psizr_snowPlus,  levels = fluxOpenPlus)
fOpenMinus    = ax3.contour(RGG129, ZGG129, psizr_snowMinus, levels = fluxOpenMinus)

plt.setp(fOpenStandard.collections, color = '0.7', linestyle = 'solid', linewidths = 0.5)
plt.setp(fOpenExact.   collections, color = '0.7', linestyle = 'solid', linewidths = 0.5)
plt.setp(fOpenPlus.    collections, color = '0.7', linestyle = 'solid', linewidths = 0.5)
plt.setp(fOpenMinus.   collections, color = '0.7', linestyle = 'solid', linewidths = 0.5)

#----------------------------------------
# Plot the X-points
#----------------------------------------

rxP_standard =  0.550951
zxP_standard = -1.320885

rxP_exactSnow =  0.556201
rxS_exactSnow =  0.560801
zxP_exactSnow = -1.485815
zxS_exactSnow = -1.545668

rxP_snowMinus =  0.534984
rxS_snowMinus =  0.595160
zxP_snowMinus = -1.452147
zxS_snowMinus = -1.589346

rxP_snowPlus =  0.581429
rxS_snowPlus =  0.515127
zxP_snowPlus = -1.457210
zxS_snowPlus = -1.586178

ax1.plot(rxP_standard,  zxP_standard,  c = 'b', marker = 'x', markersize = 10, markeredgewidth = 2.5)

ax2.plot((rxP_exactSnow + rxS_exactSnow)/2, (zxP_exactSnow + zxS_exactSnow)/2, c = 'b', marker = 'x', markersize = 10, markeredgewidth = 2.5)

ax3.plot(rxP_snowMinus, zxP_snowMinus, c = 'b', marker = 'x', markersize = 10, markeredgewidth = 2.5)
ax3.plot(rxS_snowMinus, zxS_snowMinus, c = 'b', marker = 'x', markersize = 10, markeredgewidth = 2.5)

ax4.plot(rxP_snowPlus,  zxP_snowPlus,  c = 'b', marker = 'x', markersize = 10, markeredgewidth = 2.5)
ax4.plot(rxS_snowPlus,  zxS_snowPlus,  c = 'b', marker = 'x', markersize = 10, markeredgewidth = 2.5)

#----------------------------------------
# Plot the coils
#----------------------------------------

PF1aL = Rectangle((0.2933190, -1.8223240), 0.0624840, 0.4634480)
PF1bL = Rectangle((0.3833366, -1.8949035), 0.0338328, 0.1814070)
PF1cL = Rectangle((0.5316224, -1.8968105), 0.0375412, 0.1664210)
PF2La = Rectangle((0.7178500, -1.9675000), 0.1627000, 0.0680000)
PF2Lb = Rectangle((0.7178500, -1.8866000), 0.1627000, 0.0680000)

patches = []
patches.append(PF1aL)
patches.append(PF1bL)
patches.append(PF1cL)
patches.append(PF2La)
patches.append(PF2Lb)

p1 = PatchCollection(patches, edgecolor = 'k', facecolor = [1.0, 0.75, 0.0], zorder = 10)
p2 = PatchCollection(patches, edgecolor = 'k', facecolor = [1.0, 0.75, 0.0], zorder = 10)
p3 = PatchCollection(patches, edgecolor = 'k', facecolor = [1.0, 0.75, 0.0], zorder = 10)
p4 = PatchCollection(patches, edgecolor = 'k', facecolor = [1.0, 0.75, 0.0], zorder = 10)

ax1.add_collection(p1)
ax2.add_collection(p2)
ax3.add_collection(p3)
ax4.add_collection(p4)

#----------------------------------------
# Label the coil currents
#----------------------------------------

ax1.text(0.150,  -1.320, '5.0 kA', fontsize = 7, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})
ax2.text(0.150,  -1.320, '2.8 kA', fontsize = 7, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})
ax3.text(0.150,  -1.320, '3.0 kA', fontsize = 7, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})
ax4.text(0.150,  -1.320, '2.8 kA', fontsize = 7, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})

ax1.text(0.190,  -1.930, '0 kA', fontsize = 7, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})
ax2.text(0.190,  -1.930, '0 kA', fontsize = 7, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})
ax3.text(0.190,  -1.930, '0 kA', fontsize = 7, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})
ax4.text(0.190,  -1.930, '0 kA', fontsize = 7, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})

ax1.text(0.430,  -2.000,      '0 kA', fontsize = 7, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})
ax2.text(0.430,  -2.000, r'--1.7 kA', fontsize = 7, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})
ax3.text(0.430,  -2.000, r'--1.6 kA', fontsize = 7, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})
ax4.text(0.430,  -2.000, r'--1.0 kA', fontsize = 7, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})

ax1.text(0.720,  -2.060, '1.3 kA', fontsize = 7, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})
ax2.text(0.720,  -2.060, '5.6 kA', fontsize = 7, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})
ax3.text(0.720,  -2.060, '5.3 kA', fontsize = 7, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})
ax4.text(0.720,  -2.060, '5.4 kA', fontsize = 7, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})

#----------------------------------------
# Annotate the subplots
#----------------------------------------

ax1.text(0.90, 0.90, r'\textbf{a}', fontsize = 10, transform = ax1.transAxes, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})
ax2.text(0.90, 0.90, r'\textbf{b}', fontsize = 10, transform = ax2.transAxes, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})
ax3.text(0.90, 0.90, r'\textbf{d}', fontsize = 10, transform = ax3.transAxes, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})
ax4.text(0.90, 0.90, r'\textbf{c}', fontsize = 10, transform = ax4.transAxes, fontweight = 'bold', bbox = {'facecolor': [1, 1, 1], 'edgecolor': 'none', 'pad': 2})

plt.savefig('fig1.eps', format = 'eps', dpi = 1000)
