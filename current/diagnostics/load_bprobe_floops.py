# USE THIS SCRIPT ON IRIS, WITHIN OMFIT, TO LOAD BPROBE AND FLUX LOOP DATA
# AND SAVE IT AS MAT FILE

import scipy.io as sio

k = OMFIT['CAKE']['OUTPUTS']['CAKE-0FWDNnr_155354_t04000.0-04000.0n001d050.0']['kineticEFIT']['kEQDSK'][4000]['IN1']

m = OMFIT['CAKE']['OUTPUTS']['CAKE-0FWDNnr_155354_t04000.0-04000.0n001d050.0']['kineticEFIT']['mEQDSK'][4000]

template = OMFIT['CAKE']['EFITtime']['TEMPLATES']['D3D_dprobe_112000']['IN3']

'''
struct for B-probe and flux-loop data. For defs, see
https://fusion.gat.com/theory/Efitefund and
https://fusion.gat.com/theory/Efitinputs
'''
bp_fl_data = {}  
mtags = ['expmpi', 'cmpr2', 'fwtmp2', 'saimpi','silopt','csilop','fwtsi']
ktags = ['EXPMP2', 'FWTMP2', 'COILS', 'FWTSI']
for tag in mtags:
    bp_fl_data[tag] = m[tag]['data'][0]

for tag in ktags:
    bp_fl_data[tag] = k[tag]

# probe and flux loop geometry
bp_fl_template = {}
ttags = ['XMP2','YMP2','AMP2','RSI','ZSI']
for tag in ttags:
    bp_fl_template[tag] = template[tag]
 
savedir = '/home/waij/d3d_snowflake_2020/current/diagnostics/'
sio.savemat(savedir + 'bp_fl_data.mat', bp_fl_data)
sio.savemat(savedir + 'bp_fl_template.mat', bp_fl_template)










