# script for obtaining coil current data 
# run on IRIS during an OMFIT session


import scipy.io as sio

# generate tagnames for the d3d coils
tagnames = []
for i  in range(1,10):
    tagnames.append('F' + str(i) + 'A')
    tagnames.append('F' + str(i) + 'B')

# initialize empty data matrix
data = np.zeros((len(tagnames), 49152))

# read data from mdstree
for i in range(len(tagnames)):
    ob = OMFITmdsValue('DIII-D', treename=None, shot=165288, TDI=tagnames[i])
    data[i,:] = ob.data()
    t = ob.dim_of(0)

# save mat file
fn = '/home/waij/irtv/coil_currents.mat'
sio.savemat(fn, {'coil_currents':data, 't':t, 'names':tagnames})

