load('d3d_obj_mks_struct_6565.mat')
bp = tok_data_struct.bpdata;
r = bp(2,:);
z = bp(1,:);

[iNear, dist] = dsearchn([r; z]', [XMP2; YMP2]')




  