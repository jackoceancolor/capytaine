clear;clc

pe = 4;
refine_stl('plate.stl',2,11,'plate_ref.stl','.');
for i=10:-0.5:pe
    refine_stl('plate_ref.stl',1,i,'plate_ref.stl','.');
end
refine_stl('plate_ref.stl',1,i,'plate_ref.stl','.')
refine_stl('plate.stl',2,pe,'plate_test.stl','.')
% readSTL('plate_test.stl');

fe = 1.0;
refine_stl('float.stl',2,4,'float_ref.stl','.');
for j=4:-0.5:fe
    refine_stl('float_ref.stl',1,j,'float_ref.stl','.');
end
refine_stl('float_ref.stl',1,j,'float_ref.stl','.')
refine_stl('float.stl',2,fe,'float_test.stl','.')

fclose all