

sbxprocess.m
% Aligns images in fname for all indices in idx to ref
% 
% m - mean image after the alignment
% T - optimal translation for each frame
sbxalignxpar
sbxalignx_ref
sbxalign_nonrigid_ref(fname,idx,N,ref)


sbxadoptmask.m
% Will take the segmentation from fn1 and transfer to fn2.


 sbxautocaltrig
% Auto-calibrate trigger level


vol3d(varargin)
%H = VOL3D Volume render 3-D data. 
% VOL3D uses the orthogonal plane 2-D texture mapping technique for 
% volume rending 3-D data in OpenGL. Use the 'texture' option to fine 
% tune the texture mapping technique. This function is best used with
% fast OpenGL hardware.