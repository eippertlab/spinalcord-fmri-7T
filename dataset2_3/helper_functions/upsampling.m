% upsample the given image in z direction by certain factor (scaling)
%
% needs: NLM Upsample by Pierrick Coupe
% https://sites.google.com/site/pierrickcoupe/softwares/super-resolution/monomodal
%
% Ulrike Horn
% uhorn@cbs.mpg.de
% Feb 2024

function upsampling(NLMUpsample_path, nifti_tools_path, image, scaling)

addpath(NLMUpsample_path);
addpath(nifti_tools_path);

gunzip(image);
tmp = split(image, '.nii.gz');
img = tmp{1};
nii = load_untouch_nii([img '.nii']);
fima = NLMUpsample2(double(nii.img), [1, 1, scaling]);
nii.hdr.dime.dim(4) = nii.hdr.dime.dim(4) * scaling;
nii.img = fima;
save_untouch_nii(nii, [img '_upsampled.nii'])
gzip([img '_upsampled.nii'])
delete([img '_upsampled.nii'])
delete([img '.nii'])
end