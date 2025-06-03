% shift the dilated cord mask to the area where the discs are
%
% needs: nifti_tools
%
% Ulrike Horn
% uhorn@cbs.mpg.de
% Feb 2024

function shifting_mask_to_discs(nifti_tools_path, sub, dil_cord_mask)

sub_id = ['sub-sspr' sub];

addpath(nifti_tools_path);

% take the cord mask and shift it to where the discs are
gunzip(dil_cord_mask)
tmp = split(dil_cord_mask, '.nii.gz');
mask_name = tmp{1};
nii = load_untouch_nii([mask_name '.nii']);
indices = find(nii.img==1);
[ind_x, ind_y, ind_z] = ind2sub(size(nii.img), indices);
new_img = zeros(size(nii.img));
if strcmp(sub_id, 'sub-sspr02') || strcmp(sub_id, 'sub-sspr09') ||...
        strcmp(sub_id, 'sub-sspr21') || strcmp(sub_id, 'sub-sspr17') 
    shift = 23;
elseif strcmp(sub_id, 'sub-sspr10') || strcmp(sub_id, 'sub-sspr15') || ...
        strcmp(sub_id, 'sub-sspr16') || strcmp(sub_id, 'sub-sspr20') || ...
        strcmp(sub_id, 'sub-sspr36')
    shift = 25;
elseif strcmp(sub_id, 'sub-sspr12')
    shift = 26;
elseif strcmp(sub_id, 'sub-sspr06') || strcmp(sub_id, 'sub-sspr22') ||...
    strcmp(sub_id, 'sub-sspr24') || strcmp(sub_id, 'sub-sspr28')
    shift = 19;
elseif strcmp(sub_id, 'sub-sspr23') || strcmp(sub_id, 'sub-sspr25') || ...
        strcmp(sub_id, 'sub-sspr37')
    shift = 17;
elseif strcmp(sub_id, 'sub-sspr18') || strcmp(sub_id, 'sub-sspr41')
    shift = 20;
% checked for 01, 03, 04, 05, 07, 08, 11, 13, 14, 19, 26, 27, 29, 30, 31,
% 32, 33, 34, 35, 38, 39, 40
else
    shift = 22;
end
ind_y = ind_y + shift;
indices = sub2ind(size(nii.img), ind_x, ind_y, ind_z);
new_img(indices) = 1;
nii.img = new_img;
save_untouch_nii(nii, [mask_name '_shifted.nii']);
gzip([mask_name '_shifted.nii'])
% delete unzipped versions
delete([mask_name '_shifted.nii'])
delete([mask_name '.nii'])

end