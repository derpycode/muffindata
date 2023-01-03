function [OUT_MASK] = fun_make_masks(IN_K1,IN_KMAX)
%
%   ***********************************************************************
%   *** fun_make_masks ****************************************************
%   ***********************************************************************
%
%   Create masks from a .k1 file 
%
%   \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
%   /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
%
%   ***********************************************************************

% *********************************************************************** %
% *** CREATE MASKS ****************************************************** %
% *********************************************************************** %
%
% *** initialize ******************************************************** %
%
% % paths
% addpath('y:\__git\muffindata')
%
% *** load data ********************************************************* %
%
% load data
% NOTE: test for file first ... try adding extension
% NOTE: determine topo name
if exist(IN_K1,'file')
    loc_k1 = load(IN_K1);
    if (strcmp(IN_K1(end-2:end),'.k1') || strcmp(IN_K1(end-2:end),'.K1'))
        str_k1 = IN_K1(1:end-3);
    else
        str_k1 = IN_K1;        
    end
elseif exist([IN_K1 '.k1'],'file')
    loc_k1 = load([IN_K1 '.k1']);
    str_k1 = IN_K1;
else
    % ERROR
    disp([' ** ERROR: Cannot find .k1 file: ' IN_K1 ' (with or without .k1 extension)']);
    disp([' ']);
    return;
end
% strip off N and S pole data; wrap-around boundary columns
loc_k1([1, end],:) = [];
loc_k1(:,[1, end]) = [];
%
% *** create masks ****************************************************** %
%
% first ... remove and mark land
loc_mask = loc_k1;
loc_mask(find(loc_k1 >= 90)) = NaN;
% basic ('ALL') mask (NOTE: identical to k == kmax or deeper)
loc_mask_all = loc_mask;
loc_mask_all(find(loc_mask_all < 90)) = 1.0;
% create k-level mask 
loc_mask_k = loc_mask;
loc_mask_k(find(loc_mask_k <= IN_KMAX)) = 1.0;
loc_mask_k(find(loc_mask_k > IN_KMAX)) = 0.0;
% create k-level mask's compliment
loc_mask_notk = loc_mask;
loc_mask_notk(find(loc_mask_notk <= IN_KMAX)) = 0.0;
loc_mask_notk(find(loc_mask_notk > IN_KMAX)) = 1.0;
%
% *** save mask ********************************************************* %
%
% plot and save masks
% NOTE: remember to get the file in the correct, (lom,lat) orientation!
% fprint_MASK(loc_mask_all(:,:),[str_k1 '.' 'ALL' '.mask'],true);
str_mask = [str_k1 '.' 'le' num2str(IN_KMAX) '.dat'];
plot_2dgridded(flipud(loc_mask_k),99999.0,'',str_mask,['mask -- ' str_mask],[0 1]);
fprint_MASK(loc_mask_k(:,:),str_mask,true);
str_mask = [str_k1 '.' 'gt' num2str(IN_KMAX) '.dat'];
plot_2dgridded(flipud(loc_mask_notk),99999.0,'',str_mask,['mask -- ' str_mask],[0 1]);
fprint_MASK(loc_mask_notk(:,:),str_mask,true);
%
% *** END *************************************************************** %
%
close all;
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %

end