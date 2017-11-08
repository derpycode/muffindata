function [FCMAP] = make_cmap(PCNAME,PNCOL)
% make_cmap
%
%   *******************************************************************   %
%   *** make colorbar color scale *************************************   %
%   *******************************************************************   %
%
%   make_cmap( ... )
%   creates an evenly-spaced color scale of a given number of colors
%
%   PCNAME [STRING]
%   --> name of colorbar colorscale
%   PNCOL [integer]
%   --> the number of colors
%
%   The color scale options are:
%
%   (1) One of the original MATLAB colorscales:
%       'jet','hsv','hot','cool','spring','summer','autumn','winter',
%       'gray','bone','copper','pink'
%   (2) The new MATLAB (2014 release) default (the 'jet' replacement):
%       'parula'
%       (note that this color scheme is only approximated here)
%   (3) An alternative rainbow replacements:
%       'CubicYF'
%       'LinearL'
%       http://mycarta.wordpress.com/2013/02/21/perceptual-rainbow-palette-
%       the-method/
%       http://www.mathworks.com/matlabcentral/fileexchange/28982-perceptua
%       lly-improved-colormaps/content/pmkmp/pmkmp.m
%       (note that this color scheme is only approximated here)
%   (4) One of the divergent cbrewer schemes
%       'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn'
%       http://colorbrewer2.org/
%       http://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---
%       colorbrewer-schemes-for-matlab
%       (note that these color schemes are only approximated here)
%   (5) One of the built-in/internal color schemes:
%       'anom' -- is the existing anomoly color scehem (blue-white-red)
%       '_parula' -- is an approximation of the parula m-function
%
%   The default is 'parula', and if parula.m is not available,
%   the approximateion '_parula' is used instead.
%
%   Original author:
%   Andy Ridgwell <andy@seao2.org>
%
%   License:
%   CC BY-NC-SA 3.0 (http://creativecommons.org/licenses/by-nc-sa/3.0/)
%
%   *******************************************************************   %

% *********************************************************************** %
% ***** HISTORY ********************************************************* %
% *********************************************************************** %
%
%   14/12/03: CREATED
%   15/01/01: adjustment to default colotmap
%   15/01/09: re-wrote the 're-gridding' of the color map
%   15/01/11: added 'LinearL' color scheme
%   15/01/11: renamed from make_cmap5.m
%             + substituted simple linear interpolation in place of 
%              existing bizarre and pooor code ...
%   15/01/11: added wt% color scale
%   15/02/11: added missing scaling to anom and wt% color scales
%   15/03/01: changed 'cubic' -> 'PCHIP' in interp1 on MATLAB 'advice' ...
%   15/03/03: sorted out the previous 'fix' that was bugged on old versions
%             improved warning message
%   16/02/17: corrected error message
%   17/05/15: added white-to-red-black, white-to-blue-black scale options
%
% *********************************************************************** %

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% process dummy parameters
str_name = PCNAME;
col_n = PNCOL;
%
if (col_n < 2),
    disp(['ERROR: number of colors (here = ', num2str(col_n), ') must be > 1.']);
    return;
end
% initialize color map array
loc_str_cmap = [];
% determine MUTLAB version
tmp_mutlab = version('-release');
str_mutlab = tmp_mutlab(1:4);
par_mutlab = str2num(str_mutlab);
%
% *********************************************************************** %

% *********************************************************************** %
% *** CREATE COLOR SCALE ************************************************ %
% *********************************************************************** %
%
% *** SET DEFINING TIE-POINT COLORS ************************************* %
%
switch str_name
    case {'anom'}
        loc_cmap = [  0   0 255;
                    127 255 255;
                    255 255 255;
                    255 255 127;
                    255   0   0];
        loc_cmap = loc_cmap/255.0;  
    case {'wt%'}
        loc_cmap = [ 64  32  16;
                    128  64  32;
                    255 128  64;
                    255 255 255;
                    255 255 255];
        loc_cmap = loc_cmap/255.0;  
    case {'red'}
        loc_cmap = [255 255 255;
                    255 255 127;
                    255   0   0;
                      0   0   0];
        loc_cmap = loc_cmap/255.0;  
    case {'blue'}
        loc_cmap = [255 255 255;
                    127 255 255;
                      0   0 255;
                      0   0   0];
        loc_cmap = loc_cmap/255.0;  
    case {'_parula'}
        loc_cmap = [0.208100000000000   0.166300000000000   0.529200000000000;
                    0.079475000000000   0.515900000000000   0.832825000000000;
                    0.198550000000000   0.721400000000000   0.630950000000000;
                    0.826575000000000   0.732025000000000   0.346350000000000;
                    0.976300000000000   0.983100000000000   0.053800000000000];
    case {'RdYlGn'}
        loc_cmap = [                0   0.407843137254902   0.215686274509804;
                    0.525490196078432   0.796078431372549   0.401960784313726;
                    1.000000000000000   1.000000000000000   0.749019607843137;
                    0.974509803921569   0.554901960784314   0.321568627450980;
                    0.647058823529412                   0   0.149019607843137];
    case {'RdYlBu'}
        loc_cmap = [0.192156862745098   0.211764705882353   0.584313725490196;
                    0.562745098039216   0.764705882352942   0.866666666666666;
                    1.000000000000000   1.000000000000000   0.749019607843137;
                    0.974509803921569   0.554901960784314   0.321568627450980;
                    0.647058823529412                   0   0.149019607843137];
    case {'RdGy'}
        loc_cmap = [0.101960784313725   0.101960784313725   0.101960784313725;
                    0.629411764705882   0.629411764705882   0.629411764705882;
                                    1                   1                   1;
                    0.898039215686275   0.511764705882353   0.405882352941176;
                    0.403921568627451                   0   0.121568627450980];
    case {'RdBu'}
        loc_cmap = [0.019607843137255   0.188235294117647   0.380392156862745;
                    0.417647058823529   0.674509803921569   0.817647058823529;
                    0.968627450980392   0.968627450980392   0.968627450980392;
                    0.898039215686275   0.511764705882353   0.405882352941176;
                    0.403921568627451                   0   0.121568627450980];
    case {'PuOr'}
        loc_cmap = [0.176470588235294                   0   0.294117647058824;
                    0.600000000000000   0.560784313725490   0.749019607843137;
                    0.968627450980392   0.968627450980392   0.968627450980392;
                    0.935294117647059   0.615686274509804   0.233333333333333;
                    0.498039215686275   0.231372549019608   0.031372549019608];
    case {'PRGn'}
        loc_cmap = [                0   0.266666666666667   0.105882352941176;
                    0.501960784313725   0.770588235294118   0.503921568627451;
                    0.968627450980392   0.968627450980392   0.968627450980392;
                    0.680392156862745   0.543137254901961   0.741176470588235;
                    0.250980392156863                   0   0.294117647058824];
    case {'PiYG'}
        loc_cmap = [0.152941176470588   0.392156862745098   0.098039215686275;
                    0.609803921568628   0.809803921568627   0.390196078431373;
                    0.968627450980392   0.968627450980392   0.968627450980392;
                    0.907843137254902   0.590196078431373   0.768627450980392;
                    0.556862745098039   0.003921568627451   0.321568627450980];
    case {'BrBG'}
        loc_cmap = [                0   0.235294117647059   0.188235294117647;
                    0.354901960784314   0.698039215686274   0.658823529411765;
                    0.960784313725490   0.960784313725490   0.960784313725490;
                    0.811764705882353   0.633333333333333   0.333333333333333;
                    0.329411764705882   0.188235294117647   0.019607843137255];
    case {'LinearL'}
        loc_cmap =  [0.0143	0.0143	0.0143;
                     0.1413	0.0555	0.1256;
                     0.1761	0.0911	0.2782;
                     0.1710	0.1314	0.4540;
                     0.1074	0.2234	0.4984;
                     0.0686	0.3044	0.5068;
                     0.0008	0.3927	0.4267;
                     0.0000	0.4763	0.3464;
                     0.0000	0.5565	0.2469;
                     0.0000	0.6381	0.1638;
                     0.2167	0.6966	0.0000;
                     0.3898	0.7563	0.0000;
                     0.6912	0.7795	0.0000;
                     0.8548	0.8041	0.4555;
                     0.9712	0.8429	0.7287;
                     0.9692	0.9273	0.8961]; 
    case {'CubicYF'}
        loc_cmap =  [0.5151    0.0482    0.6697;
                     0.5199    0.1762    0.8083;
                     0.4884    0.2912    0.9234;
                     0.4297    0.3855    0.9921;
                     0.3893    0.4792    0.9775;
                     0.3337    0.5650    0.9056;
                     0.2795    0.6419    0.8287;
                     0.2210    0.7123    0.7258;
                     0.2468    0.7612    0.6248;
                     0.2833    0.8125    0.5069;
                     0.3198    0.8492    0.3956;
                     0.3602    0.8896    0.2919;
                     0.4568    0.9136    0.3018;
                     0.6033    0.9255    0.3295;
                     0.7066    0.9255    0.3414;
                     0.8000    0.9255    0.3529]; 
    case {'jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink'}
        if (exist(str_name) == 0),
            disp(['ERROR: colormap ', str_name, ' does not exist.']);
            return;
        else
            loc_str_cmap = str_name;
        end
    case {'parula'}
        if (exist('parula') == 2 || exist('parula') == 5),
            loc_str_cmap = 'parula';
        else
            loc_cmap = [0.208100000000000   0.166300000000000   0.529200000000000;
                        0.079475000000000   0.515900000000000   0.832825000000000;
                        0.198550000000000   0.721400000000000   0.630950000000000;
                        0.826575000000000   0.732025000000000   0.346350000000000;
                        0.976300000000000   0.983100000000000   0.053800000000000];
            disp(['WARNING: colormap ', str_name, ' cannot be found (MATLAB version: ', num2str(par_mutlab), ' too old) => using parula-like replacement.']);
        end
    otherwise
            loc_cmap = [0.208100000000000   0.166300000000000   0.529200000000000;
                        0.079475000000000   0.515900000000000   0.832825000000000;
                        0.198550000000000   0.721400000000000   0.630950000000000;
                        0.826575000000000   0.732025000000000   0.346350000000000;
                        0.976300000000000   0.983100000000000   0.053800000000000];
        disp(['WARNING: colormap ', str_name, ' cannot be found => using parula-like default.']);
end
%
% *** CREATE COLOR SCALE ************************************************ %
%
if ~isempty(loc_str_cmap),
    cmap = colormap(eval([loc_str_cmap '(' num2str(col_n) ')']));
else
    loc_n_old = [1:1:length(loc_cmap(:,1))];
    loc_n_new = [1:((length(loc_cmap(:,1))-1)/(col_n-1)):length(loc_cmap(:,1))];
    if (par_mutlab >= 2014),
        cmap = interp1(loc_n_old(:),loc_cmap(:,:),loc_n_new(:),'pchip');
    else
        cmap = interp1(loc_n_old(:),loc_cmap(:,:),loc_n_new(:),'cubic');
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
% export data
FCMAP = cmap;
%
% *********************************************************************** %
