function [] = make_regrid_WOA2GENIE(PZI,PD,PNAME,PLNAME,PUNITS,PSCALE,PIMAX,PJMAX,PKMAX,PMASK,PFNAME,PFORMAT)
% make_regrid_WOA2GENIE
%
%   ***********************************************************************
%   *** Transform a WOA format array to GENIE netCDF **********************
%   ***********************************************************************
%
%   make_regrid_WOA2GENIE(PZI,PD,PNAME,PLNAME,PUNITS,PIMAX,PJMAX,PKMAX,PMASK,PFNAME,PDEBU)
%   ... takes a WOA gridded dataset and re-grids to the GENIE grid
%   It is important to recognise that it does not try and do anything
%   clever and it completely tied to a 1 degree x 1 degree lon/lat grid
%   The original data grid (WOA) is also assumed to start @ 0E longitude
%
%  'make_regrid_WOA2GENIE.nc' and takes 11 arguments:
%
%   PZI [REAL] (e.g. whatevername)
%   --> a 3D WOA format (360 x 180 x n) array of the data
%   PD [REAL] (e.g. whatevername)
%   --> a 1D vector of the vertical centre (t) grid of the data
%   PNAME [STRING] (e.g., 'dust flux')
%   --> short data name
%   PLNAME [STRING] (e.g., 'regridded dust depositional flux')
%   --> long data name (description)
%   PUNITS [STRING] (e.g., 'kg m-2 yr-1')
%   --> data units
%   PSCALE [REAL] (e.g. 1.0E6)
%   --> scaling factor for data
%   PIMAX [INTEGER] (e.g. 36)
%   --> grid imax
%   PJMAX [INTEGER] (e.g. 36)
%   --> grid jmax
%   PKMAX [INTEGER] (e.g. 16)
%   --> grid kmax
%   PMASK [PIMAXxPJMAX INTEGER ARRAY]
%   --> optional grid mask (the grid 'k1' topo file)
%   PFNAME [STRING] (e.g., 'observed_iron_data')
%   --> filename (if blank, one will be created)
%   PFORMAT [STRING] (e.g., 'GENIE', 'WOA')
%   --> netCDF 'format' (compatability)
%   --> use 'GENIE' to plot with plot_fields_* series of plotting functions
%
%   (lon_lo,lat_lo(, depth)) == original WOA grid
%   (lon_hi,lat_hi) == new WOA grid
%   (i,j(k)) == GENIE grid
%   and similarly: n_lon_lo, n_lat_lo (,n_depth)
%                  n_lon_hi, n_lat_hi
%                  n_i, n_j (,n_k)
%
%   Example
%           make_regrid_WOA2GENIE(fe,zaxis,'fe_conc','dissolved iron concentration','nM',1.0E6,36,36,16,[],'GENIE_fe','NONE');
%           will transform the WOA format array 'fe' to GENIE netCDF format
%           and write the variable 'fe_conc' to the file 'GENIE_fe.nc'
%
%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   15/01/26: CREATED -- adapted from make_regrid_ASCII2netcdf_GENIE.m
%             and incorporating A. Yool's get_gdep.m k-grid generator
%   15/02/16: ... moar ...
%   15/02/17: ... moar ...
%   15/02/18: first (mostly) working version!
%   15/02/19: added in area-weighting to grid cells
%   15/02/22: added format option
%   15/02/26: completed 1st working version
%   15/03/02: fixed error in WOA grid boundaries definition
%             added user-defined null data value
%             commented out hi res grid plotting (too expensive!)
%             fixed issue with sparse vertical data sets
%   15/03/05: added as input, the vertical grid of the original data
%             corrected bug in creating benthic field if no mask provided
%   15/04/26: added options for omitting data averaging;
%             horizontally and/or vertically
%             added ASCII saving of 2D benthic and 3D data
%   15/05/26: cosmetic changes only
%   15/07/22: minor adjustment to ASCII output format
%   15/12/03: altered descrtiption of 'mask' to '.k1' to avoid confusion ..
%   16/02/12: fixed filter for (GENIE) null values
%             added check for mask dimensions
%   16/02/29: moved the mask dim filter to the EXAMPLE (prior to cirshift)
%   18/02/09: added some diagnostics
%
%   ***********************************************************************

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% *** user parameters *************************************************** %
%
% set debuggin options
str_debug = 'NONE';
% set regridding parameter
% NOTE: this is the factot of resolution incrase compared to 1 degrees
%       e.g. 10 == 0.1 degree resolution re-gridding
par_regrid_scalar = 10;
% set GENIE grid paramters -- output grid offset
% NOTE: the default modern GENIE grid has -260 degrees
par_grid_i_offset_out = -260;
% set GENIE grid paramters -- input grid offset compared to data
% NOTE: no adjustment necessary if data grid starts @ 0E
par_grid_i_offset_in = 0;
% set GENIE grid paramters -- misc
par_ez0 = 0.1;
par_grid_k_max = 5000.0;
% null data value
par_data_null = 9.9E19;
% modify data averaging
opt_average_hor  = true;
opt_average_vert = true;
% modify benthic data extraction
opt_ben_sum = false;
%
% *** dummy variable processing ***************************************** %
%
% set passed parameters
data_lo = PZI;
data_grid_depth = PD;
str_dataname = PNAME;
str_dataname_long = PLNAME;
str_data_units = PUNITS;
data_scale = PSCALE;
n_i = PIMAX;
n_j = PJMAX;
n_k = PKMAX;
data_mask = PMASK;
str_filename = PFNAME;
str_format = PFORMAT;
% set strings
if isempty(str_dataname), str_dataname = 'data'; end
if isempty(str_dataname_long), str_dataname_long = str_dataname; end
if isempty(str_data_units), str_data_units = 'n/a'; end
if isempty(str_filename), str_filename = str_dataname; end
% WOA grid parameters
n_lon_woa = 360;
n_lat_woa = 180;
% MISC
switch str_debug
    case 'NONE'
        opt_plot_base  = false;
        opt_plot_raw   = false;
        opt_plot_depth = false;
        opt_plot_k     = false;
    case 'LOW'
        opt_plot_base  = true;
        opt_plot_raw   = false;
        opt_plot_depth = false;
        opt_plot_k     = true;
    case'MED'
        opt_plot_base  = true;
        opt_plot_raw   = false;
        opt_plot_depth = true;
        opt_plot_k     = true;
    case 'HIGH'
        opt_plot_base  = true;
        opt_plot_raw   = true;
        opt_plot_depth = true;
        opt_plot_k     = true;
    otherwise
        opt_plot_base  = true;
        opt_plot_raw   = false;
        opt_plot_depth = false;
        opt_plot_k     = false;
end
%
% *** create GENIE grid ************************************************* %
%
% NOTE: par_grid_i_offset_in is to enable the GENIE grid to be matched
%       to the input grid (i.e. align the Prime Meridian)
% NOTE: at this point, no account is taken of whether the final GENIE
%       grid should start at e.g. -260E (set by par_grid_i_offset_out)
%
% lon (west boundary)
for i=1:n_i,
    axis_iedge(i) = (i-1)*(360.0/n_i) + par_grid_i_offset_in;
    axis_di(i)    = (360.0/n_i);
end
axis_iedge(n_i+1)   = (n_i)*(360.0/n_i) + par_grid_i_offset_in;
axis_imid           = axis_iedge(1:n_i) + 0.5*axis_di;
axis_ibnds(1:n_i,1) = axis_iedge(1:n_i);
axis_ibnds(1:n_i,2) = axis_iedge(2:n_i+1);
axis_ibnds          = axis_ibnds';
% lat (south boundary)
for j=1:n_j,
    axis_jedge(j) = (180.0/pi)*asin((2*(j-1)/n_j) - 1.0);
    axis_jmid(j)  = (180.0/pi)*asin(((1 + 2*(j-1))/n_j) - 1.0);
end
axis_jedge(n_j+1)   = (180.0/pi)*asin((2*n_j/n_j) - 1.0);
axis_jbnds(1:n_j,1) = axis_jedge(1:n_j);
axis_jbnds(1:n_j,2) = axis_jedge(2:n_j+1);
axis_jbnds          = axis_jbnds';
% depth (bottom boundary)
z1 = par_ez0*((1.0 + 1/par_ez0)^(1.0/n_k) - 1.0);
tv4 = par_ez0*((z1/par_ez0+1)^0.5-1);
tv2 = 0;
tv1 = 0;
zro(n_k) = -tv4;
zw(n_k+1) = tv2;
for k=1:1:n_k
    if par_ez0 > 0
        tv3 = par_ez0*((z1/par_ez0+1)^k-1);
        dz(n_k-k+1) = tv3 - tv2;
        tv2 = tv3;
        tv5 = par_ez0*((z1/par_ez0+1)^(k+0.5)-1);
        if k < n_k
            dza(n_k-k) = tv5 - tv4;
        end
        tv4 = tv5;
        tv1 = tv1 + dz(n_k-k+1);
    else
        dz(k) = 1d0/n_k;
        dza(k) = 1d0/n_k;
    end
end
for k=n_k:-1:1
    if k > 1
        zro(k-1) = zro(k) - dza(k-1);
    end
    zw(k) = zw(k+1) - dz(k);
end
% set depth grid bounds
% NOTE: k counts from TOP to BOTTOM; 
%       bnd #1 is top
axis_kmid(1:n_k)    = -par_grid_k_max*zro(:);
axis_kbnds(1:n_k,1) = -par_grid_k_max*zw(2:n_k+1);
axis_kbnds(1:n_k,2) = -par_grid_k_max*zw(1:n_k);
axis_kbnds          = axis_kbnds';
%
% *** create WOA depth grid ********************************************* %
%
% number of points
n_depth = length(data_grid_depth);
% edges -- single vector format
data_grid_depthedge = zeros(n_depth+1,1);
data_grid_depthedge(1) = data_grid_depth(1) - (data_grid_depth(2)-data_grid_depth(1))/2.0;
for n=2:n_depth,
    data_grid_depthedge(n) = data_grid_depth(n-1) + (data_grid_depth(n)-data_grid_depth(n-1))/2.0;
end
data_grid_depthedge(n_depth) = data_grid_depth(n_depth) + (data_grid_depth(n_depth)-data_grid_depth(n_depth-1))/2.0;
data_grid_depthedge = data_grid_depthedge';
% edges -- bnds format
axis_depthbnds = zeros(n_depth,2);
axis_depthbnds(1:n_depth,1) = data_grid_depthedge(1:n_depth);
axis_depthbnds(1:n_depth,2) = data_grid_depthedge(2:n_depth+1);
axis_depthbnds = axis_depthbnds';
%
% *** misc initialization *********************************************** %
%
% set null/fill value
loc_nullvalue = -par_data_null;
% set date
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
% plot something ...
if ~isempty(data_mask),
    if (opt_plot_base),
        figure;
        plot_2dgridded(data_mask(:,:)',0.9E19,'','mask','grid mask');
    end
end
%
% *** CONSTANTS ********************************************************* %
%
const_rEarth = 6371000.0;
%
% *********************************************************************** %

% *********************************************************************** %
% *** RE-GRID DATA ****************************************************** %
% *********************************************************************** %
%
% *** get size and orientation of lo data in array ********************** %
%
% retrieve dimensions and sort in ascending order; check dimensions
data_in_dims = size(data_lo);
if (length(data_in_dims) ~= 3),
    disp(['ERROR: data array is not 3D.']);
    return;
end
%data_in_dims = sort(data_in_dims);
n_depth = data_in_dims(3);
% reorder array if necessary
if (data_in_dims(2) > data_in_dims(1)),
    n_lon_lo = data_in_dims(2);
    n_lat_lo = data_in_dims(1);
    tmp_data_lo = zeros(n_lon_lo,n_lat_lo,n_depth);
    for depth=1:n_depth
        tmp_data_lo(:,:,depth) = transpose(data_lo(:,:,depth));
    end
    data_lo = tmp_data_lo;
else
    n_lon_lo = data_in_dims(1);
    n_lat_lo = data_in_dims(2);
end
%
% *** re-grid: horizontal *********************************************** %
%
% create temp 2D data hi array
n_lon_hi = par_regrid_scalar*n_lon_lo;
n_lat_hi = par_regrid_scalar*n_lat_lo;
loc_data_hi = zeros(n_lon_hi,n_lat_hi);
loc_area_hi = zeros(n_lon_hi,n_lat_hi);
% create output arrays
data_depth_out = zeros(n_i,n_j,n_depth);
data_depth_count = zeros(n_i,n_j,n_depth);
% process null values
data_lo(find(data_lo >= par_data_null))  = NaN;
data_lo(find(data_lo <= -par_data_null)) = NaN;
% calculate hi grid cell areas
for lat=1:n_lat_hi,
    loc_area_hi(:,lat) = 2.0*pi*(const_rEarth^2)*( sin(pi*((-90+lat/par_regrid_scalar)/180)) - sin(pi*((-90+lat/par_regrid_scalar-1)/180)) )/(par_regrid_scalar*n_lon_lo);
end
% % plot something
% if (opt_plot_base),
%     figure;
%     plot_2dgridded(loc_area_hi(:,:)',0.9E19,'','areahi','hi resolution area grid');
% end
%
% >>> DEPTH LOOP START
for depth=1:n_depth,
    % convert lo -> hi resolution data and populate hi array
    for lon=1:n_lon_lo,
        for lat=1:n_lat_lo,
            loc_data_hi((par_regrid_scalar*(lon-1)+1):par_regrid_scalar*lon,(par_regrid_scalar*(lat-1)+1):par_regrid_scalar*lat) = data_lo(lon,lat,depth);
        end
    end
    % loop through GENIE grid
    % populate hi array
    % NOTE: 'count' should ideally be weighted by orignal data grid cell fractional area
    %       (not done here => small error induced)
    % NOTE: commented out are non-area-weighted alternatives
    for i=1:n_i,
        for j=1:n_j,
            loc_lon_min = par_regrid_scalar*int16(axis_ibnds(1,i) + 0.5);
            loc_lon_max = par_regrid_scalar*int16(axis_ibnds(2,i) - 0.5);
            loc_lat_min = par_regrid_scalar*int16(axis_jbnds(1,j) + 90.0 + 0.5);
            loc_lat_max = par_regrid_scalar*int16(axis_jbnds(2,j) + 90.0 - 0.5);
            loc_data = loc_data_hi(loc_lon_min:loc_lon_max,loc_lat_min:loc_lat_max);
            loc_area = loc_area_hi(loc_lon_min:loc_lon_max,loc_lat_min:loc_lat_max);
            loc_data_nan = find(isnan(loc_data));
            if (length(loc_data_nan) > 0),
                loc_data(loc_data_nan) = 0.0;
                loc_area(loc_data_nan) = 0.0;
            end
            loc_count = sum(sum(loc_area));
            if (loc_count > 0),
                % NOTE: if no horizontal averaging selected:
                %       sum data only
                if opt_average_hor,
                    data_depth_out(i,j,depth) = sum(sum(loc_data.*loc_area))/double(loc_count);
                else
                    data_depth_out(i,j,depth) = sum(sum(loc_data));
                end
                % NOTE: if no vertical averaging selected:
                %       maintain count at unity
                %       (to indicate valid data, no more)
                if opt_average_vert,
                    data_depth_count(i,j,depth) = data_depth_count(i,j,depth) + loc_count;
                else
                    data_depth_count(i,j,depth) = 1;
                end
            end
        end
    end
    % plot something ...
    if (opt_plot_depth),
        loc_data_depth_out   = data_depth_out(:,:,depth);
        loc_data_depth_count = data_depth_count(:,:,depth);
        loc_data_depth_out(find(loc_data_depth_count == 0)) = NaN;
        close;
        if (opt_plot_raw),
            figure;
            plot_2dgridded(data_lo(:,:,depth)',0.9E19,'',['dataraw_D' num2str(depth)],['raw data for depth = ' num2str(depth)]);
        end
        figure;
        plot_2dgridded(loc_data_depth_out',0.9E19,'',['data_D' num2str(depth)],['data for depth = ' num2str(depth)]);
        figure;
        plot_2dgridded(loc_data_depth_count',0.9E19,'',['count_D' num2str(depth)],['count for depth = ' num2str(depth)]);
    end
end
%
% <<< DEPTH LOOP END
%
% *** re-grid: vertical ************************************************* %
%
% clear local arrays
loc_data = [];
loc_count = [];
% create output arrays
data_out = zeros(n_i,n_j,n_k);
data_count = zeros(n_i,n_j,n_k);
% set default mask
if isempty(data_mask),
    data_mask = zeros(n_i,n_j);
    disp(['WARNING: 2D benthic field cannot be populated without a depth level mask being supplied.']);
end
%
% >>> K LOOP START
for k=n_k:-1:1,
    if (~strcmp(str_debug,'NONE')), disp(['k = ' num2str(k)]); end
    % find:
    %       all WOA levels with lower depth deeper (or equal to) than GENIE upper depth
    %       AND
    %       all WOA levels with upper depth shallower than GENIE lower depth
    % NOTE: axis_depthbnds(1,:) == upper; axis_depthbnds(2,:) == lower bnd
    loc_ndepth = intersect(find(axis_depthbnds(2,:)>=axis_kbnds(1,k)),find(axis_depthbnds(1,:)<axis_kbnds(2,k)));
    % determine whether depth level is wholly, or partially (and which 'end') within GENIE depth level
    if (~isempty(loc_ndepth)),
        loc_n = length(loc_ndepth);
        loc_data = zeros(n_i,n_j);
        loc_count = zeros(n_i,n_j);
        for n = loc_ndepth(1):loc_ndepth(loc_n)
            if ( (axis_depthbnds(2,n) <= axis_kbnds(2,k)) && (axis_depthbnds(1,n) >= axis_kbnds(1,k)) ),
                % CASE #1: WOA level wholly within GENIE level,
                % i.e. bottom depth shallower than GENIE bottom and top deeper than GENIE top
                loc_dD = axis_depthbnds(2,n) - axis_depthbnds(1,n);
                loc_dD_frac = 1.0;
            elseif ( (axis_depthbnds(1,n) <= axis_kbnds(1,k)) ),
                % CASE #2: WOA level top shallower than GENIE top
                loc_dD = axis_depthbnds(2,n) - axis_kbnds(1,k);
                loc_dD_frac = loc_dD/(axis_depthbnds(2,n)-axis_depthbnds(1,n));
            elseif ( (axis_depthbnds(2,n) >= axis_kbnds(2,k)) ),
                % CASE #2: WOA level bottom deeper than GENIE bottom
                loc_dD = axis_kbnds(2,k) - axis_depthbnds(1,n);
                loc_dD_frac = loc_dD/(axis_depthbnds(2,n)-axis_depthbnds(1,n));
            else
                disp(['ERROR: impossible! [something up with grid meshing']);
                return;
            end
            if ((loc_dD < 0) || (loc_dD_frac < 0.0)),
                disp(['ERROR: impossible! [negative depth interval and/or fraction]']);
                return;
            end
            % NOTE: if no vertical averaging selected:
            %       (1) add the fractional intersection of original data 
            %           depth layer with GENIE layer
            %       (2) ensure 'loc_count' is always only '1'
            %           (becomes domininator in normalization later below)
            if opt_average_vert,
                loc_data(:,:)  = loc_data(:,:)  + loc_dD*data_depth_count(:,:,n).*data_depth_out(:,:,n);
                loc_count(:,:) = loc_count(:,:) + loc_dD*data_depth_count(:,:,n);
            else
                loc_data(:,:)  = loc_data(:,:)  + loc_dD_frac*data_depth_out(:,:,n);
                loc_count(:,:) = max(loc_count(:,:),data_depth_count(:,:,n));
            end
        end
        % filter out non-value points
        % also: scale data
        for i=1:n_i,
            for j=1:n_j,
                if (loc_count(i,j) > 0),
                    data_out(i,j,k) = data_scale*loc_data(i,j)/loc_count(i,j);
                else
                    data_out(i,j,k) = NaN;
                end
            end
        end
    end
    if (opt_plot_k),
        close;
        figure;
        plot_2dgridded(data_out(:,:,k)',0.9E19,'',['data_k' num2str(k)],['data for k = ' num2str(k)]);
    end
end
% <<< K LOOP END
%
% *** apply mask ******************************************************** %
%
for k=n_k:-1:1,
    for i=1:n_i,
        for j=1:n_j,
            if (data_mask(i,j) > k), data_out(i,j,k) = NaN; end
        end
    end
end
%
% *** extract benthic data ********************************************** %
%
% create output array
data_out_ben = zeros(n_i,n_j);
for i=1:n_i,
    for j=1:n_j,
        loc_k = data_mask(i,j);
        if ((loc_k >= 1) && (loc_k <= n_k)),
            % create water-column sum if requested
            if opt_ben_sum,
                data_out(i,j,find(isnan(data_out(i,j,:)))) = 0.0;
                data_out_ben(i,j) = sum(data_out(i,j,loc_k:n_k));
            else
                data_out_ben(i,j) = data_out(i,j,loc_k);
            end
        else
            data_out_ben(i,j) = NaN;
        end
    end
end
%
% *** adjust longitude origin ******************************************* %
%
if (par_grid_i_offset_out ~= 0),
    % calculate rotation of arrays
    loc_rot = -int16(par_grid_i_offset_out/10);
    % align start of axes at 0E
    axis_imid  = circshift(axis_imid,loc_rot,2);
    axis_ibnds = circshift(axis_ibnds,loc_rot,2);
    % create negative W
    axis_imid(1:loc_rot)    = axis_imid(1:loc_rot)-360.0;
    axis_ibnds(1,1:loc_rot) = axis_ibnds(1,1:loc_rot)-360.0;
    axis_ibnds(2,1:loc_rot) = axis_ibnds(2,1:loc_rot)-360.0;
    % align start of grid at 0E
    data_out     = circshift(data_out,loc_rot,1);
    data_mask    = circshift(data_mask,loc_rot,1);
    data_out_ben = circshift(data_out_ben,loc_rot,1);
end
%
% *** GENIE format changes ********************************************** %
%
if strcmp(str_format,'GENIE'),
    % invert depth dim
    data_out = flip(data_out,3);
    axis_kmid = flip(axis_kmid,2);
    axis_kbnds = flip(axis_kbnds,2);
    % create grid edge vectors
    % NOTE: axis_lon_edges already circ shifted
    axis_zt_edges  = [axis_kbnds(1,:) axis_kbnds(2,end)];
    axis_lat_edges = [axis_jbnds(1,:) axis_jbnds(2,end)];
    axis_lon_edges = [axis_ibnds(1,:) axis_ibnds(2,end)];
end
%
% *** DIAGNOSTICS ******************************************************* %
%
% calculate volume-weighted average
% NOTE: axis_depthbnds(1,:) == upper; axis_depthbnds(2,:) == lower bnd
loc_n    = 0;
loc_wsum = 0.0; % weighted sum
loc_sumw = 0.0; % sum of weights ...
loc_sum  = 0.0;
for k=n_k:-1:1,
    for i=1:n_i,
        for j=1:n_j,
            if (~isnan(data_out(i,j,k))),
                loc_n    = loc_n + 1;
                loc_wsum = loc_wsum + data_out(i,j,k)*(axis_kbnds(2,k) - axis_kbnds(1,k));
                loc_sumw = loc_sumw + (axis_kbnds(2,k) - axis_kbnds(1,k));
                loc_sum  = loc_sum  + data_out(i,j,k);
            end
        end
    end
end
%
loc_ben_n    = 0;
loc_ben_wsum = 0.0; % weighted sum
loc_ben_sumw = 0.0; % sum of weights ...
loc_ben_sum  = 0.0;
for i=1:n_i,
    for j=1:n_j,
        if (~isnan(data_out_ben(i,j))),
            loc_ben_n    = loc_ben_n + 1;
            loc_ben_wsum = loc_ben_wsum + data_out_ben(i,j);
            loc_ben_sumw = loc_ben_sumw + 1.0;
            loc_ben_sum  = loc_ben_sum  + data_out_ben(i,j);
        end
    end
end
loc_ben_av = loc_ben_wsum/loc_ben_sumw;
% write to file
fid = fopen([str_filename '.DIAG.' str_date '.txt'],'w');
fprintf(fid,'%s\n','##################################################################################');
fprintf(fid,'\n');
fprintf(fid,'%s\n',['# of ocean data points    = ',num2str(loc_n)]);
fprintf(fid,'%s\n',['weighted mean of data     = ',num2str(loc_wsum/loc_sumw)]);
fprintf(fid,'%s\n',['mean of data              = ',num2str(loc_sum/loc_n)]);
fprintf(fid,'\n');
fprintf(fid,'%s\n',['# of benthic data points  = ',num2str(loc_ben_n)]);
fprintf(fid,'%s\n',['weighted mean of benthic  = ',num2str(loc_ben_wsum/loc_ben_sumw)]);
fprintf(fid,'%s\n',['mean of benthic data      = ',num2str(loc_ben_sum/loc_ben_n)]);
fprintf(fid,'\n');
fprintf(fid,'%s\n','##################################################################################');
fclose(fid);
%
% *********************************************************************** %

% *********************************************************************** %
% *** WRITE netCDF FILE ************************************************* %
% *********************************************************************** %
%
% *** create netCDF file ************************************************ %
%
ncid = netcdf.create([str_filename, '.', str_date, '.nc'],'NC_WRITE');
NC_GLOBAL = netcdf.getConstant('NC_GLOBAL');
%
% format compatability
if strcmp(str_format,'GENIE'),
    % *** GENIE-compatable format ******************************************* %
    % define global attributes
    netcdf.putAtt(ncid,NC_GLOBAL,'title','Gridded data');
    netcdf.putAtt(ncid,NC_GLOBAL,'long_title','Regridded on a regular lon-lat, irregular in depth, grid');
    netcdf.putAtt(ncid,NC_GLOBAL,'comments','conversion of GENIE ascii format 2D data by make_regrid_ASCII2netcdf_GENIE.m');
    netcdf.putAtt(ncid,NC_GLOBAL,'history','version as of 15/02/18');
    netcdf.putAtt(ncid,NC_GLOBAL,'Conventions','CF-1.6 / GENIE modified');
    netcdf.putAtt(ncid,NC_GLOBAL,'CreationDate',datestr(now,'yyyy/mm/dd HH:MM:SS'));
    netcdf.putAtt(ncid,NC_GLOBAL,'CreatedBy',[getenv('username'), '@', getenv('computername')]);
    netcdf.putAtt(ncid,NC_GLOBAL,'MatlabSource','make_regrid_WOA2GENIE.m');
    % define dimensions
    dimid_time = netcdf.defDim(ncid,'time',1);
    dimid_zt = netcdf.defDim(ncid,'zt',n_k);
    dimid_lat = netcdf.defDim(ncid,'lat',n_i);
    dimid_lon = netcdf.defDim(ncid,'lon',n_j);
    dimid_bnds = netcdf.defDim(ncid,'nbounds',2);
    dimid_zt_edges = netcdf.defDim(ncid,'zt_edges',n_k+1);
    dimid_lat_edges = netcdf.defDim(ncid,'lat_edges',n_i+1);
    dimid_lon_edges = netcdf.defDim(ncid,'lon_edges',n_j+1);
    % define axes -- time
    varid = netcdf.defVar(ncid,'time','double',dimid_time);
    netcdf.putAtt(ncid,varid,'standard_name','time');
    netcdf.putAtt(ncid,varid,'long_name','Year');
    netcdf.putAtt(ncid,varid,'units','years');
    netcdf.putAtt(ncid,varid,'Axis','T');
    netcdf.putAtt(ncid,varid,'_CoordinateAxisType','Time');
    varid_time = varid;
    % define axes -- D
    varid = netcdf.defVar(ncid,'zt','double',dimid_zt);
    netcdf.putAtt(ncid,varid,'standard_name','depth');
    netcdf.putAtt(ncid,varid,'long_name','Vertical distance below the surface');
    netcdf.putAtt(ncid,varid,'units','m');
    netcdf.putAtt(ncid,varid,'point_spacing','uneven');
    netcdf.putAtt(ncid,varid,'Axis','Z');
    netcdf.putAtt(ncid,varid,'positive','down');
    netcdf.putAtt(ncid,varid,'_CoordinateAxisType','Depth');
    netcdf.putAtt(ncid,varid,'bounds','depth_bnds');
    varid_zt = varid;
    % define axes -- Z bnds
    varid = netcdf.defVar(ncid,'zt_bnds','double',[dimid_bnds, dimid_zt]);
    varid_ztbnds = varid;
    % define axes -- Z edges
    varid = netcdf.defVar(ncid,'zt_edges','double',dimid_zt_edges);
    varid_zt_edges = varid;
    % define axes -- Y
    varid = netcdf.defVar(ncid,'lat','double',dimid_lat);
    netcdf.putAtt(ncid,varid,'standard_name','latitude');
    netcdf.putAtt(ncid,varid,'long_name','latitude of grid centre');
    netcdf.putAtt(ncid,varid,'units','degrees_north');
    netcdf.putAtt(ncid,varid,'point_spacing','uneven');
    netcdf.putAtt(ncid,varid,'Axis','Y');
    netcdf.putAtt(ncid,varid,'_CoordinateAxisType','Lat');
    netcdf.putAtt(ncid,varid,'bounds','lat_bnds');
    varid_lat = varid;
    % define axes -- Y bnds
    varid = netcdf.defVar(ncid,'lat_bnds','double',[dimid_bnds, dimid_lat]);
    varid_latbnds = varid;
    % define axes -- lat edges
    varid = netcdf.defVar(ncid,'lat_edges','double',dimid_lat_edges);
    varid_lat_edges = varid;
    % define axes -- X
    varid = netcdf.defVar(ncid,'lon','double',dimid_lon);
    netcdf.putAtt(ncid,varid,'standard_name','longitude');
    netcdf.putAtt(ncid,varid,'long_name','longitude of grid centre');
    netcdf.putAtt(ncid,varid,'units','degrees_east');
    netcdf.putAtt(ncid,varid,'point_spacing','even');
    netcdf.putAtt(ncid,varid,'Axis','X');
    netcdf.putAtt(ncid,varid,'_CoordinateAxisType','Lon');
    netcdf.putAtt(ncid,varid,'bounds','lon_bnds');
    varid_lon = varid;
    % define axes -- X bnds
    varid = netcdf.defVar(ncid,'lon_bnds','double',[dimid_bnds, dimid_lon]);
    varid_lonbnds = varid;
    % define axes -- lon edges
    varid = netcdf.defVar(ncid,'lon_edges','double',dimid_lon_edges);
    varid_lon_edges = varid;
    % define data variable -- data
    varid = netcdf.defVar(ncid,str_dataname,'double',[dimid_lon, dimid_lat, dimid_zt, dimid_time]);
    netcdf.putAtt(ncid,varid,'name',str_dataname);
    netcdf.putAtt(ncid,varid,'long_name',str_dataname_long);
    netcdf.putAtt(ncid,varid,'missing_value',loc_nullvalue);
    netcdf.putAtt(ncid,varid,'fillValue',loc_nullvalue);
    netcdf.putAtt(ncid,varid,'units',str_data_units);
    varid_data = varid;
    % define data variable -- benthic data
    varid = netcdf.defVar(ncid,['ben_' str_dataname],'double',[dimid_lon, dimid_lat, dimid_time]);
    netcdf.putAtt(ncid,varid,'name',['ben_' str_dataname]);
    netcdf.putAtt(ncid,varid,'long_name',['benthic ' str_dataname_long]);
    netcdf.putAtt(ncid,varid,'missing_value',loc_nullvalue);
    netcdf.putAtt(ncid,varid,'fillValue',loc_nullvalue);
    netcdf.putAtt(ncid,varid,'units',str_data_units);
    varid_data_ben = varid;
    % define data variable -- grid mask
    varid = netcdf.defVar(ncid,'grid_level','double',[dimid_lon, dimid_lat]);
    netcdf.putAtt(ncid,varid,'name','grid_level');
    netcdf.putAtt(ncid,varid,'long_name','grid level');
    netcdf.putAtt(ncid,varid,'missing_value',loc_nullvalue);
    netcdf.putAtt(ncid,varid,'fillValue',loc_nullvalue);
    netcdf.putAtt(ncid,varid,'units',str_data_units);
    varid_data_mask = varid;
    % end definition
    netcdf.endDef(ncid);
    % write axes
    netcdf.putVar(ncid,varid_zt,axis_kmid);
    netcdf.putVar(ncid,varid_lat,axis_jmid);
    netcdf.putVar(ncid,varid_lon,axis_imid);
    netcdf.putVar(ncid,varid_ztbnds,axis_kbnds);
    netcdf.putVar(ncid,varid_latbnds,axis_jbnds);
    netcdf.putVar(ncid,varid_lonbnds,axis_ibnds);
    netcdf.putVar(ncid,varid_time,[1.0]);
    netcdf.putVar(ncid,varid_zt_edges,axis_zt_edges);
    netcdf.putVar(ncid,varid_lat_edges,axis_lat_edges);
    netcdf.putVar(ncid,varid_lon_edges,axis_lon_edges);
    % write data
    netcdf.putVar(ncid,varid_data,data_out);
    % write benthic data
    netcdf.putVar(ncid,varid_data_ben,data_out_ben);
    % write level grid
    netcdf.putVar(ncid,varid_data_mask,data_mask);
else
    % *** STANDARD format *************************************************** %
    % define global attributes
    netcdf.putAtt(ncid,NC_GLOBAL,'title','Gridded data');
    netcdf.putAtt(ncid,NC_GLOBAL,'long_title','Regridded on a regular lon-lat, irregular in depth, grid');
    netcdf.putAtt(ncid,NC_GLOBAL,'comments','conversion of GENIE ascii format 2D data by make_regrid_ASCII2netcdf_GENIE.m');
    netcdf.putAtt(ncid,NC_GLOBAL,'history','version as of 15/02/18');
    netcdf.putAtt(ncid,NC_GLOBAL,'Conventions','CF-1.6');
    netcdf.putAtt(ncid,NC_GLOBAL,'CreationDate',datestr(now,'yyyy/mm/dd HH:MM:SS'));
    netcdf.putAtt(ncid,NC_GLOBAL,'CreatedBy',[getenv('username'), '@', getenv('computername')]);
    netcdf.putAtt(ncid,NC_GLOBAL,'MatlabSource','make_regrid_WOA2GENIE.m');
    % define dimensions
    dimid_lat = netcdf.defDim(ncid,'lat',n_i);
    dimid_lon = netcdf.defDim(ncid,'lon',n_j);
    dimid_D = netcdf.defDim(ncid,'depth',n_k);
    dimid_bnds = netcdf.defDim(ncid,'nbounds',2);
    % define axes -- Y
    varid = netcdf.defVar(ncid,'lat','double',dimid_lat);
    netcdf.putAtt(ncid,varid,'standard_name','latitude');
    netcdf.putAtt(ncid,varid,'long_name','latitude of grid centre');
    netcdf.putAtt(ncid,varid,'units','degrees_north');
    netcdf.putAtt(ncid,varid,'point_spacing','uneven');
    netcdf.putAtt(ncid,varid,'Axis','Y');
    netcdf.putAtt(ncid,varid,'_CoordinateAxisType','Lat');
    netcdf.putAtt(ncid,varid,'bounds','lat_bnds');
    varid_lat = varid;
    % define axes -- Y bnds
    varid = netcdf.defVar(ncid,'lat_bnds','double',[dimid_bnds, dimid_lat]);
    varid_latbnds = varid;
    % define axes -- X
    varid = netcdf.defVar(ncid,'lon','double',dimid_lon);
    netcdf.putAtt(ncid,varid,'standard_name','longitude');
    netcdf.putAtt(ncid,varid,'long_name','longitude of grid centre');
    netcdf.putAtt(ncid,varid,'units','degrees_east');
    netcdf.putAtt(ncid,varid,'point_spacing','even');
    netcdf.putAtt(ncid,varid,'Axis','X');
    netcdf.putAtt(ncid,varid,'_CoordinateAxisType','Lon');
    netcdf.putAtt(ncid,varid,'bounds','lon_bnds');
    varid_lon = varid;
    % define axes -- X bnds
    varid = netcdf.defVar(ncid,'lon_bnds','double',[dimid_bnds, dimid_lon]);
    varid_lonbnds = varid;
    % define axes -- D
    varid = netcdf.defVar(ncid,'depth','double',dimid_D);
    netcdf.putAtt(ncid,varid,'standard_name','depth');
    netcdf.putAtt(ncid,varid,'long_name','Vertical distance below the surface');
    netcdf.putAtt(ncid,varid,'units','m');
    netcdf.putAtt(ncid,varid,'point_spacing','uneven');
    netcdf.putAtt(ncid,varid,'Axis','Z');
    netcdf.putAtt(ncid,varid,'positive','down');
    netcdf.putAtt(ncid,varid,'_CoordinateAxisType','Depth');
    netcdf.putAtt(ncid,varid,'bounds','depth_bnds');
    varid_D = varid;
    % define axes -- D bnds
    varid = netcdf.defVar(ncid,'depth_bnds','double',[dimid_bnds, dimid_D]);
    varid_Dbnds = varid;
    % define data variable -- data
    varid = netcdf.defVar(ncid,str_dataname,'double',[dimid_lon, dimid_lat, dimid_D]);
    netcdf.putAtt(ncid,varid,'name',str_dataname);
    netcdf.putAtt(ncid,varid,'long_name',str_dataname_long);
    netcdf.putAtt(ncid,varid,'missing_value',loc_nullvalue);
    netcdf.putAtt(ncid,varid,'fillValue',loc_nullvalue);
    netcdf.putAtt(ncid,varid,'units',str_data_units);
    varid_data = varid;
    % end definition
    netcdf.endDef(ncid);
    % write axes
    netcdf.putVar(ncid,varid_D,axis_kmid);
    netcdf.putVar(ncid,varid_lat,axis_jmid);
    netcdf.putVar(ncid,varid_lon,axis_imid);
    netcdf.putVar(ncid,varid_Dbnds,axis_kbnds);
    netcdf.putVar(ncid,varid_latbnds,axis_jbnds);
    netcdf.putVar(ncid,varid_lonbnds,axis_ibnds);
    % write data
    netcdf.putVar(ncid,varid_data,data_out);
end
%
% *** close netCDF file ************************************************* %
%
netcdf.close(ncid);
%
% *********************************************************************** %

% *********************************************************************** %
% *** WRITE ASCII FILE ************************************************** %
% *********************************************************************** %
%
% NOTE: k counts from TOP (1) to BOTTOM (n_k)
% NOTE: also remember that fprint_2D optionally flips the entire array u-d
% save 3D data
tmp_data = [];
tmp_data = zeros(n_k*n_i,n_j);
for k=n_k:-1:1,
    tmp_data((n_k-k)*n_i+1:(n_k-k)*n_i+n_i,:) = flipud(rot90(data_out(:,:,k),1));
end
tmp_string = [str_filename '.3D.' str_date '.dat'];
fprint_2D(tmp_data,tmp_string,'%11.3e','%11.1f',false,false);
% save 2D benthic data
tmp_data = [];
tmp_data = zeros(n_i,n_j);
tmp_data = rot90(data_out_ben(1:n_i,1:n_j),1);
tmp_string = [str_filename '.2Dbenthic.' str_date '.dat'];
fprint_2D(tmp_data,tmp_string,'%11.3e','%11.1f',true,true);
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
% END
disp(['END ...'])
%
% *********************************************************************** %
