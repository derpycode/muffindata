function [] = make_regrid_WOA2GENIE_EXAMPLE_Fe()
% make_regrid_WOA2GENIE_EXAMPLE
%
%   An example of how to use the make_regrid_WOA2GENIE netCDF regridding
%   function
%   In this case -- Fe
%   Could easily be made into a flexible function by passing in the
%   'dum_*' variables
%
%   ***********************************************************************

% *** INITIALIZE ******************************************************** %
%
loc_str_nc = 'dfe_data_unpub_19032015_annWOAgrid'; % netCDF filename omitting '.nc'
loc_str_var = 'IRON'; % netCDF variable name
loc_str_longrid = 'LONGITUDE'; % netCDF depth (t) grid variable name
loc_str_depthgrid = 'DEPTH'; % netCDF depth (t) grid name (ALT: 'deptht')
loc_str_mask = 'mask_worjh2_ALL'; % GENIE k1 filename omitting '.k1'
loc_mask_offset = -260; % GENIE mask grid E longitude origin
loc_n_i = 36;
loc_n_j = 36;
loc_n_k = 16;
%
% *** LOAD DATA ********************************************************* %
%
% OPEN netCDF
ncid = netcdf.open([loc_str_nc '.nc'],'nowrite');
% load vertical grid
varid  = netcdf.inqVarID(ncid,loc_str_depthgrid);
grid_depth = netcdf.getVar(ncid,varid);
% load longitude grid
varid  = netcdf.inqVarID(ncid,loc_str_longrid);
grid_lon = netcdf.getVar(ncid,varid);
% load data
varid  = netcdf.inqVarID(ncid,loc_str_var);
[varname,xtype,dimids,natts] = netcdf.inqVar(ncid,varid);
data = double(netcdf.getVar(ncid,varid));
% check data orientation and transpose to: 360 x 180 x n (where n<180)
loc_size = size(data);
data = permute(data,[find(loc_size == 360) find(loc_size == 180) find(loc_size < 180)]);
% if longitude grid does not start @ 0E: the grid will have to be shifted
loc_data_offset = int16(grid_lon(1) - 0.5);
data = circshift(data,loc_data_offset,1);
% CLOSE netCDF
netcdf.close(ncid);
%
% *** LOAD MASK ********************************************************* %
%
% load mask
mask = load([loc_str_mask '.k1'],'-ascii');
% check mask dimensions and truncate if necessary
loc_mask_size = size(mask);
if (loc_mask_size(1) > loc_n_i),
    mask(1,:) = [];
    mask(end,:) = [];
end
if (loc_mask_size(2) > loc_n_j),
    mask(:,1) = [];
    mask(:,end) = [];
end
% align start of grid to 0E
loc_rot = int16(loc_mask_offset/10);
mask = circshift(mask,loc_rot,2);
% reorient
mask = rot90(mask,-1);
%
% *** REGRID ************************************************************ %
%
% call regridding function
% (Fe)
make_regrid_WOA2GENIE(data(:,:,:),grid_depth,'dFe','Dissolved Fe concentration','nmol kg-1',1.0,loc_n_i,loc_n_j,loc_n_k,mask,'dfe_data_unpub_19032015_annWOAgrid','GENIE');
%
% *********************************************************************** %
