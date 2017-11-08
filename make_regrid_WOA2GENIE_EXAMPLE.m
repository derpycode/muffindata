function [] = make_regrid_WOA2GENIE_EXAMPLE()
% make_regrid_WOA2GENIE_EXAMPLE
%
%   An example of how to use the make_regrid_WOA2GENIE netCDF regridding
%   function
%   In this case (uncommented) -- regridding the WOA2013 1 degree annual  
%   average phosphate climatology
%   Could easily be made into a flexible function by passing in the
%   'dum_*' variables
%
%   ***********************************************************************

% *** INITIALIZE ******************************************************** %
%
loc_str_nc = 'woa13_all_p00_01'; % netCDF filename omitting '.nc'
loc_str_var = 'p_an'; % netCDF variable name
loc_str_longrid = 'lon'; % netCDF depth (t) grid variable name
loc_str_depthgrid = 'depth'; % netCDF depth (t) grid name (ALT: 'deptht')
loc_str_mask = 'worjh2'; % GENIE k1 filename omitting '.k1'
loc_mask_offset = -260; % GENIE mask grid E longitude origin
loc_null = -999.0; % missing or fill value (get from netCDF file info)
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
% remove missing / fill values
data(find(data == loc_null))  = NaN;
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
% call regridding function [examples and alternatives given, commented-out]
% %(1) Fe
% make_regrid_WOA2GENIE(data(:,:,:),grid_depth(:),'fFe_sed','Sediment Fe flux','mol/m2/s',1.0,loc_n_i,loc_n_j,loc_n_k,mask,'Fe_flux_sed','GENIE');
% %(2a) PO4 -- NOTE: no units change (umol l-1)
% make_regrid_WOA2GENIE(data(:,:,:),'PO4','Phosphate concentration','mol kg-1',1.0),loc_n_i,loc_n_j,loc_n_k,mask,'WOA13_PO4_umoll','GENIE');
% %(2b) PO4 -- NOTE: scaling factor of 1/1.024 to convert (umol l-1) -> (umol kg-1)
% make_regrid_WOA2GENIE(data(:,:,:),'PO4','Phosphate concentration','mol kg-1',(1/1.024),loc_n_i,loc_n_j,loc_n_k,mask,'WOA13_PO4_umoll','GENIE');
% (2c) PO4 -- NOTE: scaling factor of 1/1.024E6 to convert (umol l-1) -> (mol kg-1)
make_regrid_WOA2GENIE(data(:,:,:),grid_depth,'PO4','Phosphate concentration','mol kg-1',(1/1.024E6),loc_n_i,loc_n_j,loc_n_k,mask,'WOA13_PO4_molkg','GENIE');
%
% *********************************************************************** %
