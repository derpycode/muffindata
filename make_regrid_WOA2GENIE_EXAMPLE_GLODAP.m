function [] = make_regrid_WOA2GENIE_EXAMPLE_GLODAP(PK1)
% make_regrid_WOA2GENIE_EXAMPLE_GLODAP
%
%   An example of how to use the make_regrid_WOA2GENIE netCDF regridding
%   function
%   In this case -- the GLODAP data
%   Simply pass in the name of the topography 'k1' file 
%   (in the same directory)to define which GENIE grid to re-grid to
%
%   ***********************************************************************
%%
%
close all;
%
% *** INITIALIZE re-grid from WOA to GENIE ****************************** %
%
% set re-gridding parameters
loc_str_longrid = 'longitude'; % netCDF depth (t) grid variable name
loc_str_depthgrid = 'depth'; % netCDF depth (t) grid name (ALT: 'deptht')
loc_str_mask = PK1; % GENIE k1 filename omitting '.k1'
loc_mask_offset = -260; % GENIE mask grid E longitude origin
loc_null = -999.0; % missing or fill value (get from netCDF file info)
loc_n_i = 36;
loc_n_j = 36;
loc_n_k = 16;
% set mask
mask = load([loc_str_mask '.k1'],'-ascii');
loc_mask_size = size(mask);
if (loc_mask_size(1) > loc_n_i),
    mask(1,:) = [];
    mask(end,:) = [];
end
if (loc_mask_size(2) > loc_n_j),
    mask(:,1) = [];
    mask(:,end) = [];
end
loc_rot = int16(loc_mask_offset/10);
mask = circshift(mask,loc_rot,2);
mask = rot90(mask,-1);
%
% *** DO re-grid from WOA to GENIE -- AnthCO2 *************************** %
%
% set netCFD file
loc_str_var = 'AnthCO2'; % netCDF variable name
loc_str_nc = 'AnthCO2'; % netCDF filename omitting '.nc'
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
% call regridding function
disp('... PROCESSING AnthCO2');
make_regrid_WOA2GENIE(data(:,:,:),grid_depth,'AnthCO2','Anthropogenic_CO2','mol kg-1',1.0E-6,loc_n_i,loc_n_j,loc_n_k,mask,[loc_str_mask '.AnthCO2'],'GENIE');
%
% *** DO re-grid from WOA to GENIE -- C14 ******************************* %
%
% set netCFD file
loc_str_var = 'C14'; % netCDF variable name
loc_str_nc = 'C14'; % netCDF filename omitting '.nc'
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
% call regridding function
disp('... PROCESSING C14');
make_regrid_WOA2GENIE(data(:,:,:),grid_depth,'C14','C14','permil',1.0,loc_n_i,loc_n_j,loc_n_k,mask,[loc_str_mask '.C14'],'GENIE');
%
% *** DO re-grid from WOA to GENIE -- CFC11 ***************************** %
%
% set netCFD file
loc_str_var = 'CFC11'; % netCDF variable name
loc_str_nc = 'CFC11'; % netCDF filename omitting '.nc'
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
% call regridding function
disp('... PROCESSING CFC11');
make_regrid_WOA2GENIE(data(:,:,:),grid_depth,'CFC11','CFC-11','mol kg-1',1.0E-12,loc_n_i,loc_n_j,loc_n_k,mask,[loc_str_mask '.CFC11'],'GENIE');
%
% *** DO re-grid from WOA to GENIE -- DIC ******************************* %
%
% set netCFD file
loc_str_var = 'TCO2'; % netCDF variable name
loc_str_nc = 'TCO2'; % netCDF filename omitting '.nc'
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
% call regridding function
disp('... PROCESSING TCO2');
make_regrid_WOA2GENIE(data(:,:,:),grid_depth,'TCO2','Total_CO2','mol kg-1',1.0E-6,loc_n_i,loc_n_j,loc_n_k,mask,[loc_str_mask '.TCO2'],'GENIE');
%
% *** DO re-grid from WOA to GENIE -- ALK ******************************* %
%
% set netCFD file
loc_str_var = 'PALK'; % netCDF variable name
loc_str_nc = 'PALK'; % netCDF filename omitting '.nc'
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
% call regridding function
disp('... PROCESSING PALK');
make_regrid_WOA2GENIE(data(:,:,:),grid_depth,'PALK','Potential_Alkalinity','mol kg-1',1.0E-6,loc_n_i,loc_n_j,loc_n_k,mask,[loc_str_mask '.PALK'],'GENIE');
%
% *********************************************************************** %
