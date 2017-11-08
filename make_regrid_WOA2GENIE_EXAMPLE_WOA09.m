function [] = make_regrid_WOA2GENIE_EXAMPLE_WOA09(PK1)
% make_regrid_WOA2GENIE_EXAMPLE_WOA09
%
%   An example of how to use the make_regrid_WOA2GENIE netCDF regridding
%   function
%   In this case -- the WOA2009 data
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
loc_str_longrid = 'lon'; % netCDF depth (t) grid variable name
loc_str_depthgrid = 'depth'; % netCDF depth (t) grid name (ALT: 'deptht')
loc_str_mask = PK1; % GENIE k1 filename omitting '.k1'
loc_mask_offset = -260; % GENIE mask grid E longitude origin
loc_null = 9.96921E36; % missing or fill value (get from netCDF file info)
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
% *** DO re-grid from WOA to GENIE -- T ********************************* %
%
% set netCFD file
loc_str_var = 't_an'; % netCDF variable name
loc_str_nc = 'temperature_annual_1deg'; % netCDF filename omitting '.nc'
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
disp('... PROCESSING temperature');
make_regrid_WOA2GENIE(data(:,:,:),grid_depth,'t_an','sea_water_temperature','degrees Celsius',1.0,loc_n_i,loc_n_j,loc_n_k,mask,[loc_str_mask '.t_an'],'GENIE');
%
% *** DO re-grid from WOA to GENIE -- S ********************************* %
%
% set netCFD file
loc_str_var = 's_an'; % netCDF variable name
loc_str_nc = 'salinity_annual_1deg'; % netCDF filename omitting '.nc'
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
disp('... PROCESSING salinity');
make_regrid_WOA2GENIE(data(:,:,:),grid_depth,'s_an','sea_water_salinity','1',1.0,loc_n_i,loc_n_j,loc_n_k,mask,[loc_str_mask '.s_an'],'GENIE');
%
% *** DO re-grid from WOA to GENIE -- P ********************************* %
%
% set netCFD file
loc_str_var = 'p_an'; % netCDF variable name
loc_str_nc = 'phosphate_annual_1deg'; % netCDF filename omitting '.nc'
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
disp('... PROCESSING PO4');
make_regrid_WOA2GENIE(data(:,:,:),grid_depth,'p_an','mole_concentration_of_phosphate_in_sea_water','mol kg-1',(1.0/1.024E6),loc_n_i,loc_n_j,loc_n_k,mask,[loc_str_mask '.p_an'],'GENIE');
%
% *** DO re-grid from WOA to GENIE -- O2 ******************************** %
%
% set netCFD file
loc_str_var = 'o_an'; % netCDF variable name
loc_str_nc = 'dissolved_oxygen_annual_1deg'; % netCDF filename omitting '.nc'
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
disp('... PROCESSING O2');
make_regrid_WOA2GENIE(data(:,:,:),grid_depth,'o_an','mole_concentration_of_dissolved_oxygen_in_sea_water','mol kg-1',(1.4276/32.0E3/1.024),loc_n_i,loc_n_j,loc_n_k,mask,[loc_str_mask '.o_an'],'GENIE');
% NOTE: DO solubility is calculated in milliliters per liter (ml/L) which must be multiplied by the constant 1.4276 to convert to milligrams per liter (mg/L).
%       For mass (mg) -> mol => divide by 1.0E3 and divide by 32
%       Then divide by 1.024 (~density) for l-1 -> kg-1
%
% *** DO re-grid from WOA to GENIE -- N ********************************* %
%
% set netCFD file
loc_str_var = 'n_an'; % netCDF variable name
loc_str_nc = 'nitrate_annual_1deg'; % netCDF filename omitting '.nc'
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
disp('... PROCESSING NO3');
make_regrid_WOA2GENIE(data(:,:,:),grid_depth,'n_an','mole_concentration_of_nitrate_in_sea_water','mol kg-1',(1.0/1.024E6),loc_n_i,loc_n_j,loc_n_k,mask,[loc_str_mask '.n_an'],'GENIE');
%
% *** DO re-grid from WOA to GENIE -- Si ******************************** %
%
% set netCFD file
loc_str_var = 'i_an'; % netCDF variable name
loc_str_nc = 'silicate_annual_1deg'; % netCDF filename omitting '.nc'
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
disp('... PROCESSING Si');
make_regrid_WOA2GENIE(data(:,:,:),grid_depth,'i_an','mole_concentration_of_silicate_in_sea_water','mol kg-1',(1.0/1.024E6),loc_n_i,loc_n_j,loc_n_k,mask,[loc_str_mask '.i_an'],'GENIE');
%
% *********************************************************************** %
