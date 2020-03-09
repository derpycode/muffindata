function [] = make_regrid_data2netcdf_GENIE(PDATA,PNAME,PLNAME,PUNITS,PLONO,PLATO,PZO,POFF,PMASK)
% make_regrid_data2netcdf_GENIE
%
%   ***********************************************************************
%   *** Regrid discrete data to GENIE netCDF format, regular grid *********
%   *** (3D) + ASCII forcing compatable outputs ***************************
%   ***********************************************************************
%
%   make_regrid_data2netcdf_GENIE(PDATA,PNAME,PLNAME,PUNITS,PNLON,PNLAT,PZO,POFF,PMASK)
%   (equal lat increment grid)
%   'make_regrid_data2netcdf_GENIE.nc' and takes 7 arguments:
%
%   PDATA [STRING] (e.g. 'data_d30Si.dat')
%   --> the dataset name (including extension)
%   PNAME [STRING] (e.g., 'd30Si')
%   --> short data name
%   PLNAME [STRING] (e.g., 'observed silicate isotope value in units of per mil')
%   --> long data name (description)
%   PUNITS [STRING] (e.g., '1')
%   --> data units
%       NOTE: for isotopes, the units name needs to be '1'
%   PIMAX [INTEGER] (e.g. 360)
%   --> the number of increments in longitude
%   --> e.g. 36
%   PJMAX [INTEGER] (e.g. 180)
%   --> the number of increments in latitude
%   --> e.g. 36
%   PKMAX [INTEGER] (e.g. 16)
%   --> vertical grid resolution
%   POFF [INTEGER] (e.g. -260)
%   --> the longitude grid offset (degrees)
%   PMASK [PIMAXxPJMAX INTEGER ARRAY]
%   --> optional grid mask (the grid 'k1' topo file)
%
%   Example
%           make_regrid_data2netcdf_GENIE('data_d30Si.dat','d30Si','observed water column silicon isotope value in units of per mil','1',360,180,'WOA');
%           will grid the d30Si data in file: 'data_d30Si.dat'
%           to a WOA-compatable grid (1 degree lon, lat)
%
%   NOTE: opt_benthic controls an important behavior:
%         false => the data is re-gridded straight, and any data lying
%                  deeper than the deepest depth at a particular grid point
%                  is ignore
%         true  => deeper data is enfolded into a benthic mean 
%         BUT, data deeper than the deepest grid level (anywhere)
%              is regrided regardless
%
%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   17/09/19: CREATED -- copied from make_regrid_data2netcdf
%   17/09/25: completed initial revision
%   17/09/26: edited output file string
%   17/09/28: minor edits to how masks are applied and frocing ASCII saved
%   17/12/29: changed default behavior to using data deeper than the
%             maximum depth grid point depth
%   17/12/30: various benthic bug-fixes
%
%   ***********************************************************************

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% *** START ************************************************************* %
% 
disp(['START [make_regrid_data2netcdf_GENIE] >>>'])
%
% *** dummy variable processing ***************************************** %
%
% set passed parameters
str_data_filenamein = PDATA;
str_dataname = PNAME;
str_dataname_long = PLNAME;
str_data_units = PUNITS;
n_lon = PLONO;
n_lat = PLATO;
n_D = PZO;
par_grid_lon_offset = POFF;
str_mask_filename = PMASK;
% set data strings
if isempty(str_dataname), str_dataname = str_data_filename; end
if isempty(str_dataname_long), str_dataname_long = str_data_filename; end
if isempty(str_data_units), str_data_units = 'n/a'; end
if ~isempty(str_mask_filename),
    str_data_filenameout = [str_data_filenamein '_GENIE.' str_mask_filename];
else
    str_data_filenameout = [str_data_filenamein '_GENIE'];
end
%
% *** misc (local) parameters ******************************************* %
%
% options & parameters
opt_equalarea = true;  % use equal area grid
opt_benthic   = true; % average in data deeper than bottom of depth grid
max_D = 5000.0;
% set date
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
% set null/fill value
loc_nullvalue = -9.9E19;
%
% *** create grid ******************************************************* %
%
% create GENIE grid
[go_lonm,go_lone,go_latm,go_late,go_dm,go_de] = make_genie_grid(n_lon,n_lat,n_D,max_D,par_grid_lon_offset,opt_equalarea);
% convert to naming convention consistent with original function
% longitude
loc_axis_lonmid  = go_lonm;
loc_axis_lonedge = go_lone;
loc_axis_lonbnds(1:n_lon,1) = loc_axis_lonedge(1:n_lon);
loc_axis_lonbnds(1:n_lon,2) = loc_axis_lonedge(2:n_lon+1);
loc_axis_lonbnds = loc_axis_lonbnds';
% latitude
loc_axis_latmid  = go_latm;
loc_axis_latedge = go_late;
loc_axis_latbnds(1:n_lat,1) = loc_axis_latedge(1:n_lat);
loc_axis_latbnds(1:n_lat,2) = loc_axis_latedge(2:n_lat+1);
loc_axis_latbnds = loc_axis_latbnds';
% depth
% NOTE: depth must be in reverse order (shallow to deep)
loc_axis_Dmid  = fliplr(go_dm);
loc_axis_Dedge = fliplr(go_de);
loc_axis_Dbnds(1:n_D,1) = loc_axis_Dedge(1:n_D);
loc_axis_Dbnds(1:n_D,2) = loc_axis_Dedge(2:n_D+1);
loc_axis_Dbnds = loc_axis_Dbnds';
%
% *** load mask ********************************************************* %
%
% NOTE: mask and grid in 'visual' array format 
%      (i.e. latitude as columns, counting from top to bottom)
%      [remember in MATLAB :: (rows,columns) !]
% NOTE: process k1 format
if ~isempty(str_mask_filename),
    loc_data = load(str_mask_filename,'-ascii');
    % test for extended k1 format
    if (size(loc_data) == [n_lon+2 n_lat+2]),
        loc_data(1,:)   = [];
        loc_data(end,:) = [];
        loc_data(:,1)   = [];
        loc_data(:,end) = [];
    end
    data_grid = loc_data;
    % determined whether k1 or mask
    if (max(max(loc_data)) > 1.0),
        data_mask = loc_data;
        data_mask(find(data_mask <= n_D)) = 1.0;
        data_mask(find(data_mask >= 90))  = 0.0;
    end
else
    loc_data  = zeros(n_lon,n_lat) + 1;
    data_grid = loc_data;
    data_mask = loc_data;
end
data_mask_nan = data_mask;
data_mask_nan(find(data_mask==0.0)) = NaN;
% plot something ...
str_filename = [str_data_filenameout '.' str_date '.GRID'];
plot_2dgridded(flipud(data_grid),n_D,'',str_filename,['grid']);
str_filename = [str_data_filenameout '.' str_date '.MASK'];
plot_2dgridded(flipud(data_mask),2.0,'',str_filename,['mask']);
%
% *********************************************************************** %

% *********************************************************************** %
% *** READ IN & PROCESS DATA ******************************************** %
% *********************************************************************** %
%
% *** load data ********************************************************* %
%
if (exist(str_data_filenamein,'file') == 2),
    fid = fopen(str_data_filenamein);
    C = textscan(fid, '%f %f %f %f %s', 'CommentStyle', '%');
    data_raw = cell2mat(C(1:4));
    fclose(fid);
    data_size = size(data_raw(:,:));
    nmax=data_size(1);
    if (length(find(isnan(data_raw))) > 0),
        disp(['WARNING: ', 'The data file: ', str_data_filenamein, ' does not contain numeric data in a consistent 4-column (lon,lat,depth,value), space-seperated format, or it contains NaNs. Recommend is to remove NaN containing lines (data=rmmissing(data);).']);
    end
else
    disp(['ERROR: ', 'File: ', str_data_filenamein, ' does not exist anywhere in the MATLAB directory path.']);
    return;
end
% save raw data count
diag_data_n_raw = nmax;
%
% *** filter data ******************************************************* %
%
% % option:
% % remove any data deeper that the bottom of the depth grid and update nmax
% % *or*
% % set to max_D
% n = 1;
% while (n <= nmax),
%     if (data_raw(n,3) > max_D),
%         if opt_benthic,
%             data_raw(n,3) = max_D;
%             n = n + 1;           
%         else 
%             data_raw(n,:) = [];
%             nmax = nmax - 1;
%         end
%     else
%         n = n + 1;
%     end
% end
% % move anything exactly at the bottom, slightly shallower ...
% data_raw(find(data_raw(:,3)==max_D),3) = max_D - 0.001;
% filter lon,lat values
for n = 1:nmax,
    if (data_raw(n,1) < par_grid_lon_offset), data_raw(n,1) = data_raw(n,1) + 360.0; end
    if (data_raw(n,1) >= (par_grid_lon_offset + 360.0)), data_raw(n,1) = data_raw(n,1) - 360.0; end
    if ( (data_raw(n,2) > 90.0) || (data_raw(n,2) < -90.0) ),
        disp(['ERROR: ', 'Something odd with the latitude values in file: ', str_data_filenamein, ' -- maybe latitude is not the 2nd data column as it should be?']);
        return;
    end
    if (data_raw(n,3) < 0.0),
        disp(['ERROR: ', 'Something odd with the depth values in file: ', str_data_filenamein, ' -- maybe depth is given as a negative height or depth is is not the 3rd data column as it should be?']);
        return;
    end
end
%
% *** re-grid data ****************************************************** %
%
% create empty data vectors
data_vector = zeros(nmax,5);
% create empty 3D data grids
data_gridded = zeros(n_lon,n_lat,n_D);
data_gridded = data_gridded + loc_nullvalue;
data_gridded_rawcount = zeros(n_lon,n_lat,n_D);
% create empty 2D data grid (for saving data density))
data_gridded_2Dcount = zeros(n_lon,n_lat);
% determine grid locations
% NOTE: subtract 1 from upper bound values of n for lon and lat
%       i.e. a point that lies in the cell with lon index n, fits criteria:
%       lat >= lat_edge(n)
%       lat < lat_edge(n+1)
% NOTE: deal with special cases of lon or lat grid boundary values
% NOTE: remember, the k grid vector is counted with 1 == surface
for n = 1:nmax,
    loc_lon_n = intersect(find(data_raw(n,1)>=loc_axis_lonedge(:)),find(data_raw(n,1)<loc_axis_lonedge(:))-1);
    if (abs(data_raw(n,1)) == abs(par_grid_lon_offset)), loc_lon_n = 1; end
    loc_lat_n = intersect(find(data_raw(n,2)>=loc_axis_latedge(:)),find(data_raw(n,2)<loc_axis_latedge(:))-1);
    if (data_raw(n,2) == -90),   loc_lat_n = 1; end
    if (data_raw(n,2) == 90),    loc_lat_n = n_lat; end
    loc_D_n   = intersect(find(data_raw(n,3)>=loc_axis_Dedge(:)),find(data_raw(n,3)<loc_axis_Dedge(:))-1);
    if (data_raw(n,3) >= max_D), loc_D_n = n_D; end
    if (isempty(loc_lon_n*loc_lat_n*loc_D_n)),
        disp(['ERROR: ', 'Failed to regrid ... check lon,lat values and/or grid resolution choice. Error info:']);
        disp(['n = ',num2str(n),' : ',num2str(data_raw(n,4)),' @ ','(',num2str(data_raw(n,1)),',',num2str(data_raw(n,2)),',',num2str(data_raw(n,3)),')',' == ','(',num2str(loc_lon_n),',',num2str(loc_lat_n),',',num2str(loc_D_n),')']);
        return;
    end
    data_vector(n,:) = [loc_lon_n; loc_lat_n; loc_D_n; data_raw(n,4); 1]';
    data_gridded_rawcount(loc_lon_n,loc_lat_n,loc_D_n) = data_gridded_rawcount(loc_lon_n,loc_lat_n,loc_D_n) + 1;
end
%
% make conform to topo:
% move deeper depth levels that local deepest level, to local deepest level
% NOTE: k=1 is the surface in the data array indexing ...
loc_grid = fliplr(data_grid');
n = 1;
while (n <= nmax),
    % extract data from data vector
    loc_lon_n = data_vector(n,1);
    loc_lat_n = data_vector(n,2);
    loc_D_n   = data_vector(n,3);
    % convert from normal k grid counting with 1 == deep,
    % to 1 == surface
    loc_D_n = n_D - loc_D_n + 1;
    % extract local grid depth level
    loc_D_nmax = loc_grid(loc_lon_n,loc_lat_n);
    % test for land ELSE cell depth index value < grid value 
    if (loc_D_nmax >= 90),
        % remove data point
        data_vector(n,:) = [];
        nmax = nmax - 1;
    elseif (loc_D_n < loc_D_nmax),
        if opt_benthic,
            % set depth level equal to deepest level at that grid point
            data_vector(n,3) = n_D - loc_D_nmax + 1;
            n = n + 1;
        else
            % remove data point
            data_vector(n,:) = [];
            nmax = nmax - 1;        
        end
    else
        n = n + 1;
    end
end
% save valid data count
diag_data_n_valid = nmax;
%
% average data
% NOTE: remove duplicate locations (having added the value from there to the running average)
n = 1;
while (n <= nmax),
    loc_coordinate = data_vector(n,1:3);
    m = n + 1;
    n_dup = 0;
    while (m <= nmax),
        %disp(['n = ',num2str(n),' : ','m = ',num2str(m)])
        if (data_vector(m,1:3) == loc_coordinate),
            %disp(['Duplicate @ ','n = ',num2str(n),', ','m = ',num2str(m)])
            n_dup = n_dup + 1;
            data_vector(n,4) = (n_dup*data_vector(n,4) + data_vector(m,4))/(n_dup + 1);
            data_vector(n,5) = data_vector(n,5) + 1;
            data_vector(m,:) = [];
            nmax = nmax - 1;
        else
            m = m + 1;
        end
    end
    n = n + 1;
end
% save final data count
diag_data_n = nmax;
%
% populate 3D grid
for n = 1:nmax,
    loc_lon_n = data_vector(n,1);
    loc_lat_n = data_vector(n,2);
    loc_D_n   = data_vector(n,3);
    data_gridded(loc_lon_n,loc_lat_n,loc_D_n) = data_vector(n,4);
    data_gridded_2Dcount(loc_lon_n,loc_lat_n) = data_gridded_2Dcount(loc_lon_n,loc_lat_n) + 1;
end
%
% *** SAVE ASCII ******************************************************** %
%
% NOTE: re-orientate to (j,i) (from (lon,lat))
% NOTE: data starts in array top LH corner as: (par_grid_lon_offset, -90S)
% surface
loc_data = data_gridded(:,:,1);
if ( max(max(loc_data)) > loc_nullvalue)
    % save surface data
    str_filename = [str_data_filenameout '.SUR'];
    loc_data = flipud(loc_data');
    plot_2dgridded(flipud(loc_data),loc_nullvalue,'',[str_filename '.NOmask'],['data']);
    loc_data = data_mask.*loc_data;
    fprint_2DM(loc_data,data_mask,[str_filename '.dat'],'%14.6e','%14.1f',true,false);
    loc_data(find(loc_data == 0.0)) = loc_nullvalue;
    fprint_2DM(loc_data,[],[str_filename '.NULLmask.dat'],'%14.6e','%14.1f',true,false);
    loc_data = data_mask_nan.*loc_data;
    plot_2dgridded(flipud(loc_data),loc_nullvalue,'',[str_filename '.NaNmask'],['data']);
end
% benthic surface
% NOTE: k=1 is the surface in the data array indexing ...
% extract benthic data surface
loc_data  = zeros(n_lon,n_lat) + loc_nullvalue;
loc_grid = fliplr(data_grid');
grid_data_count = zeros(n_lon,n_lat);
n_sur = 0;
n_ben = 0;
for n = 1:nmax,
    % extract data from data vector
    loc_lon_n = data_vector(n,1);
    loc_lat_n = data_vector(n,2);
    loc_D_n   = data_vector(n,3);
    loc_n_n   = data_vector(n,5);
    % set local count
    loc_data_count = grid_data_count(loc_lon_n,loc_lat_n);
    % convert from normal k grid counting with 1 == deep,
    % to 1 == surface
    % => test needs to be of cell depth index value >= grid value
    loc_D_nmax = loc_grid(loc_lon_n,loc_lat_n);
    loc_D_nmax = n_D - loc_D_nmax + 1;
    % count surface grid points
    if (loc_D_n == 1),
        n_sur = n_sur+1;
    end
    % benthic
    if (loc_D_n == loc_D_nmax),
        loc_data(loc_lon_n,loc_lat_n) = data_gridded(loc_lon_n,loc_lat_n,loc_D_n);
        grid_data_count(loc_lon_n,loc_lat_n) = loc_n_n;
        n_ben = n_ben+1;
    end
end
if ( max(max(loc_data)) > loc_nullvalue)
    % save benthic data
    str_filename = [str_data_filenameout '.BEN'];
    loc_data = flipud(loc_data');
    plot_2dgridded(flipud(loc_data),loc_nullvalue,'',[str_filename '.NOmask'],['data']);
    loc_data = data_mask.*loc_data;
    fprint_2DM(loc_data,data_mask,[str_filename '.dat'],'%14.6e','%14.1f',true,false);
    loc_data(find(loc_data == 0.0)) = loc_nullvalue;
    fprint_2DM(loc_data,[],[str_filename '.NULLmask.dat'],'%14.6e','%14.1f',true,false);
    loc_data = data_mask_nan.*loc_data;
    plot_2dgridded(flipud(loc_data),loc_nullvalue,'',[str_filename '.NaNmask'],['data']);
    % save benthic count data
    loc_data = flipud(grid_data_count');
    loc_data = data_mask_nan.*loc_data;
    plot_2dgridded(flipud(loc_data),loc_nullvalue,'',[str_filename '.data_cnt'],['data']);
end
%
% *** DIAGNOSTICS ******************************************************* %
%
fid = fopen([str_data_filenameout '.' str_date '.txt'],'w');
fprintf(fid,'%s\n','##################################################################################');
fprintf(fid,'\n');
fprintf(fid,'%s\n',['# of raw data points     = ',num2str(diag_data_n_raw)]);
% fprintf(fid,'%s\n',['sum of raw data          = ',num2str(diag_raw_sum)]);
% fprintf(fid,'%s\n',['mean of raw data         = ',num2str(diag_raw_mean)]);
fprintf(fid,'%s\n',['# of (grid) valid data   = ',num2str(diag_data_n_valid)]);
fprintf(fid,'\n');
fprintf(fid,'%s\n',['# of gridded data points = ',num2str(diag_data_n)]);
fprintf(fid,'%s\n',['# (check)                = ',num2str(sum(sum(data_gridded_2Dcount)))]);
fprintf(fid,'%s\n',['sum of gridded data      = ',num2str(sum(data_vector(:,4)))]);
fprintf(fid,'%s\n',['mean of gridded data     = ',num2str(sum(data_vector(:,4))/diag_data_n)]);
fprintf(fid,'\n');
fprintf(fid,'%s\n',['# of gridded sur points  = ',num2str(n_sur)]);
fprintf(fid,'%s\n',['# of gridded ben points  = ',num2str(n_ben)]);
fprintf(fid,'\n');
fprintf(fid,'%s\n','##################################################################################');
fclose(fid);
%
% *********************************************************************** %

% *********************************************************************** %
% *** WRITE netCDF FILE ************************************************* %
% *********************************************************************** %
%
% create netCDF file
ncid = netcdf.create([str_data_filenameout, '.', str_date, '.nc'],'NC_WRITE');
NC_GLOBAL = netcdf.getConstant('NC_GLOBAL');
% define global attributes
netcdf.putAtt(ncid,NC_GLOBAL,'title','Gridded data');
netcdf.putAtt(ncid,NC_GLOBAL,'long_title','Regridded on a regular lon-lat, irregular in depth, grid');
netcdf.putAtt(ncid,NC_GLOBAL,'comments','conversion of ascii format data by make_regrid_data2netcdf_GENIE.m');
netcdf.putAtt(ncid,NC_GLOBAL,'history','version as of 14/08/14');
netcdf.putAtt(ncid,NC_GLOBAL,'Conventions','CF-1.6');
netcdf.putAtt(ncid,NC_GLOBAL,'CreationDate',datestr(now,'yyyy/mm/dd HH:MM:SS'));
netcdf.putAtt(ncid,NC_GLOBAL,'CreatedBy',[getenv('username'), '@', getenv('computername')]);
netcdf.putAtt(ncid,NC_GLOBAL,'MatlabSource','make_regrid_data2netcdf_GENIE.m');
% define dimensions
dimid_D = netcdf.defDim(ncid,'depth',n_D);
dimid_lat = netcdf.defDim(ncid,'lat',n_lat);
dimid_lon = netcdf.defDim(ncid,'lon',n_lon);
di_bnds = netcdf.defDim(ncid,'nbounds',2);
di_time = netcdf.defDim(ncid,'time',1);
% define axes -- Z
varid = netcdf.defVar(ncid,'depth','double',dimid_D);
netcdf.putAtt(ncid,varid,'standard_name','depth');
netcdf.putAtt(ncid,varid,'long_name','Vertical distance below the surface');
netcdf.putAtt(ncid,varid,'units','m');
netcdf.putAtt(ncid,varid,'point_spacing','uneven');
netcdf.putAtt(ncid,varid,'Axis','Z');
netcdf.putAtt(ncid,varid,'positive','down');
netcdf.putAtt(ncid,varid,'_CoordinateAxisType','Depth');
netcdf.putAtt(ncid,varid,'bounds','D_bnds');
varid_D = varid;
% define axes -- D bnds
varid = netcdf.defVar(ncid,'D_bnds','double',[di_bnds, dimid_D]);
varid_Dbnds = varid;
% define axes -- Y
varid = netcdf.defVar(ncid,'lat','double',dimid_lat);
netcdf.putAtt(ncid,varid,'standard_name','latitude');
netcdf.putAtt(ncid,varid,'long_name','latitude of grid centre');
netcdf.putAtt(ncid,varid,'units','degrees_north');
netcdf.putAtt(ncid,varid,'point_spacing','even');
netcdf.putAtt(ncid,varid,'Axis','Y');
netcdf.putAtt(ncid,varid,'_CoordinateAxisType','Lat');
netcdf.putAtt(ncid,varid,'bounds','lat_bnds');
varid_lat = varid;
% define axes -- Y bnds
varid = netcdf.defVar(ncid,'lat_bnds','double',[di_bnds, dimid_lat]);
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
varid = netcdf.defVar(ncid,'lon_bnds','double',[di_bnds, dimid_lon]);
varid_lonbnds = varid;
% define data variable -- averaged data
varid = netcdf.defVar(ncid,str_dataname,'double',[dimid_lon, dimid_lat, dimid_D]);
netcdf.putAtt(ncid,varid,'name',str_dataname);
netcdf.putAtt(ncid,varid,'long_name',str_dataname_long);
netcdf.putAtt(ncid,varid,'missing_value',loc_nullvalue);
netcdf.putAtt(ncid,varid,'fillValue',loc_nullvalue);
netcdf.putAtt(ncid,varid,'units',str_data_units);
varid_data = varid;
% define data variable -- raw data density
varid = netcdf.defVar(ncid,'data_density','double',[dimid_lon, dimid_lat, dimid_D]);
netcdf.putAtt(ncid,varid,'name','data_density');
netcdf.putAtt(ncid,varid,'long_name','Density of data per grid cell prior to averaging');
netcdf.putAtt(ncid,varid,'missing_value',loc_nullvalue);
netcdf.putAtt(ncid,varid,'fillValue',loc_nullvalue);
netcdf.putAtt(ncid,varid,'units','m^-3');
varid_data_rawcount = varid;
% define data variable -- raw data density
varid = netcdf.defVar(ncid,'data_inventory','double',[dimid_lon, dimid_lat]);
netcdf.putAtt(ncid,varid,'name','data_inventory');
netcdf.putAtt(ncid,varid,'long_name','Inventory of processed data per grid point');
netcdf.putAtt(ncid,varid,'missing_value',loc_nullvalue);
netcdf.putAtt(ncid,varid,'fillValue',loc_nullvalue);
netcdf.putAtt(ncid,varid,'units','m^-2');
varid_data_2Dcount = varid;
% end definition
netcdf.endDef(ncid);
% write axes
netcdf.putVar(ncid,varid_D,loc_axis_Dmid);
netcdf.putVar(ncid,varid_lat,loc_axis_latmid);
netcdf.putVar(ncid,varid_lon,loc_axis_lonmid);
netcdf.putVar(ncid,varid_Dbnds,loc_axis_Dbnds);
netcdf.putVar(ncid,varid_latbnds,loc_axis_latbnds);
netcdf.putVar(ncid,varid_lonbnds,loc_axis_lonbnds);
% write data
netcdf.putVar(ncid,varid_data,data_gridded);
netcdf.putVar(ncid,varid_data_rawcount,data_gridded_rawcount);
netcdf.putVar(ncid,varid_data_2Dcount,data_gridded_2Dcount);
% close netCDF file
netcdf.close(ncid);
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
% END
close all;
disp(['<<< END [make_regrid_data2netcdf_GENIE]']);
disp([' ']);
%
% *********************************************************************** %
