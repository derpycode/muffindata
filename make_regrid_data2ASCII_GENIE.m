function [FNDATA] = make_regrid_data2ASCII_GENIE(PDATA,PIMAX,PJMAX,POFF,PMASK,PNAME,PCOAST)
% make_regrid_data2ASCII_GENIE
%
%   ***********************************************************************
%   *** Regrid discrete data to (GENIE) gridded ASCII *********************
%   *** 2D! ***************************************************************
%   ***********************************************************************
%
%   make_regrid_data2ASCII_GENIE(PDATA,PNAME,PNi,PNj)
%   (equal j increment grid)
%   'make_regrid_data2netcdf.nc' and takes 7 arguments:
%
%   PDATA [STRING] (e.g. 'data_d30Si.dat')
%   --> the dataset name (including extension)
%   PIMAX [INTEGER] (e.g. 36)
%   --> the number of increments in igitude
%   --> e.g. 360 for a 1 degree grid
%   PJMAX [INTEGER] (e.g. 36)
%   --> the number of increments in jitude
%   --> e.g. 180 for a 1 degree grid
%   POFF [INTEGER] (e.g. -260)
%   --> the longitude grid offset (degrees)
%   PMASK [PIMAXxPJMAX INTEGER ARRAY]
%   --> optional grid mask (the grid 'k1' topo file)
%   PNAME [STRING]
%   --> filename root
%   PCOAST [LOGICAL] (OPTIONAL)
%   --> whether to move masked near-coastal data to coastal cells
%
%   Example
%
%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   16/11/16: CREATED [from make_regrid_data2netcdf]
%   16/11/17: fixed to saved grid orientation
%   16/11/20: added passing of filename root
%   17/08/17: added simple diagnostics output
%   17/09/07: add parameter to allow near-coastal data to be assigned 
%             to the coast
%   17/09/11: finished 'coast' inclusion algorithm
%             made new passed parameter optional
%   17/09/24: fixed minor bug
%
%   ***********************************************************************

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% *** START ************************************************************* %
% 
disp(['START [make_regrid_data2ASCII_GENIE] >>>'])
%
% *** dummy variable processing ***************************************** %
%
% set passed parameters
str_data_filename = PDATA;
n_i = PIMAX;
n_j = PJMAX;
par_grid_i_offset = POFF;
str_mask_filename = PMASK;
str_name = PNAME;
if ~exist('PCOAST','var'), PCOAST = false; end
opt_coast = PCOAST;
%
% *** options *********************************************************** %
%
% no duplication
% => omit move if there are multiple adjacent coastal cells
opt_coast_nodupl = false; 
%
% *** create grid ******************************************************* %
%
% i (west boundary)
for i=1:n_i,
    axis_iedge(i) = (i-1)*(360.0/n_i) + par_grid_i_offset;
    axis_di(i)    = (360.0/n_i);
end
axis_iedge(n_i+1)   = (n_i)*(360.0/n_i) + par_grid_i_offset;
axis_imid           = axis_iedge(1:n_i) + 0.5*axis_di;
axis_ibnds(1:n_i,1) = axis_iedge(1:n_i);
axis_ibnds(1:n_i,2) = axis_iedge(2:n_i+1);
axis_ibnds          = axis_ibnds';
% j (south boundary)
for j=1:n_j,
    axis_jedge(j) = (180.0/pi)*asin((2*(j-1)/n_j) - 1.0);
    axis_jmid(j)  = (180.0/pi)*asin(((1 + 2*(j-1))/n_j) - 1.0);
end
axis_jedge(n_j+1)   = (180.0/pi)*asin((2*n_j/n_j) - 1.0);
axis_jbnds(1:n_j,1) = axis_jedge(1:n_j);
axis_jbnds(1:n_j,2) = axis_jedge(2:n_j+1);
axis_jbnds          = axis_jbnds';
%
% *** misc initialization *********************************************** %
%
% set null/fill value
loc_nullvalue = -9.9E19;
% set date
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
%
% *********************************************************************** %

% *********************************************************************** %
% *** READ IN & PROCESS DATA ******************************************** %
% *********************************************************************** %
%
% *** load data ********************************************************* %
%
if (exist(str_data_filename,'file') == 2),
    fid = fopen(str_data_filename);
    C = textscan(fid, '%f %f %f %s', 'CommentStyle', '%');
    data_raw = cell2mat(C(1:3));
    fclose(fid);
else
    disp(['ERROR: ', 'File: ', str_data_filename, ' does not exist anywhere in the MATLAB directory path.']);
    return;
end
% remove NaNs
data_raw(find(isnan(data_raw(:,3))),:) = [];
% determine data size
data_size = size(data_raw(:,:));
nmax = data_size(1);
nmax_raw = nmax;
%
% load mask & reorient
if ~isempty(str_mask_filename),
    data_mask = load(str_mask_filename,'-ascii');
    data_mask = rot90(data_mask,-1);
else
    data_mask = zeros(n_i,n_j) + 1;
end
%
% *** fliter data ******************************************************* %
%
% filter lon,lat values
for n = 1:nmax,
    if (data_raw(n,1) < par_grid_i_offset), data_raw(n,1) = data_raw(n,1) + 360.0; end
    if (data_raw(n,1) > (par_grid_i_offset + 360.0)), data_raw(n,1) = data_raw(n,1) - 360.0; end
    if ( (data_raw(n,2) > 90.0) || (data_raw(n,2) < -90.0) ),
        disp(['ERROR: ', 'Something odd with the latitude values in file: ', str_data_filename, ' -- maybe latitude is not the 2nd data column as it should be?']);
        return;
    end
end
%
% *** re-grid data ****************************************************** %
%
% create empty data vectors
data_vector = zeros(nmax,3);
count_vector = zeros(nmax,1);
% create empty 3D data grids
data_gridded = NaN*zeros(n_i,n_j);
%%%data_gridded = data_gridded + loc_nullvalue;
data_gridded_rawcount = zeros(n_i,n_j);
% create empty 2D data grid (for saving data density))
data_gridded_count = zeros(n_i,n_j);
% determine grid locations
% NOTE: subtract 1 from upper bound values of n for i and j
%       i.e. a point that lies in the cell with i index n, fits criteria:
%       j >= j_edge(n)
%       j < j_edge(n+1)
% NOTE: deal with special cases of i or j grid boundary values
for n = 1:nmax,
    loc_i_n = intersect(find(data_raw(n,1)>=axis_iedge(:)),find(data_raw(n,1)<axis_iedge(:))-1);
    if (abs(data_raw(n,1)) == abs(par_grid_i_offset)), loc_i_n = 1; end
    if (abs(data_raw(n,1)) == (360-abs(par_grid_i_offset))), loc_i_n = n_i; end
    loc_j_n = intersect(find(data_raw(n,2)>=axis_jedge(:)),find(data_raw(n,2)<axis_jedge(:))-1);
    if (data_raw(n,2) == -90), loc_j_n = 1; end
    if (data_raw(n,2) == 90), loc_j_n = n_j; end
    if (isempty(loc_i_n*loc_j_n)),
        disp(['ERROR: ', 'Failed to regrid ... check i,j values and/or grid resolution choice. Error info:']);
        disp(['n = ',num2str(n),' : ',num2str(data_raw(n,3)),' @ ','(',num2str(data_raw(n,1)),',',num2str(data_raw(n,2)),')',' == ','(',num2str(loc_i_n),',',num2str(loc_j_n),')']);
        return;
    end
    %disp(['n = ',num2str(n),' : ',num2str(data_raw(n,3)),' @ ','(',num2str(data_raw(n,1)),',',num2str(data_raw(n,2)),')',' == ','(',num2str(loc_i_n),',',num2str(loc_j_n),')']);
    data_vector(n,:) = [loc_i_n; loc_j_n; data_raw(n,3)]';
    data_gridded_rawcount(loc_i_n,loc_j_n) = data_gridded_rawcount(loc_i_n,loc_j_n) + 1;
end
% store raw data diagnostics
diag_raw_n    = nmax;
diag_raw_sum  = sum(data_vector(:,3));
diag_raw_mean = diag_raw_sum/diag_raw_n;
%
% *** deal with near coastal data *************************************** %
%
% move land (near coastal) points to ocean
% NOTE: allow data to be duplicated in multiple coastal cells if
%       there are multiple adjacent ocean cells, IF SELECTED
% NOTE: <data_mask> must be transposed when
%       fun_grid_cell_neighbour_search_v is called
% NOTE: fun_grid_cell_neighbour_search_v returns (j,i) format
% NOTE: last parameter in fun_grid_cell_neighbour_search_v sets
%       whether or not to include diagnoals in the search
if opt_coast,
    n = 1;
    while (n <= nmax),
        loc_i = data_vector(n,1);
        loc_j = data_vector(n,2);
        if (data_mask(loc_i,loc_j) == 0),
            [vloc] = fun_grid_cell_neighbour_search_v(loc_j,loc_i,data_mask',true);
            if opt_coast_nodupl,
                loc_size = 0;
            else
                loc_size = size(vloc,1);
            end
            if (loc_size > 0),
                for o = 1:loc_size,
                    data_vector = [data_vector; [vloc(o,2) vloc(o,1) data_vector(n,3)]];
                end
                data_vector(n,:) = [];
                nmax = nmax-1;
            else
                n = n + 1;
            end
        else
            n = n + 1;
        end
    end
end
% update nmax
nmax = size(data_vector,1);
%
% *** average data ****************************************************** %
%
% NOTE: remove duplicate locations 
%       (having added the value from there to the running average)
n = 1;
while (n <= nmax),
    loc_coordinate  = data_vector(n,1:2);
    count_vector(n) = 1;
    m = n + 1;
    n_dup = 0;
    while (m <= nmax),
        %disp(['n = ',num2str(n),' : ','m = ',num2str(m)])
        if (data_vector(m,1:2) == loc_coordinate),
            %disp(['Duplicate @ ','n = ',num2str(n),', ','m = ',num2str(m)])
            n_dup = n_dup + 1;
            data_vector(n,3) = (n_dup*data_vector(n,3) + data_vector(m,3))/(n_dup + 1);
            count_vector(n) = count_vector(n) + 1;
            data_vector(m,:) = [];
            nmax = nmax - 1;
        else
            m = m + 1;
        end
    end
    n = n + 1;
end
%
% *** prepare gridded data ********************************************** %
%
% populate 2D grid
loc_sum = 0;
for n = 1:nmax,
    loc_lon_n = data_vector(n,1);
    loc_lat_n = data_vector(n,2);
    data_gridded(loc_lon_n,loc_lat_n) = data_vector(n,3);
    data_gridded_count(loc_lon_n,loc_lat_n) = count_vector(n);
    loc_sum = loc_sum + data_vector(n,3);
end
% store re-gridded data diagnostics
diag_grid_n    = nmax;
diag_grid_sum  = loc_sum;
diag_grid_mean = loc_sum/diag_grid_n;
%
% *********************************************************************** %

% *********************************************************************** %
% *** WRITE DATA ******************************************************** %
% *********************************************************************** %
%
% *** 2D **************************************************************** %
%
% plot mask
if ~isempty(data_mask),
    figure;
    plot_2dgridded(data_mask(:,:)',1+1,'',[str_name, '.MASK'],'grid mask');
end
% plot data
figure;
plot_2dgridded(data_gridded(:,:)',abs(loc_nullvalue),'',[str_name, '.DATA'],'grided data');
% plot data count
figure;
plot_2dgridded(data_gridded_count(:,:)',nmax_raw+1,'',[str_name, '.DATACOUNT'],'data count');
% plot data + mask
if ~isempty(data_mask),
    data_tmp = data_gridded;
    for i = 1:n_i,
        for j = 1:n_j,
            if (data_mask(i,j) == 0),
                data_tmp(i,j) = abs(loc_nullvalue);
        end
    end
    figure;
    plot_2dgridded(data_tmp(:,:)',abs(loc_nullvalue),'',[str_name, '.DATAMASK'],'masked gridded data');
end
% create filename
% NOTE: re-orientate to (j,i) (from (lon,lat)
str_filename_data = [str_name, '.GRIDDED_DATA', '.', str_date, '.dat'];
str_filename_datacount = [str_name, '.GRIDDED_DATACOUNT', '.', str_date, '.dat'];
% save gridded data
fprint_2DM(data_gridded(:,:)',data_mask(:,:)',str_filename_data,'%14.6e','%14.1f',false,false);
fprint_2DM(data_gridded_count(:,:)',data_mask(:,:)',str_filename_datacount,'%14.6e','%14.1f',false,false);
%
% *** 1D **************************************************************** %
%
% write out data in 1D format
% filted masked data from vector
for n = 1:nmax
    if (data_mask(data_vector(n,1),data_vector(n,2)) == 0)
        data_vector(n,:) = [];
        nmax = nmax - 1;
    end
end
% write data vector [lon,lat,value]
str_filename_datavector = [str_name, '.ijMEANVALUE', '.', str_date, '.dat'];
fprint_1Dn_d([data_vector(:,1) data_vector(:,2) data_vector(:,3)],str_filename_datavector);
%
% *** DIAGNOSTICS ******************************************************* %
%
fid = fopen([str_name '.' str_date '.txt'],'w');
fprintf(fid,'%s\n','##################################################################################');
fprintf(fid,'\n');
fprintf(fid,'%s\n',['# of raw data points     = ',num2str(diag_raw_n)]);
fprintf(fid,'%s\n',['sum of raw data          = ',num2str(diag_raw_sum)]);
fprintf(fid,'%s\n',['mean of raw data         = ',num2str(diag_raw_mean)]);
fprintf(fid,'\n');
fprintf(fid,'%s\n',['# of gridded data points = ',num2str(diag_grid_n)]);
fprintf(fid,'%s\n',['sum of gridded data      = ',num2str(diag_grid_sum)]);
fprintf(fid,'%s\n',['mean of gridded data     = ',num2str(diag_grid_mean)]);
fprintf(fid,'\n');
fprintf(fid,'%s\n','##################################################################################');
fclose(fid);
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
% RETURN FUNCTION VALUE
FNDATA = data_gridded;
%
% END
disp(['<<< END [make_regrid_data2ASCII_GENIE]']);
disp([' ']);
%
% *********************************************************************** %
