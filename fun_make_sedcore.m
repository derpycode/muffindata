function [OUT_SEDCORE] = fun_make_sedcore(DUM_NI,DUM_NJ,DUM_ORIGIN,DUM_EA,IN_PATH,IN_FILENAME)
%
%   **********************************************************************
%   *** fun_make_sedcore **************************************************
%   ***********************************************************************
%
%   Create sedcore input file 
%
%   \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
%   /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
%
%   ***********************************************************************

% *********************************************************************** %
% *** CREATE SEDCORE **************************************************** %
% *********************************************************************** %
%
% *** initialize ******************************************************** %
%
% dummy variables
n_i = DUM_NI;
n_j = DUM_NJ;
grid_lon_origin = DUM_ORIGIN;  % grid origin (degrees from zero)
grid_equalarea  = DUM_EA;      % equal area grid? [logical]
str_inpath      = IN_PATH;     % [string]
str_infilename  = IN_FILENAME; % [string]
%
% *** load sedcore location data **************************************** %
%
% read data(!)
[file_data] = fun_read_file([str_inpath str_infilename]);
% extract data cell array, and vector format
cdata    = file_data.cdata;
v_format = file_data.vformat;
% determine number of rows and columns
n_size    = size(cdata);
n_rows    = n_size(1);
n_columns = length(v_format);
% initialize valid format flag
flag_format = false;
% parse call array data
% NOTE: valid formats are:
%       3: lon,lad,depth
%       4: lon,lat,depth,value
%       4: lon,lad,depth,label
%       5: lon,lat,depth,value,label
% NOTE: labels will not be parsed and saved in output (for now)
switch n_columns
    case 3
        if (sum(v_format(1:3)) == 0)
            % lon, lat, value
            data_raw = cell2mat(cdata(:,1:3));
            % flag for a valid format
            flag_format = true;
        end
    case 4
        if (sum(v_format(1:4)) == 0)
            % lon, lat, depth, value
            data_raw  = cell2mat(cdata(:,1:4));
            % flag for a valid format
            flag_format      = true;
        elseif (sum(v_format(1:3)) == 0)
            % lon, lat, value, LABEL
            data_raw  = cell2mat(cdata(:,1:3));
            % flag for a valid format
            flag_format      = true;
        end
    case 5
        if (sum(v_format(1:4)) == 0)
            % lon, lat, depth, value, LABEL
            data_raw  = cell2mat(cdata(:,1:4));
            % flag for a valid format
            flag_format      = true;
        end
    otherwise
        % NOTHING
end
% report invalid format and exit
if (~flag_format)
    disp([' ']);
    disp([' * ERROR: Data format not recognized:']);
    disp(['   Number of data columns found == ' num2str(n_columns)']);
    disp(['   Columns must be: comma/tab/space-separated.']);
    disp([' ']);
    fclose(fid);
    return;   
end
% determine final data size:
% (1) set number of data rows
% (2) pad data value if necessary
data_size = size(data_raw(:,:));
nmax = data_size(1);
cmax = data_size(2);
% check for incomplete file read
if (nmax < n_rows)
    disp([' ']);
    disp([' * ERROR: Problem in data format (e.g. string must not contain spaces):']);
    disp(['   Number of data rows read == ' num2str(nmax)]);
    disp(['   Number of lines in file  == ' num2str(n_rows)]);
    disp([' ']);
    return;
end
% check for mixed up lon-lat ... i.e. not (LON, LAT) format
if ( (min(data_raw(:,2)) < -90) || (max(data_raw(:,2)) > 90) )
    loc_tmp = overlaydata_raw(:,1);
    data_raw(:,1) = data_raw(:,2);
    data_raw(:,2) = loc_tmp;
    disp([' ']);
    disp([' * WARNING: lon-lat is not in (LON, LAT) column order format:']);
    disp(['   Swapping ...'])
    disp([' ']);
end
% create (i,j) from (lon,lat)
% precondition lon
for n = 1:nmax
    if (data_raw(n,1) >= (360.0 + grid_lon_origin))
        data_raw(n,1) = data_raw(n,1) - 360.0;
    end
    if (data_raw(n,1) < grid_lon_origin)
        data_raw(n,1) = data_raw(n,1) + 360.0;
    end
end
% copy data_raw
data_ij = data_raw;
% convert (lon,lat) overlay data to (i,j)
% NOTE: function 'calc_find_ij' takes input in order: (lon,lat)
%       i.e., the same as the raw overlay data, which is (lon,lat) (i.e., (i,j)) format
for n = 1:nmax
    data_ij(n,1:2) = calc_find_ij(data_raw(n,1),data_raw(n,2),grid_lon_origin,n_i,n_j);
end
% grid (and average per cell) data
% NOTE: data vector length is re-calculated and the value of nmax reset
m=0;
for i = 1:n_i
    for j = 1:n_j
        samecell_locations = find((data_ij(:,1)==i)&(data_ij(:,2)==j));
        samecell_n = length(samecell_locations);
        if (samecell_n > 0)
            m=m+1;
            if (cmax == 3)
                samecell_mean = [mean(data_ij(samecell_locations,3))];
            elseif (cmax == 4)
                samecell_mean = [mean(data_ij(samecell_locations,3)), mean(data_ij(samecell_locations,4))];
            end
            % record re-gridded data
            data(m,1)      = i;
            data(m,2)      = j;
            data(m,3:cmax) = samecell_mean;
            % record averaging (samecell_n)
            data(m,cmax+1) = samecell_n;
        end
    end
end
%
% *** RETURN ************************************************************ %
%
OUT_SEDCORE = data;
%
% *** END *************************************************************** %
%
close all;
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %

end