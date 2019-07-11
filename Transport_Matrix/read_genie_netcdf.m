%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read_genie_netcdf
% Description: reads netcdf files from GENIE and processes for use with 
% matrices.  Data is last time slice of run, NaNs for missing values and
% dimensions are permuted (k,j,i).  
%
% Inputs:
% path                 - path to folder i.e. exp1/biogem/ 
%                        (leave empty if netcdf is sitting in current folder)
% vec_flag             - flag for 3D/vector output, 1: vector, 0: 3D 
% vec_info             - indexing for vector output as structure
%                        (i.e. v_index)
% varargin             - list of variable names to extract from netcdf
%
% Outputs:
% varargin             - names of variable outputs
%
% Author: J.D.Wilson 26/03/2014
%
% Example:              
%[PO4,SAL,LON,LAT]=read_genie_netcdf('',1,v_index,'ocn_PO4','ocn_sal','lon','lat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [ varargout ] = read_genie_netcdf( path , vec_flag , vec_info , varargin )

narg=size(varargin,2);    % number of files to get from netcdf

%newpath=cat(2,path,'fields_biogem_3d.nc');
newpath=path;

% process and output data
for n=1:narg
    
    temp_data=extract_data(newpath,varargin{n});
    
    temp_data=temp_data(:,:,:,end);             % choose last time slice
    temp_data(temp_data>1e30)=NaN;              % missing values -> NaNs
        
    if ndims(temp_data)==3
        temp_data=permute(temp_data,[3 2 1]);   % permute dimensions for vectorising
    end
    
    % output vector data if selected
    if vec_flag==1
        v_index=vec_info;
        varargout{n}=f2v(temp_data,v_index.i,v_index.j,v_index.rk);
    else
        varargout{n}=temp_data;
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%% SUB FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basic function to read in data from netcdf file
function [ file_out ] = extract_data ( netcdf_path , variable_name )


% netcdf file
fname=netcdf_path;

% open netcdf file
ncid=netcdf.open(fname,'NC_NOWRITE');

% retrieve variable data
file_out=netcdf.getVar(ncid,netcdf.inqVarID(ncid,variable_name));

% close netcdf file
netcdf.close(ncid);

% end of nested function
end



end

