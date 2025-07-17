% -------------------------------------------------------------------------
% load_genie_matrix
% 
% Description: loads matrix data in coordinate format as netcdf file
%              produced in genie and writes out matrices and index to
%              matlab file
% Author: Jamie D. Wilson 08/07/2019
%
% Inputs: 
% - netcdf_filename: filename of matrix data 
%                    (usually: 'transport_matrix_COO.nc')
% - save_TM_name: name of output file
% -------------------------------------------------------------------------

function [ ] = load_genie_matrix ( netcdf_filename , save_TM_name )


% load data from file
m_row=ncread('transport_matrix_COO.nc','coo_row');
m_col=ncread('transport_matrix_COO.nc','coo_col');
m_season=ncread('transport_matrix_COO.nc','coo_avg_n');
m_data=ncread('transport_matrix_COO.nc','coo_val');

% get wet grid-box number
nb=numel(unique(m_row));
if numel(unique(m_col))~=nb
    disp('Number of grid-boxes does not add up')
    return
end

% detect whether seasonal or annual and make TMs from file
n_seasons=numel(unique(m_season));
if n_seasons>1
    varnames_exclude={'A'};
    for n=1:n_seasons
        varname = genvarname({'A'},varnames_exclude);
        eval([...
            varname{1}...
            '=sparse(double(m_row(m_season==n)),double(m_col(m_season==n)),m_data(m_season==n),nb,nb,numel(m_data(m_season==n)));'...
            ]);
        varnames_exclude(n+1)=varname;
        clear varname
    end
elseif n_seasons==1
    A=sparse(m_row,m_col,m_data,nb,nb,numel(m_data));
end
    
% create associated index file
i=ncread('transport_matrix_COO.nc','index_i');
j=ncread('transport_matrix_COO.nc','index_j');
k=ncread('transport_matrix_COO.nc','index_k');
    
v_index=struct(...
'i',ones(nb,1),...
'j',ones(nb,1),...
'k',ones(nb,1),...
'rk',ones(nb,1)...
);

n_k=max(k);

v_index.i=i;
v_index.j=j;
v_index.k=k;
v_index.rk=n_k-k+1;


% create extra indexing variables
Ii=find(v_index.rk>min(v_index.rk));
Ib=find(v_index.rk==min(v_index.rk));

% save files
save(save_TM_name,'A*')
save('matrix_vars','v_index','Ii','Ib','nb')

end