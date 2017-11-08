function []  = make_grid_border_mask(PK1,PDIAG,PDEPTH)
%
%%

% *********************************************************************** %
% *** CREATE BORDER MASK ************************************************ %
% *********************************************************************** %
%
% *** initialize ******************************************************** %
%
% copy passed parameters
str_k1   = PK1;
opt_diag = PDIAG;
n_depth  = PDEPTH;
% set default paths
str_dirin  = pwd;
str_dirout = pwd;
% set filename
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
str_nameout = [str_k1];
% load k1 and process file
[grid_k1,grid_mask,imax,jmax] = fun_read_k1(pwd,str_k1,'k1');
% icreate & nitialzie borders array
grid_bdrs = zeros(jmax,imax);
%
% *** iterate through search depth ************************************** %
%
for n = 1:n_depth
    % search across grid
    for i = 1:imax
        for j = 1:jmax
            % test for cell being land
            if ~grid_mask(j,i),
                % look for ocean cells and mark borders
                [grid_bdrs] = fun_grid_cell_neighbour_search(j,i,grid_mask,grid_bdrs,opt_diag);
            end
        end
    end
    % update mask
    % => make border cell as 'land'
    grid_mask = grid_mask & ~grid_bdrs;
end
% convert from logical to real and replace mask with borders
grid_mask = double(grid_bdrs);
%
% *** plot & save mask ************************************************** %
%
% plot orignal k1 grid
plot_2dgridded(flipud(grid_k1),9999,'',[str_nameout '.originalk1'],['original grid']);
% plot borders mask
plot_2dgridded(flipud(grid_mask),9999,'',[str_nameout '.bordersmask'],['mask']);
% save borders mask
fprint_2DM(grid_mask(:,:),[],[str_nameout, '.', str_date, '.bordersmask.dat'],'%4.1f','%4i',true,false);
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
%%%
%
% *********************************************************************** %

