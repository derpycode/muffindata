function [] = fun_process_ensemble_cycle(DUM_ENSEMBLE,DUM_YEAR,DUM_NAME)
%
%   ***********************************************************************
%   *** fun_process_ensemble_2d *******************************************
%   ***********************************************************************
%
%   DUM_ENSEMBLE == ensemble name 
%   (omitting the numerical code of individual ensemble members)
%   e.g. for the 2x3 (6 member) ensemble:
%   run00,run01,run02,run10,run11,run12
%   >> fun_process_ensemble_2d('run',9999.5,'myensemble')
%   DUM_YEAR == year of model data to extract
%   DUM_NAME == name to assign to the output files
%               (which could be different to and more meaningful than
%                e.g. DUM_ENSEMBLE)
%
%   NOTE: the ensmeble member number consists of 2 digits, 
%         each between 0 and 9
%         this code is represented in the comments by ##
%
%   For the selection of variables to extract and plot/analyse, 
%   some commond examples are given below
%   NOTE: these must be (uncommented and) copy-pasted in the appropriate 
%         section of the code
% 
%   time-series:
%
%   \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
%   % atmopsheric pCO2
%   m = m+1;
%   data(m).dataname = 'atm_pCO2';
%   data(m).datacol  = 3;
%   data(m).scale    = 1.0E6;
%   data(m).dataunit = 'pCO_2 (ppm)';
%   data(m).minmax   = [180 280];
%   % atmopsheric pCO2 d13C
%   m = m+1;
%   data(m).dataname = 'atm_pCO2_13C';
%   data(m).datacol  = 3;
%   data(m).scale    = 1.0;
%   data(m).dataunit = 'pCO_2 \delta^{13}C (o/oo)';
%   data(m).minmax   = [-7.0 -6.0];
%   % global POC export
%   m = m+1;
%   data(m).dataname = 'fexport_POC';
%   data(m).datacol  = 2;
%   data(m).scale    = 12.0E-15;
%   data(m).dataunit = 'POC export (PgC yr^{-1})';
%   data(m).minmax   = [6.0 12.0];
%   % global mean [O2]
%   m = m+1;
%   data(m).dataname = 'ocn_O2';
%   data(m).datacol  = 3;
%   data(m).scale    = 1.0E+6;
%   data(m).dataunit = '[O_2] (\mumol kg^{-1})';
%   data(m).minmax   = [120 180];
%   % AMOC
%   m = m+1;
%   data(m).dataname = 'AMOC strength';
%   data(m).datacol  = 0;
%   data(m).scale    = 1.0;
%   data(m).dataunit = 'AMOC (Sv)';
%   data(m).minmax   = [0 20];
%   % model skill score of global ocean salinty
%   m = m+1;
%   data(m).dataname = 'ocn_sal';
%   data(m).datacol  = 0;
%   data(m).scale    = 1.0;
%   data(m).dataunit = 'MSS (n/a)';
%   data(m).minmax   = [0.45 0.55];
%   % model skill score of global ocean PO4
%   m = m+1;
%   data(m).dataname = 'ocn_PO4';
%   data(m).datacol  = 0;
%   data(m).scale    = 1.0;
%   data(m).dataunit = 'MSS (n/a)';
%   data(m).minmax   = [0.6 0.7];
%   % model skill score of global ocean O2
%   m = m+1;
%   data(m).dataname = 'ocn_O2';
%   data(m).datacol  = 0;
%   data(m).scale    = 1.0;
%   data(m).dataunit = 'MSS (n/a)';
%   data(m).minmax   = [0.4 0.6];
%   /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
% 
%   time-slices:
%
% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
%   % extract AMOC
%   n = n + 1;
%   str = plot_fields_biogem_3d_i(loc_str_exp,'','ocn_temp','',loc_year,-1,0,['mask_worlg4_Atlantic.dat'],1.0,0.0,30.0,30,'','plot_fields_settings_OPSI_ATL',[loc_str_exp '.MOC.ATL']);
%   loc_data = str.moc_max;
%   data(n).array(y,x) = data(n).scale*loc_data;
%   % extract and compare ocean temperature to observations
%   n = n + 1;
%   str = plot_fields_biogem_3d_k(loc_str_exp,'worlg4.WOA13_Temperature_degC_worlg4.151207.nc','ocn_temp','T',loc_year,1,0,'',1.0,-4.0,4.0,40,'','plot_fields_SETTINGS_ANOM',[loc_str_exp '.temp.AMON']);
%   loc_data = str.statm_m;
%   data(n).array(y,x) = data(n).scale*loc_data;
%   % extract and compare ocean salinity to observations
%   n = n + 1;
%   str = plot_fields_biogem_3d_k(loc_str_exp,'worlg4.WOA13_Salinity_worlg4.151207.nc','ocn_sal','Sa',loc_year,1,0,'',1.0,-2.0,2.0,40,'','plot_fields_SETTINGS_ANOM',[loc_str_exp '.sal.AMON']);
%   loc_data = str.statm_m;
%   data(n).array(y,x) = data(n).scale*loc_data;
%   % extract and compare ocean [PO4] to observations
%   n = n + 1;
%   n_BESTS(1) = n;
%   str = plot_fields_biogem_3d_k(loc_str_exp,'worlg4.WOA13_PO4_molkg-1_worlg4.151208.nc','ocn_PO4','PO4',loc_year,1,0,'',1.0E-6,-0.5,0.5,40,'','plot_fields_SETTINGS_ANOM',[loc_str_exp '.PO4.AMON']);
%   loc_data = str.statm_m;
%   data(n).array(y,x) = data(n).scale*loc_data;
%   % extract and compare ocean [O2] to observations
%   n = n + 1;
%   n_BESTS(2) = n;
%   str = plot_fields_biogem_3d_k(loc_str_exp,'worlg4.WOA13_O2_molkg-1_worlg4.151208.nc','ocn_O2','O2',loc_year,1,0,'',1.0E-6,-20.0,20.0,40,'','plot_fields_SETTINGS_ANOM',[loc_str_exp '.O2.AMON']);
%   loc_data = str.statm_m;
%   data(n).array(y,x) = data(n).scale*loc_data;
% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
%
%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   20/11/25: updated 'help' :o)
%             added ability to pick up missing ensemble members
%             + pick out specific time rather than just end of time-series
%   21/03/25: modified warning error behaviour
%             examples and commenting, tidy-up
%
%   ***********************************************************************

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% *** initialize ******************************************************** %
%
% set ensemble file-string
str_ensemble = DUM_ENSEMBLE; % ensemble name
loc_year     = DUM_YEAR;     % time-slice
str_name     = DUM_NAME;
% set results directory parameters
str_dir     = 'cgenie_output';
str_archive = '.tar.gz';
str_ts_root = 'biogem_series';
str_ts_ext  = '.res';
% get directory listing of ensemble member folders
struct_dir = dir([str_dir '/' str_ensemble '*']);
% set date
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
% create empty 'best' vector
n_BESTS = [];
%
% *** initialize -- user options **************************************** %
%
% user-options for ensemble
str_sep  = '.';   % seperator string (if any) between str_ensemble and ##
str_data = '';    % data file name (if any)
%
% --- STEP #1 ----------------------------------------------------------- %
% define ensemble x axis
% NOTE: the length of the tick label structure defines the size of the 
%       first dimension of the ensemble, e.g.
%       struct_plot.xticks = {'DEFAULT'; '+seaice'; '+MLDexport'; '+MLDexport+seaice'; '+diagMLD'; '+diagMLD+seaice'};
% NOTE: for a 1-D (vertical) ensemble, set a single dummy label, e.g.
%       struct_plot.xticks = {'1D'};
% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
% define x axis label
struct_plot.xlabel = 'XLABEL';
% define x axis tick label angle
struct_plot.xtickangle = 45.0;
%
% --- STEP #2 ----------------------------------------------------------- %
% define ensemble y axis
% NOTE: the length of the tick label structure defines the size of the 
%       second dimension of the ensemble, e.g.
%       struct_plot.yticks = {'1.00'; '0.95'; '0.90'; '0.85'; '0.80'; '0.75'; '0.70'};
% NOTE: for a 1-D (vertical) ensemble, set a single dummy label, e.g.
%       struct_plot.xticks = {'1D'};
% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
% define y axis label
struct_plot.ylabel = 'YLABEL';
%
% initialize count
m = 0;
%
% --- STEP #3 ----------------------------------------------------------- %
% define time-series variables to extract and plot
% -- see help for examples
% NOTE: remember that m must be incremented by 1 for each added variable
% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
% 
% set number of time-series data
n_par_ts = m;
%
% --- STEP #4a ---------------------------------------------------------- %
% define netCDF variables to extract and plot/analyse
% -- see help for examples
% NOTE: remember that m must be incremented by 1 for each added variable
% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
%
% *** initialize -- process user options ******************************** %
%
% set total number of data
n_par = m;
% double-check (defined) number of requested variables
if (n_par ~= length(data))
    % ERROR
    disp([' ** ERROR: Inconsistency in the number of defined variables to process']);
    disp([' ']);
    return;
end
% define number of x- (column) and y- (row) axis parameters in ensemble
xmax = length(struct_plot.xticks);
ymax = length(struct_plot.yticks);
% set year to extract
loc_years = num2cell(repmat(loc_year,1,n_par));
[data.year] = loc_years{:};
%
% *** set up results arrays ********************************************* %
%
% create an array of zeros associated with each parameter
for n=1:n_par
    data(n).array = zeros(ymax,xmax);
end
%
% *********************************************************************** %
%
% *********************************************************************** %
% *** EXTRACT RESULTS AND PROCESS ENSEMBLE ****************************** %
% *********************************************************************** %
%
% loop through each (x,y) ensemble pair
% NOTE: although x and y start in the loop at a value of 1, 
%       the ensemble member numbers count up from zero
% >>>>
for x=1:xmax
    % >>>>
    for y=1:ymax
        %
        % *** prepare results directory ********************************* %
        %
        % re-create experiment name
        loc_str_exp = [str_ensemble '.' num2str(x-1) num2str(y-1)];
        disp([' >> exp == ' loc_str_exp]);
        % test for occurrence of tar.gz extension
        if exist([str_dir '/' loc_str_exp],'dir')
            loc_flag_unpack      = false;
            loc_flag_exptmissing = false;
        elseif exist([str_dir '/' [loc_str_exp str_archive]],'file')
            disp(['    UN-PACKING ...']);
            untar([str_dir '/' [loc_str_exp str_archive]],str_dir);
            loc_flag_unpack      = true;
            loc_flag_exptmissing = false;
        else
            % ERROR (report as 'warning' and keep going)
            disp([' ** WARNING: Cannot find either results directory or archive file of experiment: ' loc_str_exp]);
            disp([' ']);
            loc_flag_unpack      = false;
            loc_flag_exptmissing = true;
        end
        %
        % *** idenfity any cycles *************************************** %
        %
        if ~loc_flag_exptmissing
            % identify cycles
            [output] = fun_make_analysis_cycle(loc_str_exp);
            % set cycles present flag
            if isnan(output.n_max)
                loc_flag_cycles = false;
            else
                loc_flag_cycles = true;
                disp([' ** WARNING: Cycle identified @ ensemble member: ' num2str(x-1) num2str(y-1)]);
            end
        end
        %
        %
        % *** extract specific variables ******************************** %
        %
        % NOTE: the specified (time-slice) year is looked for
        %       and incomplete/crashed/missing runs are assigned NaN   
        % NOTE: loc_flag_exptmissing does not have to be used explicitly
        %       (becasue missing files and time-slices are tested for)
        for n=1:n_par_ts
            % set filename
            loc_str_file = [str_ts_root '_' data(n).dataname str_ts_ext];
            % test for file
            if exist([str_dir '/' loc_str_exp '/biogem/' loc_str_file],'file')
                % read data
                loc_array = load([str_dir '/' loc_str_exp '/biogem/' loc_str_file],'ascii');
                loc_data_i = find(loc_array(:,1) == data(n).year);
                if ~isempty(loc_data_i)
                    loc_data = loc_array(loc_data_i,data(n).datacol);
                    % test for data difference request
                    if (isfield(data,'datacolD'))
                        if (~isempty(data(n).datacolD))
                            loc_data = loc_data - loc_array(loc_data_i,data(n).datacolD);
                        end
                    end
                    % write data
                    data(n).array(y,x) = data(n).scale*loc_data;
                else
                    data(n).array(y,x) = NaN;
                    disp([' ** WARNING: Cannot find time-point: ' num2str(data(n).year)]);
                end
                % replace data if cycles are present
                if loc_flag_cycles
                    data(n).array(y,x) = data(n).scale*mean(loc_array(output.n_min:output.n_max,data(n).datacol)); 
                    disp([' **          Replacing point data: ' data(n).dataname ' with cycle mean.']);                   
                end
            else
                data(n).array(y,x) = NaN;
                disp([' ** WARNING: Cannot find file: ' [str_dir '/' loc_str_exp '/biogem/' loc_str_file]]);
            end
        end
        %
        % *** process 3D netCDF ***************************************** %
        %
        if (n_par > n_par_ts)
            n = n_par_ts;
            if loc_flag_exptmissing
                % missing experiments => set all NaNs for results
                for n=n_par_ts+1:n_par
                    data(n).array(y,x) = NaN;
                end
            else
                % --- STEP #4b ------------------------------------------ %
                % extract and plot/analyse netCDF variables 
                % \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
                % /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
                if (n ~= n_par)
                    % ERROR (fatal)
                    disp([' ** ERROR: Inconsistency in the number of processed variables']);
                    disp([' ']);
                    return;
                end
            end
        end
        %
        % *** clean up ************************************************** %
        %
        % remote unpacked dir
        if loc_flag_unpack
            disp(['    REMOVE DIR']);
            rmdir([str_dir '/' loc_str_exp],'s');
        end
        
    end
    % <<<<
end
% <<<<
%
% *********************************************************************** %
%
% *********************************************************************** %
% *** PLOT DATA ARRAYS ************************************************** %
% *********************************************************************** %
%
% create gridded plots of ensemble
for n=1:n_par
    % test for data difference request (and modify filename)
    if (isfield(data,'datacolD'))
        if (~isempty(data(n).datacolD))
            data(n).dataname = ['D' data(n).dataname];
        end
    end
    % construct filename and create plot
    struct_plot.filename = [str_name '.' data(n).dataname '.' num2str(data(n).datacol)];
    % PRESCRIBED SCALE
    struct_plot.unitslabel = data(n).dataunit;
    plot_2dgridded2(data(n).array,data(n).minmax,'',struct_plot);
%     % AUTOSCALE
%     plot_2dgridded2(data(n).array,[],'',struct_plot);
    % save data
    % NOTE: y-axis is opposite to as displayed in the plot 
    %      (counting rows down)
    fprint_2DM(data(n).array,[],[struct_plot.filename '.dat'],'%10.4f','%10.4f',true,false);
end
%
% *********************************************************************** %
%
% *********************************************************************** %
% *** ANALYSE BEST(S) *************************************************** %
% *********************************************************************** %
%
for n=1:length(n_BESTS)
    %
    n_BEST = n_BESTS(n);
    %
    % *** find best ensemble member ************************************* %
    %
    mss_BEST = max(max(data(n_BEST).array)); % best model skill score
    I = find(data(n_BEST).array == mss_BEST);
    [n_y,n_x] = ind2sub([ymax xmax],I);
    %
    % *** print best stats ********************************************** %
    %
    fid = fopen([str_name '.' data(n_BEST).dataname '.' num2str(n_x-1) num2str(n_y-1) '.STATS.txt'], 'wt');
    fprintf(fid, '\n');
    fprintf(fid, '=== MSS STATS SUMMARY === \n');
    fprintf(fid, '\n');
    fprintf(fid, 'Best (x,y) : %d %d \n', n_x,n_y);
    fprintf(fid, '(ensemble notation: .%d%d ) \n', n_x-1,n_y-1);
    fprintf(fid, '\n');
    fprintf(fid, '------------------------- \n');
    for n=1:n_par
        fprintf(fid, [data(n).dataname ' = %8.4f \n'], data(n).array(n_y,n_x));
    end
    fprintf(fid, '------------------------- \n');
    fprintf(fid, 'BEST: \n');
    fprintf(fid, [data(n_BEST).dataname ' = %8.4f \n'], data(n_BEST).array(n_y,n_x));
    fprintf(fid, '------------------------- \n');
    fprintf(fid, '\n');
    fprintf(fid, '========================= \n');
    fprintf(fid, '\n');
    fclose(fid);
    %
    % *** prepare results directory ************************************* %
    %
    % NOTE: we already know that the 'best' experiment must exist!
    % re-create experiment name
    loc_str_exp = [str_ensemble '.' num2str(n_x-1) num2str(n_y-1)];
    % test for occurrence of tar.gz extension
    if exist([str_dir '/' loc_str_exp],'dir')
        loc_flag_unpack = false;
    elseif exist([str_dir '/' [loc_str_exp str_archive]],'file')
        disp(['    UN-PACKING ...']);
        untar([str_dir '/' [loc_str_exp str_archive]],str_dir);
        loc_flag_unpack = true;
    end
    %
    % *** analyse ******************************************************* %
    %
    loc_str_name = [str_name '.' num2str(n_x-1) num2str(n_y-1)];
%     % PO4, O2
%     plot_fields_biogem_3d_k(loc_str_exp,'worjh2.p_an.200709.nc','ocn_PO4','p_an',loc_year,1,16,'',1.0E-6,-0.5,0.5,40,'','plot_fields_SETTINGS_ANOM',[loc_str_name '.k.PO4.ANOM.SUR']);
%     plot_fields_biogem_3d_i(loc_str_exp,'worjh2.p_an.200709.nc','ocn_PO4','p_an',loc_year,1,0,'mask_worjh2_AtlanticALL.dat',1.0E-6,-0.5,0.5,40,'','plot_fields_SETTINGS_ANOM',[loc_str_name '.i.PO4.ANOM.ATLALL']);
%     plot_fields_biogem_3d_i(loc_str_exp,'worjh2.p_an.200709.nc','ocn_PO4','p_an',loc_year,1,0,'mask_worjh2_PacificALL.dat',1.0E-6,-0.5,0.5,40,'','plot_fields_SETTINGS_ANOM',[loc_str_name '.i.PO4.ANOM.PACALL']);
%     plot_fields_biogem_3d_i(loc_str_exp,'worjh2.p_an.200709.nc','ocn_PO4','p_an',loc_year,1,0,'mask_worjh2_IndianALL.dat',1.0E-6,-0.5,0.5,40,'','plot_fields_SETTINGS_ANOM',[loc_str_name '.i.PO4.ANOM.INDALL']);
%     plot_fields_biogem_3d_k(loc_str_exp,'worjh2.o_an.200709.nc','ocn_O2','o_an',loc_year,1,16,'',1.0E-6,-20.0,20.0,40,'','plot_fields_SETTINGS_ANOM',[loc_loc_str_name '.k.O2.ANOM.SUR']);
%     plot_fields_biogem_3d_i(loc_str_exp,'worjh2.o_an.200709.nc','ocn_O2','o_an',loc_year,1,0,'mask_worjh2_AtlanticALL.dat',1.0E-6,-50.0,50.0,40,'','plot_fields_SETTINGS_ANOM',[loc_str_name '.i.O2.ANOM.ATLALL']);
%     plot_fields_biogem_3d_i(loc_str_exp,'worjh2.o_an.200709.nc','ocn_O2','o_an',loc_year,1,0,'mask_worjh2_PacificALL.dat',1.0E-6,-50.0,50.0,40,'','plot_fields_SETTINGS_ANOM',[loc_str_name '.i.O2.ANOM.PACALL']);
%     plot_fields_biogem_3d_i(loc_str_exp,'worjh2.o_an.200709.nc','ocn_O2','o_an',loc_year,1,0,'mask_worjh2_IndianALL.dat',1.0E-6,-50.0,50.0,40,'','plot_fields_SETTINGS_ANOM',[loc_str_name '.i.O2.ANOM.INDALL']);
    %
    % *** clean up ****************************************************** %
    %
    % (optional) remove unpacked dir
    if loc_flag_unpack
        disp(['    KEEP UNPACKED BEST RUN!']);
        %     disp(['    REMOVE DIR']);
        %     rmdir([str_dir '/' loc_str_exp],'s');
    end
    %
    % ******************************************************************* %
    %
end
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
%
close all;
%
% *********************************************************************** %
