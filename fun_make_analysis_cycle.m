function [DUM_OUTPUT] = fun_make_analysis_cycle(DUM_EXP)
%
%   ***********************************************************************
%   *** fun_make_analysis_cycle *******************************************
%   ***********************************************************************
%
%
%
%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   21/03/26: committed
%
%   ***********************************************************************


% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% *** USER OPTIONS ****************************************************** %
%
% set reference data
% --> data from which cycles will be identified from
ref_data_str = 'ocn_temp';
% reference data column
ref_col_n    = 4;
% set reference data tolerance
% --> minimum variability to indicate the existance of a cycle
ref_data_tol = 0.2;
% set reference data tolerance interval
% --> time interval at end to test for a critical amplitude of variability
tol_time_min = 8000.0;
% set spinup time
% --> time interval at start of record to exclude
ref_time_min = 4000.0;
% set additional (optional) time-series to plot
opt1_data_str = 'misc_opsi';
opt2_data_str = 'misc_opsi';
opt1_col_n    = 2;
opt2_col_n    = 3;
%
% *** initialize ******************************************************** %
%
% set experiment file-string
str_experiment = DUM_EXP;
% set resutls directory parameters
str_dir     = 'cgenie_output';
str_archive = '.tar.gz';
str_ts_root = 'biogem_series';
str_ts_ext  = '.res';
% set date
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
%
% *** load and process reference time-series **************************** %
%
ref_data = load([str_dir '/' str_experiment '/biogem/biogem_series_' ref_data_str str_ts_ext],'-ascii');
%
time_min = 0.0;
time_max = (ref_data(end,1) + 0.5)/1000.0;
n_max    = length(ref_data);
% determine tolorance bound
n_min_tol = min(find(ref_data(:,1) >= tol_time_min));
% determine spinup bound
n_min = min(find(ref_data(:,1) >= ref_time_min));
%
% *** plot time-series ************************************************** %
%
% plot full time-series
loc_str = [str_experiment '.' ref_data_str '_FULL'];
plot_timeseries_biogem(str_experiment,'',time_min,time_max,ref_data_str,ref_col_n,opt1_data_str,opt1_col_n,opt2_data_str,opt2_col_n,'plot_timeseries_SETTINGS_cycles',loc_str);
%
% *** load secondary data *********************************************** %
%
if ~isempty(opt1_data_str)
    opt1_data = load([str_dir '/' str_experiment '/biogem/biogem_series_' opt1_data_str str_ts_ext],'-ascii');
end
if ~isempty(opt2_data_str)
    opt2_data = load([str_dir '/' str_experiment '/biogem/biogem_series_' opt2_data_str str_ts_ext],'-ascii');
end
%
% *********************************************************************** %
%
% *********************************************************************** %
% *** IDENTIFY CYCLES *************************************************** %
% *********************************************************************** %
%
% *** test for a cycle potentially existing ***************************** %
%
% test for minimum variability (cycles) existing
% NOTE: within only the tolorance testing internal
ref_min = min(ref_data(n_min_tol:n_max,ref_col_n));
ref_max = max(ref_data(n_min_tol:n_max,ref_col_n));
if ((ref_max-ref_min) > ref_data_tol)
    flag_cycle = true;
else
    flag_cycle = false;
end
%
if flag_cycle
    %
    % *** identify cycles *********************************************** %
    %
    % set reference (initial) value
    ref_init = ref_data(end,ref_col_n);
    % determine initial sign of change
    if (ref_data(end-1,ref_col_n) > ref_init)
        flag_up = true;
    else
        flag_up = false;
    end
    % initiailze loc min,max
    loc_data_old = ref_init;
    % initialize storage vectors for mins and maxs
    data_min_n = [];
    data_max_n = [];
    % work back through record -- find first and then second min or max
    % (assume wavelength based on that)
    % NOTE: assume record not starting on a min,max (???)
    for n = n_max-1:-1:n_min
        loc_data = ref_data(n,ref_col_n);
        if (loc_data < loc_data_old)
            if (flag_up)
                % record maximum
                data_max_n = [data_max_n n+1];
                % reset
                flag_up = false;
            end
            % update min,max
            loc_data_old = loc_data;
        elseif (loc_data > loc_data_old)
            if (~flag_up)
                % record minimum
                data_min_n = [data_min_n n+1];
                % reset
                flag_up = true;
            end
            % update min,max
            loc_data_old = loc_data;
        else
            %%% nothing
        end
    end
    % determine whether a min or max occurs first
    % and set results array accordingly
    if (data_max_n(1) > data_min_n(1))
        data_peak_n = data_max_n(1:2);
    else
        data_peak_n = data_min_n(1:2);
    end
    % determine number of sample points inbetween the first pair of mins/maxs
    loc_dn    = (data_peak_n(1)-data_peak_n(2));
    % determine the the position, one wavelength prior to the record end
    loc_n_min = n_max - loc_dn;
    %
    % *** plot time-series of single complete cycle ********************* %
    %
    % last full cycle
    time_min = ref_data(loc_n_min,1)/1000.0;
    time_max = ref_data(n_max,1)/1000.0;
    % plot last complete cycle
    loc_str = [str_experiment '.' ref_data_str '_END'];
    plot_timeseries_biogem(str_experiment,'',time_min,time_max,ref_data_str,ref_col_n,opt1_data_str,opt1_col_n,opt2_data_str,opt2_col_n,'plot_timeseries_SETTINGS_cycles',loc_str);
    %
    % *** calculate and return stats ************************************ %
    %
    % NOTE: period in years
    % reference data
    loc_amplitude = ref_data(data_max_n(1),ref_col_n) - ref_data(data_min_n(1),ref_col_n);
    loc_period    = ref_data(data_peak_n(1),1) - ref_data(data_peak_n(2),1);
    loc_average   = mean(ref_data(data_peak_n(2):data_peak_n(1),ref_col_n));
    output.amplitude = loc_amplitude;
    output.period    = loc_period;
    output.amplitude = loc_average;
    % cycle limits
    output.n_min     = loc_n_min;
    output.n_max     = n_max;
    output.t_min     = ref_data(loc_n_min,1);
    output.t_max     = ref_data(n_max,1);
    %
else
    %
    % *** return null values ******************************************** %
    %
    output.amplitude = NaN;
    output.period    = NaN;
    output.amplitude = NaN;
    output.n_min     = NaN;
    output.n_max     = NaN;
    output.t_min     = NaN;
    output.t_max     = NaN;
    %
end
%
% *********************************************************************** %
%
DUM_OUTPUT = output;
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
