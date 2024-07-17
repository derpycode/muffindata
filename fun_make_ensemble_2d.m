function [] = fun_make_ensemble_2d(STR_TEMPLATE,STR_PARAMS)
%
%   ***********************************************************************
%   *** fun_make_ensemble_2d **********************************************
%   ***********************************************************************
%
%   fun_make_ensemble_2d
%   creates all the individual user-config files needed for a
%   2D parameter ensemble
%   also creates an ensemble description file: '*.KEY.txt'
%   The resulting ensemble can be submitted by modifying the BASH script:
%   sub_ens.muffin.sh
%
%   STR_TEMPLATE == full filename of the 'template' (default) user-config
%                   NOTE: this filename can be anything you want and
%                         will not form part of the ensemble name
%   STR_PARAMS   == full filename of the parameter configuration file
%                   NOTE: this filename can be anything you want, but:
%                         (i)  it will form the ensemble name
%                         (ii) if it ends in '.dat' or '.txt', the
%                              extension will be stripped off 
%
%  fun_make_ensemble_2d must be run form the same directory that contains:
%  (a) the 'template' (default) user-config file
%  (b) the file containing the parameter configuration
%
%   The format of the file containing the parameter configuration is
%   (by column number):
%   (1) parameter change number (ID) [integer]
%   (2) parameter name [string]
%   (3) default parameter value
%   (4-n) any number of parameter modifiers [real]
%   (n+1) line termination:
%         -99 == end of data line
%        -999 == end of LAST data line
%   NOTE: where more than 1 parameter shares the same ID 
%         (and hence are modified together in the same ensemble member)
%         adding the (duplicate) vector is not necessary
% 
%   EXAMPLE 1
%
%   The parameter configuration file contents:
%
%   1 go_0 1.0E-3 0.0 10.0 100.0 -99
%   2 ea_0 1.0E6  1.0 0.8 0.6 0.4 0.2 0.0 -999
%
%   would create a 3 x 6 ensemble,
%   modifying param setting go_0=1.0E-3 by vector [0.0 10.0 100.0]
%   and parameter setting ea_0=1.0E6 by vector [1.0 0.8 0.6 0.4 0.2 0.0]
%
%   EXAMPLE 2
%
%   The parameter configuration file contents:
%
%   1 go_0 1.0E-3 0.0 10.0 100.0 -99
%   1 go_9 1.0E-9 1.0 2.0 3.0 -99
%   2 ea_0 1.0E6  1.0 0.8 0.6 0.4 0.2 0.0 -999
%
%   would create a 3 x 6 ensemble as example #1, 
%   but also modifying go_9 at the same time as go_0
%
%   EXAMPLE 3
%
%   The parameter configuration file contents:
%
%   1 go_0 1.0E-3 0.0 10.0 100.0 -999
%
%   creates a 6 x 1 (1D) ensemble in which only the 1st parameter varies.
%
%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   20/02/29: set count of ensemble # to start from 0 (rather than 1)
%   20/11/25: updated to allow independent additional 1st or 2nd
%             axis parameters
%   24/07/17: improved commenting
%
%   ***********************************************************************

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% *** initialize ******************************************************** %
%
% set date
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
% set template user-config filename
str_template = STR_TEMPLATE;
% set parameter definition filename string
str_file = [STR_PARAMS];
% define a structure for later
s = struct;
%
% *** load ensemble parameter file ************************************** %
%
% open parameter configuration file ...
fid = fopen(str_file);
% initialize the loop
nl = 0;
go_lines = true;
loc_id = 0;
% scan through the file
while go_lines
    % read fixed format data portion (excluding vector)
    C = textscan(fid, '%d %s %f', 1, 'CommentStyle', '%');
    % initialize vector read
    go_columns = true;
    v = [];
    % read vector elements ...
    % ... but only if the parameter ID is unique (not the same a previous)
    if (cell2mat(C(1)) > loc_id)
        while go_columns
            loc_C = textscan(fid, '%f', 1, 'CommentStyle', '%');
            if (cell2mat(loc_C) == -99)
                go_columns = false;
                nl = nl + 1;
            elseif (cell2mat(loc_C) == -999)
                go_columns = false;
                go_lines = false;
                nl = nl + 1;
            else
                v = [v cell2mat(loc_C)];
            end
        end
        loc_unique = true;
    else
        while go_columns
            loc_C = textscan(fid, '%f', 1, 'CommentStyle', '%');
            if (cell2mat(loc_C) == -99)
                go_columns = false;
                nl = nl + 1;
            else
                v = [v cell2mat(loc_C)];
            end
        end
        loc_unique = false;
    end
    % update structure data
    s(nl).id = cell2mat(C(1));
    s(nl).paramname = char(C{2}(:));
    s(nl).defaultvalue = cell2mat(C(3));
    s(nl).vector = v;
    s(nl).unique = loc_unique;
    % update current ID
    loc_id = s(nl).id;
end
% ... close
fclose(fid);
% determine total number of parameters
par_pmax = length([s(:).id]);
% determine number of parameter elements to modify
% (not necessarily the same as the total number of parameters)
par_nmax = max([s(:).id]);
% determine number of parameter modifications to make (for each parameter)
for n=1:par_nmax
    loc_v = find([s(:).id] == n);
    loc_n = loc_v(1);
    par_vmax(n) = length(s(loc_n).vector);
end
% check for a 'duplicate' parameter / 'extended 1D' enesmble
flag_ext1d = false;
if (par_nmax == 2)
    if (strcmp(s(1).paramname,s(2).paramname))
        flag_ext1d = true;
    end
end
%
% *** create ensemble parameter summary file **************************** %
%
% try and shorten parameter filename string
str_name = str_file;
if (strcmp(str_file(end-3:end),'.dat')), str_name = str_file(1:end-4); end
if (strcmp(str_file(end-3:end),'.txt')), str_name = str_file(1:end-4); end
%
str_keyfile = [str_date '.' str_name '.KEY.txt'];
fid0 = fopen(str_keyfile, 'w');
loc_str = [' ***********************************************************'];
fprintf(fid0, '%s\n', loc_str);
loc_str = ['     ensemble name                    : ' str_date '.' str_name];
fprintf(fid0, '%s\n', loc_str);
loc_str = ['     template user-config filename    : ' str_template];
fprintf(fid0, '%s\n', loc_str);
loc_str = ['     parameter specification filename : ' str_file];
fprintf(fid0, '%s\n', loc_str);
for n=1:par_pmax
    if (s(n).unique)
        loc_str = [' -----------------------------------------------------------'];
        fprintf(fid0, '%s\n', loc_str);
        loc_str = ['     ensemble axis #' num2str(s(n).id)];
        fprintf(fid0, '%s\n', loc_str);
    end
    loc_str = ['     parameter                        : ' s(n).paramname];
    fprintf(fid0, '%s\n', loc_str);    
    loc_str = ['     default value                    : ' num2str(s(n).defaultvalue)];
    fprintf(fid0, '%s\n', loc_str);    
    loc_str = ['     parameter modifiers              : ' num2str(s(n).vector)];
    fprintf(fid0, '%s\n', loc_str);    
end
loc_str = [' ***********************************************************'];
fprintf(fid0, '%s\n', loc_str);
loc_str = ['     individual ensemble members:'];
fprintf(fid0, '%s\n', loc_str);
%
% *********************************************************************** %
%
% *********************************************************************** %
% *** CREATE ENSEMBLE CONFIG FILES ************************************** %
% *********************************************************************** %
%
% *** DOUBLE PARAMETER CHANGES ****************************************** %
%
% NOTE: count ensemble # from 0 (rather than 1)
% loop through all unique parameter modifications
for n=1:par_vmax(1)
    for m=1:par_vmax(2)
        % create a vector of where we are in the nested loop
        loc_vn = [n m];
        % copy template
        str_templatefilein  = [str_template];
        str_templatefileout = [str_date '.' str_file '.' num2str(n-1) num2str(m-1)];
        copyfile(str_templatefilein,str_templatefileout,'f');
        % open sesame! (file pipe of ensemble member user-config)
        fid = fopen(str_templatefileout, 'a+');
        % write parameter file header info
        loc_str = '# ';
        fprintf(fid, '%s\n', loc_str);
        loc_str = '# --- generated by MATLAB with love :) --------------------------------';
        fprintf(fid, '%s\n', loc_str);
        loc_str = '# ';
        fprintf(fid, '%s\n', loc_str);
        % write parameters -- loop through all parameters
        if (flag_ext1d)
            % special case of 'extended 1D' ensemble
            loc_str = ['# parameter: ' s(1).paramname ' -- default value (' num2str(s(1).defaultvalue) ') modified by factor: ' num2str(s(1).vector(loc_vn(s(1).id))+s(2).vector(loc_vn(s(2).id)))];
            fprintf(fid, '%s\n', loc_str);
            loc_str = [s(1).paramname '=' num2str((s(1).vector(loc_vn(s(1).id))+s(2).vector(loc_vn(s(2).id)))*s(1).defaultvalue)];
            fprintf(fid, '%s\n', loc_str);
        else
            for p=1:par_pmax
                % add comment line
                loc_str = ['# parameter: ' s(p).paramname ' -- default value (' num2str(s(p).defaultvalue) ') modified by factor: ' num2str(s(p).vector(loc_vn(s(p).id)))];
                fprintf(fid, '%s\n', loc_str);
                % add parameter line
                loc_str = [s(p).paramname '=' num2str(s(p).vector(loc_vn(s(p).id))*s(p).defaultvalue)];
                fprintf(fid, '%s\n', loc_str);
            end
        end
        % add parameter file end marker
        loc_str = '# ';
        fprintf(fid, '%s\n', loc_str);
        % close file pipe of ensemble member user-config
        fclose(fid);
        % add ensemble key file info
        loc_str = ['member #' num2str(n-1) num2str(m-1) ':'];
        fprintf(fid0, '%s\n', loc_str);
        for p=1:par_pmax
            % add comment line to ensemble key file
            loc_str = ['           ' s(p).paramname '=' num2str(s(p).vector(loc_vn(s(p).id))*s(p).defaultvalue)];
            fprintf(fid0, '%s\n', loc_str);
        end
    end
end
%
% *********************************************************************** %
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
% add terminal comment line to ensemble key file
loc_str = [' ***********************************************************'];
fprintf(fid0, '%s\n', loc_str);
% close file pipe of ensemble key file
fprintf(fid0, '\n', loc_str);
fclose(fid0);
%
% *********************************************************************** %