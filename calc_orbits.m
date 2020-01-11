% function calc_orbits(solutionname, t1,t2)
% Calculate orbital parameters for cGENIE
%
% INPUT
%   solutionname: filename of an astronomical solution
%       e.g., 'La2004.dat'
%   t1: selected time interval start, must be a non-negative integer
%       e.g., 55000. Default: 0
%   t2: selected time interval end, must be a non-negative integer
%       e.g., 56000. Default: 100000
%
% OUTPUT
%   A *.dat file used for cGENIE transient model.
%   should be saved in the folder:
%       cgenie.muffin/genie-embm/data/input/
%
%   This file 5 column without header
%       col #1: age
%       col #2: eccentricity
%       col #3: Sine of obliquity
%       col #4: mean longitude of perihelion
%       col #5: angular time of year
%
% Calls for
%   La2004.dat  % solution data file
%       must be a 4 column data file without header
%       col #1: age (unit is ka, non-negetive)
%       col #2: eccentricity
%       col #3: longitude of perihelion, in degree
%       col #4: obliquity, in degree
%
% Example #1
%
% solutionname = 'La2004.dat';
% t1 = 55000;
% t2 = 56000;
% calc_orbits(solutionname,t1,t2)
%
% Example #2
% solutionname = 'La2004.dat';
% calc_orbits(solutionname)
%
% Example #3
% calc_orbits()
%
% By Sandra Kirtland Turner (UCR)
%
% Revised by Mingsong Li (Penn State), Jan, 2020
%
% Reference:
%   Berger, A., 1978. Long-term variations of daily insolation and 
%   Quaternary climatic changes. Journal of the atmospheric sciences 
%   35, 2362-2367.
%
function [] = calc_orbits(solutionname,t1,t2)

if nargin < 1; solutionname = 'La2004.dat'; end
% read data
[~,filename,~] = fileparts(solutionname);
orbit_vals=load(solutionname);

if nargin < 3; t2 = orbit_vals(end,1);end
if nargin < 2; t1 = orbit_vals(1,1);end

if t1 < 0; error('Error, t1 must be >= 0'); end
if t2 < 1; error('Error, t2 must be > 0');  end

t1 = min(round(t1), round(t2));
t2 = max(round(t1), round(t2));
% First load the appropriate file formatted as yr(kyr),ecc,long. peri(deg),

%% Calculate orbital parameters
orbit_vals = orbit_vals(t1+1:t2+1,:);
nyrs=length(orbit_vals);
yr=orbit_vals(:,1);
ecc=orbit_vals(:,2);
zw=orbit_vals(:,3);
obl=orbit_vals(:,4);

for i=1:nyrs
    zw(i)=180+zw(i);
    if zw(i)>360
        zw(i)=zw(i)-360;
    end
    zv(i)=pi/180*(360-zw(i));
    zm(i)=zv(i);
    for ii=1:10
        zm(i)=zv(i)-(2*ecc(i)-0.25*ecc(i)*ecc(i)*ecc(i))*sin(zm(i))-1.25*ecc(i)*ecc(i)*sin(2*zm(i))-(13/12)*ecc(i)*ecc(i)*ecc(i)*sin(3*zm(i));
    end
    dn(i)=zm(i)*365.24/(2*pi);
    zprff(i)=dn(i)-79.265;
    zgamma(i)=-zprff(i);
    zw2(i)=zw(i)-zw(1);
    if zgamma(i)<0
        zgamma(i)=365.24+zgamma(i);
    end
end

% Gamma = mean longitude of perihelion
gamma=0.5*(zv+zm)';

% Sin obl
sin_obl=sin(pi/180*obl);

% Tau = zgamma
tau=zgamma';

% Output data
orbits_out=horzcat(yr,ecc,sin_obl,gamma,tau);

% Write output to file in format for cGENIE
% yr(kyr), ecc, sin(obl), mean longitude, tau(gamma)
filename1 = [filename,'_',num2str(t1),'_',num2str(t2)];
dlmwrite([filename1 '_cGENIE.dat'],orbits_out,'delimiter',' ','precision',9)
end