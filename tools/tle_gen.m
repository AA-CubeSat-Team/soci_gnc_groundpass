function [orbit_tle] = tle_gen(YMDHMS, INC, RAAN, ECC, AOP, MNA, SMA, varargin)
%TLE_GEN
%
% Script that generates a formatted TLE.
%
% Inputs:
%   - YMDHMS : the Year-Month-Day-Hour-Min-Sec at epoch as a [6x1] vector
%   - INC : orbit inclination in degrees
%   - RAAN : right ascension of ascending node in degrees
%   - ECC : orbit eccentricity
%   - AOP : argument of perigee in degrees
%   - MNA : mean anomaly in degrees
%   - SMA : semimajor axis in km
%   - varargin :
%       - MNM_dot : derivative of mean motion in revs/day/day
%       - B_star : base of the B_star drag term
%       - B_star_ex : exponent of the B_star drag term
%
% T. Reynolds

% Constants
twopi = 2*pi;
KM2M  = 1e3;            % convert km to m
rE    = 6378.137e3;     % radius of Earth in m
muE   = 3.986004418e14; % Earth gravitational parameter in m3/s2
J2000 = 2451545;        % Julian date of J2000 epoch

% Run input checks
if( nargin < 7 )
    error('TLE_gen: not enough inputs')
end

if( nargin > 7 )
    MNM_dot = varargin{1}; %#ok
    % pad this with zeros like ECC when needed
    if( nargin > 8 )
        if( nargin == 8 )
            error('TLE_gen: Must give both B_star and B_star_ex')
        else
            B_star = varargin{2};
            B_star_ex = varargin{3};
        end
    end
else
    % Assumed values from SWISSCUBE -- WAG for now
    B_star      = 0.32923;
    B_star_ex   = -4;
    sMNM_dot     = strcat('+.','00000',num2str(145));
end

if( numel(YMDHMS) < 6 )
    error('TLE_gen: Not enough entries in YMDHMS; must be 6')
end

% map SMA to mean motion for use in TLE
a   = KM2M * SMA;
T   = twopi/sqrt(muE) * (a^(3/2));
MNM = 86400/T; % gives the answer in revs/day

% extract epoch datetime terms
year     = YMDHMS(1);
fullyear = year;
if( year > 2000 )
    year = year - 2000;
else
    year = year - 1900;
end
month   = YMDHMS(2);
day     = YMDHMS(3);
hour    = YMDHMS(4);
min     = YMDHMS(5);
sec     = YMDHMS(6);

% Find the decimal day of year
Lmonth  = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
if( mod(year,4) == 0 )
    Lmonth(2) = 29;
end
day_dec = 0.0;
for m = 1:month-1
    day_dec = day_dec + Lmonth(m);
end
day_dec = day_dec + day + hour/24 + min/1440 + sec/86400;

% Julian date at Jan 1, epoch_year, 00:00:00
JD_begin_of_year    = 367*fullyear - floor((7/4)*(fullyear + floor(10/12))) + ...
                        floor(275/9) + 1721013.5;
JD_UTC_epoch_J2000  = day_dec + JD_begin_of_year - J2000;

% Reformat the eccentricity
nz  = abs(round(log10(ECC)));
ECC = ECC * 10^7;
switch nz
    case 0
        sECC = num2str(ECC,'%7.0f');
    case 1
        sECC = strcat('0',num2str(ECC,'%6.0f'));
    case 2
        sECC = strcat('00',num2str(ECC,'%5.0f'));
    case 3 
        sECC = strcat('000',num2str(ECC,'%4.0f'));
    case 4 
        sECC = strcat('0000',num2str(ECC,'%3.0f'));
    case 5
        sECC = strcat('00000',num2str(ECC,'%2.0f'));
    case 6
        sECC = strcat('000000',num2str(ECC,'%1.0f'));
    case 7
        sECC = '0000000';
end

% Output a text file with this TLE
fID = fopen('TLEs/test_TLE.txt','w');
fprintf(fID,'1 00000U 20001Z   %02.0f%012.8f %s  00000-0 %+5.0f%2.0f 0  0017\n',...
            year,day_dec,sMNM_dot,B_star*10^5,B_star_ex);
fprintf(fID,'2 00000 %08.4f %08.4f %s %08.4f %08.4f %10.8f563537',...
            INC,RAAN,sECC,AOP,MNA,MNM);
fclose(fID);

if (nargout>1)
    % Output a TLE in the format used for sim
    orbit_tle = [ year;
                  JD_UTC_epoch_J2000;
                  B_star*10^(B_star_ex);
                  INC;
                  RAAN;
                  ECC/1e7;
                  AOP;
                  MNA;
                  MNM ];
end

end