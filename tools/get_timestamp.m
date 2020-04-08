function time = get_timestamp(JD)
% This algorithm is based on Vallado 4e algorithm 22, pp 202.
%
% Inputs
%   - JD : a scalar absolute Julian Date
%
% Outputs
%   - time : a (6,1) vector of [Y,M,D,h,m,s]

% Constants
JD_1900     = 2415019.5;
Lmonth      = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

% Convert JD to date-time
T1900   = (JD - JD_1900)/365.25; % base epoch is 1900
year    = 1900 + floor(T1900);

leapyrs = floor((year - 1900 - 1)*0.25);
days    = (JD - JD_1900) - ((year - 1900)*365 + leapyrs );

if days < 1
    year    = year - 1;
    leapyrs = floor((year - 1900 - 1)*0.25);
    days    = (JD - JD_1900) - ((year - 1900)*365 + leapyrs );
end

if mod(year,4) == 0
    Lmonth(2) = 29;
end

dayofyear   = floor(days);

day     = 0;
month   = 0;

while day < dayofyear
    month   = month + 1;
    day = day + Lmonth(month);
end

%dayofmonth = dayofyear - (day - Lmonth(month));

tau     = 24*( days-dayofyear );
hour    = floor( tau );
min     = floor( 60*(tau - hour) );
sec     = 3600*( tau - hour - (min/60) );

time    = [ year, month, dayofyear, hour, min, sec ];
end
