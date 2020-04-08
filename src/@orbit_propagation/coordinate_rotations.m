function [ecef_to_eci, ppef_to_veci, mod_to_eci, teme_to_eci] = coordinate_rotations(obj)
% COORDINATE_ROTATIONS excludes the polar motion transformation. Hence it
% includes the Precession, Nutation, and Sidereal Time rotation matrices.

% UPDATED: Excluding polar motion since it would need to be updated ~1 per
% week to keep accurate. For our mission duration and requirements,
% accounting for polar motion is not necessary. -T. Reynolds, 9.11.17

% Constants
% Everything should be precomputed in radians
JC_J2000_ut1 = obj.JC_J2000_ut1;
JC_J2000_TT  = obj.JC_J2000_TT;
deg2rad      = pi/180;
asec2rad     = 1/3600*deg2rad;
w_prec       = obj.params.w_prec;

% Evaluate Delaunay parameters
epsb_1980 = deg2rad * ( 23.439291 ...
                        - 0.0130042 * JC_J2000_TT ...
                        - 1.64e-7 * JC_J2000_TT^2 ...
                        + 5.04e-7 * JC_J2000_TT^3);

% Evaluate Delaunay parameters (correct equations from 4e errata)
M_moon	=   deg2rad * ( 134.96298139 ...
                        + (1325 * 360 + 198.8673981) * JC_J2000_TT ...
                      	+ 0.0086972 * JC_J2000_TT^2 ...
                       	+ 1.78e-5 * JC_J2000_TT^3);
               
M_sun	=   deg2rad * ( 357.52772333 ...
                        + (99 * 360 + 359.0503400) * JC_J2000_TT ...
                      	- 0.0001603 * JC_J2000_TT^2 ...
                     	- 3.3e-6 * JC_J2000_TT^3);
                
u_moon	=   deg2rad * ( 93.27191028 ...
                        + (1342 * 360 + 82.0175381) * JC_J2000_TT ...
                      	- 0.0036825 * JC_J2000_TT^2 ...
                        + 3.1e-6 * JC_J2000_TT^3);
                  
D_sun	=   deg2rad * ( 297.85036306 ...
                        + (1236 * 360 + 307.1114800) * JC_J2000_TT ...
                        - 0.0019142 * JC_J2000_TT^2 ...
                     	+ 5.3e-6 * JC_J2000_TT^3);
                   
O_moon	=   deg2rad * ( 125.04452222 ...
                        - (5 * 360 + 134.1362608) * JC_J2000_TT ...
                    	+ 0.0020708 * JC_J2000_TT^2 ...
                        + 2.2e-6 * JC_J2000_TT^3);              
                                        
% Largest Nutation coefficents (out of 106)                    
nc = [ ...
1   0   0   0   0   1   -171996 -174.2  92025   8.9
9   0   0   2   -2  2   -13187  -1.6    5736    -3.1
31  0   0   2   0   2   -2274   -0.2    977     -0.5
2   0   0   0   0   2   2062    0.2     -895    0.5
10  0   1   0  	0   0   1426    -3.4    54      -0.1
32  1   0   0   0   0   712     0.1     -7      0.0
11  0   1   2   -2  2   -517    1.2     224     -0.6
33  0   0   2   0   1   -386    -0.4    200     0.0
34  1   0   2   0   2   -301    0.0     129     -0.1
12  0   -1  2   -2  2   217     -0.5    -95     0.3
35  1   0   0   -2  0   -158    0.0     -1      0.0
13  0   0   2   -2  1   129     0.1     -70     0.0 
36  -1  0   2   0   2   123     0.0     -53     0.0 
38  1   0   0   0   1   63      0.1     -33     0.0
37  0   0   0   2   0   63      0.0     -2      0.0
40  -1  0   2   2   2   -59     0.0     26      0.0
39  -1  0   0   0   1   -58     -0.1    32      0.0
41  1   0   2   0   1   -51     0.0     27      0.0
14  2   0   0   -2  0   48      0.0     1       0.0
3   -2  0   2   0   1   46      0.0     -24     0.0
42  0   0   2   2   2   -38     0.0     16      0.0
45  2   0   2   0   2   -31     0.0     13      0.0
43  2   0   0   0   0   29      0.0     -1      0.0
44  1   0   2   -2  2   29      0.0     -12     0.0
46  0   0   2   0   0   26      0.0     -1      0.0
15  0   0   2   -2  0   -22     0.0     0       0.0
47  -1  0   2   0   1   21      0.0     -10     0.0
16  0   2   0   0   0   17      -0.1    0       0.0
18  0   2   2   -2  2   -16     0.1     7       0.0
48  -1  0   0   2   1   16      0.0     -8      0.0 ...
];

% Map coefficients to radians (this should be embedded in the data)
% nc(:,end-3:end) = nc(:,end-3:end) * 0.0001 * asec2rad;

% Compute nutation in longitude (psi), and latitude (eps)
dpsi_1980 = 0; deps_1980 = 0;
for i = 1:size(nc,1)
    api         =   nc(i,2) * M_moon + nc(i,3) * M_sun + nc(i,4) * u_moon + nc(i,5) * D_sun + nc(i,6) * O_moon;
    dpsi_1980   =   0.0001 * (nc(i,7) + nc(i,8) * JC_J2000_TT) * sin(api) + dpsi_1980;
    deps_1980   =   0.0001 * (nc(i,9) + nc(i,10) * JC_J2000_TT) * cos(api) + deps_1980;
end

eps_1980 = epsb_1980 + asec2rad * deps_1980;
tod_to_mod = rot1(-epsb_1980) * rot3(asec2rad * dpsi_1980) * rot1(eps_1980);

% Generate Precession matrix
zeta    =   asec2rad * (2306.2181 * JC_J2000_TT + 0.30188 * JC_J2000_TT^2 + 0.017998 * JC_J2000_TT^3);
theta   =   asec2rad * (2004.3109 * JC_J2000_TT - 0.42665 * JC_J2000_TT^2 - 0.041833 * JC_J2000_TT^3);
z       =   asec2rad * (2306.2181 * JC_J2000_TT + 1.09468 * JC_J2000_TT^2 + 0.018203 * JC_J2000_TT^3);

mod_to_eci = rot3(zeta) * rot2(-theta) * rot3(z);

% Compute equation of Equinoxes
eqe_1980 = asec2rad * (dpsi_1980  * cos(epsb_1980) + 0.00264 * sin(O_moon) + 0.000063 * sin(2 * O_moon));
teme_to_tod = rot3(-eqe_1980);

% % Compute sidereal rotation
% GMST_1982 = asec2rad * (67310.54841 ...
%                         + (876600 * 3600 + 8640184.812866) * jd_ut1_j2000_century ...
%                         + 0.093104 * jd_ut1_j2000_century^2 ...
%                         - 6.2e-6 * jd_ut1_j2000_century^3);
GMST_1982 = (67310.54841 ...
                        + (876600 * 3600 + 8640184.812866) * JC_J2000_ut1 ...
                        + 0.093104 * JC_J2000_ut1^2 ...
                        - 6.2e-6 * JC_J2000_ut1^3);
GMST_1982 = mod(GMST_1982,sign(GMST_1982)*86400) / 240;                 
GMST_1982 = deg2rad * mod(GMST_1982, 360);
GAST_1982 = GMST_1982 + eqe_1980;

% Compute sidereal rotation
% GMST_1982 = asec2rad * (67310.54841 ...
%                         + (876600 * 3600 + 8640184.812866) * jd_ut1_j2000_century ...
%                         + 0.093104 * jd_ut1_j2000_century^2 ...
%                         - 6.2e-6 * jd_ut1_j2000_century^3);
% GMST_1982 = mod(GMST_1982, 2*pi);
% GAST_1982 = GMST_1982 + eqe_1980;

% Compute rotation from pseudo-Earth fixed frame (does not include polar
% motion)
pef_to_tod = rot3(-GAST_1982);

% % Compute polar motion (small angle approximation is assumed)
% itrf_to_pef = polar(polar_motion_rad(1), polar_motion_rad(2));

% Compute composite rotations
teme_to_eci = mod_to_eci * tod_to_mod * teme_to_tod;    
ecef_to_eci = mod_to_eci * tod_to_mod * pef_to_tod;% * itrf_to_pef;   % for vecef
ppef_to_veci = ecef_to_eci * skew([0; 0; w_prec]);     % for w x recef

end

function y = rot1(u)
y = [1 0 0; 0 cos(u) sin(u); 0 -sin(u) cos(u)];
end

function y = rot2(u)
y = [cos(u) 0 -sin(u); 0 1 0; sin(u) 0 cos(u)];
end

function y = rot3(u)
y = [cos(u) sin(u) 0; -sin(u) cos(u) 0; 0 0 1];
end

function y = skew(u)
y = [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];
end

function w = polar(x,y)
w = [1 0 -x; 0 1 y; x -y 1];
end


