% This script propagates the current SOCi orbit that is stored in
%
%   fsw_params.env_estimation.orb_estimation.sgp4.orbit_tle
%
% in order to plot the ground track and predict UW ground station passes.
% To change the orbit, you must change what TLE is loaded by `sim_init.m`
% 
% For fast prediction, turn off the plotting function by setting PLOT=false
% and DEBUG=false. This will propagate the orbit using the SGP-4 method,
% which is AN APPROXIMATION ONLY. The longer you propagate the orbit for,
% the less you should trust the results. Short term predictions (<4 days)
% should be fine.
%
% The plot that is generate has the camera's field of view outlined in red,
% the ground track plotted in red, and the UW ground station in cyan. After
% 30 minutes of propagation the tail of the ground track will start to
% disappear - this is just so that the plot doesn't get too big and start
% to really slow down.
%
% A ground station pass is defined by the elevation above the horizon, and
% can be set by changing the variable `elev`. For example, SOCi
% communication hardware should be able to communicate with the satellite
% whenever it has an elevation of at least 15 degress above the local
% horizon at the UW ground station.
%
% Printing format is YYYY:MM:DD - HH:MM:SS.S
%
% YOU MUST RUN SIM_INIT.M FIRST BEFORE RUNNING THIS FILE
%
% T. Reynolds -- RAIN Lab

close all; clc;
addpath(genpath('src/'))
addpath('tools/')
addpath('TLEs/')

% simulation options
MET     = 0.0;       % initial MET in [s]
MET_end  = 86400.0;       % end MET in [s]
theta    = 56/2;         % camera FoV half angle
elev_deg = 10;           % elevation above ground station where s/c is visible
lat_T    = 47.655548;    % ground station latitude
lon_T    = -122.303200;  % ground station longitude
alt_T    = 0.0;          % ground station altitude

% instantiate an orbit propagation object
propagator = orbit_propagation();
propagator.gs.init(lat_T,lon_T,alt_T);

% set custom options
propagator.dT = 1;

kk        = 1;   
dR        = 30;
kk_tot    = 200/propagator.dT;
psc       = cell(kk_tot,1);
lens      = NaN(100,1);
LLDATA    = NaN(2,MET_end/(30*propagator.dT));

% initialize figure
if (propagator.make_plot)
    fig = figure;
    ax  = axes('Parent',fig);
    axis([-180 180 -90 90]);
    Earth_im = imread('Flat_earth.jpg');
    imagesc([-180 180],[-90 90],flipud(Earth_im));
    set(gca,'YDir','normal');
    hold on
    pc = plot(1,1,'Visible','off');  
    plot(lon_T,lat_T,'c*','MarkerSize',3,'MarkerFaceColor','c')
%     hold off
    set(fig,'Units','Normalized','Position',[0.46736,0.42667,0.53264,0.51111]);
    set(gca,'XTick',-180:30:180);
    set(gca,'YTick',-90:15:90);
    set(gca,'XGrid','on','YGrid','on','GridAlpha',0.2);
end

% main loop
while(true)
   
    % current UTC time
    JD_J2000_utc = propagator.JD_J2000_utc;
    
    % map current time to other frames
    JD_J2000_ut1 = propagator.JD_J2000_ut1;
    JD_J2000_TT  = propagator.JD_J2000_TT;
    T_J2000_ut1  = propagator.JC_J2000_ut1;
    T_J2000_TT   = propagator.JC_J2000_TT;
        
    % coordinate transformations
    [ecef_to_eci,~,~,teme_to_eci] = propagator.coordinate_rotations();
    eci_to_ecef = ecef_to_eci';
    
    % propagate position with SGP4
    flag = propagator.sgp4();
    
    % map TEME position to ECEF frame in meters
    pos_ecef_m = propagator.pos_ecef_m;
    
    % map to LLA
    lla = ecef2lla( reshape(pos_ecef_m,1,3) );
    
    % get camera fov
    [ll_circle] = get_ground_circle(lla,propagator.params.R_earth_m,theta);
    
    % plot point
    if (propagator.make_plot)
        if (mod(MET,dR)==0)
            kk_ = mod(kk,kk_tot)+1;
            if (kk>kk_tot)
                delete(psc{kk_});
            end
            figure(fig)
            delete(pc);
            psc{kk_} = plot(lla(2),lla(1),'ro','MarkerSize',3,'MarkerFaceColor','r');
    %         pc = plot(ll_circle(:,2),ll_circle(:,1),'r');
            kk = kk + 1;
        end
    else
        if (mod(MET,dR)==0)
            LLDATA(:,kk) = [ lla(2); lla(1) ];
            kk = kk + 1;
        end
    end
    
    if (propagator.verbose)
        fprintf('MET=%2.2f s  lat=%2.2f deg  long=%2.2f deg  alt=%2.2f km  FLAG=%d',...
                MET,lla(1),lla(2),lla(3)*m2km,flag);
    end
    
    % check for a ground pass
    propagator.check_groundpass(elev_deg);
    
    % increment time
    propagator.increment_time();
    
    % exit loop when we reach MET_end
    if (propagator.MET_s>MET_end)
        break;
    end
end

lens    = rmmissing(lens);
max_len = max(lens);
min_len = min(lens);
mu_len  = mean(lens);
std_len = std(lens);
fprintf('==================================================\n');
fprintf('TOTAL passes: %d, TIME overhead: %2.2f +/- %2.2fs \n\n',...
            passes,mu_len,std_len);
        
% Plot ground track afterwards if not done during sim
if (~propagator.make_plot)
    fig = figure;
    ax  = axes('Parent',fig);
    axis([-180 180 -90 90]);
    Earth_im = imread('Flat_earth.jpg');
    imagesc([-180 180],[-90 90],flipud(Earth_im));
    set(gca,'YDir','normal');
    hold on
    pc = plot(1,1,'Visible','off');  
    plot(lon_T,lat_T,'c*','MarkerSize',3,'MarkerFaceColor','c')
    plot(LLDATA(1,:),LLDATA(2,:),'ro','MarkerSize',3,'MarkerFaceColor','r')
    set(fig,'Units','Normalized','Position',[0.46736,0.42667,0.53264,0.51111]);
    set(gca,'XTick',-180:30:180);
    set(gca,'YTick',-90:15:90);
    set(gca,'XGrid','on','YGrid','on','GridAlpha',0.2);
end

edges = linspace(min_len,max_len,10);
quants_len = quantile(lens,[0.25,0.5,0.75]);
figure, hold on, grid on, box on
histogram(lens,edges)
for k = 1:numel(quants_len)
    plot([quants_len(k) quants_len(k)],get(gca,'Ylim'),'r--','LineWidth',1)
end