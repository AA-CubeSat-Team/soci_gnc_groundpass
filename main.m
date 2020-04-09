% This script propagates an orbit, displays the ground track, and identifies
% ground station passes.
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
% T. Reynolds -- RAIN Lab

close all; clc;
addpath(genpath('src/'))
addpath('tools/')
addpath('TLEs/')

% simulation options
MET_end  = 86400.0;      % end MET in [s]
theta    = 56/2;         % camera FoV half angle
elev_deg = 10;           % elevation above ground station where s/c is visible
lat_T    = 47.655548;    % ground station latitude
lon_T    = -122.303200;  % ground station longitude
alt_T    = 0.0;          % ground station altitude

% instantiate an orbit propagation object
propagator = orbit_propagation();
% initialize the ground station
propagator.gs.init(lat_T,lon_T,alt_T);

% set custom propagator options
propagator.dT = 1;
propagator.make_plot = false;
propagator.params.plot.plot_density_s = 60;
propagator.params.plot.show_camera = true;

% initialize figure
propagator.initialize_plot();

% main loop
while(true)
                   
    % propagate position with SGP4
    flag = propagator.sgp4();
        
    % update plot
    propagator.update_plot(theta);
    
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

lens = rmmissing(propagator.len_passes_s);
max_len = max(lens);
min_len = min(lens);
mu_len  = mean(lens);
std_len = std(lens);
fprintf('==================================================\n');
fprintf('TOTAL passes: %d, TIME overhead: %2.2f +/- %2.2fs \n\n',...
            propagator.passes,mu_len,std_len);
        
% Plot ground track afterwards if not done during sim
plot(propagator.lla(2,:),propagator.lla(1,:),'ro','MarkerSize',3,'MarkerFaceColor','r')

if (propagator.passes>1)
    edges = linspace(min_len,max_len,10);
    quants_len = quantile(lens,[0.25,0.5,0.75]);
    figure, hold on, grid on, box on
    histogram(lens,edges)
    for k = 1:numel(quants_len)
        plot([quants_len(k) quants_len(k)],get(gca,'Ylim'),'r--','LineWidth',1)
    end
end