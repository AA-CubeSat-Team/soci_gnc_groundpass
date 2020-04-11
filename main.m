% This script propagates an orbit, displays the ground track, and identifies
% ground station passes. Function currently only supports *one* ground
% station
% 
% You can enable/disable the plotting and logging functionality by setting 
%   propagator.make_plot = true/false
%   propagator.enable_logging(true/false)
% Logged data will be sent to the `ROOT/logs/` folder
% 
% In general, the right way to change any parameters is to use the syntax
%   propagator.<paramter> = <value>
%   propagator.params.<parameter> = <value>
% A list of parameters (and their defaults) can be found in the parameters.m 
% class file 
%  
% This tool propagates the s/c orbit using the SGP-4 method,
% which is AN APPROXIMATION ONLY. The longer you propagate the orbit for,
% the less you should trust the results. Short term predictions (<4 days)
% should be fine.
%
% The plot that is generated can have the camera's field of view outlined in red,
% the ground track plotted in red, and the UW ground station in cyan. The
% camera's FoV can be turned on/off using
%   propagator.params.plot.show_camera = true/false
%
% After 30 minutes of propagation the tail of the ground track will start to
% disappear - this is just so that the plot doesn't get too big and start
% to really slow down.
%
% A ground station pass is defined by the elevation above the horizon, and
% can be set by changing the variable `elev_deg`. For example, SOCi
% communication hardware should be able to communicate with the satellite
% whenever it has an elevation of at least 10-15 degress above the local
% horizon at the UW ground station.
%
% Printing format for passes is YYYY:MM:DD - HH:MM:SS.S
%
% At the end of the run, if any passes were recorded there is a figure that
% pops up to show the distribution of the length of the passes. 
%
% T. Reynolds -- RAIN Lab

clearvars; close all; clc;
addpath(genpath('src/'))
addpath('tools/')
addpath('TLEs/')

% simulation options
MET_end  = 2400.0;      % end MET in [s]
theta    = 56/2;         % camera FoV half angle
elev_deg = 10;           % elevation above ground station where s/c is visible
lat_T    = 47.655548;    % ground station latitude
lon_T    = -122.303200;  % ground station longitude
alt_T    = 0.0;          % ground station altitude

% instantiate an orbit propagation object
propagator = orbit_propagation('ISS_TLE.txt');

% initialize the ground station
propagator.gs.init(lat_T,lon_T,alt_T);

% set custom propagator options
propagator.dT = 1;
propagator.make_plot = true;
propagator.params.plot.plot_density_s = 60;
propagator.params.plot.show_camera = false;
propagator.params.plot.view_num_points = 30;
propagator.enable_logging(true);

% initialize figure
propagator.initialize_plot(elev_deg);

% main loop
while(true)
                   
    % propagate position with SGP4
    propagator.sgp4();
        
    % update plot
    propagator.update_plot(theta);
        
    % check for a ground pass
    propagator.check_groundpass(elev_deg);
    
    % increment time
    propagator.increment_time();
    
    % exit loop when we reach MET_end
    if (propagator.MET_s>MET_end)
        break;
    end
end

% post process the results; updates some plots and computes the statistics
% of the ground pass data
propagator.post_process()