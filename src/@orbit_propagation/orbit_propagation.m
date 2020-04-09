classdef orbit_propagation < handle
    
    properties
        MET_s(1,1) double
        params(1,1) parameters
        gs(1,1) ground_station
        orbit_tle(9,1) double
        make_plot(1,1) logical = false
        verbose(1,1) logical = false
        dT(1,1) double = 1.0
    end
    
    properties (SetAccess = private)
       passes(1,1) uint32 = 0
       passing(1,1) logical = false
       len_passes_s(:,1) double
       pos_teme_km(3,1) double
       vel_teme_kmps(3,1) double
       lla(3,:) double
    end
    
    properties (Dependent)
        JD_J2000_utc(1,1) double
        JD_J2000_ut1(1,1) double
        JD_J2000_TT(1,1) double
        JC_J2000_ut1(1,1) double
        JC_J2000_TT(1,1) double
        pos_ecef_m(3,1) double
        vel_ecef_mps(3,1) double
        view_num_points(1,1) double
    end
    
    methods (Access=public)
        function obj = orbit_propagation(TLE_filename)
            if ( nargin<1 || isempty(TLE_filename) )
                TLE_filename = 'ISS_TLE.txt';
            end
            obj.orbit_tle = get_tle(TLE_filename);
        end
        % function to increment the internal Mission Elapsed Time (MET)
        function increment_time(obj)
            obj.MET_s = obj.MET_s + obj.dT;
        end
        
        % function for coordinate rotations
        [ecef_to_eci, ppef_to_veci, mod_to_eci, teme_to_eci] = coordinate_rotations(obj);
        % function for orbit propagation using SGP4
        flag = sgp4(obj);
        
        % function to check for a ground pass
        check_groundpass(obj,min_elevation)
        % function to initialize the plot
        function initialize_plot(obj)
            figure(1), clf, hold on
            axis([-180 180 -90 90]);
            Earth_im = imread('figs/Flat_earth.jpg');
            imagesc([-180 180],[-90 90],flipud(Earth_im));
            plot(obj.gs.lon_deg,obj.gs.lat_deg,'c*',...
                    'MarkerSize',3,'MarkerFaceColor','c')
            set(gcf,'Units','Normalized','Position',[0.5,0.4,0.5,0.5]);
            set(gca,'YDir','normal',...
                    'XTick',-180:30:180,...
                    'YTick',-90:15:90,...
                    'XGrid','on','YGrid','on',...
                    'GridAlpha',0.2);
        end
        % function to update the plot
        function update_plot(obj,theta)
            opts = obj.params.plot;
            % compute to LLA
            lla_now = ecef2lla( obj.pos_ecef_m' );
            if (obj.make_plot)
                if (mod(obj.MET_s,opts.plot_density_s)==0)
                    count = opts.counter; 
                    count_mod = mod(count,obj.view_num_points)+1;
                    if (count>obj.view_num_points)
                        delete(opts.handles{count_mod});
                    end
                    figure(1)
                    % get camera fov
                    if (opts.show_camera)
                        % get lat/lon of intersection of camera FoV and 
                        % spherical Earth
                        [ll_circle] = get_ground_circle(lla_now,...
                                                obj.params.R_earth_m,theta);
                        try
                            delete(opts.camera_handle)
                        catch
                            % no worries
                        end
                        % plot camera FoV
                        opts.camera_handle = ...
                                        plot(ll_circle(:,2),ll_circle(:,1),'r');
                    end
                    % plot current point
                    opts.handles{count_mod} = plot(lla_now(2),lla_now(1),'ro',...
                                        'MarkerSize',3,'MarkerFaceColor','r');
                    % increment counter
                    opts.counter = count + 1; 
                end
            else
                if (mod(obj.MET_s,obj.params.plot.plot_density_s)==0)
                    count = obj.params.plot.counter;
                    obj.lla(:,count) = lla_now;
                    obj.params.plot.counter = count + 1; 
                end
            end
        end
    end
    
    % get/set methods
    methods
        function JD_J2000_utc = get.JD_J2000_utc(obj)
           JD_J2000_utc = obj.orbit_tle(2) + obj.params.sec2day * obj.MET_s; 
        end
        function JD_J2000_ut1 = get.JD_J2000_ut1(obj)
            JD_J2000_ut1 = obj.JD_J2000_utc + obj.params.dut1;
        end
        function JD_J2000_TT = get.JD_J2000_TT(obj)
            temp = obj.params.sec2day * (obj.params.utc2gps + obj.params.gps2tt);
            JD_J2000_TT = obj.JD_J2000_utc + temp;
        end
        function JC_J2000_ut1 = get.JC_J2000_ut1(obj)
            JC_J2000_ut1 = obj.JD_J2000_ut1 * obj.params.JD2cent;
        end
        function JC_J2000_TT = get.JC_J2000_TT(obj)
            JC_J2000_TT = obj.JD_J2000_TT * obj.params.JD2cent;
        end
        function pos_ecef_m = get.pos_ecef_m(obj)
            [ecef_2_eci,~,~,teme_2_eci] = obj.coordinate_rotations();
            pos_teme_m = obj.params.km2m * obj.pos_teme_km;
            pos_ecef_m = ecef_2_eci' * teme_2_eci * pos_teme_m;
        end
        function vel_ecef_mps = get.vel_ecef_mps(obj)
            [ecef_2_eci,ppef_2_veci,~,teme_2_eci] = obj.coordinate_rotations();
            vel_teme_mps = obj.params.km2m * obj.vel_teme_kmps;
            vel_ecef_mps = ecef_2_eci' * ( teme_2_eci * vel_teme_mps ...
                                            - ppef_2_veci * obj.pos_ecef_m );
        end
        function view_num_points = get.view_num_points(obj)
           view_num_points = ceil(obj.params.plot.view_last_s/obj.dT);
        end
    end
end