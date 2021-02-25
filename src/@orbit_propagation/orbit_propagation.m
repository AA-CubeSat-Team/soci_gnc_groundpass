classdef orbit_propagation < handle
    
    properties
        TLE(1,1) string
        MET_s(1,1) double
        params(1,1) parameters
        gs(1,1) ground_station
        orbit_tle(9,1) double
        make_plot(1,1) logical = false
        dT(1,1) double = 1.0
    end
    
    properties (SetAccess = private)
       passing(1,1) logical = false
       len_passes_s(:,1) double
       pos_teme_km(3,1) double
       vel_teme_kmps(3,1) double
       lla(3,:) double
       slr(1,:) double
       logging(1,1) logical
       log_file(1,1) string
       sgp4_flag(1,1) int32
    end
    
    properties(Access=private)
        log_file_id(1,1) double
        log_file_head(1,1) string
        log_file_spec(1,4) string
    end
    
    properties (Dependent,Hidden)
        passes(1,1) uint32
        JD_J2000_utc(1,1) double
        JD_J2000_ut1(1,1) double
        JD_J2000_TT(1,1) double
        JC_J2000_ut1(1,1) double
        JC_J2000_TT(1,1) double
        pos_ecef_m(3,1) double
        vel_ecef_mps(3,1) double
    end
    
    methods (Access=public)
        % basic constructor
        function obj = orbit_propagation(TLE_filename)
            if ( nargin<1 || isempty(TLE_filename) )
                TLE_filename = 'ISS_TLE.txt';
            end
            [~,filename,~] = fileparts(TLE_filename);
            obj.TLE = filename;
            obj.orbit_tle = get_tle(TLE_filename);
        end
        % function to increment the internal Mission Elapsed Time (MET)
        function increment_time(obj)
            obj.MET_s = obj.MET_s + obj.dT;
        end
        % function to enable data logging, create and open the log file
        function enable_logging(obj,enable)
            obj.logging = enable;
            rn = datetime;
            datetime_fmt = 'yy_mm_dd_HH_MM_SS';
            if ~exist('logs', 'dir')
                mkdir('logs')
            end
            obj.log_file = strcat('logs/',obj.TLE,'_',...
                                        datestr(rn,datetime_fmt),'.txt');
            obj.log_file_head = strcat(' MET [s]  flag  lat [dg]  lon [dg]',...
                                       '  alt [km]  elev [dg]  rho_z [km]',...
                                       '  slant_range [km]');
            if (enable)
                obj.log_file_id = fopen( obj.log_file, 'wt' );
                % write the header to the file
                fprintf(obj.log_file_id,obj.log_file_head);
                fprintf(obj.log_file_id,'\n');
            end
            obj.log_file_spec(1,1) = '%08.2f   %d ';
            obj.log_file_spec(1,2) = '   %+06.2f    %+07.2f     %06.2f ';
            obj.log_file_spec(1,3) = '   %+06.2f     %+09.2f ';
            obj.log_file_spec(1,4) = '   %07.2f\n';
        end
        % function for coordinate rotations
        [ecef_to_eci, ppef_to_veci, mod_to_eci, teme_to_eci] = coordinate_rotations(obj);
        % function for orbit propagation using SGP4
        sgp4(obj);
        % function to check for a ground pass
        check_groundpass(obj,min_elevation)
        % function to post process the results
        post_process(obj)
        % function to initialize the plot
        function initialize_plot(obj,elev_deg)
            % initialize printout
            fprintf('--------------------------------------------------\n');
            fprintf('   _       _     ___  _____ \n')
            fprintf('  / \\     / \\   |       |   \n')
            fprintf(' /___\\   /___\\  |       |   \n')
            fprintf('/     \\ /     \\ |___    |   \n\n')
            fprintf('Orbital Pass Predictor Tool: ')
            fprintf(' All times are UTC%+2d\n',obj.gs.TZ)
            fprintf('By: T. P. Reynolds\n')
            % get data for field of view circle
            gs_pos_ecef_m = obj.gs.pos_ecef_m;
            ortho_dir = gs_pos_ecef_m./norm(gs_pos_ecef_m);
            % compute orbit semi-major axis to use to get radius of FoV
            mean_motion_radps = obj.orbit_tle(9) * (2*pi/obj.params.day2sec); 
            a_m = (obj.params.mu_earth_m3ps2/(mean_motion_radps^2))^(1/3);
            % altitude to compute FoV cone
            h_m = a_m - obj.params.R_earth_m;
            % radius of FoV cone at this altitude
            radius_m = h_m * tand(90-elev_deg);
            % center of the FoV circle in SEZ frame (origin at ground station)
            center_sez_m = h_m.*ortho_dir;
            angles  = linspace(0,2*pi);
            % coordinates of the FoV circle in the SEZ frame
            fov_sez_m = center_sez_m + radius_m .* [ cos(angles); 
                                                     sin(angles); 
                                                     zeros(size(angles)) ];
            % map FoV coordinates to ECEF frame
            R_sez2ecef = rot(obj.gs.lon_deg,'z')' * rot(90-obj.gs.lat_deg,'y')';
            fov_ecef_m = gs_pos_ecef_m + R_sez2ecef * fov_sez_m; 
            % compute the LLA of the FoV coordinates
            fov_lla = ecef2lla(fov_ecef_m');
            figure(1), clf, hold on
            axis([-180 180 -90 90]);
            Earth_im = imread('figs/Flat_earth.jpg');
            imagesc([-180 180],[-90 90],flipud(Earth_im));
            plot(obj.gs.lon_deg,obj.gs.lat_deg,'c*',...
                    'MarkerSize',3,'MarkerFaceColor','c')
            plot(fov_lla(:,2),fov_lla(:,1),'c')
            set(gcf,'Units','Normalized','Position',[0.5,0.4,0.5,0.5]);
            set(gca,'YDir','normal',...
                    'XTick',-180:30:180,...
                    'YTick',-90:15:90,...
                    'XGrid','on','YGrid','on',...
                    'GridAlpha',0.2,...
                    'TickLabelInterpreter','latex',...
                    'FontSize',15);
            xlabel('Longitude [deg]','FontSize',16,'Interpreter','latex')
            ylabel('Latitude [deg]','FontSize',16,'Interpreter','latex')
        end
        % function to update the plot
        function update_plot(obj,theta)
            opts = obj.params.plot;
            % compute current LLA
            lla_now = ecef2lla( obj.pos_ecef_m' );
            if (obj.make_plot)
                if (mod(obj.MET_s,opts.plot_density_s)==0)
                    count = opts.counter; 
                    count_mod = mod(count,opts.view_num_points)+1;
                    if (count>opts.view_num_points)
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
                obj.params.plot = opts;
            else
                if (mod(obj.MET_s,obj.params.plot.plot_density_s)==0)
                    count = obj.params.plot.counter;
                    obj.lla(:,count) = lla_now;
                    obj.params.plot.counter = count + 1; 
                end
            end
            if (obj.logging)
                fprintf(obj.log_file_id,obj.log_file_spec(1,2),...
                    lla_now(1),lla_now(2),lla_now(3)*obj.params.m2km);
            end
        end
        % function to save the time at the start of a ground pass
        function start_pass(obj,time)
            % input time can be any frame/value, so long as it is
            % consistent with the input to end_pass 
            obj.len_passes_s = orbit_propagation.push(obj.len_passes_s,time);
        end
        % function to save the time at the end of a ground pass
        function end_pass(obj,time)
            % input time can be any frame/value, so long as it is
            % consistent with the input to start_pass 
            obj.len_passes_s(end) = time - obj.len_passes_s(end);
        end
    end
        
    methods (Static)
        function vec_out = push(vec_in,input)
            vec_out = [vec_in;input];
        end
        function check_sgp4_flag(flag)
            switch flag
                case 0
                    err = false;
                case 1
                    err = true;
                    msg = 'Warning: low altitude (<220km)';
                case {-1,-2,-3,-4}
                    err = true;
                    msg = sprintf('Warning: SGP4 reports error code %d',...
                                    abs(flag));
            end
            if (err)
                fprintf(msg);
            end
        end
    end
    
    % get/set methods
    methods
        function passes = get.passes(obj)
           passes = numel(obj.len_passes_s);
        end
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
    end
end