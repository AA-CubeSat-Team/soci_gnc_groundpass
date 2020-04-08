classdef orbit_propagation < handle
    
    properties
        MET_s(1,1) double
        params(1,1) parameters
        gs(1,1) ground_station
        orbit_tle(9,1) double
        pos_teme_km(3,1) double
        vel_teme_kmps(3,1) double
        make_plot(1,1) logical = false
        verbose(1,1) logical = false
        dT(1,1) double = 1.0
    end
    
    properties (SetAccess = private)
       passes(1,1) uint32 = 0
       passing(1,1) logical = false
       len_passes_s(:,1) double
    end
    
    properties (Dependent)
        JD_J2000_utc(1,1) double
        JD_J2000_ut1(1,1) double
        JD_J2000_TT(1,1) double
        JC_J2000_ut1(1,1) double
        JC_J2000_TT(1,1) double
        pos_ecef_m(3,1) double
        vel_ecef_mps(3,1) double
    end
    
    methods (Access=public)
        function obj = orbit_propagation(TLE_filename)
            if ( nargin<1 || isempty(TLE_filename) )
                TLE_filename = 'ISS_TLE.txt';
            end
            obj.orbit_tle = get_tle(TLE_filename);
        end
        
        function increment_time(obj)
            obj.MET_s = obj.MET_s + obj.dT;
        end
        
        % function for coordinate rotations
        [ecef_to_eci, ppef_to_veci, mod_to_eci, teme_to_eci] = coordinate_rotations(obj);
        % main function for orbit propagation using SGP4
        flag = sgp4(obj);
        % function to check for a ground pass
        check_groundpass(obj,min_elevation)
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
    end
end