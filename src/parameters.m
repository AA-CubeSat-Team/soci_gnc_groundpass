classdef parameters
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % time parameters
        dut1(1,1) double = 0
        sec2day(1,1) double = 1/86400
        day2sec(1,1) double = 86400
        JDJ2000(1,1) double = 2451545.0
        JD2cent(1,1) double = 1/36525
        cent2JD(1,1) double = 36525
        gps2utc(1,1) double = -18
        utc2gps(1,1) double = 18
        gps2tai(1,1) double = 19.0
        dAT(1,1) double = 37
        gps2tt(1,1) double
        tt2gps(1,1) double
        % length parameters
        km2m(1,1) double = 1e3
        m2km(1,1) double = 1e-3
        R_earth_m(1,1) double = 6378.137e3
        % sim parameters
        dT(1,1) double = 1
        % misc
        w_prec(1,1) double = 7.292115146706979e-5 
    end
    
    properties (Dependent)
        gps2ut1(1,1) double
    end
    
    methods 
        function obj = parameters()
            obj.gps2utc = obj.gps2tai - obj.dAT;
            obj.utc2gps = -obj.gps2utc;
            obj.gps2tt  = 32.184 + obj.gps2tai;
            obj.tt2gps  = - obj.gps2tt;
        end
    end
    
    methods
        function gps2ut1 = get.gps2ut1(obj)
            gps2ut1 = obj.gps2utc + obj.dut1;
        end
    end
end

