classdef ground_station < handle
    
    properties
        lat_deg(1,1) double
        lon_deg(1,1) double
        alt_m(1,1) double
        pos_ecef_m(3,1) double
    end
    
    methods
        function init(obj,lat_deg,lon_deg,alt_m)
            obj.lat_deg     = lat_deg;
            obj.lon_deg     = lon_deg;
            obj.alt_m       = alt_m;
            obj.pos_ecef_m  = lla2ecef([lat_deg,lon_deg,alt_m])';
        end
    end
end

