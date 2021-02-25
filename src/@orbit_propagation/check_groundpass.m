function check_groundpass(obj,min_elevation)

err_ecef_m  = obj.pos_ecef_m - obj.gs.pos_ecef_m;
err_sez_m   = rot(90-obj.gs.lat_deg,'y') * rot(obj.gs.lon_deg,'z') * err_ecef_m;
err_mag_m   = norm(err_sez_m);
s_ell       = err_sez_m(3)/err_mag_m;
% c_ell       = norm(err_sez_m(1:2))/err_mag_m;
ell         = atan2d( err_sez_m(3), norm(err_sez_m(1:2)) );

if (obj.logging)
    fprintf(obj.log_file_id,obj.log_file_spec(1,3),...
                                         ell,err_sez_m(3)*obj.params.m2km);
end

% check if we are overhead
if (s_ell>=0)
    % print current slant range to log file
    if (obj.logging)
        fprintf(obj.log_file_id,obj.log_file_spec(1,4),...
                                        err_mag_m*obj.params.m2km);
    end
    % starting or ending a pass
    if (abs(ell)>(min_elevation) && ~obj.passing)
        % get current time stamp
        JD      = obj.JD_J2000_utc + obj.params.JDJ2000;
        % correct for local time relative to UTC
        JD = JD + obj.gs.TZ * obj.params.hr2day;
        time    = get_timestamp(JD);
        % set passing == true and record the current MET at pass begin
        obj.passing = true;
        obj.start_pass(obj.MET_s);
        % print lines
        fprintf('--------------------------------------------------\n');
        fprintf('START ground pass: #%02d at %04d:%02d:%02d - %02d:%02d:%04.1f\n',...
            obj.passes,time(1),time(2),time(3),time(4),time(5),time(6));
        
    elseif ( (abs(ell)<min_elevation) && obj.passing )
        
        % get current time stamp
        JD      = obj.JD_J2000_utc + obj.params.JDJ2000;
        % correct for local time relative to UTC
        JD = JD + obj.gs.TZ * obj.params.hr2day;
        time    = get_timestamp(JD);
        % set passing == false and record the duration of pass
        obj.passing = false;
        obj.end_pass(obj.MET_s);
        % print lines
        fprintf('END ground pass:   #%02d at %04d:%02d:%02d - %02d:%02d:%04.1f\n',...
            obj.passes,time(1),time(2),time(3),time(4),time(5),time(6));
    end
else
    if (obj.logging)
        fprintf(obj.log_file_id,'\n');
    end
end

end

