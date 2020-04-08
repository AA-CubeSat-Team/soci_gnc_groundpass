function check_groundpass(obj,min_elevation)

err_ecef_m  = obj.pos_ecef_m - obj.gs.pos_ecef_m;
err_sez_m   = rot(90-obj.gs.lat_deg,'y') * rot(obj.gs.lon_deg,'z') * err_ecef_m;
err_mag     = norm(err_sez_m);
s_ell       = err_sez_m(3)/err_mag;
c_ell       = norm(err_sez_m(1:2))/err_mag;
ell         = acosd(c_ell);

if (obj.verbose)
    fprintf('elev=%2.2f deg  rho_z=%2.2f\n',ell,err_sez_m(3));
end

% check if we are overhead
if (s_ell>=0)
    % starting or ending a pass
    if (abs(ell)>(min_elevation) && ~obj.passing)
        obj.passes  = obj.passes + 1;
        JD      = obj.JD_J2000_utc + obj.params.JDJ2000;
        time    = get_timestamp(JD);
        fprintf('--------------------------------------------------\n');
        fprintf('START ground pass: #%02d at %04d:%02d:%02d - %02d:%02d:%04.1f\n',...
            obj.passes,time(1),time(2),time(3),time(4),time(5),time(6));
        obj.passing = true;
        obj.len_passes_s(obj.passes) = obj.MET_s;
    elseif ( (abs(ell)<min_elevation) && obj.passing )
        JD      = obj.JD_J2000_utc + obj.params.JDJ2000;
        time    = get_timestamp(JD);
        fprintf('END ground pass:   #%02d at %04d:%02d:%02d - %02d:%02d:%04.1f\n',...
            obj.passes,time(1),time(2),time(3),time(4),time(5),time(6));
        obj.passing = false;
        obj.len_passes_s(obj.passes) = obj.MET_s - obj.len_passes_s(obj.passes);
    end
end

end

