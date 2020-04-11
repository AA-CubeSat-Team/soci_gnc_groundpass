function [ll_circle] = get_ground_circle(lla,Re,theta)

c   = tand(theta);
z0  = lla(3) + Re;
z   = (c^2*z0 + sqrt(c^2*(Re^2-z0^2) + Re^2))/(c^2+1);
a   = sqrt(Re^2-z^2);
xi  = linspace(0,2*pi);

lambda  = lla(2);       % longitude
phi     = 90 - lla(1);  % co-lattitude
coordT  = [ a.*cos(xi); a.*sin(xi); z.*ones(size(xi)) ];
ECEF    = rot(-lambda,'z') * rot(-phi,'y') * coordT;
LLA     = ecef2lla(ECEF');

ll_circle = LLA(:,1:2);

end

