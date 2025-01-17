function meqoe = classical_to_equinoctial(coe)

    % Converts classical orbital elements into Modified Equinoctial Orbital
    % elements (MEqOE)

    format long;
    % Classical Orbital Elements
    a = coe(1);  % Semi-major axis (km)
    e = coe(2);  % Eccentricity
    i = coe(3);  % Inclination (deg.)
    capomega = coe(4); % RAAN (deg.)
    omega = coe(5);  % Argument of Perigee (deg.)
    nu = coe(6); % True Anomaly (deg.)

    % Compute MEqOE
    p = a * (1 - e^2); % Semi-Latus Rectum (km)
    f = e * cosd(omega + capomega);
    g = e * sind(omega + capomega);
    h = tand(i / 2) * cosd(capomega);
    k = tand(i / 2) * sind(capomega);
    L = mod(deg2rad(nu + omega + capomega), 2 * pi); % Normalize to [0, 2*pi] (rad)

    % MEqOE vector
    meqoe = [p, f, g, h, k, L];

end
