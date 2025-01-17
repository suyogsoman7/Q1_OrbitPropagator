% Main Script for Converting Equinoctial States to Cartesian and Plotting
clear; close all; clc;
format long;

% Constants
mu = 398600.4418; % Earth's gravitational parameter (km^3/s^2)

% File names
input_file = 'trajectory.csv';           % Input file with equinoctial elements
output_file = 'trajectory_cartesian.csv'; % Output file for Cartesian states

% Read equinoctial data from the input file
data = readmatrix(input_file);          % Read CSV
times = data(:, 1);                     % Extract time column
equinoctial_states = data(:, 2:end);    % Extract equinoctial states (p, f, g, h, k, L)


cartesian_states = zeros(size(equinoctial_states, 1), 6);

% Perform conversion for each row of states
for i = 1:size(equinoctial_states, 1)
    cartesian_states(i, :) = equinoctial_to_cartesian(mu, equinoctial_states(i, :));
end

% Combine time and Cartesian data
output_data = [times, cartesian_states];

% Open the output file for writing
fileID = fopen(output_file, 'w');

% Write header
fprintf(fileID, 'Time [s], Rx [km], Ry [km], Rz [km], Vx [km/s], Vy [km/s], Vz [km/s]\n');

% Write the data 
for i = 1:size(output_data, 1)
    fprintf(fileID, '%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f\n', output_data(i, :));
end

% Close the file
fclose(fileID);
fprintf('Data written to %s\n', output_file);

% Plot the trajectory
plot_trajectory(output_data);

% Animate the trajectory
animate_trajectory(output_data);

% Function to convert a single equinoctial state to Cartesian state
function rv = equinoctial_to_cartesian(mu, evec)
    % Extract equinoctial elements
    p = evec(1); f = evec(2); g = evec(3);
    h = evec(4); k = evec(5); L = evec(6);

    % Compute required intermediate values
    k2 = k^2; h2 = h^2;
    tkh = 2 * k * h;
    s2 = 1 + h2 + k2;

    cL = cos(L); sL = sin(L);
    w = 1 + f * cL + g * sL;
    r = p / w;

    smp = sqrt(mu / p);

    % Unit vectors in orbital plane
    fhat = [(1 - k2 + h2), tkh, -2 * k] / s2;
    ghat = [tkh, (1 + k2 - h2), 2 * h] / s2;

    % Orbital position and velocity in-plane
    x = r * cL; y = r * sL;
    xdot = -smp * (g + sL);
    ydot = smp * (f + cL);

    % Combine into Cartesian state vectors
    rv(1:3) = x * fhat + y * ghat;       % Position (rx, ry, rz)
    rv(4:6) = xdot * fhat + ydot * ghat; % Velocity (vx, vy, vz)
end

% Function to plot the trajectory in 3D
function plot_trajectory(data)
    % Extract Cartesian positions
    rx = data(:, 2);
    ry = data(:, 3);
    rz = data(:, 4);

    % Plot trajectory in 3D
    figure(1);
    plot3(rx, ry, rz, 'b-', 'LineWidth', 1.5); hold on;
    
    % Highlight start and end points
    plot3(rx(1), ry(1), rz(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g'); % Start (green)
    plot3(rx(end), ry(end), rz(end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % End (red)

    % Add labels and grid
    xlabel('x [km]');
    ylabel('y [km]');
    zlabel('z [km]');
    grid on;
    axis equal;
    title('3D Trajectory');
    legend('Trajectory', 'Start Point', 'End Point');
end

% Function to animate the trajectory in 3D
function animate_trajectory(data)
    % Extract Cartesian positions and time
    rx = data(:, 2);
    ry = data(:, 3);
    rz = data(:, 4);
    times = data(:, 1);

    
    figure(2);
    hold on;
    grid on;
    axis equal;
    xlabel('x [km]');
    ylabel('y [km]');
    zlabel('z [km]');
    title('Trajectory Animation');

    % Plot full trajectory for reference
    plot3(rx, ry, rz, 'b--', 'LineWidth', 1); 

    % Add the animated point
    h = plot3(rx(1), ry(1), rz(1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

    % Set default 3D view
    view(45, 30); % Azimuth = 45°, Elevation = 30°
    rotateAngle = 0; % Initial rotation angle

    % Animation loop
    for i = 1:length(times)
        % Update point position
        set(h, 'XData', rx(i), 'YData', ry(i), 'ZData', rz(i));
        
        % Smoothly rotate the view for interactivity
        rotateAngle = rotateAngle + 0.2; % Increment rotation angle
        view(45 + rotateAngle, 30); % Update view dynamically

        % Pause to create smooth animation
        pause(0.05); % Adjust to control animation speed
    end

    % Highlight the end point after animation
    plot3(rx(end), ry(end), rz(end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
end
