function double()

%INPUT
m1 = input('Enter mass m1 (kg): ');
m2 = input('Enter mass m2 (kg): ');
L1 = input('Enter length L1 (m): ');
L2 = input('Enter length L2 (m): ');
theta1_0 = input('Enter initial angular position of bob 1 (rad): ');
theta2_0 = input('Enter initial angular position of bob 2 (rad): ');
omega1_0 = input('Enter initial angular velocity of bob 1 (rad/s): ');
omega2_0 = input('Enter initial angular velocity of bob 2 (rad/s): ');

%CONSTANTS
g = 9.81;
y0 = [theta1_0; omega1_0; theta2_0; omega2_0];

%TIME SETTINGS
t_final = 10;
dt = 0.1;
tspan = 0:dt:t_final;

%SOLVE SYSTEM
opts = odeset('RelTol',1e-10, 'AbsTol',1e-10);
[t, y] = ode45(@(t,y) equations(t, y, m1, m2, L1, L2, g), tspan, y0, opts);

%COORDINATES
x1 = L1 * sin(y(:,1));
y1 = -L1 * cos(y(:,1));
x2 = x1 + L2 * sin(y(:,3));
y2 = y1 - L2 * cos(y(:,3));

%ENERGY, VELOCITY, ACCELERATION CALCULATION
theta1 = y(:,1); omega1 = y(:,2);
theta2 = y(:,3); omega2 = y(:,4);

% Velocities
v1 = L1 * omega1;
v2 = sqrt((L1*omega1).^2 + (L2*omega2).^2 + ...
    2*L1*L2.*omega1.*omega2.*cos(theta1 - theta2));

% Accelerations
alpha1 = zeros(size(t));
alpha2 = zeros(size(t));
for i = 1:length(t)
    dydt_i = equations([], y(i,:)', m1, m2, L1, L2, g);
    alpha1(i) = dydt_i(2);
    alpha2(i) = dydt_i(4);
end
a1 = L1 * alpha1;
a2 = sqrt( ...
    (L1*alpha1).^2 + (L2*alpha2).^2 + ...
    2*L1*L2.*alpha1.*alpha2.*cos(theta1 - theta2));

% Kinetic Energy
T = 0.5 * m1 .* v1.^2 + 0.5 * m2 .* v2.^2;

% Potential Energy
y1_pos = -L1 * cos(theta1);
y2_pos = y1_pos - L2 * cos(theta2);
V = m1 * g * y1_pos + m2 * g * y2_pos;

%ANIMATION
figure('Color', 'w');
for i = 1:length(t)
    clf;

    % Plot rods
    plot([0 x1(i)], [0 y1(i)], 'r', 'LineWidth', 2); hold on;
    plot([x1(i) x2(i)], [y1(i) y2(i)], 'b', 'LineWidth', 2);

    % Plot masses
    plot(x1(i), y1(i), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    plot(x2(i), y2(i), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');

    % Trail
    trail_duration = 0.5;
    trail_start_idx = find(t >= t(i) - trail_duration, 1);
    if ~isempty(trail_start_idx) && trail_start_idx < i
        trail_len = i - trail_start_idx + 1;
        cmap = jet(trail_len);
        for j = 1:trail_len-1
            idx = trail_start_idx + j - 1;
            plot([x2(idx) x2(idx+1)], [y2(idx) y2(idx+1)], '-', ...
                'Color', cmap(j,:), 'LineWidth', 1.5);
            plot([x1(idx) x1(idx+1)], [y1(idx) y1(idx+1)], '-', ...
                'Color', cmap(j,:), 'LineWidth', 1);
        end
    end

    axis equal;
    axis([-L1-L2 L1+L2 -L1-L2 L1+L2]);
    title(sprintf('Double Pendulum at t = %.2f s', t(i)), 'Color', 'k');
    xlabel('X'); ylabel('Y');
    drawnow;
end

%PLOTTING
figure;

subplot(3,2,1);
plot(t, T, 'r');
title('Total Kinetic Energy vs Time'); xlabel('Time (s)'); ylabel('Kinetic Energy (J)');

subplot(3,2,2);
plot(t, V, 'b');
title('Total Potential Energy vs Time'); xlabel('Time (s)'); ylabel('Potential Energy (J)');

subplot(3,2,3);
plot(t, v1, 'g');
title('Velocity of Bob 1 (v1) vs Time'); xlabel('Time (s)'); ylabel('v1 (m/s)');

subplot(3,2,4);
plot(t, v2, 'm');
title('Velocity of Bob 2 (v2) vs Time'); xlabel('Time (s)'); ylabel('v2 (m/s)');

subplot(3,2,5);
plot(t, a1, 'c');
title('Acceleration of Bob 1 (a1) vs Time'); xlabel('Time (s)'); ylabel('a1 (m/s^2)');

subplot(3,2,6);
plot(t, a2, 'k');
title('Acceleration of Bob 2 (a2) vs Time'); xlabel('Time (s)'); ylabel('a2 (m/s^2)');

end

%SYSTEM EQUATIONS
function dydt = equations(~, y, m1, m2, L1, L2, g)
theta1 = y(1); omega1 = y(2);
theta2 = y(3); omega2 = y(4);
delta = theta2 - theta1;

den1 = (m1 + m2) * L1 - m2 * L1 * cos(delta)^2;
den2 = (L2 / L1) * den1;

domega1 = (m2 * L1 * omega1^2 * sin(delta) * cos(delta) + ...
    m2 * g * sin(theta2) * cos(delta) + ...
    m2 * L2 * omega2^2 * sin(delta) - ...
    (m1 + m2) * g * sin(theta1)) / den1;

domega2 = (-m2 * L2 * omega2^2 * sin(delta) * cos(delta) + ...
    (m1 + m2) * g * sin(theta1) * cos(delta) - ...
    (m1 + m2) * L1 * omega1^2 * sin(delta) - ...
    (m1 + m2) * g * sin(theta2)) / den2;

dydt = [omega1; domega1; omega2; domega2];
end
