%% Earth’s Magnetic Field Vector through the Orbit of the Satellite Using Dipole Model
%%
% Created by Anýl Can Ateþ

%% Inputs

Me = 7.943 * 10^15; % Wb.m, the magnetic dipole moment of the Earth
mu = 3.98601 * 10^14; % m^3/s^2 the Earth Gravitational constant

n = 53; % The student number in class list is n

inc = deg2rad(80 + (0.5*n)); % deg to rad, the orbit inclination

we = 7.29 * 10^-5; % rad/s, the spin rate of the Earth

% Inclination of Geomagnetic axis relative to Earth’s daily rotating axis
e = deg2rad(11.7); % the magnetic dipole tilt

% the distance between the center of mass of the satellite and the Earth
r0 = (6378.140 + 500 + 2*n) * 1000; % m

% the angular velocity of the orbit with respect to the inertial frame
w0 = (mu/r0^3)^(1/2);

% The iteration number
N_i = 54000;

% The sample time
dt = 0.1; % seconds

t = 0; % starting time
j = 1; % initial index

%% Loop

while j <= N_i
    time(j) = t;
    
    H1(j) = (Me/r0^3) * (cos(w0*t)*(cos(e)*sin(inc)-sin(e)*cos(inc)*cos(we*t)) ...
    - sin(w0*t)*sin(e)*sin(we*t));
    
    H2(j) = (-Me/r0^3) * (cos(e)*cos(inc)+sin(e)*sin(inc)*cos(we*t));
    
    H3(j) = (2*Me/r0^3) * (sin(w0*t)*(cos(e)*sin(inc)-sin(e)*cos(inc)*cos(we*t)) ...
    - 2*sin(w0*t)*sin(e)*sin(we*t));  

    j = j + 1;
    t = t + dt;
end

% Direction cosine elements of the magnetic field vector
Hx0 = (1./sqrt(H1.^2+H2.^2+H3.^2)).*H1;
Hy0 = (1./sqrt(H1.^2+H2.^2+H3.^2)).*H2;
Hz0 = (1./sqrt(H1.^2+H2.^2+H3.^2)).*H3;

%% Plotting

% Plot Hx, Hy and Hz with time
subplot(3,1,1)
plot(time,H1)
grid on 
xlim([0 5400]);
xlabel('Time(s)')
ylabel('H_x (Tesla)')
title('Time versus H_x')

subplot(3,1,2)
plot(time,H2)
grid on
xlim([0 5400]);
xlabel('Time (s)')
ylabel('H_y (Tesla)')
title('Time versus H_y')

subplot(3,1,3)
plot(time,H3)
grid on
xlim([0 5400]);
xlabel('Time (s)')
ylabel('H_z (Tesla)')
title('Time versus H_z')

sgtitle('The elements of the Earth’s Magnetic Field vector H')

% Direction cosine elements of the magnetic field vector, H0
% Plot Hx0, Hy0 and Hz0 for direction cosines

figure;
subplot(3,1,1)
plot(time,Hx0)
grid on
xlim([0 5400]);
xlabel('Time(s)')
ylabel('H_{x0}')

subplot(3,1,2)
plot(time,Hy0)
grid on
xlim([0 5400]);
xlabel('Time(s)')
ylabel('H_{y0}')

subplot(3,1,3)
plot(time,Hz0)
grid on
xlim([0 5400]);
xlabel('Time(s)')
ylabel('H_{z0}')

sgtitle('Direction cosine elements of the magnetic field vector')