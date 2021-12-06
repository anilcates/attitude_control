%% Simulation of Magnetometer Measurements
%% 
% Created by Anýl Can Ateþ

%%

n = 53;

% Initial data of quaternions
q1 = 0.002 * n;
q2 = 0.001 * n;
q3 = 0.005 * n;
q4 = sqrt(1-q1^2-q2^2-q3^2);

% The initial data of the satellite’s angular velocities
wx0 = 0.0002 + 0.0001 * n;
wy0 = 0.0003 + 0.0001 * n;
wz0 = 0.0004 + 0.0001 * n;

% The initial moments of inertia of the satellite
Jx = 2.1 * 10^-3; % m^4
Jy = 2 * 10^-3; % m^4
Jz = 1.9 * 10^-3; % m^4

% The angular orbit velocity of satellite
w_orbit = 0.0011; % rad/s

% The disturbance torque acting on the satellite
N_T = 3.6 * 10^-10; % N.m

% The iteration number
N_i = 54000;

% The sample time
dt = 0.1 % seconds

i = 1; % initial index

%%

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


%%

%% CONTINUATION

% Inputs

o_mc = 0.008; 
o_m = 1.66 % muT, the standard deviation of each magnetometer error

%  the magnetometer bias vector components
bx = 3;
by = 5;
bz = 6;

% the magnetometer bias vector in terms of direction cosines components
bxc = 0.04;
byc = 0.06;
bzc = 0.08;

% identity matrix 3x3
id_m = eye(3);

% the zero mean Gaussian white noise with the characteristic of
ni = id_m*o_m^2;

t = 0; % starting time
i = 1; % initial index

% iterations
while i <= N_i
    % save values of indices for time and angular velocities
    time(i) = t;
    wx(i) = wx0;
    wy(i) = wy0;
    wz(i) = wz0;

    % Euler method
    wx_new = wx0 + (dt/Jx)*(Jy-Jz)*wz0*wy0 + (dt/Jx)*N_T;
    wy_new = wy0 + (dt/Jy)*(Jz-Jx)*wx0*wz0 + (dt/Jy)*N_T;
    wz_new = wz0 + (dt/Jz)*(Jx-Jy)*wx0*wy0 + (dt/Jz)*N_T;
    wx0 = wx_new;
    wy0 = wy_new;
    wz0 = wz_new;
    
    % save value of indices for quaternions
    q_1(i) = q1;
    q_2(i) = q2;
    q_3(i) = q3;
    q_4(i) = q4;
    
    q1_new = q1 - 0.5*dt*(q2*wx0 + q3*wy0 + q4*wz0);
    q2_new = q2 + 0.5*dt*(q1*wx0 - q4*wy0 + q3*wz0);
    q3_new = q3 + 0.5*dt*(q4*wx0 + q1*wy0 - q2*wz0);
    q4_new = q4 - 0.5*dt*(q3*wx0 - q2*wy0 - q1*wz0);
    q1 = q1_new;
    q2 = q2_new;
    q3 = q3_new;
    q4 = q4_new;
       
    
    % Transformation matrix
    C_T(i).a = [q1^2-q2^2-q3^2+q4^2 2*(q1*q2+q3*q4) 2*(q1*q3-q2*q4)
        2*(q1*q2-q3*q4) -q1^2+q2^2-q3^2+q4^2 2*(q2*q3-q1*q4)
        2*(q1*q3+q2*q4) 2*(q2*q3-q1*q4) -q1^2-q2^2+q3^2+q4^2];
    
    i = i + 1;
    t = t + dt;
end

j = 1;

% iterations
while j <= N_i
    
    H1(j) = (Me/r0^3) * (cos(w0*t)*(cos(e)*sin(inc)-sin(e)*cos(inc)*cos(we*t)) ...
    - sin(w0*t)*sin(e)*sin(we*t));
    
    H2(j) = (-Me/r0^3) * (cos(e)*cos(inc)+sin(e)*sin(inc)*cos(we*t));
    
    H3(j) = (2*Me/r0^3) * (sin(w0*t)*(cos(e)*sin(inc)-sin(e)*cos(inc)*cos(we*t)) ...
    - 2*sin(w0*t)*sin(e)*sin(we*t));  
    
    % Direction cosine elements of the magnetic field vector
    Hx0(j) = (1./sqrt(H1(j).^2+H2(j).^2+H3(j).^2)).*H1(j);
    Hy0(j) = (1./sqrt(H1(j).^2+H2(j).^2+H3(j).^2)).*H2(j);
    Hz0(j) = (1./sqrt(H1(j).^2+H2(j).^2+H3(j).^2)).*H3(j);
    
    % The elements of the Earth's magnetic field vector
    H(j).a = [H1(j);H2(j);H3(j)]; 
    % Direction cosine elements of the magnetic field vector
    H0(j).a = [Hx0(j); Hy0(j); Hz0(j)];
    
    % Measuring magnetic field vector (body frame)
    Bm = C_T(j).a * H(j).a * 10^6 + [bx;by;bz] +o_m*randn; % (muT)
    % True value of magnetic Field vector (body frame)
    Bms = C_T(j).a * H(j).a * 10^6; % (muT)
    
    % Measuring direction cosine of magnetic field vector (body frame)
    Bmc = C_T(j).a * H0(j).a + [bxc;byc;bzc] +o_mc*randn;
    % True value of direction cosine of magnetic field vector (body frame)
    Bmcs = C_T(j).a * H0(j).a;
    
    % x,y,z components of Bm
    Bmx(j) = Bm(1);
    Bmy(j) = Bm(2);
    Bmz(j) = Bm(3);
    
    % x,y,z components of Bm*
    Bmsx(j) = Bms(1); 
    Bmsy(j) = Bms(2); 
    Bmsz(j) = Bms(3);
    
    % x,y,z components of Bmc
    Bmcx(j) = Bmc(1);
    Bmcy(j) = Bmc(2);
    Bmcz(j) = Bmc(3);
    
    % x,y,z components of Bmc*
    Bmcsx(j) = Bmcs(1);
    Bmcsy(j) = Bmcs(2);
    Bmcsz(j) = Bmcs(3);
    
    j = j + 1;
    t = t + dt;
end

%% Plot section

subplot(3,1,1)
plot(time, Bmx)
grid on
xlabel('Time (s)')
ylabel('$B_{m_{x}}$ ($\mu T$)','Interpreter', 'latex ')
xlim([0 5400])
title('Time versus x component of B_m')

subplot(3,1,2)
plot(time, Bmy)
grid on
xlabel('Time (s)')
ylabel('$B_{m_{y}}$ ($\mu T$)','Interpreter', 'latex ')
xlim([0 5400])
title('Time versus y component of B_m')

subplot(3,1,3)
plot(time, Bmz)
grid on
xlabel('Time (s)')
ylabel('$B_{m_{z}}$ ($\mu T$)','Interpreter', 'latex ')
xlim([0 5400])
title('Time versus z component of B_m')

sgtitle('Time versus Measured Magnetic Field Vector, B_m')

figure;

subplot(3,1,1)
plot(time, Bmsx)
grid on
xlabel('Time (s)')
ylabel('$B_{m_{x}}^*$ ($\mu T$)','Interpreter', 'latex ')
xlim([0 5400])
title('Time versus x component of B_m^*')

subplot(3,1,2)
plot(time, Bmsy)
grid on
xlabel('Time (s)')
ylabel('$B_{m_{y}}^*$ ($\mu T$)','Interpreter', 'latex ')
xlim([0 5400])
title('Time versus y component of B_m^*')

subplot(3,1,3)
plot(time, Bmsz)
grid on
xlabel('Time (s)')
ylabel('$B_{m_{z}}^*$ ($\mu T$)','Interpreter', 'latex ')
xlim([0 5400])
title('Time versus z component of B_m^*')

sgtitle('Time versus True Magnetic Field Vector, B_{m}^*')

figure;

subplot(3,1,1)
plot(time, Bmcx)
grid on
xlabel('Time (s)')
ylabel('$B_{mc_{x}}$','Interpreter', 'latex ')
xlim([0 5400])
title('Time versus x component of B_{mc}')

subplot(3,1,2)
plot(time, Bmcy)
grid on
xlabel('Time (s)')
ylabel('$B_{mc_{y}}$','Interpreter', 'latex ')
xlim([0 5400])
title('Time versus y component of B_{mc}')

subplot(3,1,3)
plot(time, Bmcz)
grid on
xlabel('Time (s)')
ylabel('$B_{mc_{z}}$','Interpreter', 'latex ')
xlim([0 5400])
title('Time versus z component of B_{mc}')

sgtitle('Time versus Measured Direction Cosine of Magnetic Field Vector, B_{mc}')

figure;

subplot(3,1,1)
plot(time, Bmcsx)
grid on
xlabel('Time (s)')
ylabel('$B_{mc_{x}}^*$','Interpreter', 'latex ')
xlim([0 5400])
title('Time versus x component of B_{mc}^*')

subplot(3,1,2)
plot(time, Bmcsy)
grid on
xlabel('Time (s)')
ylabel('$B_{mc_{y}}^*$','Interpreter', 'latex ')
xlim([0 5400])
title('Time versus y component of B_{mc}^*')

subplot(3,1,3)
plot(time, Bmcsz)
grid on
xlabel('Time (s)')
ylabel('$B_{mc_{z}}^*$','Interpreter', 'latex ')
xlim([0 5400])
title('Time versus z component of B_{mc}^*')

sgtitle('Time versus True Direction Cosine of Magnetic Field Vector, B_{mc}^*')


figure;

subplot(3,1,1)
plot(time,Bmx)
grid on
hold on
plot(time,Bmsx,'LineWidth',1.3)
hold off
xlabel('Time (s)')
ylabel('($\mu T$)','Interpreter','latex')
xlim([0 5400])
title('Time versus x component of Magnetic Field Vector')
legend('Measured','True')

subplot(3,1,2)
plot(time,Bmy)
grid on
hold on
plot(time,Bmsy,'LineWidth',1.3)
hold off
xlabel('Time (s)')
ylabel('($\mu T$)','Interpreter','latex')
xlim([0 5400])
title('Time versus y component of Magnetic Field Vector')
legend('Measured','True')

subplot(3,1,3)
plot(time,Bmz)
grid on
hold on
plot(time,Bmsz,'LineWidth',1.3)
hold off
xlabel('Time (s)')
ylabel('($\mu T$)','Interpreter','latex')
xlim([0 5400])
title('Time versus z component of Magnetic Field Vector')
legend('Measured','True')

sgtitle('Time versus Magnetic Field Vector of the Earth in Body Frame')


figure;

subplot(3,1,1)
plot(time,Bmcx)
grid on
hold on
plot(time,Bmcsx)
hold off
xlabel('Time (s)')
xlim([0 5400])
title('Time versus x component of Direction Cosine of Magnetic Field Vector')
legend('Measured','True')

subplot(3,1,2)
plot(time,Bmcy)
grid on
hold on
plot(time,Bmcsy)
hold off
xlabel('Time (s)')
xlim([0 5400])
title('Time versus y component of Direction Cosine of Magnetic Field Vector')
legend('Measured','True')

subplot(3,1,3)
plot(time,Bmcz)
grid on
hold on
plot(time,Bmcsz)
hold off
xlabel('Time (s)')
xlim([0 5400])
title('Time versus z component of Direction Cosine of Magnetic Field Vector')
legend('Measured','True')

sgtitle('Time versus Direction Cosine Elements of the Magnetic Field Vector of the Earth in Body Frame')
