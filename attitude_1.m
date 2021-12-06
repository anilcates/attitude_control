%% Simulation of The Low Earth Orbit Satellite’s Attitude Dynamics
%%
% Created by Anýl Can Ateþ

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

t_num = 0; % starting time
i = 1; % initial index

% The mathematical model of the satellite’s rotational motion 
% about its center of mass is given as:
while i <= N_i
    % save values of indices for time and angular velocities
    time(i) = t_num;
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
    
    % Euler method
    q1_new = q1 - 0.5*dt*(q2*wx0 + q3*wy0 + q4*wz0);
    q2_new = q2 + 0.5*dt*(q1*wx0 - q4*wy0 + q3*wz0);
    q3_new = q3 + 0.5*dt*(q4*wx0 + q1*wy0 - q2*wz0);
    q4_new = q4 - 0.5*dt*(q3*wx0 - q2*wy0 - q1*wz0);
    q1 = q1_new;
    q2 = q2_new;
    q3 = q3_new;
    q4 = q4_new;
    
    % Propagation of quaternion with time
    % Quaternion rates: qdot = [q1dot; q2dot; q3dot; q4dot]
    qdot = 0.5*[q1 -q2 -q3 -q4; q2 q1 -q4 q3; q3 q4 q1 -q2; q4 -q3 q2 q1]...
    *[0; wx0; wy0; wz0];
    
    % separate the qdot vector
    q1dot = qdot(1,:);
    q2dot = qdot(2,:);
    q3dot = qdot(3,:);
    q4dot = qdot(4,:);
    
    % save value of indices for quaternion rates
    q1_rate(i) = q1dot;
    q2_rate(i) = q2dot;
    q3_rate(i) = q3dot;
    q4_rate(i) = q4dot;
    
    % Transformation matrix 
    C_T(i).a = [q1^2-q2^2-q3^2+q4^2 2*(q1*q2+q3*q4) 2*(q1*q3-q2*q4)
        2*(q1*q2-q3*q4) -q1^2+q2^2-q3^2+q4^2 2*(q2*q3-q1*q4)
        2*(q1*q3+q2*q4) 2*(q2*q3-q1*q4) -q1^2-q2^2+q3^2+q4^2];
    
    % save value of indices for tranformation matrix 
   % C{i} = C_T;
 
  
    i = i + 1;
    t_num = t_num + dt;
end

% quaternion plots
subplot(2,2,1)
plot(time,q_1)
xlabel('Time(s)')
ylabel('q_1')

subplot(2,2,2)
plot(time,q_2)
xlabel('Time(s)')
ylabel('q_2')

subplot(2,2,3)
plot(time,q_3)
xlabel('Time(s)')
ylabel('q_3')

subplot(2,2,4)
plot(time,q_4)
xlabel('Time(s)')
ylabel('q_4')

sgtitle('Quaternions') 

figure;
% angular velocity plots
subplot(2,2,1)
plot(time,wx)
xlabel('Time(s)')
ylabel('\omega_x(rad/s)')

subplot(2,2,2)
plot(time,wy)
xlabel('Time(s)')
ylabel('\omega_y(rad/s)')

subplot(2,2,3)
plot(time,wz)
xlabel('Time(s)')
ylabel('\omega_z(rad/s)')

sgtitle('Angular Velocities')

figure;
% quaternion rate plots
subplot(2,2,1)
plot(time,q1_rate)
xlabel('Time(s)')
ylabel('$\dot{q_1}$', 'Interpreter','latex')

subplot(2,2,2)
plot(time,q2_rate)
xlabel('Time(s)')
ylabel('$\dot{q_2}$', 'Interpreter','latex')

subplot(2,2,3)
plot(time,q3_rate)
xlabel('Time(s)')
ylabel('$\dot{q_3}$', 'Interpreter','latex')

subplot(2,2,4)
plot(time,q4_rate)
xlabel('Time(s)')
ylabel('$\dot{q_4}$', 'Interpreter','latex')

sgtitle('Quaternion Rates, $\dot{q}$', 'Interpreter','latex')