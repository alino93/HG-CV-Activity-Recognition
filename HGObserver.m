function [x_hat, L] = HGObserver(a_body, W, pos, th) % acc_body (3xN [m/s^2]), omega_body(3xN [rad/s]), position measurement(3xN [m]) in global coordiante, theta for HG
%% Observer gain
la= 3.33;
%th = 30;

A = zeros(12);
B = zeros(12,6);
C = zeros(6,12);
temp = [1,3,5,7,9,11; 2,4,6,8,10,12; 1,2,3,4,5,6]';
for i=1:1:length(temp)
    A(temp(i,1),temp(i,2)) = 1;
    B(temp(i,2),temp(i,3)) = 1;
    C(temp(i,3),temp(i,1)) = 1;
end

I = eye(12);
nsys = 12;
nout = 6;
% Define the decision variables P and R
P = sdpvar(nsys, nsys);
R = sdpvar(nsys, nout);

% Set LMI constraints
lmi = A'*P + P*A - C'*R'- R * C + la * I;

% Construct LMIs for solver
F = [P >= 0];
F = F + [lmi <= 0];

% Choose solver
ops = sdpsettings('solver','sedumi'); % need to install
ops.verbose = 0;
% Solve the LMI
diagnostics = optimize(F, [ ], ops);
% optimize(Constraints, Objective, options)

%Check the results of optimization
if diagnostics.problem == 0
    disp('Feasible from Solver')
elseif diagnostics.problem == 1
    disp('lnfeasible from Solver')
else
    disp('Something else happened')
end
    
% Display the results
P = double(P);
R = double(R);
K = P\R;
L = diag([th,th^2,th,th^2,th,th^2,th,th^2,th,th^2,th,th^2]) * K;

%% Nonlinear observer quat estimation
T = 0.01;
g = 9.81;
ge = norm(a_body(:,1));
a_n = a_body(:,1)/ge ; % normalized body acc

%time
t = 0:T:T*(length(W(1,:))-1);

% high pass a_body to keep only linear acc
a_lin = highpass(a_body,1,100);
%% Observer
% initialize
x_hat = zeros(9,length(t));
z_hat = zeros(nsys,length(t));

x_hat(:,1) = [0;0;0;0;0;0;atan2(a_n(2),a_n(3));-asin(a_n(1));0.01];

x_h = x_hat(7:9,1);
x_dot = [1 sin(x_h(1))*tan(x_h(2)) cos((x_h(1)))*tan(x_h(2));
         0 cos(x_h(1)) -sin(x_h(1));
         0 sin(x_h(1))/cos(x_h(2)) cos(x_h(1))/cos(x_h(2))] * W(:,1);

z_hat(:,1) = [0;0;0;0;0;0;
    x_h(1);x_dot(1);x_h(2);x_dot(2);x_h(3);x_dot(3)];

% global linear acc

% iterate
for i=2:1:length(t)
    ge = norm(a_body(:,i-1));
    a_n = a_body(:,i-1)/ge;
    
    w = W(:,i-1);
    x_h = x_hat(7:9,i-1);
    a_global = eul2rotm([x_h(3), x_h(2), x_h(1)], 'ZYX') * a_lin(:,i-1);
    x_dot = [1 sin(x_h(1))*tan(x_h(2)) cos((x_h(1)))*tan(x_h(2));
             0 cos(x_h(1)) -sin(x_h(1));
             0 sin(x_h(1))/cos(x_h(2)) cos(x_h(1))/cos(x_h(2))] * w;
    f_4 = x_dot(1) * tan(x_h(2)) * (w(2)*cos(x_h(1)) - w(3)*sin(x_h(1))) + ...
        x_dot(2)/cos(x_h(2))^2 * (w(3) * cos(x_h(1)) + w(2) * sin(x_h(1))); 
    f_5 = -x_dot(1) * (w(3)*cos(x_h(1)) + w(2)*sin(x_h(1)));
    f_6 = x_dot(1) * (w(2)*cos(x_h(1)) - w(3)*sin(x_h(1))) / cos(x_h(2)) + ...
        x_dot(2) * tan(x_h(2)) * (w(3) * cos(x_h(1)) + w(2) * sin(x_h(1)))/cos(x_h(2));
   
    f_hat = [a_global;f_4;f_5;f_6];
    
    % Measurement eq: pos_x;pos_y;pos_z;roll;pich;yaw
    % we can neglect yaw since dont have measuremnt for it
    y = [pos(:,i-1);1*atan2(a_n(2),a_n(3));-asin(a_n(1));z_hat(5,i-1)];
    
    z_hat(:,i) = z_hat(:,i-1) + T * (A * z_hat(:,i-1) + ...
    B * f_hat + L * (y - C * z_hat(:,i-1)));
    
    % suppress results outside -pi,pi
    if z_hat(7,i)> pi
        z_hat(7,i) = pi;
    elseif z_hat(7,i)<-pi
        z_hat(7,i) = -pi;
    end
    x_hat(1,i) = z_hat(1,i);
    x_hat(2,i) = z_hat(3,i);
    x_hat(3,i) = z_hat(5,i);
    x_hat(4,i) = z_hat(2,i);
    x_hat(5,i) = z_hat(4,i);
    x_hat(6,i) = z_hat(6,i);
    x_hat(7,i) = z_hat(7,i);
    x_hat(8,i) = z_hat(9,i);
    % we can neglect yaw since dont have measuremnt for it
    % x_hat(9,i) = z_hat(11,i);
    x_hat(9,i) = x_hat(9,i-1);

%     if x_dot(3)<0.4
%         x_dot(3) = 0; 
%     end
    % x_hat(3,i) = x_hat(3,i-1) + T * (x_dot(3));
       
end
end
