%% ------------------------------------------------------------------------------------ %%
%   author: Giorgio Simonini
%   date:   9/12/2022
%   info: This file created an elastic model (SEA) for a given robot
%   convention:
%       - Links start from 0, first bring the base frame to the first
%           joint, last bring the last joint to the end effector
%       - Frames start from 0, there are as many frames as links
%       - Frames have the origins on the phisical joints (0 and ee
%           excluded)
%       - Joints start from 1, frames and joints have same origins and
%           indexes
%       - Link Li brings the Joint {Ji} to the Joint {Ji+1} with a 
%           rototraslational matrix, defined by a traslation and a 
%           rotation using the XYZ parametrization
%       - Center of mass of link Li is expressed in frame {Si}
%       - Inertia tensor of link Li is defined on the center of mass in the
%           frame {Si}
%       - the joint variables are included in the trasformation matrices
% ---------------------------------------------------------------------------------------%
% TODO: 
% - default configuration
% - add FF_SEA
% - regressor kernel is experimental
% - working numeric regressor kernel
% - do not recalculate matrix on linear parameters
% - generalize on parameters choice

function out = create_model_symb(conf)

%% ----- INIT ----- %%
mkdir functions
N_joints = conf.N_joints;        % number of joints
Q_MIN = conf.Q_MIN;
Q_MAX = conf.Q_MAX;
par_nl = conf.par_nl;		% nonlinear parameters, non used in adaptive
CREATE_FUNCTIONS = conf.CREATE_FUNCTIONS;
SAVE_MODEL = conf.SAVE_MODEL;
CREATE_REGRESSOR_KERNEL = 0;% conf.CREATE_REGRESSOR_KERNEL;
REGRESSOR_KERNEL_SYMB = conf.REGRESSOR_KERNEL_SYMB;
% - hyperparameters - %
N_frames = N_joints + 1;    % joints + base + ee
joints_type = conf.joints_type;  % rotoidal:'R' or prismatic:'P'
% - known parameters - %
% syms tau a1 a2 g real
g = sym('g', 'real');
par_nl = [par_nl, g];		% nonlinear parameters, non used in adaptive
par_nl = reshape(par_nl, [numel(par_nl),1]);		% reshape in 1-D vector
% - Frames - %
frames_tr = conf.frames_tr;
frames_or = conf.frames_or;

%% VARIABLES AND PARAMETERS
disp('Variables and parameters initialization ...')
tic
% - variables - %
q = sym('q%d', [N_joints,1], 'real');
dq = sym('dq%d', [N_joints,1], 'real');
ddq = sym('ddq%d', [N_joints,1], 'real');
d3q = sym('d3q', [N_joints,1], 'real');
d4q = sym('d4q', [N_joints,1], 'real');
x = sym('x%d', [N_joints,1], 'real');
dx = sym('dx%d', [N_joints,1], 'real');
ddx = sym('ddx%d', [N_joints,1], 'real');
% - parameters - %
masses = sym('m%d', [N_joints,1], 'real');                  % Masses
center_of_masses = sym('l%d_%d', [N_joints,3], 'real');     % Center of Masses
link_inertia = sym('i%d_%d', [N_joints,6], 'real');         % Inertia Tensors
link_damp = sym('d%d', [N_joints,1], 'real');               % Link dampings
stiffness = sym('k%d', [N_joints,1], 'real');               % Coupling Stiffnesses
motor_inertia = sym('bm%d', [N_joints,1], 'real');          % Motor Inertia
motor_damp = sym('dm%d', [N_joints,1], 'real');             % Motor Damping
par = [masses, center_of_masses, link_inertia, link_damp, stiffness, motor_inertia, motor_damp];  % parameters
par = reshape(par,[N_joints*size(par,2),1]);
p = sym('p%d_%d', [N_joints, 14], 'real');
disp('done!')
toc

%% FRAMES CHAIN
disp('Building frames chain ...')
tic
T_0_i = sym(zeros(4,4,N_joints+1));
for index = 1 : N_joints+1    % joints + end-effector
    if index <= N_joints
		if joints_type{index} == 'P'
			frames_tr(index, 3) = frames_tr(index, 3) + q(index);
		elseif joints_type{index} == 'R'
			frames_or(index, 3) = frames_or(index, 3) + q(index);
		end
        R_tmp = R_x(frames_or(index,1)) * R_y(frames_or(index,2)) * R_z(frames_or(index,3));
        d_tmp = frames_tr(index, :)';
        if index == 1
            T_0_i(:,:,index) = [R_tmp, d_tmp; 0,0,0, 1];
        else
            T_0_i(:,:,index) = T_0_i(:,:,index-1)*[R_tmp, d_tmp; 0,0,0, 1];
        end
    else
        R_tmp = R_x(frames_or(index,1)) * R_y(frames_or(index,2)) * R_z(frames_or(index,3));
        d_tmp = frames_tr(index, :)';
        T_0_i(:,:,index) = T_0_i(:,:,index-1)*[R_tmp, d_tmp; 0,0,0, 1];
    end
end
T_0_i = combine(T_0_i,'sincos');
T_0_i = simplify(T_0_i);
disp('done!')
toc
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(T_0_i,'file','functions\get_T_0_i', 'Vars', {q, par_nl});
	disp('done!')
    toc
end

%% JACOBIANS
disp('Computing Jacobians ...')
tic
% The position jacobian is the differentiation of the position
J_i = sym(zeros(6, N_joints, size(T_0_i, 3)));
J_i_or = sym(zeros(3,N_joints));
for index = 1 : size(T_0_i, 3)
    d_0_i = T_0_i(1:3, 4, index);
    J_i_pos = simplify(jacobian(d_0_i, q));
    J_i_pos = combine(J_i_pos,'sincos');
    J_i_pos = simplify(J_i_pos);
    for index2 = 1 : index
		if index <= N_joints
        	J_i_or(:,index2) = T_0_i(1:3, 3, index2);
		end
    end
    J_i(:,:,index) = [J_i_pos; J_i_or];
end
J_ee = J_i(:,:,end);
disp('done!')
toc
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(J_ee,'file','functions\get_J_ee', 'Vars', {q, par, par_nl});
    matlabFunction(J_i,'file','functions\get_J_i', 'Vars', {q, par, par_nl});
    disp('done!')
    toc
end
%% INERTIA MATRIX
disp('Computing Inertia Matrix ...')
tic
% - center of mass jacobians - %
J_cm_i = sym(zeros(6, N_joints, N_joints));
J_i_or = sym(zeros(3,N_joints));
for index = 1 : N_joints
    d_0_i = T_0_i(1:3,4,index) + T_0_i(1:3,1:3,index)*(center_of_masses(index,:)');
	d_0_i = simplify(d_0_i);
    J_i_pos = simplify(jacobian(d_0_i, q));
    J_i_pos = combine(J_i_pos,'sincos');
    J_i_pos = simplify(J_i_pos);
	for index2 = 1 : index
		if index <= N_joints
        	J_i_or(:,index2) = T_0_i(1:3, 3, index2);
		end
    end
    J_cm_i(:,:,index) = [J_i_pos; J_i_or];
end
if CREATE_FUNCTIONS
    disp(' Creating Jacobian function ...')
    matlabFunction(J_cm_i,'file','functions\get_J_cm_i', 'Vars', {q, par, par_nl});
    disp('done!')
    toc
end
% - inertia matrix - %
B = sym(zeros(N_joints));
for index = 1 : N_joints
    I_tmp = link_inertia(index,:);
    I_cm_i = [I_tmp(1), I_tmp(4), I_tmp(5); I_tmp(4), I_tmp(2), I_tmp(6); I_tmp(5), I_tmp(6), I_tmp(3)];
    I_cm_i = simplify(T_0_i(1:3,1:3,index) * I_cm_i * (T_0_i(1:3,1:3,index)'));
    I_cm_i = combine(I_cm_i,'sincos');
    I_cm_i = simplify(I_cm_i);
    M_cm_i = [masses(index)*eye(3), zeros(3); zeros(3), I_cm_i];
    B = B + J_cm_i(:,:,index)' * M_cm_i * J_cm_i(:,:,index);
end
toc
disp('simplifying ...')
B = simplify(B);
B = combine(B,'sincos');
B = simplify(B);
disp('done!')
toc
if CREATE_FUNCTIONS 
    disp(' Creating Inertia function ...')
    matlabFunction(B,'file','functions\get_B', 'Vars', {q, par, par_nl});
    disp('done!')
    toc
end

%% CORIOLIS MATRIX
disp('Computing Coriolis Matrix ...')
tic
C = sym(zeros(N_joints));
for index = 1 : 1 : N_joints
    for index2 = 1 : 1 : N_joints
        gamma_sum = sym(0);
        for index3 = 1 : 1 : N_joints
            gamma = 1/2 * (diff(B(index,index2), q(index3)) + diff(B(index,index3), q(index2)) - diff(B(index2,index3), q(index)));
            gamma = simplify(gamma);
            gamma_sum = ( gamma_sum + gamma ) * dq(index3);
        end
        clear gamma;
        C(index,index2) = gamma_sum;
    end
end
toc
disp('simplifying ...')
C = simplify(C);
C = combine(C,'sincos');
C = simplify(C);
disp('done!')
toc
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(C,'file','functions\get_C', 'Vars', {q, dq, par, par_nl});
    disp('done!')
    toc
end

%% LINK DAMPING MATRIX
disp('Computing Link damping Matrix ...')
D = diag(link_damp);
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(D,'file','functions\get_D', 'Vars', {par, par_nl});
end
disp('done!')

%% GRAVITY TERMS
disp('Computing Gravity terms ...')
tic
G = sym(zeros(N_joints,1));
grav = [0 0 g]';
for index = 1 : N_joints
    G = G - J_cm_i(:,:,index)'*[grav*masses(index);zeros(3,1)];
end
G = simplify(G);
G = combine(G,'sincos');
G = simplify(G);
disp('done!')
toc
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(G,'file','functions\get_G', 'Vars', {q, par, par_nl});
    disp('done!')
    toc
end

%% STIFFNESS MATRIX
disp('Computing Stiffness Matrix ...')
K = diag(stiffness);
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(K,'file','functions\get_K', 'Vars', {par, par_nl});
end
disp('done!')

%% MOTOR INERTIA MATRIX
disp('Computing Motor Inertia Matrix ...')
Bm = diag(motor_inertia);
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(Bm,'file','functions\get_Bm', 'Vars', {par, par_nl});
end
disp('done!')

%% MOTOR DAMPING MATRIX
disp('Computing Motor Damping Matrix ...')
Dm = diag(motor_damp);
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(Dm,'file','functions\get_Dm', 'Vars', {par, par_nl});
end
disp('done!')

%% MODEL LINEAR IN PARAMETERS
disp('Creating linear parameters ...')
% p1 = m;
% p(2:4) = m*l;
% p(5:7) = I + m*l^2;
% p(8:10) = ixy, iyz, ixz;
% p11 = d;
% p12 = k;
% p13 = bm;
% p14 = dm;
% % ----- INVERSE -----
% m = p1;
% l = p(2:4)/p1;
% I = p(5:7) - (p(2:4)^2)/p1;
masses_lin = p(:, 1);                  % Masses
center_of_masses_lin = p(:,2:4)./p(:,1); 	% Center of Masses
link_inertia_lin = sym(zeros(N_joints,6));	% Link inertia
link_inertia_lin(:,1) = p(:,5) - (p(:,3).^2)./p(:,1) - (p(:,4).^2)./p(:,1);
link_inertia_lin(:,2) = p(:,6) - (p(:,2).^2)./p(:,1) - (p(:,4).^2)./p(:,1);
link_inertia_lin(:,3) = p(:,7) - (p(:,2).^2)./p(:,1) - (p(:,3).^2)./p(:,1);
link_inertia_lin(:,4) = p(:,8) + (p(:,2).*p(:,3))./p(:,1);
link_inertia_lin(:,5) = p(:,9) + (p(:,2).*p(:,4))./p(:,1);
link_inertia_lin(:,6) = p(:,10) + (p(:,3).*p(:,4))./p(:,1);
link_damp_lin = p(:,11);               % Link dampings
stiffness_lin = p(:,12);               % Coupling Stiffnesses
motor_inertia_lin = p(:,13);           % Motor Inertia
motor_damp_lin = p(:,14);              % Motor Damping
p = reshape(p,[N_joints*14,1]);
disp('done!')

%% LINEAR INERTIA MATRIX
disp('computing linear Inertia Matrix ...')
tic
% - center of mass jacobians - %
J_cm_i_p = sym(zeros(6, N_joints, N_joints));
J_i_or = sym(zeros(3,N_joints));
for index = 1 : N_joints
    d_0_i = T_0_i(1:3,4,index) + T_0_i(1:3,1:3,index)*(center_of_masses_lin(index,:)');
	d_0_i = simplify(d_0_i);
    d_0_i = combine(d_0_i,'sincos');
    d_0_i = simplify(d_0_i);
    J_i_pos = simplify(jacobian(d_0_i, q));
	for index2 = 1 : index
		if index <= N_joints
        	J_i_or(:,index2) = T_0_i(1:3, 3, index2);
        end
    end
    J_i_or = simplify(J_i_or);
    J_cm_i_p(:,:,index) = [J_i_pos; J_i_or];
end
if CREATE_FUNCTIONS 
    disp(' Creating Jacobian function ...')
    matlabFunction(J_cm_i_p,'file','functions\get_J_cm_i_p', 'Vars', {q, p, par_nl});
    disp('done!')
    toc
end
Bp = sym(zeros(N_joints));
for index = 1 : N_joints
    I_tmp = link_inertia_lin(index,:);
    I_cm_i = [I_tmp(1), I_tmp(4), I_tmp(5); I_tmp(4), I_tmp(2), I_tmp(6); I_tmp(5), I_tmp(6), I_tmp(3)];
    I_cm_i = simplify(T_0_i(1:3,1:3,index) * I_cm_i * (T_0_i(1:3,1:3,index)'));
    M_cm_i = simplify([masses_lin(index)*eye(3), zeros(3); zeros(3), I_cm_i]);
    Bp = Bp + J_cm_i_p(:,:,index)' * M_cm_i * J_cm_i_p(:,:,index);
end
toc
disp('simplifying ...')
Bp = simplify(Bp);
Bp = combine(Bp,'sincos');
Bp = simplify(Bp);
disp('done!')
toc
if CREATE_FUNCTIONS 
    disp(' Creating Inertia function ...')
    matlabFunction(Bp,'file','functions\get_Bp', 'Vars', {q, p, par_nl});
    disp('done!')
    toc
end

%% LINEAR CORIOLIS MATRIX
disp('Computing linear Coriolis Matrix ...')
tic
Cp = sym(zeros(N_joints));
for index = 1 : 1 : N_joints
    for index2 = 1 : 1 : N_joints
        gamma_sum = sym(0);
        for index3 = 1 : 1 : N_joints
            gamma = 1/2 * (diff(Bp(index,index2), q(index3)) + diff(Bp(index,index3), q(index2)) - diff(Bp(index2,index3), q(index)));
            gamma = simplify(gamma);
            gamma_sum = ( gamma_sum + gamma ) * dq(index3);
        end
        clear gamma;
        Cp(index,index2) = gamma_sum;
    end
end
toc
disp('simplifying ...')
Cp = simplify(Cp);
Cp = combine(Cp,'sincos');
Cp = simplify(Cp);
disp('done!')
toc
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(Cp,'file','functions\get_Cp', 'Vars', {q, dq, p, par_nl});
    disp('done!')
    toc
end

%% LINEAR LINK DAMPING MATRIX
disp('Computing linear Link damping Matrix ...')
Dp = diag(link_damp_lin);
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(Dp,'file','functions\get_Dp', 'Vars', {p, par_nl});
end
disp('done!')

%% LINEAR GRAVITY TERMS
disp('Computing linear Gravity terms ...')
tic
Gp = sym(zeros(N_joints,1));
grav = [0 0 g]';
for index = 1 : N_joints
    Gp = Gp - J_cm_i_p(:,:,index)'*[grav*masses_lin(index);zeros(3,1)];
end
Gp = simplify(Gp);
Gp = combine(Gp,'sincos');
Gp = simplify(Gp);
disp('done!')
toc
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(Gp,'file','functions\get_Gp', 'Vars', {q, p, par_nl});
    disp('done!')
    toc
end

%% LINEAR STIFFNESS MATRIX
disp('Computing linear Stiffness Matrix ...')
Kp = diag(stiffness_lin);
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(Kp,'file','functions\get_Kp', 'Vars', {p, par_nl});
end
disp('done!')

%% LINEAR MOTOR INERTIA MATRIX
disp('Computing linear Motor Inertia Matrix ...')
Bmp = diag(motor_inertia_lin);
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(Bmp,'file','functions\get_Bmp', 'Vars', {p, par_nl});
end
disp('done!')

%% LINEAR MOTOR DAMPING MATRIX
disp('Computing linear Motor Damping Matrix ...')
Dmp = diag(motor_damp_lin);
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(Dmp,'file','functions\get_Dmp', 'Vars', {p, par_nl});
end
disp('done!')

%% - Matrix derivatives - %%
% Inertia %
disp(' Inertia matrix derivatives ...')
tic
Bp_dot = Bp;
Bp_ddot = Bp;
for index = 1:size(Bp,2)
    Bp_dot(:,index) = jacobian(Bp(:,index),q)*dq;
    Bp_ddot(:,index) = jacobian(Bp_dot(:,index),q)*dq + jacobian(Bp_dot(:,index),dq)*ddq;
end
toc
disp('simplifying ...')
Bp_dot = simplify(Bp_dot);
Bp_ddot = simplify(Bp_ddot);
disp('done!')
toc
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(Bp_dot,'file','functions\get_Bp_dot', 'Vars', {q, dq, p, par_nl});
    matlabFunction(Bp_ddot,'file','functions\get_Bp_ddot', 'Vars', {q, dq, ddq, p, par_nl});
	disp('done!')
    toc
end

% Coriolis %
disp(' Coriolis matrix derivatives ...')
tic
Cp_dot = Cp;
Cp_ddot = Cp;
for index = 1:size(Cp,2)
    Cp_dot(:,index) = jacobian(Cp(:,index),q)*dq + jacobian(Cp_dot(:,index),dq)*ddq;
    Cp_ddot(:,index) = jacobian(Cp_dot(:,index),q)*dq + jacobian(Cp_dot(:,index),dq)*ddq + jacobian(Cp_dot(:,index),ddq)*d3q;
end
toc
disp('simplifying ...')
Cp_dot = simplify(Cp_dot);
Cp_ddot = simplify(Cp_ddot);
disp('done!')
toc
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(Cp_dot,'file','functions\get_Cp_dot', 'Vars', {q, dq, ddq, p, par_nl});
    matlabFunction(Cp_ddot,'file','functions\get_Cp_ddot', 'Vars', {q, dq, ddq, d3q, p, par_nl});
	disp('done!')
    toc
end

% Gravity %
disp(' Gravity derivatives ...')
tic
Gp_dot = jacobian(Gp,q)*dq;
Gp_ddot = jacobian(Gp_dot,q)*dq + jacobian(Gp_dot,dq)*ddq;
toc
disp('simplifying ...')
Gp_dot = simplify(Gp_dot);
Gp_ddot = simplify(Gp_ddot);
disp('done!')
toc
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(Gp_dot,'file','functions\get_Gp_dot', 'Vars', {q, dq, p, par_nl});
    matlabFunction(Gp_ddot,'file','functions\get_Gp_ddot', 'Vars', {q, dq, ddq, p, par_nl});
	disp('done!')
    toc
end
disp('done!')

%% MODELS AND REGRESSOR
disp('Computing regressor ...')
% - MODELs - %
% f1 = B*ddq + (C+D)*dq + G - K*(x-q);
% f2 = Bm*ddx + Dm*dx + K*(x-q);
% f = [f1; f2];
fp1 = Bp*ddq + (Cp+Dp)*dq + Gp - Kp*(x-q);
fp2 = Bmp*ddx + Dmp*dx + Kp*(x-q);
fp = [fp1; fp2];
% - Newton Laws - %
% tau_sea = K*(x-q); % + D_sea*(qm_dot - q_dot);
% F_l = K*(x-q) - (C+D)*dq - G;
% ddq_l = B\F_l;
% F_m = tau - Dm*dx - K*(x-q);
% ddq_m = Bm\F_m;
% - Regressor - %
Y = jacobian(fp, p);
disp('done!')
toc
disp('simplifying ...')
Y = simplify(Y);
Y = combine(Y,'sincos');
Y = simplify(Y);
Y = collect(expand(Y));     % p/p=1
disp('done!')
toc
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(Y,'file','functions\get_Y', 'Vars', {q, dq, ddq, x, dx, ddx, par_nl});
    disp('done!')
    toc
end

%% Build regressor kernel and projectors
addpath functions\
disp('Creating regressor kernel ...')
if CREATE_REGRESSOR_KERNEL == 1
    N_it = 100;
    
    if REGRESSOR_KERNEL_SYMB
        N_it = ceil(size(Y,2)/size(Y,1));
        q_ = sym('q_',[N_joints, N_it], 'real');
        dq_ = sym('dq_',[N_joints, N_it], 'real');
        ddq_ = sym('ddq_',[N_joints, N_it], 'real');
        x_ = sym('x_',[N_joints, N_it], 'real');
        dx_ = sym('dx_',[N_joints, N_it], 'real');
        ddx_ = sym('ddx_',[N_joints, N_it], 'real');
        Y = sym(zeros(1,14*N_joints));
    else
        % qmin = repmat(Q_MIN, 1, N_it);
        % qmax = repmat(Q_MAX, 1, N_it);
        % rng(1);
        % q_ = qmin + rand(N_joints, N_it).*(qmax-qmin);
        % dq_ = qmin + rand(N_joints, N_it).*(qmax-qmin);
        % ddq_ = qmin + rand(N_joints, N_it).*(qmax-qmin);
        % x_ = qmin + rand(N_joints, N_it).*(qmax-qmin);
        % dx_ = qmin + rand(N_joints, N_it).*(qmax-qmin);
        % ddx_ = qmin + rand(N_joints, N_it).*(qmax-qmin);
        % Y = zeros(1,14*N_joints);
    end
    
    for it = 1:N_it
        Y_tmp = get_Y(q_(:,it), dq_(:,it), ddq_(:,it), x_(:,it), dx_(:,it), ddx_(:,it), par_nl);
        Y = [Y; Y_tmp];
    end
    
    if REGRESSOR_KERNEL_SYMB
        n = null(Y);
        n = simplify(n);
        n = combine(n, 'sincos');
        n = simplify(n);
    else
        n = null(Y, 1e-6);
    end
    
    %% - normalized null - %%
    for index = 1:size(n,2)
        norm_tmp = norm(n(:,index));
        n(:,index) = n(:,index)/norm_tmp;
    end
    if REGRESSOR_KERNEL_SYMB
        n = simplify(n);
        n = combine(n, 'sincos');
        n = simplify(n);
        null_Y_sym = n;
    end
    null_Y = n;
    
    if CREATE_FUNCTIONS 
        disp(' Creating function ...')
        matlabFunction(null_Y_sym,'file','..\functions\get_null_Y_sym', 'Vars', {par_nl});
        disp('done!')
        toc
    end
end

%% Save model
if SAVE_MODEL
    disp(' Saving model ...')
    tic
    model.N_joints = N_joints;
    model.joints_type = joints_type;
    model.frames_tr = frames_tr;
    model.frames_or = frames_or;
    model.B = B;
    model.Bp = Bp;
    model.Bp_dot = Bp_dot;
    model.Bp_ddot = Bp_ddot;
    model.C = C;
    model.Cp = Cp;
    model.Cp_dot = Cp_dot;
    model.Cp_ddot = Cp_ddot;
    model.D = D;
    model.Dp = Dp;
    model.G = G;
    model.Gp = Gp;
    model.Gp_dot = Gp_dot;
    model.Gp_ddot = Gp_ddot;
    model.K = K;
    model.Kp = Kp;
    model.Bm = Bm;
    model.Bmp = Bmp;
    model.Dm = Dm;
    model.Dmp = Dmp;
    model.J_i = J_i;
    model.J_cm_i = J_cm_i;
    model.J_cm_i_p = J_cm_i_p;
    model.J_ee = J_ee;
    model.T_0_i = T_0_i;
    model.Y = Y;
    if CREATE_REGRESSOR_KERNEL == 1
        model.null_Y = null_Y;
    end
    save ('model.mat', 'model');
    disp('done!')
    toc
end

%% OTHER
function R = R_x(phi)
    R = [1, 0, 0; 0, cos(phi), -sin(phi); 0, sin(phi), cos(phi)];
end

function R = R_y(theta)
    R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
end

function R = R_z(psi)
    R = [cos(psi), -sin(psi), 0; sin(psi), cos(psi), 0; 0, 0, 1];
end

function S = skew(v)
S = [0 -v(3) v(2);...
    v(3)  0  -v(1);...
    -v(2) v(1) 0 ];
end 

function v = vect(S)
    
%     if S(1,1)~=0 || S(2,2)~=0 || S(3,3)~=0
%         error('Input is not a skew matrix: diagonal element not equal to 0!')
%     elseif ~isequal(S(3,2),-S(2,3)) || ~isequal(S(1,3),-S(3,1)) || ~isequal(S(2,1),-S(1,2))
%         error('Input is not a skew matrix')
%     end
    
    v = [S(3,2); S(1,3); S(2,1)];
end
out = 1;
end

