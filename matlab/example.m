%% ----- INIT ----- %%
clear all
conf.N_joints = 2;               % number of joints
conf.Q_MIN = [-pi; -3*pi/4];
conf.Q_MAX = [pi; 3*pi/4];
% - defines - %
conf.CREATE_FUNCTIONS = 1;
conf.SAVE_MODEL = 1;
conf.CREATE_REGRESSOR_KERNEL = 1;   % esperimental
conf.REGRESSOR_KERNEL_SYMB = 1;
% - hyperparameters - %
conf.joints_type = {'R', 'R'};  % rotoidal:'R' or prismatic:'P'
% - known parameters - %
syms nopar real
conf.par_nl = [nopar];		% nonlinear parameters, non used in adaptive
% - Frames - %
conf.frames_tr = sym([0,0,-0.1;   0,0,0.2;      0,-0.2,0]);%;     0,0,0.2]);
conf.frames_or = sym([-pi,0,0;   -pi/2,0,0;    pi/2,0,0]);%;     -pi/2,0,0]);
% frames_tr = sym('t%d_%d', [N_frames, 3], 'real');         % traslation XYZ
% frames_or = sym('o%d_%d', [N_frames, 3], 'real');         % rotation   Euler XYZ

create_model_symb(conf)