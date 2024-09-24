%% ----- INIT ----- %%
clear all
conf.N_joints = 3;               % number of joints
conf.Q_MIN = [-pi; -pi; -pi];
conf.Q_MAX = [pi; pi; pi];
% - defines - %
conf.CREATE_FUNCTIONS = 1;
conf.SAVE_MODEL = 1;
conf.CREATE_REGRESSOR_KERNEL = 0;   % esperimental
conf.REGRESSOR_KERNEL_SYMB = 0;
% - hyperparameters - %
conf.joints_type = {'R', 'R', 'R'};  % rotoidal:'R' or prismatic:'P'
% - known parameters - %
syms nopar real
conf.par_nl = [nopar];		% nonlinear parameters, non used in adaptive
% - Frames - %
conf.frames_tr = sym([0,0,1;   0,0,0;      1,0,0;  1,0,0]);%;     0,0,0.2]);
conf.frames_or = sym([0,0,0;   pi/2,0,0;   0,0,0;  -pi/2,0,0]);
% frames_tr = sym('t%d_%d', [N_frames, 3], 'real');         % traslation XYZ
% frames_or = sym('o%d_%d', [N_frames, 3], 'real');         % rotation   Euler XYZ

create_model_symb(conf)



