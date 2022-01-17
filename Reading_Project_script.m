%% Ahmet Burak KOC - kocab@hotmail.com

%% Initialization
clear; clc; close all
addpath('funs')
% 3x3 Nominal moment of inertia [kg.m^2] eq(49)
    J0 = [  800.027   0           0         ;
            0         839.930     0         ;
            0         0           289.930   ];
    % 3x3 Uncertain inertia matrix [kg.m^2]
    %Delta_J = diag([50 30 20]');
    Delta_J = diag([0 0 0]');
    % 3x3 Actual inertia matrix
    J       = J0 + Delta_J;
    

% Flexible parameters
    % Vibration-mode frequency matrix i=1,2,3,4
    omega_n = [1.0973, 1.2761, 1.6538, 2.2893]';
    % Ratio of vibration-mode damping i=1,2,3,4
    zeta    = [0.05, 0.06, 0.08, 0.025]';
    % Coupling matrix (delta) 4x3 matrix
    coupling_matrix = [     6.45637     1.27814     2.15629     ;
                            -1.25619    0.91756     -1.67264    ;
                            1.11687     2.48901     -0.83674    ;
                            1.23637     -2.65810    -1.12503    ];
                        
    % K, L -> stiffness and damping matrices
    K   =   diag(omega_n.^2);
    L   =   diag(2*omega_n.*zeta);

    % 4x1 Initial flexible vibration
    eta_init    = [0.01242, 0.01584, -0.01749, 0.01125]';
    
    % 3x1 Initial Angular velocity (omega)
    omega_init  = [0, 0, 0]';
    
    % 4x1 Initial attitude in quaternion
    q_init = [0.3, -0.2, -0.3, 0.8832]' ;

    %qd_init = zeros(4,1) ;
    qd_init = [0, 0, 0, 1]' ;
    
    % disturbances
    d = zeros(3,1);

    % sliding_mode
    epsilon = 0.1   ;
    iota1   = 3     ;
    iota2   = 5     ;
    q_p     = 3/5   ;
    k11 = 0.5 ;
    k12 = 0.5 ;
    k13 = 0.5 ;
    k21 = 0.5 ;
    k22 = 0.5 ;
    k23 = 0.5 ;
    k1  = [k11, k12, k13]';
    k2  = [k21, k22, k23]';
    K1 = diag(k1);
    K2 = diag(k2);
    
    % Chebyshev Polynomial
    % order_num = 2
    sigma1  = 100   ;
    sigma2  = 0.01  ;
    weight_init = rand(3,25) ;
    
    % Neural Network Controller
    % J0 ;
    tau     = 10    ;
    chi     = 25    ;
    gamma   = 0.35  ;
    lambda  = 0.6   ;
    
    % Actuator parameters
    %%%%%%%%%%% to be changed later %%%%%%%%%%
    um = 68; % maximum torque Nm
    
    
    save('Reading_Project.mat')