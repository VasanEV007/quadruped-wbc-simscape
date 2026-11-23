function p = params()
%PARAMS Robot inertial, geometric, contact, VMC, gait, navigation & 6-DOF arm parameters.
%   p = PARAMS() returns a struct containing parameters for the 82.6 kg
%   quadruped robot (8 actuated leg joints: hip pitch and knee pitch per
%   leg, no hip abduction) equipped with a 6-axis robotic arm and vacuum
%   gripper for MATLAB Simscape Multibody whole-body control.

    %% 1. Mass Budget (Total Unloaded Robot = 82.6 kg)
    p.m_trunk = 50.0;                           % Main torso/chassis mass [kg]
    p.m_thigh = 3.2;                            % Single thigh link mass [kg]
    p.m_shank = 1.6;                            % Single shank link mass [kg]
    p.m_leg = p.m_thigh + p.m_shank;            % Single leg mass (4.8 kg) [kg]
    p.num_legs = 4;                             % Number of legs (FL, FR, RL, RR)

    % 6-Axis Robotic Arm + Vacuum Gripper Mass Breakdown
    p.m_arm_links = [3.0, 3.2, 2.4, 1.4, 1.2, 0.8]; % 6-DOF arm link masses [kg]
    p.m_arm = sum(p.m_arm_links);               % Total 6-axis arm mass (12.0 kg) [kg]
    p.m_gripper = 1.4;                          % Vacuum gripper assembly mass [kg]
    p.m_payload = 0.0;                          % External payload held by vacuum gripper [kg]
    p.m_box = 0.50;                             % Mass of pick-and-place payload box [kg]

    % Total System Mass (50.0 + 4*4.8 + 12.0 + 1.4 = 82.6 kg)
    p.m_total = p.m_trunk + p.num_legs * p.m_leg + p.m_arm + p.m_gripper + p.m_payload;
    p.g = 9.81;                                 % Gravitational acceleration [m/s^2]

    %% 2. Inertial Properties (Composite Torso + Arm at Nominal Stow)
    p.I_body = diag([3.85, 8.40, 9.15]);
    p.I_thigh_yy = 0.032;                       % Thigh pitch inertia about CoM [kg*m^2]
    p.I_shank_yy = 0.016;                       % Shank pitch inertia about CoM [kg*m^2]

    % Nominal CoM offset in body geometric frame [x; y; z] [m]
    p.r_com_offset = [0.0; 0.0; 0.0];
    p.trunk_com_offset = [-0.104158; 0.0; 0.0]; % Trunk internal CoM counterbalance [m]

    %% 3. Kinematic & Geometric Parameters
    % Body frame convention: +X forward, +Y left, +Z up
    p.L_body = 0.70;                            % Longitudinal hip-to-hip distance [m]
    p.W_body = 0.44;                            % Lateral hip-to-hip distance [m]
    p.L_thigh = 0.34;                           % Thigh link length (hip to knee) [m]
    p.L_shank = 0.34;                           % Shank link length (knee to foot) [m]
    p.r_foot = 0.025;                           % Spherical foot contact radius [m]

    p.leg_names = {'FL', 'FR', 'RL', 'RR'};
    p.hip_offsets = [ ...
         p.L_body/2,  p.L_body/2, -p.L_body/2, -p.L_body/2;  % X [m]
         p.W_body/2, -p.W_body/2,  p.W_body/2, -p.W_body/2;  % Y [m]
         0.0,         0.0,         0.0,         0.0          % Z [m]
    ];

    % Nominal standing configuration (q_hip = -pi/4 rad, q_knee = +pi/2 rad)
    p.q_hip_nom = -pi / 4;                      % Nominal hip pitch angle [rad]
    p.q_knee_nom = pi / 2;                      % Nominal knee pitch angle [rad]
    p.q_leg_nom = [p.q_hip_nom; p.q_knee_nom];
    p.z_nom = p.L_thigh * cos(p.q_hip_nom) + p.L_shank * cos(p.q_hip_nom + p.q_knee_nom);

    % Joint limits (Hip Pitch and Knee Pitch)
    p.q_hip_lim = [-1.40, 1.40];                % Hip pitch limits [rad]
    p.q_knee_lim = [0.25, 2.50];                % Knee pitch limits [rad]
    p.tau_hip_max = 180.0;                      % Max hip pitch torque [Nm]
    p.tau_knee_max = 220.0;                     % Max knee pitch torque [Nm]

    %% 4. Ground Contact Stiffness & Damping (Simscape Spatial Contact)
    p.k_contact = 1.5e5;                        % Normal ground contact stiffness [N/m]
    p.b_contact = 2.5e3;                        % Normal ground contact damping [N*s/m]
    p.mu_static = 0.85;                         % Static friction coefficient
    p.mu_kinetic = 0.70;                        % Dynamic/kinetic friction coefficient
    p.mu_vmc = 0.65;                            % Conservative friction pyramid coeff for VMC
    p.v_transition = 0.01;                      % Friction Stribeck transition velocity [m/s]

    %% 5. Virtual Model Control (VMC) & WLS Force Distribution Parameters
    p.Kp_com = diag([2200.0, 1800.0, 4800.0]);  % [Fx, Fy, Fz] proportional gains [N/m]
    p.Kd_com = diag([350.0,  280.0,  650.0]);   % [Fx, Fy, Fz] derivative gains [N*s/m]

    p.Kp_rot = diag([3600.0, 4800.0, 1400.0]);  % [Tx, Ty, Tz] proportional gains [Nm/rad]
    p.Kd_rot = diag([360.0,  480.0,  180.0]);   % [Tx, Ty, Tz] derivative gains [Nm*s/rad]

    p.S_wrench = diag([25.0, 5.0, 80.0, 90.0, 90.0, 35.0]);
    p.alpha_reg = 1.0e-3;                       % Tikhonov force regularization weight [1/N^2]
    p.W_force = diag(repmat([p.alpha_reg, 50.0, p.alpha_reg], 1, 4));

    p.Fz_min = 15.0;                            % Minimum unilateral normal force [N]
    p.Fz_max = 600.0;                           % Maximum vertical ground reaction force [N]
    p.wls_max_iter = 60;                        % Max iterations for projected WLS solver
    p.wls_tol = 1.0e-7;                         % Convergence tolerance for projected WLS

    %% 6. Stage 3: Diagonal Trot Gait Scheduler & Swing-Leg PD Parameters
    p.T_gait = 0.40;                            % Full trot cycle period [s]
    p.duty_factor = 0.50;                       % Stance duty factor beta in (0, 1)
    p.t_stand_init = 0.20;                      % Initial 4-leg settle window before trot [s]
    p.h_swing = 0.06;                           % Nominal flat-ground swing apex clearance [m]
    p.k_raibert = 0.045;                        % Raibert step placement velocity gain [s]

    p.Kp_swing_cart = diag([2800.0, 0.0, 3200.0]); % [N/m]
    p.Kd_swing_cart = diag([160.0,  0.0, 180.0]);  % [N*s/m]

    p.Kp_swing_joint = [260.0; 280.0];          % [Nm/rad]
    p.Kd_swing_joint = [22.0;  22.0];           % [Nm*s/rad]

    p.Kp_stance_joint = [320.0; 360.0];         % [Nm/rad]
    p.Kd_stance_joint = [28.0;  30.0];          % [Nm*s/rad]

    %% 7. Stage 4: Step Climbing, Foothold Edge Avoidance & Pure-Pursuit Navigation
    % 10 cm Step Platform Geometry (front edge at x = 1.20 m, back edge at x = 1.85 m)
    p.x_step_start = 1.20;                      % Rising step edge X coordinate [m]
    p.x_step_end   = 1.85;                      % Falling step edge X coordinate [m]
    p.y_step_min   = -0.80;                     % Step platform right lateral boundary [m]
    p.y_step_max   = 2.20;                      % Step platform left lateral boundary [m]
    p.h_step       = 0.10;                      % Step platform height (10 cm) [m]

    p.d_edge_margin = 0.055;                    % Minimum foothold clearance from step edge (>= 0.05 m) [m]
    p.h_swing_step  = 0.145;                    % Swing apex clearance during step climb/descent (>= 0.14 m) [m]
    p.Kv_raibert    = 0.045;                    % Raibert velocity tracking feedback gain [s]
    p.c_turn_scale  = 0.35;                     % Step-length down-scaling factor with |omega_z| [s/rad]
    p.max_step_x    = 0.13;                     % Maximum sagittal foothold offset from hip [m]

    p.R_uturn        = 0.70;                    % Semicircular U-turn radius [m]
    p.x_uturn_center = 2.15;                    % U-turn arc entry/exit X coordinate [m]
    p.y_uturn_center = 0.70;                    % U-turn semicircle center Y coordinate [m]
    p.x_return_goal  = 0.55;                    % Final return X coordinate in front of table [m]
    p.L_lookahead    = 0.40;                    % Pure-pursuit lookahead distance (0.35-0.45 m) [m]
    p.v_cruise       = 0.42;                    % Nominal forward cruise velocity [m/s]
    p.v_step_cross   = 0.34;                    % Reduced forward speed across step transitions [m/s]
    p.omega_z_max    = 0.85;                    % Maximum yaw rate command [rad/s]

    %% 8. Stage 5: 6-Axis Arm Kinematics, Vacuum Gripper & 22.4 s Mission Schedule
    p.arm_base_offset  = [0.22; 0.0; 0.14];     % Arm base mount in body frame [x; y; z] [m]
    p.arm_link_lengths = [0.10, 0.28, 0.26, 0.08, 0.04, 0.04]; % 6 link lengths [m]
    p.arm_ee_offset    = 0.035;                 % Vacuum suction cup tip offset along tool axis [m]
    p.q_arm_stow       = [0.0; 0.45; 1.85; 0.8416; 0.0; 0.0]; % Compact stowed arm joint angles [rad]

    % Arm joint PD gains & torque limits
    p.Kp_arm = [420.0; 550.0; 480.0; 260.0; 160.0; 120.0]; % [Nm/rad]
    p.Kd_arm = [38.0;  48.0;  42.0;  22.0;  14.0;  10.0];  % [Nm*s/rad]
    p.tau_arm_max = [120.0; 150.0; 120.0; 60.0; 40.0; 30.0]; % [Nm]

    % Pick-and-Place Table & Box Payload Geometry
    % Table sits in median between outbound lane (y=0) and return lane (y=1.40 m)
    p.table_center_xy = [0.30; 0.70];           % Table center [x; y] in world frame [m]
    p.table_dims      = [0.50, 0.84, 0.45];     % Table [length_x, width_y, height_z] [m]
    p.table_height    = 0.45;                   % Table top surface elevation [m]
    p.box_side        = 0.10;                   % Payload box cubic edge length [m]
    p.box_init_pos    = [0.38; 0.34; 0.50];     % Initial box center on table [x; y; z] [m]
    p.box_place_pos   = [0.22; 1.06; 0.50];     % Target place box center on table [x; y; z] [m]

    % Key End-Effector Poses in Body Frame [x; y; z] relative to Torso Center [m]
    % At nominal standing height z_torso = z_nom + r_foot = 0.505833 m,
    % box top surface (z_world = 0.55 m) is at z_body = 0.55 - 0.505833 = +0.044167 m.
    p.p_ee_home_b      = [0.36; 0.00; 0.24];    % Home / retracted EE pose in body frame [m]
    p.p_ee_carry_b     = [0.35; 0.00; 0.26];    % Tucked carry EE pose in body frame [m]
    p.p_ee_prepick_b   = [0.38; 0.34; 0.11];    % Pre-pick hover above box in body frame [m]
    p.p_ee_pick_b      = [0.38; 0.34; 0.0442];  % Pick contact at box top surface in body frame [m]
    p.p_ee_preplace_b  = [0.33; 0.34; 0.11];    % Pre-place hover above table (at return pose) [m]
    p.p_ee_place_b     = [0.33; 0.34; 0.0442];  % Place contact on table surface in body frame [m]

    % Vacuum Suction Gripper Parameters
    p.suction_prox_thresh   = 0.015;            % Proximity threshold (< 0.015 m) to engage suction [m]
    p.vacuum_suction_max_force = 120.0;         % Max vacuum suction holding force [N]
    p.kp_suction            = 8500.0;           % Suction weld attraction stiffness [N/m]
    p.kd_suction            = 280.0;            % Suction weld attraction damping [N*s/m]

    % 22.4 s Autonomous Mission Phase Boundary Timestamps [s]
    p.t_phase1_end = 3.5;                       % Phase 1: Pick (0.0 - 3.5 s)
    p.t_phase2_end = 9.0;                       % Phase 2: Outbound & Step Climb (3.5 - 9.0 s)
    p.t_phase3_end = 14.5;                      % Phase 3: Wide-Arc U-Turn (9.0 - 14.5 s)
    p.t_phase4_end = 19.5;                      % Phase 4: Return Climb & Table Approach (14.5 - 19.5 s)
    p.t_mission_end = 22.4;                     % Phase 5: Place & Retract (19.5 - 22.4 s)

    if nargout == 0
        assignin('caller', 'p', p);
    end
end
