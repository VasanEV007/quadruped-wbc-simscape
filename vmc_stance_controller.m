function [tau_legs, F_grf, wrench_cmd, info] = vmc_stance_controller(state, des_state, contact_mask, q_legs, p)
%VMC_STANCE_CONTROLLER Weighted least-squares GRF distributor & VMC controller.
%   [tau_legs, F_grf, wrench_cmd, info] = VMC_STANCE_CONTROLLER(state, des_state, contact_mask, q_legs, p)
%   computes the desired 6-DOF centroidal wrench using Virtual Model Control
%   (VMC) with full gravity compensation for the 82.6 kg quadruped + 6-axis
%   arm + vacuum gripper, distributes ground reaction forces (GRFs) across
%   the stance feet via friction-pyramid-constrained Weighted Least-Squares
%   (WLS), and maps the resulting forces to hip-pitch and knee-pitch joint
%   torques via the 3D leg Jacobians.
%
%   Inputs:
%     state        - Struct with current torso state:
%                      .p_com : 3x1 CoM position in world frame [m]
%                      .v_com : 3x1 CoM linear velocity in world frame [m/s]
%                      .rpy   : 3x1 Euler angles [roll; pitch; yaw] [rad]
%                      .omega : 3x1 Angular velocity in world frame [rad/s]
%     des_state    - Struct with desired torso references:
%                      .p_com         : 3x1 desired CoM position [m]
%                      .v_com         : 3x1 desired CoM velocity [m/s]
%                      .rpy           : 3x1 desired Euler angles [rad]
%                      .omega         : 3x1 desired angular velocity [rad/s]
%                      .a_com_ff      : (Optional) 3x1 feedforward accel [m/s^2]
%                      .wrench_arm_ff : (Optional) 6x1 arm disturbance wrench
%     contact_mask - 4x1 binary/logical stance mask [c_FL; c_FR; c_RL; c_RR]
%                    (1 = stance, 0 = swing)
%     q_legs       - 2x4 matrix (or 8x1 vector) of joint angles [q_hip; q_knee]
%                    for legs [FL, FR, RL, RR] [rad]
%     p            - (Optional) Parameter struct from params()
%
%   Outputs:
%     tau_legs   - 2x4 matrix of joint torque commands [tau_hip; tau_knee]
%                  for [FL, FR, RL, RR] [Nm]
%     F_grf      - 3x4 matrix of distributed world-frame ground reaction
%                  forces [Fx; Fy; Fz] for [FL, FR, RL, RR] [N]
%     wrench_cmd - 6x1 desired centroidal wrench [Fx; Fy; Fz; Tx; Ty; Tz]
%                  in world frame [N; Nm]
%     info       - Diagnostic struct (A matrix, actual wrench, WLS residual,
%                  friction utilization, foot positions)

    if nargin < 5 || isempty(p)
        p = params();
    end
    if nargin < 4 || isempty(q_legs)
        q_legs = repmat(p.q_leg_nom, 1, 4);
    end
    if nargin < 3 || isempty(contact_mask)
        contact_mask = [1; 1; 1; 1];
    end
    if nargin < 2 || isempty(des_state)
        des_state = default_torso_state(p);
    end
    if nargin < 1 || isempty(state)
        state = default_torso_state(p);
    end

    % Reshape q_legs to 2x4 if passed as 8x1 column vector
    if numel(q_legs) == 8
        q_legs = reshape(double(q_legs), 2, 4);
    else
        error('vmc_stance_controller:InvalidQLegs', ...
            'q_legs must be a 2x4 matrix or 8x1 vector.');
    end
    contact_mask = logical(contact_mask(:));
    if numel(contact_mask) ~= 4
        error('vmc_stance_controller:InvalidContactMask', ...
            'contact_mask must have 4 elements [FL; FR; RL; RR].');
    end

    % Fill optional fields in state and des_state
    state = ensure_state_fields(state, p);
    des_state = ensure_state_fields(des_state, p);

    %% 1. Orientation Matrices & SO(3) Orientation Error
    R_wb = rpy_to_rotm(state.rpy);          % Body-to-world rotation matrix
    R_des = rpy_to_rotm(des_state.rpy);     % Desired body-to-world rotation

    % Orientation error in world frame: e_R = 0.5 * vee(R_des * R_wb' - R_wb * R_des')
    E_rot = 0.5 * (R_des * R_wb.' - R_wb * R_des.');
    e_R = [E_rot(3, 2); E_rot(1, 3); E_rot(2, 1)];

    %% 2. VMC Centroidal Wrench Command + Gravity Compensation
    e_pos = des_state.p_com - state.p_com;
    e_vel = des_state.v_com - state.v_com;
    e_omega = des_state.omega - state.omega;

    % Gravity compensation vector in world frame (+Z up)
    F_gravity_comp = [0.0; 0.0; p.m_total * p.g];

    % Linear force command [Fx; Fy; Fz]
    F_cmd = p.Kp_com * e_pos + p.Kd_com * e_vel + ...
            F_gravity_comp + p.m_total * des_state.a_com_ff;

    % Rotate rotational PD gains from body frame to world frame
    Kp_rot_w = R_wb * p.Kp_rot * R_wb.';
    Kd_rot_w = R_wb * p.Kd_rot * R_wb.';

    % Gravity moment compensation if CoM is offset from geometric center
    r_com_w = R_wb * p.r_com_offset;
    tau_gravity_comp = cross(r_com_w, F_gravity_comp);

    % Angular torque command [Tx; Ty; Tz]
    tau_cmd = Kp_rot_w * e_R + Kd_rot_w * e_omega + tau_gravity_comp;

    % Net desired centroidal wrench (compensating 6-axis arm reaction wrench)
    wrench_cmd = [F_cmd; tau_cmd] - des_state.wrench_arm_ff;

    %% 3. Centroidal Moment Matrix A (6x12) & Foot Kinematics
    A = zeros(6, 12);
    p_feet_body = zeros(3, 4);
    p_feet_world = zeros(3, 4);
    J_legs_body = zeros(3, 2, 4);

    for i = 1:4
        [p_foot_b, J_3d_b] = leg_kinematics(q_legs(:, i), i, p);
        p_feet_body(:, i) = p_foot_b;
        J_legs_body(:, :, i) = J_3d_b;

        % Lever arm from CoM to foot i in world frame
        r_foot_w = R_wb * (p_foot_b - p.r_com_offset);
        p_feet_world(:, i) = state.p_com + r_foot_w;

        if contact_mask(i)
            cols = (3*i - 2):(3*i);
            A(1:3, cols) = eye(3);
            A(4:6, cols) = skew3(r_foot_w);
        end
    end

    %% 4. Weighted Least-Squares (WLS) with Friction Pyramid Constraints
    % Objective: min_F (A*F - wrench_cmd)' * S * (A*F - wrench_cmd) + F' * W * F
    % Equivalent QP: min_F 0.5 * F' * H * F + g_qp' * F
    S = p.S_wrench;
    W = p.W_force;
    H = A.' * S * A + W;
    H = 0.5 * (H + H.'); % Ensure strict symmetry
    g_qp = -A.' * S * wrench_cmd;

    F_vec = zeros(12, 1);
    active_legs = find(contact_mask).';

    if ~isempty(active_legs)
        active_cols = zeros(1, 3 * numel(active_legs));
        for idx = 1:numel(active_legs)
            leg_i = active_legs(idx);
            active_cols((3*idx - 2):(3*idx)) = (3*leg_i - 2):(3*leg_i);
        end

        H_act = H(active_cols, active_cols);
        g_act = g_qp(active_cols);

        % Step 4a: Unconstrained analytical WLS solution
        F_unc_act = -H_act \ g_act;
        F_vec(active_cols) = F_unc_act;

        % Step 4b: Check friction pyramid & normal force feasibility
        F_proj = project_all_legs(F_vec, contact_mask, p);
        is_unconstrained_feasible = (norm(F_proj - F_vec) <= 1.0e-9);
        num_iters = 0;

        if ~is_unconstrained_feasible
            % Step 4c: Accelerated Projected Gradient (FISTA) with Adaptive Restart
            L_lip = max(eig(H_act));
            step_size = 1.0 / max(L_lip, 1.0e-6);

            x_curr = F_proj;
            y_curr = x_curr;
            t_fista = 1.0;

            for iter = 1:p.wls_max_iter
                num_iters = iter;
                grad_y = zeros(12, 1);
                grad_y(active_cols) = H_act * y_curr(active_cols) + g_act;

                x_next = project_all_legs(y_curr - step_size * grad_y, contact_mask, p);

                if norm(x_next - x_curr) <= p.wls_tol
                    x_curr = x_next;
                    break;
                end

                t_next = 0.5 * (1.0 + sqrt(1.0 + 4.0 * t_fista^2));
                y_next = x_next + ((t_fista - 1.0) / t_next) * (x_next - x_curr);

                % Adaptive gradient restart for monotonic convergence
                if (y_curr - x_next).' * (x_next - x_curr) > 0
                    t_next = 1.0;
                    y_next = x_next;
                end

                x_curr = x_next;
                y_curr = y_next;
                t_fista = t_next;
            end
            F_vec = x_curr;
        end
    else
        is_unconstrained_feasible = true;
        num_iters = 0;
    end

    F_grf = reshape(F_vec, 3, 4);

    %% 5. Map World GRFs to Joint Torques: tau_i = -J_3d_i' * (R_wb' * F_grf_i)
    tau_legs = zeros(2, 4);
    friction_ratio = zeros(1, 4);

    for i = 1:4
        if contact_mask(i)
            F_world_i = F_grf(:, i);
            F_body_i = R_wb.' * F_world_i;
            J_i = J_legs_body(:, :, i);

            % Stance joint torque: tau = -J^T * F_grf_body
            tau_raw = -J_i.' * F_body_i;

            % Saturate at actuator torque limits
            tau_legs(1, i) = min(max(tau_raw(1), -p.tau_hip_max),  p.tau_hip_max);
            tau_legs(2, i) = min(max(tau_raw(2), -p.tau_knee_max), p.tau_knee_max);

            friction_ratio(i) = max(abs(F_world_i(1)), abs(F_world_i(2))) / ...
                                max(p.mu_vmc * F_world_i(3), 1.0e-9);
        end
    end

    %% 6. Populate Diagnostic Info Struct
    if nargout >= 4
        info.A = A;
        info.wrench_actual = A * F_vec;
        info.wrench_error = info.wrench_actual - wrench_cmd;
        info.wls_cost = (info.wrench_error.' * S * info.wrench_error) + (F_vec.' * W * F_vec);
        info.p_feet_body = p_feet_body;
        info.p_feet_world = p_feet_world;
        info.friction_ratio = friction_ratio;
        info.is_unconstrained_feasible = is_unconstrained_feasible;
        info.solver_iterations = num_iters;
    end
end

%% ---------------- Helper Functions ----------------

function s = default_torso_state(p)
    s.p_com = [0.0; 0.0; p.z_nom];
    s.v_com = [0.0; 0.0; 0.0];
    s.rpy = [0.0; 0.0; 0.0];
    s.omega = [0.0; 0.0; 0.0];
    s.a_com_ff = [0.0; 0.0; 0.0];
    s.wrench_arm_ff = zeros(6, 1);
end

function s = ensure_state_fields(s, p)
    if ~isfield(s, 'p_com') || isempty(s.p_com)
        s.p_com = [0.0; 0.0; p.z_nom];
    end
    if ~isfield(s, 'v_com') || isempty(s.v_com)
        s.v_com = [0.0; 0.0; 0.0];
    end
    if ~isfield(s, 'rpy') || isempty(s.rpy)
        s.rpy = [0.0; 0.0; 0.0];
    end
    if ~isfield(s, 'omega') || isempty(s.omega)
        s.omega = [0.0; 0.0; 0.0];
    end
    if ~isfield(s, 'a_com_ff') || isempty(s.a_com_ff)
        s.a_com_ff = [0.0; 0.0; 0.0];
    end
    if ~isfield(s, 'wrench_arm_ff') || isempty(s.wrench_arm_ff)
        s.wrench_arm_ff = zeros(6, 1);
    end
    s.p_com = double(s.p_com(:));
    s.v_com = double(s.v_com(:));
    s.rpy = double(s.rpy(:));
    s.omega = double(s.omega(:));
    s.a_com_ff = double(s.a_com_ff(:));
    s.wrench_arm_ff = double(s.wrench_arm_ff(:));
end

function R = rpy_to_rotm(rpy)
%RPY_TO_ROTM Convert [roll; pitch; yaw] (ZYX convention) to 3x3 rotation matrix.
    cr = cos(rpy(1)); sr = sin(rpy(1));
    cp = cos(rpy(2)); sp = sin(rpy(2));
    cy = cos(rpy(3)); sy = sin(rpy(3));

    R = [ ...
        cy*cp, cy*sp*sr - sy*cr, cy*sp*cr + sy*sr; ...
        sy*cp, sy*sp*sr + cy*cr, sy*sp*cr - cy*sr; ...
        -sp,   cp*sr,            cp*cr             ...
    ];
end

function S = skew3(v)
%SKEW3 Return 3x3 skew-symmetric cross-product matrix [v]_x.
    S = [ ...
         0.0,  -v(3),  v(2); ...
         v(3),  0.0,  -v(1); ...
        -v(2),  v(1),  0.0   ...
    ];
end

function F_proj = project_all_legs(F_vec, contact_mask, p)
%PROJECT_ALL_LEGS Project 12x1 force vector onto per-leg friction pyramids.
    F_proj = zeros(12, 1);
    for i = 1:4
        if contact_mask(i)
            cols = (3*i - 2):(3*i);
            F_proj(cols) = project_friction_pyramid(F_vec(cols), p.mu_vmc, p.Fz_min, p.Fz_max);
        end
    end
end

function u_best = project_friction_pyramid(v, mu, fz_min, fz_max)
%PROJECT_FRICTION_PYRAMID Exact Euclidean projection onto 3D truncated friction pyramid.
%   Solves: min_u ||u - v||_2^2
%           s.t.  fz_min <= u_z <= fz_max
%                 -mu * u_z <= u_x <= mu * u_z
%                 -mu * u_z <= u_y <= mu * u_z
    vx = v(1);
    vy = v(2);
    vz = v(3);

    % Candidate 1D optimal u_z values for the 4 active-face combinations
    % plus box endpoints
    uz_candidates = [ ...
        vz; ...
        (vz + mu * abs(vx)) / (1.0 + mu^2); ...
        (vz + mu * abs(vy)) / (1.0 + mu^2); ...
        (vz + mu * (abs(vx) + abs(vy))) / (1.0 + 2.0 * mu^2); ...
        fz_min; ...
        fz_max  ...
    ];

    best_dist_sq = inf;
    u_best = [0.0; 0.0; fz_min];

    for k = 1:numel(uz_candidates)
        uz = min(max(uz_candidates(k), fz_min), fz_max);
        ux = min(max(vx, -mu * uz), mu * uz);
        uy = min(max(vy, -mu * uz), mu * uz);

        dist_sq = (ux - vx)^2 + (uy - vy)^2 + (uz - vz)^2;
        if dist_sq < best_dist_sq
            best_dist_sq = dist_sq;
            u_best = [ux; uy; uz];
        end
    end
end
