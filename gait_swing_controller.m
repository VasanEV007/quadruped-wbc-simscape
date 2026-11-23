function [contact_mask, tau_swing, p_foot_des, v_foot_des, phase_legs, q_des_legs] = ...
    gait_swing_controller(t, q_legs, dq_legs, v_com_body, v_cmd_body, gait_mode, p, p_td_body_ext, h_swing_legs_ext)
%GAIT_SWING_CONTROLLER Diagonal trot scheduler and operational-space swing PD controller.
%   [contact_mask, tau_swing, p_foot_des, v_foot_des, phase_legs, q_des_legs] = ...
%       GAIT_SWING_CONTROLLER(t, q_legs, dq_legs, v_com_body, v_cmd_body, gait_mode, p, p_td_body_ext, h_swing_legs_ext)
%   computes the diagonal trot phase schedule, contact mask, smooth cycloidal
%   swing-foot trajectories (with adaptive step-climb apex clearance h_swing),
%   and operational-space PD + gravity compensation joint torques for swing legs.

    if nargin < 7 || isempty(p)
        p = params();
    end
    if nargin < 5 || isempty(v_cmd_body)
        v_cmd_body = [0.0; 0.0; 0.0];
    end
    if nargin < 4 || isempty(v_com_body)
        v_com_body = [0.0; 0.0; 0.0];
    end
    if nargin < 6 || isempty(gait_mode)
        if norm(v_cmd_body(1:2)) < 1.0e-4
            gait_mode = 0;
        else
            gait_mode = 1;
        end
    end
    if nargin < 3 || isempty(dq_legs)
        dq_legs = zeros(2, 4);
    end
    if nargin < 2 || isempty(q_legs)
        q_legs = repmat(p.q_leg_nom, 1, 4);
    end

    has_ext_td = (nargin >= 8) && ~isempty(p_td_body_ext);
    if nargin >= 9 && ~isempty(h_swing_legs_ext)
        h_swing_vec = double(h_swing_legs_ext(:));
    else
        h_swing_vec = p.h_swing * ones(4, 1);
    end

    q_mat   = reshape(double(q_legs), 2, 4);
    dq_mat  = reshape(double(dq_legs), 2, 4);
    v_com_b = double(v_com_body(:));
    v_cmd_b = double(v_cmd_body(:));

    %% 1. Diagonal Trot Phase Scheduler
    % Diagonal Pair A: FL (1) & RR (4) at phase offset 0.0
    % Diagonal Pair B: FR (2) & RL (3) at phase offset 0.5
    phase_offsets = [0.0; 0.5; 0.5; 0.0];
    beta   = p.duty_factor;
    T_gait = p.T_gait;

    contact_mask = ones(4, 1);
    phase_legs   = zeros(4, 1);

    if gait_mode >= 1 && t >= p.t_stand_init
        t_active = t - p.t_stand_init;
        for i = 1:4
            phi_i = mod(t_active / T_gait + phase_offsets(i), 1.0);
            phase_legs(i) = phi_i;
            if phi_i < beta
                contact_mask(i) = 1; % Stance phase
            else
                contact_mask(i) = 0; % Swing phase
            end
        end
    end

    %% 2. Swing-Leg Trajectory & Operational-Space PD + Gravity Compensation
    tau_swing  = zeros(2, 4);
    p_foot_des = zeros(3, 4);
    v_foot_des = zeros(3, 4);
    q_des_legs = repmat(p.q_leg_nom, 1, 4);

    L1 = p.L_thigh;
    L2 = p.L_shank;
    T_stance = beta * T_gait;
    T_swing  = max((1.0 - beta) * T_gait, 1.0e-3);

    for i = 1:4
        hip_i = p.hip_offsets(:, i);
        q_i   = q_mat(:, i);
        dq_i  = dq_mat(:, i);
        h_sw  = h_swing_vec(i);

        [p_foot_3d, J_3d, ~, J_2d] = leg_kinematics(q_i, i, p);
        v_foot_3d = J_3d * dq_i;

        if has_ext_td
            x_touch_rel = p_td_body_ext(1, i) - hip_i(1);
            z_touch_rel = p_td_body_ext(3, i);
            x_lift_rel  = -x_touch_rel;
            z_lift_rel  = z_touch_rel;
        else
            x_lift_rel  = min(max(-0.5 * v_cmd_b(1) * T_stance, -0.12), 0.12);
            x_touch_rel = min(max( 0.5 * v_cmd_b(1) * T_stance + ...
                p.k_raibert * (v_com_b(1) - v_cmd_b(1)), -0.12), 0.12);
            z_lift_rel  = -p.z_nom;
            z_touch_rel = -p.z_nom;
        end

        if contact_mask(i) == 1
            % Stance leg: desired terrain-adjusted foothold under hip
            if has_ext_td
                s_st = phase_legs(i) / max(beta, 1.0e-3);
                x_st_rel = (1.0 - s_st) * x_touch_rel + s_st * x_lift_rel;
                z_st_rel = z_touch_rel;
            else
                x_st_rel = 0.0;
                z_st_rel = -p.z_nom;
            end
            p_foot_des(:, i) = hip_i + [x_st_rel; 0.0; z_st_rel];
            v_foot_des(:, i) = [0.0; 0.0; 0.0];
            q_des_legs(:, i) = solve_leg_ik_2d(x_st_rel, z_st_rel, L1, L2);
            tau_swing(:, i)  = [0.0; 0.0];
        else
            % Swing leg: normalized swing progress s in [0, 1]
            s = (phase_legs(i) - beta) / (1.0 - beta);
            s = min(max(s, 0.0), 1.0);
            ds_dt = 1.0 / T_swing;

            % Cycloidal X progression (zero velocity & acceleration at s=0, s=1)
            sigma_x  = s - sin(2.0 * pi * s) / (2.0 * pi);
            dsigma_x = 1.0 - cos(2.0 * pi * s);
            x_rel_des  = x_lift_rel + (x_touch_rel - x_lift_rel) * sigma_x;
            vx_rel_des = (x_touch_rel - x_lift_rel) * dsigma_x * ds_dt;

            % Smooth cycloidal Z clearance trajectory with apex h_sw at s=0.5
            sigma_z  = 0.5 * (1.0 - cos(2.0 * pi * s));
            dsigma_z = pi * sin(2.0 * pi * s);
            z_base     = (1.0 - sigma_x) * z_lift_rel + sigma_x * z_touch_rel;
            dz_base_dt = (z_touch_rel - z_lift_rel) * dsigma_x * ds_dt;

            z_rel_des  = z_base + h_sw * sigma_z;
            vz_rel_des = dz_base_dt + h_sw * dsigma_z * ds_dt;

            p_des_i = hip_i + [x_rel_des; 0.0; z_rel_des];
            v_des_i = [vx_rel_des; 0.0; vz_rel_des];
            p_foot_des(:, i) = p_des_i;
            v_foot_des(:, i) = v_des_i;

            % Analytical 2-DOF Inverse Kinematics for desired swing posture
            q_des_i = solve_leg_ik_2d(x_rel_des, z_rel_des, L1, L2);
            q_des_legs(:, i) = q_des_i;

            % Desired joint velocity from 2D Jacobian inverse (damped against singularity)
            det_J2 = J_2d(1,1)*J_2d(2,2) - J_2d(1,2)*J_2d(2,1);
            if abs(det_J2) > 1.0e-3
                dq_des_i = J_2d \ [vx_rel_des; vz_rel_des];
            else
                dq_des_i = [0.0; 0.0];
            end

            % Operational-space Cartesian PD force: F_sw = Kp_sw * e_p + Kd_sw * e_v
            e_p = p_des_i - p_foot_3d;
            e_v = v_des_i - v_foot_3d;
            F_sw = p.Kp_swing_cart * e_p + p.Kd_swing_cart * e_v;
            tau_op_pd = J_3d.' * F_sw;

            % Auxiliary joint-space IK PD torque
            tau_joint_pd = p.Kp_swing_joint .* (q_des_i - q_i) + ...
                           p.Kd_swing_joint .* (dq_des_i - dq_i);

            % Exact 2-link leg gravity compensation torque: dV_leg / dq
            q1  = q_i(1);
            q12 = q_i(1) + q_i(2);
            tau_grav_1 = (0.5 * p.m_thigh * L1 + p.m_shank * L1) * p.g * sin(q1) + ...
                          0.5 * p.m_shank * L2 * p.g * sin(q12);
            tau_grav_2 =  0.5 * p.m_shank * L2 * p.g * sin(q12);
            tau_grav = [tau_grav_1; tau_grav_2];

            tau_raw = tau_op_pd + tau_joint_pd + tau_grav;
            tau_swing(1, i) = min(max(tau_raw(1), -p.tau_hip_max),  p.tau_hip_max);
            tau_swing(2, i) = min(max(tau_raw(2), -p.tau_knee_max), p.tau_knee_max);
        end
    end
end

function q_ik = solve_leg_ik_2d(x_rel, z_rel, L1, L2)
%SOLVE_LEG_IK_2D Analytical 2-DOF sagittal inverse kinematics (positive knee branch).
    r_sq = x_rel^2 + z_rel^2;
    cos_q2 = (r_sq - L1^2 - L2^2) / (2.0 * L1 * L2);
    cos_q2 = min(max(cos_q2, -0.995), 0.995);
    q2_des = acos(cos_q2);
    q1_des = atan2(-x_rel, -z_rel) - atan2(L2 * sin(q2_des), L1 + L2 * cos(q2_des));
    q_ik = [q1_des; q2_des];
end
