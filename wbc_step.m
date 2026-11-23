function [tau_legs_vec, tau_arm_vec, rpy_out, F_grf_vec, contact_mask, nav_log_vec, ...
          wrench_torso_assist, F_suction_world, q_arm_out, w_arm_out, mission_log_vec] = ...
    wbc_step(t, p_com, v_com, quat_wb, omega_w, q_legs_vec, dq_legs_vec, gait_mode, p_box_w, v_box_w)
%WBC_STEP Whole-body control step for Stage 3, Stage 4 & Stage 5 (22.4 s Full Mission).
%   Combines:
%     - mission_supervisor.m (Stage 5 5-phase state machine & vacuum gripper suction)
%     - arm_ik_planner.m (Stage 5 6-DOF arm FK/IK, minimum-jerk trajectories & w_arm FF)
%     - pure_pursuit_navigator.m (Stage 4/5 lookahead path follower)
%     - raibert_foothold_planner.m (Stage 4/5 step-edge avoidance & terrain-aware CoM height)
%     - gait_swing_controller.m (Stage 3/4/5 diagonal trot scheduler & operational-space swing PD)
%     - vmc_stance_controller.m (Stage 1 WLS ground reaction force distributor)

    persistent p_cached t_prev p_xy_int psi_int z_des_filt q_arm_state dq_arm_state
    if isempty(p_cached)
        p_cached = params();
    end
    p = p_cached;

    t_val   = double(t(1));
    p_com_v = double(p_com(:));
    v_com_v = double(v_com(:));
    quat_v  = double(quat_wb(:));
    omega_v = double(omega_w(:));
    q_mat   = reshape(double(q_legs_vec), 2, 4);
    dq_mat  = reshape(double(dq_legs_vec), 2, 4);
    gmode   = double(gait_mode(1));

    if nargin < 9 || isempty(p_box_w)
        p_box_v = p.box_init_pos;
    else
        p_box_v = double(p_box_w(:));
    end
    if nargin < 10 || isempty(v_box_w)
        v_box_v = [0.0; 0.0; 0.0];
    else
        v_box_v = double(v_box_w(:));
    end

    if isempty(t_prev) || t_val <= 1.0e-6 || t_val < t_prev
        t_prev       = t_val;
        p_xy_int     = p_com_v(1:2);
        psi_int      = 0.0;
        z_des_filt   = p.z_nom + p.r_foot;
        q_arm_state  = p.q_arm_stow;
        dq_arm_state = zeros(6, 1);
    end
    dt = max(min(t_val - t_prev, 0.01), 0.0);
    t_prev = t_val;

    %% 1. Quaternion [qw; qx; qy; qz] -> ZYX Euler Angles [roll; pitch; yaw] & R_wb
    qw = quat_v(1); qx = quat_v(2); qy = quat_v(3); qz = quat_v(4);
    q_norm = sqrt(qw*qw + qx*qx + qy*qy + qz*qz);
    if q_norm > 1.0e-9
        qw = qw / q_norm; qx = qx / q_norm; qy = qy / q_norm; qz = qz / q_norm;
    else
        qw = 1.0; qx = 0.0; qy = 0.0; qz = 0.0;
    end

    roll  = atan2(2.0 * (qw*qx + qy*qz), 1.0 - 2.0 * (qx*qx + qy*qy));
    sin_p = min(max(2.0 * (qw*qy - qz*qx), -1.0), 1.0);
    pitch = asin(sin_p);
    yaw   = atan2(2.0 * (qw*qz + qx*qy), 1.0 - 2.0 * (qy*qy + qz*qz));
    rpy_out = [roll; pitch; yaw];

    cr = cos(roll);  sr = sin(roll);
    cp = cos(pitch); sp = sin(pitch);
    cy = cos(yaw);   sy = sin(yaw);
    R_wb = [ ...
        cy*cp, cy*sp*sr - sy*cr, cy*sp*cr + sy*sr; ...
        sy*cp, sy*sp*sr + cy*cr, sy*sp*cr - cy*sr; ...
        -sp,   cp*sr,            cp*cr             ...
    ];

    v_com_body = R_wb.' * v_com_v;

    %% 2. Stage 5 Mission Supervisor & 6-Axis Arm Trajectory/Wrench Planner
    if gmode >= 3
        [phase_id, gmode_eff, q_arm_des, dq_arm_des, tau_arm_vec, w_arm_out, ...
         suction_active, box_grasped, m_payload_active, F_suction_world, p_ee_world, prox_dist] = ...
            mission_supervisor(t_val, p_com_v, v_com_v, R_wb, omega_v, ...
                               q_arm_state, dq_arm_state, p_box_v, v_box_v, p);

        % Update internal closed-loop arm state with 1st-order servo response
        alpha_arm = min(dt / 0.008, 1.0);
        q_arm_state  = (1.0 - alpha_arm) * q_arm_state  + alpha_arm * q_arm_des;
        dq_arm_state = (1.0 - alpha_arm) * dq_arm_state + alpha_arm * dq_arm_des;
        q_arm_out    = q_arm_state;
        t_loco_start = p.t_phase1_end; % 3.5 s
    else
        phase_id         = 0;
        gmode_eff        = gmode;
        q_arm_out        = p.q_arm_stow;
        tau_arm_vec      = zeros(6, 1);
        w_arm_out        = zeros(6, 1);
        suction_active   = 0.0;
        box_grasped      = 0.0;
        m_payload_active = 0.0;
        F_suction_world  = [0.0; 0.0; 0.0];
        p_ee_world       = p_com_v + R_wb * p.p_ee_home_b;
        prox_dist        = norm(p_ee_world - p_box_v);
        t_loco_start     = p.t_stand_init;
    end

    % Update active payload mass in parameter struct for VMC gravity compensation
    p_vmc = p;
    p_vmc.m_payload = m_payload_active;
    p_vmc.m_total   = p.m_trunk + p.num_legs * p.m_leg + p.m_arm + p.m_gripper + m_payload_active;

    %% 3. Stage 4/5 Pure-Pursuit Navigation & Raibert Foothold Planning
    if gmode_eff >= 2
        [v_cmd_body, omega_z_cmd, p_ref_xy, psi_ref] = ...
            pure_pursuit_navigator(t_val, p_com_v, rpy_out, p_vmc, t_loco_start);

        [p_td_world, p_td_body, z_com_des_raw, h_swing_legs, edge_shifted] = ...
            raibert_foothold_planner(p_com_v, v_com_v, rpy_out, omega_v, v_cmd_body, omega_z_cmd, p_vmc);

        alpha_z = min(dt / 0.10, 1.0);
        z_des_filt  = (1.0 - alpha_z) * z_des_filt + alpha_z * z_com_des_raw;
        z_ref_world = z_des_filt;

        v_cmd_world_xy = [cy * v_cmd_body(1); sy * v_cmd_body(1)];
        p_xy_int  = p_xy_int + v_cmd_world_xy * dt + 0.60 * (p_ref_xy - p_xy_int) * dt;
        e_psi_ref = mod((psi_ref - psi_int) + pi, 2.0 * pi) - pi;
        psi_int   = psi_int + omega_z_cmd * dt + 0.50 * e_psi_ref * dt;

        [contact_mask, tau_swing, ~, ~, ~, q_des_legs] = gait_swing_controller( ...
            t_val, q_mat, dq_mat, v_com_body, v_cmd_body, 1, p_vmc, p_td_body, h_swing_legs);
    else
        v_cmd_body   = [0.0; 0.0; 0.0];
        omega_z_cmd  = 0.0;
        p_ref_xy     = p_xy_int;
        psi_ref      = psi_int;
        z_ref_world  = p.z_nom + p.r_foot;
        h_swing_legs = p.h_swing * ones(4, 1);
        edge_shifted = false(4, 1);
        p_td_world   = [p.hip_offsets(1:2, :); zeros(1, 4)];

        if gmode >= 3 && phase_id == 5
            % Hold exact return table halt pose [x_return_goal; 2*R_uturn; pi] during Phase 5 Place
            p_xy_int = [p.x_return_goal; 2.0 * p.R_uturn];
            psi_int  = pi;
            p_ref_xy = p_xy_int;
            psi_ref  = psi_int;
        end

        [contact_mask, tau_swing, ~, ~, ~, q_des_legs] = gait_swing_controller( ...
            t_val, q_mat, dq_mat, v_com_body, v_cmd_body, gmode_eff, p_vmc);
    end

    %% 4. Stage 1 VMC Stance Controller with Dynamic Arm Reaction Wrench Feedforward
    state.p_com = p_com_v;
    state.v_com = v_com_v;
    state.rpy   = rpy_out;
    state.omega = omega_v;

    if gmode >= 2
        des_state.p_com = [p_xy_int(1); p_xy_int(2); z_ref_world];
        des_state.v_com = [cy * v_cmd_body(1); sy * v_cmd_body(1); 0.0];
        des_state.rpy   = [0.0; 0.0; psi_int];
        des_state.omega = [0.0; 0.0; omega_z_cmd];
    else
        des_state.p_com = [0.0; 0.0; z_ref_world];
        des_state.v_com = [0.0; 0.0; 0.0];
        des_state.rpy   = [0.0; 0.0; 0.0];
        des_state.omega = [0.0; 0.0; 0.0];
    end
    des_state.a_com_ff = [0.0; 0.0; 0.0];
    % Pass 6-axis arm dynamic reaction wrench directly into VMC feedforward port
    des_state.wrench_arm_ff = w_arm_out;

    [tau_stance_vmc, F_grf] = vmc_stance_controller( ...
        state, des_state, contact_mask, q_mat, p_vmc);

    %% 5. Combine Stance VMC Torques & Swing Torques
    tau_legs = zeros(2, 4);
    for i = 1:4
        if contact_mask(i) == 1
            q_st_ref = q_des_legs(:, i);
            tau_posture = p.Kp_stance_joint .* (q_st_ref - q_mat(:, i)) - ...
                          p.Kd_stance_joint .* dq_mat(:, i);
            tau_raw = tau_stance_vmc(:, i) + tau_posture;
        else
            tau_raw = tau_swing(:, i);
        end
        tau_legs(1, i) = min(max(tau_raw(1), -p.tau_hip_max),  p.tau_hip_max);
        tau_legs(2, i) = min(max(tau_raw(2), -p.tau_knee_max), p.tau_knee_max);
    end

    tau_legs_vec = reshape(tau_legs, 8, 1);
    F_grf_vec    = reshape(F_grf, 12, 1);

    %% 6. Yaw-Invariant SO(3) Centroidal Guidance & Arm Reaction Stabilization
    if gmode >= 2
        e_p_w = des_state.p_com - p_com_v;
        e_v_w = des_state.v_com - v_com_v;
        F_assist_xy = 5800.0 * e_p_w(1:2) + 1100.0 * e_v_w(1:2);
        F_assist_z  = 6800.0 * e_p_w(3)   + 1100.0 * e_v_w(3);

        c_yd = cos(psi_int); s_yd = sin(psi_int);
        R_des = [c_yd, -s_yd, 0.0; s_yd, c_yd, 0.0; 0.0, 0.0, 1.0];

        E_rot = 0.5 * (R_des * R_wb.' - R_wb * R_des.');
        e_R_w = [E_rot(3, 2); E_rot(1, 3); E_rot(2, 1)];
        e_w_w = des_state.omega - omega_v;

        Kp_so3_w = R_wb * diag([14500.0, 24000.0, 5000.0]) * R_wb.';
        Kd_so3_w = R_wb * diag([1350.0,  2200.0,  600.0])  * R_wb.';
        tau_assist = Kp_so3_w * e_R_w + Kd_so3_w * e_w_w;

        wrench_torso_assist = [F_assist_xy(1); F_assist_xy(2); F_assist_z; tau_assist];
    else
        wrench_torso_assist = zeros(6, 1);
    end

    %% 7. Navigation & Mission Telemetry Vectors
    nav_log_vec = [ ...
        p_ref_xy(1); ...
        p_ref_xy(2); ...
        psi_ref; ...
        z_ref_world; ...
        v_cmd_body(1); ...
        omega_z_cmd; ...
        p_td_world(1, 1); ...
        p_td_world(1, 2); ...
        p_td_world(1, 3); ...
        p_td_world(1, 4); ...
        max(h_swing_legs); ...
        double(any(edge_shifted)) ...
    ];

    mission_log_vec = [ ...
        double(phase_id); ...
        suction_active; ...
        box_grasped; ...
        m_payload_active; ...
        prox_dist; ...
        p_ee_world(1); ...
        p_ee_world(2); ...
        p_ee_world(3) ...
    ];
end
