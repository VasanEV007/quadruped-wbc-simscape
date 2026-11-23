function [phase_id, gait_mode_active, q_arm_des, dq_arm_des, tau_arm_cmd, w_arm_world, ...
          suction_active, box_grasped, m_payload_active, F_suction_world, p_ee_world, prox_dist, ee_err_world_m] = ...
    mission_supervisor(t, p_com, v_com, R_wb, omega_w, q_arm_meas, dq_arm_meas, p_box_w, v_box_w, p)
%MISSION_SUPERVISOR Stage 5/6 22.4 s autonomous mission supervisor with 3-DOF decoupled arm flow.
%   Combines seamless_flow_planner.m (C^2 Hermite quintic trajectories & overlapped
%   pre-positioning) with arm_3dof_controller.m (decoupled 3-DOF analytical IK,
%   active base-twist velocity cancellation, and dynamic reaction wrench w_arm).

    if nargin < 10 || isempty(p)
        p = params();
    end

    persistent latched_grasp
    if isempty(latched_grasp) || t <= 1.0e-5
        latched_grasp = false;
    end

    t_val    = double(t(1));
    p_com_v  = double(p_com(:));
    v_com_v  = double(v_com(:));
    omega_v  = double(omega_w(:));
    q_arm_v  = double(q_arm_meas(:));
    dq_arm_v = double(dq_arm_meas(:));
    p_box_v  = double(p_box_w(:));
    v_box_v  = double(v_box_w(:));

    % Extract [roll; pitch; yaw] from R_wb for seamless_flow_planner
    pitch = asin(min(max(-R_wb(3, 1), -1.0), 1.0));
    roll  = atan2(R_wb(3, 2), R_wb(3, 3));
    yaw   = atan2(R_wb(2, 1), R_wb(1, 1));
    rpy_v = [roll; pitch; yaw];

    %% 1. Stage 6 Seamless C^2 Trajectory Planner (with Overlapped Pre-Positioning)
    [p_ee_des_w, v_ee_des_w, a_ee_des_w, phase_id, gait_mode_active, suction_cmd] = ...
        seamless_flow_planner(t_val, p_com_v, v_com_v, rpy_v, omega_v, p);

    if latched_grasp && suction_cmd
        m_payload_curr = p.m_box;
    else
        m_payload_curr = 0.0;
    end

    %% 2. Stage 6 Decoupled 3-DOF Analytical IK, Base-Twist Compensation & Computed Torque
    [q_arm_des, dq_arm_des, ~, tau_arm_cmd, w_arm_world, p_ee_world_act, ee_err_world_m] = ...
        arm_3dof_controller(p_ee_des_w, v_ee_des_w, a_ee_des_w, ...
                            p_com_v, v_com_v, R_wb, omega_v, ...
                            q_arm_v, dq_arm_v, m_payload_curr, p);

    p_ee_world = p_ee_world_act;

    %% 3. Vacuum Gripper Proximity Check (< 0.015 m) & Suction Weld Force
    p_box_top_w = p_box_v + [0.0; 0.0; 0.5 * p.box_side];
    prox_dist   = norm(p_ee_world - p_box_top_w);

    if suction_cmd
        if ~latched_grasp && (prox_dist < p.suction_prox_thresh)
            latched_grasp = true;
        end
    else
        latched_grasp = false;
    end

    suction_active = double(suction_cmd);
    box_grasped    = double(latched_grasp);

    if latched_grasp
        m_payload_active = p.m_box;
        % Use smooth C^2 world-frame reference p_ee_des_w for suction attraction
        % to avoid 1 kHz ZOH staircase steps on the 0.5 kg box DAE solver
        p_box_target_w = p_ee_des_w - [0.0; 0.0; 0.5 * p.box_side];
        e_box_p = p_box_target_w - p_box_v;
        e_box_v = v_ee_des_w - v_box_v;

        F_raw = 2800.0 * e_box_p + 150.0 * e_box_v + ...
                [0.0; 0.0; p.m_box * p.g];
        F_mag = norm(F_raw);
        if F_mag > p.vacuum_suction_max_force
            F_raw = F_raw * (p.vacuum_suction_max_force / F_mag);
        end
        F_suction_world = F_raw;
    else
        m_payload_active = 0.0;
        F_suction_world  = [0.0; 0.0; 0.0];
    end
end
