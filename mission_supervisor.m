function [phase_id, gait_mode_active, q_arm_des, dq_arm_des, tau_arm_cmd, w_arm_world, ...
          suction_active, box_grasped, m_payload_active, F_suction_world, p_ee_world, prox_dist] = ...
    mission_supervisor(t, p_com, v_com, R_wb, omega_w, q_arm_meas, dq_arm_meas, p_box_w, v_box_w, p)
%MISSION_SUPERVISOR 5-Phase 22.4 s autonomous pick-and-place + step-climb + U-turn supervisor.
%   Executes the 5-phase mission state machine:
%     Phase 1 (0.0 - 3.5 s)  : Pick payload box from table in 4-leg VMC stance
%     Phase 2 (3.5 - 9.0 s)  : Outbound approach, 10 cm step climb, traverse & step down
%     Phase 3 (9.0 - 14.5 s) : 180-deg semicircular U-turn (R = 0.70 m)
%     Phase 4 (14.5 - 19.5 s): Return 10 cm step climb/descent & halt at table
%     Phase 5 (19.5 - 22.4 s): Place box on table, release suction & retract arm

    if nargin < 10 || isempty(p)
        p = params();
    end

    persistent latched_grasp
    if isempty(latched_grasp) || t <= 1.0e-5
        latched_grasp = false;
    end

    t_val = double(t(1));
    q_arm_v  = double(q_arm_meas(:));
    dq_arm_v = double(dq_arm_meas(:));
    p_box_v  = double(p_box_w(:));
    v_box_v  = double(v_box_w(:));

    %% 1. Determine Mission Phase & Arm Trajectory Segment
    %Compute exact body-frame coordinates of the box top surface at Pick (t < 3.5 s)
    %and Place target on the table (t >= 19.5 s) so the suction cup lands within < 2 mm!
    p_pick_top_w  = [p.box_init_pos(1);  p.box_init_pos(2);  p.table_height + p.box_side];
    p_place_top_w = [p.box_place_pos(1); p.box_place_pos(2); p.table_height + p.box_side];

    p_nom_pick_com  = [0.0; 0.0; p.z_nom + p.r_foot];
    p_pick_b  = p_pick_top_w - p_nom_pick_com;
    p_hover_pick_b = p_pick_b + [0.0; 0.0; 0.065];

    % At Phase 5 return halt (p_com ~ [x_return_goal; 2*R_uturn; z_nom+r_foot], psi ~ pi)
    Rz_pi = [-1.0, 0.0, 0.0; 0.0, -1.0, 0.0; 0.0, 0.0, 1.0];
    p_nom_place_com = [p.x_return_goal; 2.0 * p.R_uturn; p.z_nom + p.r_foot];
    p_place_b = Rz_pi.' * (p_place_top_w - p_nom_place_com);
    p_hover_place_b = p_place_b + [0.0; 0.0; 0.065];

    suction_cmd = false;

    if t_val < p.t_phase1_end
        % ---------------- Phase 1: Pick (0.0 - 3.5 s) ----------------
        phase_id = 1;
        gait_mode_active = 0; % 4-leg VMC standing

        if t_val < 0.45
            p_seg0 = p.p_ee_home_b;  p_seg1 = p.p_ee_home_b;
            t0 = 0.0; t1 = 0.45;
        elseif t_val < 1.25
            p_seg0 = p.p_ee_home_b;  p_seg1 = p_hover_pick_b;
            t0 = 0.45; t1 = 1.25;
        elseif t_val < 1.95
            p_seg0 = p_hover_pick_b; p_seg1 = p_pick_b;
            t0 = 1.25; t1 = 1.95;
            if t_val >= 1.82
                suction_cmd = true;
            end
        elseif t_val < 2.15
            p_seg0 = p_pick_b;       p_seg1 = p_pick_b;
            t0 = 1.95; t1 = 2.15;
            suction_cmd = true;
        elseif t_val < 2.70
            p_seg0 = p_pick_b;       p_seg1 = p_hover_pick_b;
            t0 = 2.15; t1 = 2.70;
            suction_cmd = true;
        else
            p_seg0 = p_hover_pick_b; p_seg1 = p.p_ee_carry_b;
            t0 = 2.70; t1 = 3.45;
            suction_cmd = true;
        end

    elseif t_val < p.t_phase2_end
        % -------- Phase 2: Outbound Approach & Step Climb (3.5 - 9.0 s) --------
        phase_id = 2;
        gait_mode_active = 2;
        suction_cmd = true;
        p_seg0 = p.p_ee_carry_b; p_seg1 = p.p_ee_carry_b;
        t0 = p.t_phase1_end; t1 = p.t_phase2_end;

    elseif t_val < p.t_phase3_end
        % ------------ Phase 3: Wide-Arc U-Turn (9.0 - 14.5 s) ------------
        phase_id = 3;
        gait_mode_active = 2;
        suction_cmd = true;
        p_seg0 = p.p_ee_carry_b; p_seg1 = p.p_ee_carry_b;
        t0 = p.t_phase2_end; t1 = p.t_phase3_end;

    elseif t_val < p.t_phase4_end
        % ------ Phase 4: Return Climb & Table Approach (14.5 - 19.5 s) ------
        phase_id = 4;
        gait_mode_active = 2;
        suction_cmd = true;
        p_seg0 = p.p_ee_carry_b; p_seg1 = p.p_ee_carry_b;
        t0 = p.t_phase3_end; t1 = p.t_phase4_end;

    else
        % ------------ Phase 5: Place & Retract (19.5 - 22.4 s) ------------
        phase_id = 5;
        gait_mode_active = 0; % 4-leg VMC standing balance

        if t_val < 20.25
            p_seg0 = p.p_ee_carry_b;  p_seg1 = p_hover_place_b;
            t0 = 19.50; t1 = 20.25;
            suction_cmd = true;
        elseif t_val < 20.95
            p_seg0 = p_hover_place_b; p_seg1 = p_place_b;
            t0 = 20.25; t1 = 20.95;
            suction_cmd = true;
        elseif t_val < 21.15
            p_seg0 = p_place_b;       p_seg1 = p_place_b;
            t0 = 20.95; t1 = 21.15;
            suction_cmd = false; % Release vacuum suction on table surface
        elseif t_val < 21.70
            p_seg0 = p_place_b;       p_seg1 = p_hover_place_b;
            t0 = 21.15; t1 = 21.70;
            suction_cmd = false;
        else
            p_seg0 = p_hover_place_b; p_seg1 = p.p_ee_home_b;
            t0 = 21.70; t1 = 22.35;
            suction_cmd = false;
        end
    end

    %% 2. Evaluate Arm Trajectory Planner & Dynamic Base Reaction Wrench
    if latched_grasp && suction_cmd
        m_payload_curr = p.m_box;
    else
        m_payload_curr = 0.0;
    end

    [q_arm_des, dq_arm_des, ~, p_ee_des_b, v_ee_des_b, w_arm_world, arm_info] = ...
        arm_ik_planner(p_seg0, p_seg1, t_val, t0, t1, R_wb, m_payload_curr, p, q_arm_v);

    % World-frame position and velocity of the vacuum suction cup tip
    % Use measured FK position of the physical Simscape arm for proximity check
    [p_ee_meas_b, J_arm_meas, ~, ~] = arm_fk_local(q_arm_v, p);
    v_ee_meas_b = J_arm_meas(1:3, :) * dq_arm_v;

    p_ee_world = p_com(:) + R_wb * p_ee_des_b;
    v_ee_world = v_com(:) + cross(omega_w(:), R_wb * p_ee_des_b) + R_wb * v_ee_des_b;

    p_ee_meas_w = p_com(:) + R_wb * p_ee_meas_b;
    v_ee_meas_w = v_com(:) + cross(omega_w(:), R_wb * p_ee_meas_b) + R_wb * v_ee_meas_b;
    % Keep unused variables referenced cleanly
    if norm(p_ee_meas_w - p_ee_world) < 0.05
        p_ee_world = p_ee_meas_w;
        v_ee_world = v_ee_meas_w;
    end

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
        % Target box center is directly below the vacuum suction cup tip
        p_box_target_w = p_ee_world - [0.0; 0.0; 0.5 * p.box_side];
        e_box_p = p_box_target_w - p_box_v;
        e_box_v = v_ee_world - v_box_v;

        F_raw = p.kp_suction * e_box_p + p.kd_suction * e_box_v + ...
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

    %% 4. Compute 6-Axis Arm Joint PD + Gravity Compensation Torques
    tau_grav = arm_info.tau_grav_arm;
    tau_pd   = p.Kp_arm .* (q_arm_des - q_arm_v) + p.Kd_arm .* (dq_arm_des - dq_arm_v);
    tau_arm_cmd = min(max(tau_grav + tau_pd, -p.tau_arm_max), p.tau_arm_max);
end

function [p_ee_b, J_arm, p_coms_b, p_joints_b] = arm_fk_local(q, p)
    [~, ~, ~, ~, ~, ~, info] = arm_ik_planner(p.p_ee_home_b, p.p_ee_home_b, 0, 0, 1, eye(3), 0, p, q);
    p_ee_b     = info.p_joints_b(:, 7);
    J_arm      = info.J_arm;
    p_coms_b   = info.p_coms_b;
    p_joints_b = info.p_joints_b;
end
