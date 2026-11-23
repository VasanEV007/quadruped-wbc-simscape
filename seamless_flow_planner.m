function [p_ee_des_w, v_ee_des_w, a_ee_des_w, phase_id, gait_mode_active, suction_cmd, flow_info] = ...
    seamless_flow_planner(t, p_com_w, v_com_w, rpy, omega_w, p)
%SEAMLESS_FLOW_PLANNER C^2-continuous overlapped world-frame trajectory planner for 3-DOF arm.
%   [p_ee_des_w, v_ee_des_w, a_ee_des_w, phase_id, gait_mode_active, suction_cmd, flow_info] = ...
%       SEAMLESS_FLOW_PLANNER(t, p_com_w, v_com_w, rpy, omega_w, p)
%   generates C^2-continuous (position, velocity, acceleration) world-frame
%   end-effector trajectories with:
%     1. Overlapped Pre-Positioning:
%        - Overlaps post-pick carry tuck with outbound trot start (t in [2.65, 3.80] s).
%        - Begins extending and aligning the 3-DOF arm toward the place target
%          during the final 1.0 s of the return approach trot (t in [18.50, 19.50] s)
%          before the legs come to a halt.
%     2. General Hermite Quintic (C^2) Velocity Blending with non-zero boundary
%        velocities across phase knots so acceleration and jerk remain smooth.
%     3. World-Frame Task-Space Stabilization during trotting, step climbing,
%        and the 0.70 m U-turn (decoupled from body heave/pitch/roll bobbing).

    if nargin < 6 || isempty(p)
        p = params();
    end

    t_val   = double(t(1));
    p_com_v = double(p_com_w(:));
    v_com_v = double(v_com_w(:));
    rpy_v   = double(rpy(:));
    omega_v = double(omega_w(:));

    psi   = rpy_v(3);
    c_psi = cos(psi);
    s_psi = sin(psi);
    % Yaw-only horizontal frame R_yaw (decoupled from body roll/pitch and vertical heave!)
    R_yaw = [c_psi, -s_psi, 0.0; s_psi, c_psi, 0.0; 0.0, 0.0, 1.0];

    % Smooth terrain-aware reference base altitude (ignores high-frequency gait bobbing!)
    z_terrain_smooth = query_smooth_terrain_local(p_com_v(1), p_com_v(2), p);
    z_base_smooth    = p.z_nom + p.r_foot + z_terrain_smooth;
    p_base_smooth_w  = [p_com_v(1); p_com_v(2); z_base_smooth];
    v_base_smooth_w  = [v_com_v(1); v_com_v(2); 0.0];
    w_yaw_only       = [0.0; 0.0; omega_v(3)];

    %% 1. Key World-Frame & Yaw-Stabilized Waypoints
    % Fixed world-frame Pick and Place coordinates on the Pick-and-Place Table
    p_pick_top_w  = [p.box_init_pos(1);  p.box_init_pos(2);  p.table_height + p.box_side];
    p_place_top_w = [p.box_place_pos(1); p.box_place_pos(2); p.table_height + p.box_side];

    p_home_0_w  = [0.0; 0.0; p.z_nom + p.r_foot] + p.p_ee_home_b;
    p_hover1_w  = p_pick_top_w + [0.0; 0.0; 0.065];
    v_glide_in  = [0.01; 0.04; -0.055]; % Non-zero boundary velocity at hover entry
    v_glide_out = [-0.02; -0.06; 0.060]; % Non-zero boundary velocity at lift-to-tuck transition

    % World-stabilized carry pose riding smoothly with horizontal yaw frame (immune to heave/roll/pitch)
    p_carry_w = p_base_smooth_w + R_yaw * p.p_ee_carry_b;
    v_carry_w = v_base_smooth_w + cross(w_yaw_only, R_yaw * p.p_ee_carry_b);

    % Fixed world-frame Place hover and Place surface targets at return table station
    p_hover2_w  = p_place_top_w + [0.0; 0.0; 0.065];
    p_home_end_w = [p.x_return_goal; 2.0 * p.R_uturn; p.z_nom + p.r_foot] + ...
                   [-1, 0, 0; 0, -1, 0; 0, 0, 1] * p.p_ee_home_b;

    %% 2. Mission Phase & Gait Mode Schedule (with Overlapped Transitions)
    if t_val < p.t_phase1_end
        phase_id = 1;
        % Overlap: begin trot warm-up at t = 3.30 s while arm finishes carry blend
        if t_val < 3.35
            gait_mode_active = 0;
        else
            gait_mode_active = 2;
        end
    elseif t_val < p.t_phase2_end
        phase_id = 2;
        gait_mode_active = 2;
    elseif t_val < p.t_phase3_end
        phase_id = 3;
        gait_mode_active = 2;
    elseif t_val < p.t_phase4_end
        phase_id = 4;
        gait_mode_active = 2;
    else
        phase_id = 5;
        gait_mode_active = 0;
    end

    %% 3. Seamless C^2 Hermite Quintic Trajectory Evaluation Across All Segments
    suction_cmd = false;
    is_overlapped_reach = false;

    if t_val < 0.35
        % Segment 1a: Initial smooth hold at home
        [p_ee_des_w, v_ee_des_w, a_ee_des_w] = hermite_quintic_3d( ...
            t_val, 0.0, 0.35, p_home_0_w, [0;0;0], [0;0;0], p_home_0_w, [0;0;0], [0;0;0]);

    elseif t_val < 1.25
        % Segment 1b: Reach toward pre-pick hover with non-zero arrival velocity v_glide_in
        [p_ee_des_w, v_ee_des_w, a_ee_des_w] = hermite_quintic_3d( ...
            t_val, 0.35, 1.25, p_home_0_w, [0;0;0], [0;0;0], p_hover1_w, v_glide_in, [0;0;0]);

    elseif t_val < 1.90
        % Segment 1c: Continuous C^2 blend from hover velocity v_glide_in down to box top surface
        [p_ee_des_w, v_ee_des_w, a_ee_des_w] = hermite_quintic_3d( ...
            t_val, 1.25, 1.90, p_hover1_w, v_glide_in, [0;0;0], p_pick_top_w, [0;0;0], [0;0;0]);
        if t_val >= 1.78
            suction_cmd = true;
        end

    elseif t_val < 2.10
        % Segment 1d: Vacuum suction dwell on box top surface
        [p_ee_des_w, v_ee_des_w, a_ee_des_w] = hermite_quintic_3d( ...
            t_val, 1.90, 2.10, p_pick_top_w, [0;0;0], [0;0;0], p_pick_top_w, [0;0;0], [0;0;0]);
        suction_cmd = true;

    elseif t_val < 2.75
        % Segment 1e: Lift box through hover with non-zero upward/tuck velocity v_glide_out
        [p_ee_des_w, v_ee_des_w, a_ee_des_w] = hermite_quintic_3d( ...
            t_val, 2.10, 2.75, p_pick_top_w, [0;0;0], [0;0;0], p_hover1_w, v_glide_out, [0;0;0]);
        suction_cmd = true;

    elseif t_val < 3.65
        % Segment 1f (Overlapped Pick->Trot Blend): C^2 blend from v_glide_out into moving carry pose
        [p_ee_des_w, v_ee_des_w, a_ee_des_w] = hermite_quintic_3d( ...
            t_val, 2.75, 3.65, p_hover1_w, v_glide_out, [0;0;0], p_carry_w, v_carry_w, [0;0;0]);
        suction_cmd = true;
        is_overlapped_reach = true;

    elseif t_val < 18.50
        % Segments 2-4a (World-Stabilized Carry + Independent 3-DOF Spatial Float Wave):
        % While the quadruped climbs the 10 cm step and executes the R=0.70 m U-turn,
        % the 3-DOF arm maintains a world-stabilized altitude and executes a smooth
        % C^2 independent 3-DOF task-space trajectory decoupled from gait bobbing!
        suction_cmd = true;
        t_float = t_val - 3.65;
        T_float = 18.50 - 3.65;
        s_env   = min(max(t_float / T_float, 0.0), 1.0);
        w_f     = 2.0 * pi / 3.0;

        % sin^2(pi*s) envelope has zero value AND zero derivative at t=3.65s and t=18.50s
        env    = sin(pi * s_env)^2;
        denv   = (pi / T_float) * sin(2.0 * pi * s_env);

        wave_p = [0.025 * sin(w_f * t_float); 0.035 * sin(0.5 * w_f * t_float); 0.020 * (1.0 - cos(w_f * t_float))];
        wave_v = [0.025 * w_f * cos(w_f * t_float); 0.035 * 0.5 * w_f * cos(0.5 * w_f * t_float); 0.020 * w_f * sin(w_f * t_float)];
        wave_a = [-0.025 * w_f^2 * sin(w_f * t_float); -0.035 * 0.25 * w_f^2 * sin(0.5 * w_f * t_float); 0.020 * w_f^2 * cos(w_f * t_float)];

        dp_wave_b = env * wave_p;
        dv_wave_b = denv * wave_p + env * wave_v;
        da_wave_b = env * wave_a;

        p_ee_des_w = p_carry_w + R_yaw * dp_wave_b;
        v_ee_des_w = v_carry_w + R_yaw * dv_wave_b + cross(w_yaw_only, R_yaw * dp_wave_b);
        a_ee_des_w = R_yaw * da_wave_b;

    elseif t_val < 19.75
        % Segment 4b/5a (OVERLAPPED PRE-POSITIONING DURING FINAL 1.0 s OF APPROACH TROT):
        % Smoothly blends from the live world-stabilized carry state (p_carry_w, v_carry_w)
        % into the fixed world-frame pre-place hover target (p_hover2_w, v_preplace_out)!
        suction_cmd = true;
        is_overlapped_reach = true;
        v_preplace_out = [-0.015; -0.030; -0.045];

        [p_ee_des_w, v_ee_des_w, a_ee_des_w] = hermite_quintic_3d( ...
            t_val, 18.50, 19.75, p_carry_w, v_carry_w, [0;0;0], p_hover2_w, v_preplace_out, [0;0;0]);

    elseif t_val < 20.85
        % Segment 5b: C^2 continuous descent from pre-place hover velocity onto table surface
        suction_cmd = true;
        v_preplace_out = [-0.015; -0.030; -0.045];
        [p_ee_des_w, v_ee_des_w, a_ee_des_w] = hermite_quintic_3d( ...
            t_val, 19.75, 20.85, p_hover2_w, v_preplace_out, [0;0;0], p_place_top_w, [0;0;0], [0;0;0]);

    elseif t_val < 21.10
        % Segment 5c: Release vacuum suction on table surface
        suction_cmd = false;
        [p_ee_des_w, v_ee_des_w, a_ee_des_w] = hermite_quintic_3d( ...
            t_val, 20.85, 21.10, p_place_top_w, [0;0;0], [0;0;0], p_place_top_w, [0;0;0], [0;0;0]);

    elseif t_val < 21.65
        % Segment 5d: Smooth lift-off after release with non-zero retract velocity
        suction_cmd = false;
        v_ret_mid = [0.02; 0.08; 0.05];
        [p_ee_des_w, v_ee_des_w, a_ee_des_w] = hermite_quintic_3d( ...
            t_val, 21.10, 21.65, p_place_top_w, [0;0;0], [0;0;0], p_hover2_w, v_ret_mid, [0;0;0]);

    else
        % Segment 5e: C^2 blend from retract velocity back to stowed home pose
        suction_cmd = false;
        v_ret_mid = [0.02; 0.08; 0.05];
        [p_ee_des_w, v_ee_des_w, a_ee_des_w] = hermite_quintic_3d( ...
            t_val, 21.65, 22.40, p_hover2_w, v_ret_mid, [0;0;0], p_home_end_w, [0;0;0], [0;0;0]);
    end

    if nargout >= 7
        flow_info.is_overlapped_reach = is_overlapped_reach;
        flow_info.z_base_smooth       = z_base_smooth;
    end
end

%% ---------------- General C^2 Hermite Quintic Spline Solver ----------------

function [p_t, v_t, a_t] = hermite_quintic_3d(t, t0, t1, p0, v0, a0, p1, v1, a1)
%HERMITE_QUINTIC_3D Evaluate 3D C^2 quintic spline with arbitrary boundary velocities & accels.
    T = max(t1 - t0, 1.0e-3);
    tau = min(max((t - t0) / T, 0.0), 1.0);

    tau2 = tau * tau;
    tau3 = tau2 * tau;
    tau4 = tau3 * tau;
    tau5 = tau4 * tau;

    c0 = p0;
    c1 = v0 * T;
    c2 = 0.5 * a0 * (T * T);
    c3 = 10.0 * (p1 - p0) - (6.0 * v0 + 4.0 * v1) * T + 0.5 * (a1 - 3.0 * a0) * (T * T);
    c4 = -15.0 * (p1 - p0) + (8.0 * v0 + 7.0 * v1) * T + 0.5 * (3.0 * a0 - 2.0 * a1) * (T * T);
    c5 = 6.0 * (p1 - p0) - 3.0 * (v0 + v1) * T + 0.5 * (a1 - a0) * (T * T);

    p_t = c0 + c1*tau + c2*tau2 + c3*tau3 + c4*tau4 + c5*tau5;
    v_t = (c1 + 2.0*c2*tau + 3.0*c3*tau2 + 4.0*c4*tau3 + 5.0*c5*tau4) / T;
    a_t = (2.0*c2 + 6.0*c3*tau + 12.0*c4*tau2 + 20.0*c5*tau3) / (T * T);
end

function z_s = query_smooth_terrain_local(x_w, y_w, p)
    if y_w < p.y_step_min || y_w > p.y_step_max
        z_s = 0.0;
        return;
    end
    w_blend = 0.24;
    s_rise = 0.5 * (1.0 + tanh((x_w - p.x_step_start) / w_blend));
    s_fall = 0.5 * (1.0 - tanh((x_w - p.x_step_end)   / w_blend));
    z_s = p.h_step * s_rise * s_fall;
end
