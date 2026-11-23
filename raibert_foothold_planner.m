function [p_td_world, p_td_body, z_com_des, h_swing_legs, edge_shifted, z_hips_terrain] = ...
    raibert_foothold_planner(p_com, v_com, rpy, omega_w, v_cmd_body, omega_z_cmd, p)
%RAIBERT_FOOTHOLD_PLANNER Raibert foothold placement, step-edge avoidance & terrain CoM height.
%   [p_td_world, p_td_body, z_com_des, h_swing_legs, edge_shifted, z_hips_terrain] = ...
%       RAIBERT_FOOTHOLD_PLANNER(p_com, v_com, rpy, omega_w, v_cmd_body, omega_z_cmd, p)
%   computes target touchdown coordinates (x_td, y_td, z_td) in both world and
%   body frames for all 4 legs [FL, FR, RL, RR], enforces a safety margin
%   (>= 0.05 m) away from the 10 cm step edges (x_step_start, x_step_end),
%   computes the terrain-aware VMC desired CoM height z_com_des, and elevates
%   swing apex clearance (h_swing >= 0.14 m) during step climb/descent.

    if nargin < 7 || isempty(p)
        p = params();
    end

    p_com_v = double(p_com(:));
    v_com_v = double(v_com(:));
    rpy_v   = double(rpy(:));
    omega_v = double(omega_w(:));
    v_cmd_b = double(v_cmd_body(:));
    wz_cmd  = double(omega_z_cmd(1));

    psi   = rpy_v(3);
    c_psi = cos(psi);
    s_psi = sin(psi);
    Rz_2d = [c_psi, -s_psi; s_psi, c_psi];

    % Convert world CoM velocity to horizontal body frame
    v_com_xy_b = Rz_2d.' * v_com_v(1:2);
    vx_body = v_com_xy_b(1);
    wz_act  = omega_v(3);

    % 1. Turn-rate scaling factor: automatically scales down step length as |omega_z| increases
    wz_eff = max(abs(wz_act), abs(wz_cmd));
    turn_scale = 1.0 / (1.0 + p.c_turn_scale * wz_eff);

    T_stance = p.duty_factor * p.T_gait;
    d_margin = max(p.d_edge_margin, 0.055); % Strictly >= 0.05 m safety margin

    p_td_world     = zeros(3, 4);
    p_td_body      = zeros(3, 4);
    h_swing_legs   = p.h_swing * ones(4, 1);
    edge_shifted   = false(4, 1);
    z_hips_terrain = zeros(4, 1);
    z_hips_smooth  = zeros(4, 1);

    %% 2. Evaluate Hip Projections & Smooth Terrain Elevation Beneath Hips
    p_hips_xy_w = zeros(2, 4);
    for i = 1:4
        hip_b = p.hip_offsets(:, i);
        p_hips_xy_w(:, i) = p_com_v(1:2) + Rz_2d * hip_b(1:2);
        z_hips_terrain(i) = query_step_terrain(p_hips_xy_w(1, i), p_hips_xy_w(2, i), p);
        z_hips_smooth(i)  = query_smooth_terrain(p_hips_xy_w(1, i), p_hips_xy_w(2, i), p);
    end

    % Terrain-aware desired CoM height in world frame: z_nom + r_foot above average terrain
    z_terrain_avg = mean(z_hips_smooth);
    z_com_des = p.z_nom + p.r_foot + z_terrain_avg;

    %% 3. Raibert Foothold Computation & Step-Edge Clearance per Leg
    edges_x = [p.x_step_start, p.x_step_end];

    for i = 1:4
        hip_b = p.hip_offsets(:, i);
        hip_w = p_hips_xy_w(:, i);

        % Differential left/right longitudinal velocity for yaw turning without abduction
        vx_leg_act = vx_body - hip_b(2) * wz_act;
        vx_leg_cmd = v_cmd_b(1) - hip_b(2) * wz_cmd;

        % Raibert velocity heuristic along sagittal axis +X_B
        dx_raibert_b = turn_scale * (0.5 * T_stance * vx_leg_act + ...
                                     p.Kv_raibert * (vx_leg_act - vx_leg_cmd));
        dx_raibert_b = min(max(dx_raibert_b, -p.max_step_x), p.max_step_x);

        % Candidate touchdown XY in world frame (constrained to leg sagittal plane y_rel_b = 0)
        p_cand_xy_w = hip_w + Rz_2d * [dx_raibert_b; 0.0];
        x_td_w = p_cand_xy_w(1);
        y_td_w = p_cand_xy_w(2);

        % Enforce strict step-edge clearance |x_td_w - x_edge| >= d_margin (>= 0.055 m)
        in_y_span = (y_td_w >= p.y_step_min) && (y_td_w <= p.y_step_max);
        if in_y_span
            for e_idx = 1:2
                x_edge = edges_x(e_idx);
                dist_to_edge = x_td_w - x_edge;

                if abs(dist_to_edge) < d_margin
                    edge_shifted(i) = true;
                    x_minus = x_edge - d_margin;
                    x_plus  = x_edge + d_margin;

                    if dist_to_edge >= 0.0
                        x_target_w = x_plus;
                    else
                        x_target_w = x_minus;
                    end

                    if abs(c_psi) > 0.20
                        dx_shifted_b = (x_target_w - hip_w(1)) / c_psi;
                        if abs(dx_shifted_b) > p.max_step_x + 0.025
                            if x_target_w == x_plus
                                x_target_w = x_minus;
                            else
                                x_target_w = x_plus;
                            end
                            dx_shifted_b = (x_target_w - hip_w(1)) / c_psi;
                        end
                        dx_raibert_b = min(max(dx_shifted_b, -p.max_step_x - 0.025), p.max_step_x + 0.025);
                        p_cand_xy_w  = hip_w + Rz_2d * [dx_raibert_b; 0.0];
                        x_td_w = p_cand_xy_w(1);
                        y_td_w = p_cand_xy_w(2);
                    end

                    % Guarantee strict world-X clearance even if reach-clamped
                    if abs(x_td_w - x_edge) < d_margin
                        if x_td_w >= x_edge
                            x_td_w = x_plus;
                        else
                            x_td_w = x_minus;
                        end
                        dp_b = Rz_2d.' * ([x_td_w; y_td_w] - hip_w);
                        dx_raibert_b = min(max(dp_b(1), -p.max_step_x - 0.03), p.max_step_x + 0.03);
                    end
                end
            end
        end

        % Target touchdown height z_td from terrain map
        z_td_w = query_step_terrain(x_td_w, y_td_w, p);
        p_td_world(:, i) = [x_td_w; y_td_w; z_td_w];

        % Target touchdown in body frame relative to CoM
        z_rel_b = (z_td_w + p.r_foot) - z_com_des;
        z_rel_b = min(max(z_rel_b, -0.56), -0.38);
        p_td_body(:, i) = [hip_b(1) + dx_raibert_b; hip_b(2); z_rel_b];

        % Elevate swing apex clearance to h_swing_step (>= 0.14 m) when near or crossing step edges
        dist_start = min(abs(x_td_w - p.x_step_start), abs(hip_w(1) - p.x_step_start));
        dist_end   = min(abs(x_td_w - p.x_step_end),   abs(hip_w(1) - p.x_step_end));
        if (dist_start < 0.35 || dist_end < 0.35 || z_td_w ~= z_hips_terrain(i)) && in_y_span
            h_swing_legs(i) = max(p.h_swing_step, 0.14);
        else
            h_swing_legs(i) = p.h_swing;
        end
    end
end

%% ---------------- Local Terrain Helper Functions ----------------

function z_t = query_step_terrain(x_w, y_w, p)
%QUERY_STEP_TERRAIN Exact height of the 10 cm step platform at world (x_w, y_w).
    if x_w >= p.x_step_start && x_w <= p.x_step_end && ...
       y_w >= p.y_step_min   && y_w <= p.y_step_max
        z_t = p.h_step;
    else
        z_t = 0.0;
    end
end

function z_s = query_smooth_terrain(x_w, y_w, p)
%QUERY_SMOOTH_TERRAIN Smoothly blended terrain elevation for C^1 CoM height reference.
    if y_w < p.y_step_min || y_w > p.y_step_max
        z_s = 0.0;
        return;
    end
    w_blend = 0.22; % Smooth ramp half-width across step edges [m]
    s_rise = 0.5 * (1.0 + tanh((x_w - p.x_step_start) / w_blend));
    s_fall = 0.5 * (1.0 - tanh((x_w - p.x_step_end)   / w_blend));
    z_s = p.h_step * s_rise * s_fall;
end
