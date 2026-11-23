function [v_cmd_body, omega_z_cmd, p_ref_xy, psi_ref, p_look_xy, path_xy, nav_info] = ...
    pure_pursuit_navigator(t, p_com, rpy, p, t_loco_start)
%PURE_PURSUIT_NAVIGATOR Lookahead pure-pursuit path follower for step climb & 0.7 m U-turn.
%   [v_cmd_body, omega_z_cmd, p_ref_xy, psi_ref, p_look_xy, path_xy, nav_info] = ...
%       PURE_PURSUIT_NAVIGATOR(t, p_com, rpy, p, t_loco_start)
%   generates smooth forward velocity v_x_cmd and yaw rate omega_z_cmd commands
%   (with zero lateral body velocity command v_y_cmd = 0 to avoid lateral skidding)
%   along the 3-stage course:
%     1. Straight approach along +X (y = 0) across the 10 cm step (x in [1.20, 1.85] m),
%     2. 180-deg semicircular U-turn arc of radius R = 0.70 m (y: 0 -> 1.40 m),
%     3. Straight return along -X (y = 1.40 m) back over the 10 cm step toward the table.

    if nargin < 4 || isempty(p)
        p = params();
    end
    if nargin < 5 || isempty(t_loco_start)
        t_loco_start = p.t_stand_init;
    end

    % Build and cache dense reference path (ds = 0.01 m)
    persistent cached_path
    if isempty(cached_path) || cached_path.R_uturn ~= p.R_uturn || ...
       cached_path.x_uturn ~= p.x_uturn_center
        cached_path = build_uturn_course(p);
    end

    path_xy   = cached_path.xy;
    path_psi  = cached_path.psi;
    path_s    = cached_path.s;
    path_kap  = cached_path.kappa;
    s_total   = cached_path.s_total;
    N_pts     = size(path_xy, 2);

    p_xy = double(p_com(1:2));
    p_xy = p_xy(:);
    psi  = double(rpy(3));

    %% 1. Closest Point Projection Along Reference Course
    best_cost = inf;
    idx_closest = 1;
    for k = 1:N_pts
        dx = path_xy(1, k) - p_xy(1);
        dy = path_xy(2, k) - p_xy(2);
        dpsi = wrap_to_pi_local(path_psi(k) - psi);
        cost_k = (dx*dx + dy*dy) + 0.15 * (dpsi*dpsi);
        if cost_k < best_cost
            best_cost = cost_k;
            idx_closest = k;
        end
    end

    s_curr   = path_s(idx_closest);
    p_ref_xy = path_xy(:, idx_closest);
    psi_ref  = path_psi(idx_closest);

    %% 2. Lookahead Target Point Selection (L_d in [0.35, 0.45] m)
    L_d = p.L_lookahead;
    s_look = min(s_curr + L_d, s_total);

    idx_look = idx_closest;
    while idx_look < N_pts && path_s(idx_look) < s_look
        idx_look = idx_look + 1;
    end
    p_look_xy = path_xy(:, idx_look);
    psi_look  = path_psi(idx_look);

    %% 3. Velocity Profile with Step Climb/Descent Slowdown & Goal Arrival
    if t < t_loco_start
        vx_cmd = 0.0;
    else
        ramp_up = min((t - t_loco_start) / 0.35, 1.0);

        % Slow down when approaching or crossing the 10 cm step edges
        dist_edge_start = abs(p_xy(1) - p.x_step_start);
        dist_edge_end   = abs(p_xy(1) - p.x_step_end);
        dist_edge = min(dist_edge_start, dist_edge_end);

        if dist_edge < 0.38
            alpha_step = dist_edge / 0.38;
            v_target = p.v_step_cross + (p.v_cruise - p.v_step_cross) * alpha_step;
        else
            v_target = p.v_cruise;
        end

        % Smooth arrival deceleration near the end of the return path
        s_rem = max(s_total - s_curr, 0.0);
        if s_rem < 0.45
            ramp_down = s_rem / 0.45;
        else
            ramp_down = 1.0;
        end
        if s_rem < 0.02
            ramp_down = 0.0;
        end

        vx_cmd = v_target * ramp_up * ramp_down;
    end

    %% 4. Lookahead Pure-Pursuit Steering Law
    c_psi = cos(psi);
    s_psi = sin(psi);
    dp_w  = p_look_xy - p_xy;

    dx_b =  c_psi * dp_w(1) + s_psi * dp_w(2);
    dy_b = -s_psi * dp_w(1) + c_psi * dp_w(2);

    L_eff_sq = max(dx_b*dx_b + dy_b*dy_b, 0.04);
    kappa_pp = (2.0 * dy_b) / L_eff_sq;

    e_psi = wrap_to_pi_local(psi_look - psi);
    if vx_cmd > 1.0e-3
        omega_raw = vx_cmd * kappa_pp + 1.35 * e_psi;
    else
        omega_raw = 0.0;
    end

    omega_z_cmd = min(max(omega_raw, -p.omega_z_max), p.omega_z_max);
    v_cmd_body  = [vx_cmd; 0.0; 0.0];

    if nargout >= 7
        nav_info.s_curr      = s_curr;
        nav_info.s_total     = s_total;
        nav_info.kappa_pp    = kappa_pp;
        nav_info.kappa_ref   = path_kap(idx_closest);
        nav_info.cross_track = hypot(p_xy(1) - p_ref_xy(1), p_xy(2) - p_ref_xy(2));
        nav_info.heading_err = wrap_to_pi_local(psi_ref - psi);
    end
end

%% ---------------- Local Course Builder & Math Helpers ----------------

function course = build_uturn_course(p)
%BUILD_UTURN_COURSE Construct outbound straight, R=0.70 m U-turn, and return straight.
    ds = 0.01;
    x_turn = p.x_uturn_center;
    R      = p.R_uturn;
    x_goal = p.x_return_goal;

    % Segment 1: Straight outbound from (0, 0) to (x_turn, 0) along +X
    x1   = 0.0 : ds : x_turn;
    y1   = zeros(size(x1));
    psi1 = zeros(size(x1));
    kap1 = zeros(size(x1));

    % Segment 2: 180-deg Semicircular U-Turn of radius R = 0.70 m
    arc_len = pi * R;
    n_arc   = max(round(arc_len / ds), 50);
    theta   = linspace(-pi/2, pi/2, n_arc);
    theta   = theta(2:end);
    x2      = x_turn + R * cos(theta);
    y2      = R + R * sin(theta);
    psi2    = theta + pi/2;
    kap2    = (1.0 / R) * ones(size(theta));

    % Segment 3: Straight return from (x_turn, 2*R) to (x_goal, 2*R) along -X
    x3   = (x_turn - ds) : -ds : x_goal;
    y3   = (2.0 * R) * ones(size(x3));
    psi3 = pi * ones(size(x3));
    kap3 = zeros(size(x3));

    xy    = [x1, x2, x3; y1, y2, y3];
    psi   = [psi1, psi2, psi3];
    kappa = [kap1, kap2, kap3];

    diffs = diff(xy, 1, 2);
    seg_d = sqrt(sum(diffs.^2, 1));
    s_vec = [0.0, cumsum(seg_d)];

    course.xy      = xy;
    course.psi     = psi;
    course.kappa   = kappa;
    course.s       = s_vec;
    course.s_total = s_vec(end);
    course.R_uturn = R;
    course.x_uturn = x_turn;
end

function ang = wrap_to_pi_local(ang_in)
    ang = mod(ang_in + pi, 2.0 * pi) - pi;
end
