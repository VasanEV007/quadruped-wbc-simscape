%RUN_STAGE4_VERIFICATION Execute Stage 4 Step Climbing, Edge Avoidance & 0.7 m U-Turn test.
cd('C:\Users\eluma\.gemini\antigravity\scratch\quadruped_wbc_simscape');
clear('params', 'leg_kinematics', 'vmc_stance_controller', ...
      'gait_swing_controller', 'raibert_foothold_planner', ...
      'pure_pursuit_navigator', 'wbc_step', 'build_quadruped_plant', 'wire_wbc_controller');
rehash;

p = params();

%% 1. Wire Stage 4 Controller (Gait Mode = 2: Full Course Navigation)
t_wire = tic;
mdl_path = wire_wbc_controller('quadruped_wbc', 2);
fprintf('Stage 4 model build & wiring time : %.2f s\n', toc(t_wire));

% Open model in Desktop GUI with Mechanics Explorer if desktop is active
has_gui = logical(usejava('jvm')) && logical(usejava('desktop'));
if has_gui
    open_system('quadruped_wbc');
    set_param('quadruped_wbc', 'SimMechanicsOpenEditorOnUpdate', 'on');
end

%% 2. Execute Full Locomotion Course Simulation (13.5 s)
sim_duration = '13.5';
set_param('quadruped_wbc', 'StopTime', sim_duration);
set_param('quadruped_wbc/Gait_Mode_Cmd', 'Value', '2');
clear('wbc_step', 'pure_pursuit_navigator');

fprintf('Running Stage 4 Full Course Simulation (%s s: Step Climb -> R=0.7m U-Turn -> Step Descent)...\n', sim_duration);
t_sim = tic;
simOut = sim('quadruped_wbc', 'StopTime', sim_duration);
wall_time = toc(t_sim);

%% 3. Extract & Analyze Telemetry
pos_ts  = simOut.get('sim_torso_pos');
vel_ts  = simOut.get('sim_torso_vel');
rpy_ts  = simOut.get('sim_torso_rpy');
q_ts    = simOut.get('sim_q_legs');
dq_ts   = simOut.get('sim_dq_legs');
grf_ts  = simOut.get('sim_F_grf');
mask_ts = simOut.get('sim_contact_mask');
nav_ts  = simOut.get('sim_nav_log');

t_vec = pos_ts.Time;
N_t   = numel(t_vec);
pos_d = fmt_mat(pos_ts.Data, 3, N_t);
vel_d = fmt_mat(vel_ts.Data, 3, N_t);
q_d   = fmt_mat(q_ts.Data,   8, N_t);
dq_d  = fmt_mat(dq_ts.Data,  8, N_t);

rpy_raw  = fmt_mat(rpy_ts.Data,  3,  numel(rpy_ts.Time));
grf_raw  = fmt_mat(grf_ts.Data,  12, numel(grf_ts.Time));
mask_raw = fmt_mat(mask_ts.Data, 4,  numel(mask_ts.Time));
nav_raw  = fmt_mat(nav_ts.Data,  12, numel(nav_ts.Time));

rpy_d  = interp1(rpy_ts.Time,  rpy_raw,  t_vec, 'previous', 'extrap');
grf_d  = interp1(grf_ts.Time,  grf_raw,  t_vec, 'previous', 'extrap');
mask_d = interp1(mask_ts.Time, mask_raw, t_vec, 'previous', 'extrap');
nav_d  = interp1(nav_ts.Time,  nav_raw,  t_vec, 'previous', 'extrap');

roll_deg  = rad2deg(rpy_d(:, 1));
pitch_deg = rad2deg(rpy_d(:, 2));
yaw_deg   = rad2deg(rpy_d(:, 3));

% Compute 3D World Foot Positions, Clearances & Lateral Slip During Stance
foot_pos_w = zeros(N_t, 3, 4);
foot_lat_slip_uturn = [];

for k = 1:N_t
    psi_k = rpy_d(k, 3);
    c_p = cos(psi_k); s_p = sin(psi_k);
    R_k = [c_p, -s_p, 0; s_p, c_p, 0; 0, 0, 1];
    omega_k = [0; 0; nav_d(k, 6)];
    v_com_b = R_k.' * vel_d(k, :).';

    for i = 1:4
        q_i  = [q_d(k, 2*i - 1);  q_d(k, 2*i)];
        dq_i = [dq_d(k, 2*i - 1); dq_d(k, 2*i)];
        [pf_b, J3_b] = leg_kinematics(q_i, i, p);
        pf_w = pos_d(k, :).' + R_k * pf_b;
        % Bottom of spherical foot in world frame
        pf_w(3) = pf_w(3) - p.r_foot;
        foot_pos_w(k, :, i) = pf_w.';

        % Lateral foot velocity in body frame during U-turn stance
        if mask_d(k, i) == 1 && abs(nav_d(k, 6)) > 0.15
            vf_b = v_com_b + cross(omega_k, pf_b) + J3_b * dq_i;
            % Lateral slip relative to pure rolling arc
            lat_slip_i = abs(vf_b(2) - omega_k(3) * pf_b(1));
            foot_lat_slip_uturn(end+1) = lat_slip_i; %#ok<SAGROW>
        end
    end
end

% Check step-edge clearance of all planned touchdowns x_td during active locomotion
idx_loco = t_vec >= p.t_stand_init;
td_x_all = nav_d(idx_loco, 7:10);
min_edge_dist_m = min(abs(td_x_all(:) - p.x_step_start));
max_foot_z_m    = squeeze(max(foot_pos_w(:, 3, :), [], 1));

% Cross-track error along course
cross_track_err = hypot(pos_d(:, 1) - nav_d(:, 1), pos_d(:, 2) - nav_d(:, 2));

peak_roll_deg  = max(abs(roll_deg));
peak_pitch_deg = max(abs(pitch_deg));
rms_roll_deg   = sqrt(mean(roll_deg.^2));
rms_pitch_deg  = sqrt(mean(pitch_deg.^2));
mean_lat_slip  = mean(foot_lat_slip_uturn);

fprintf('\n=== Stage 4: Step Climbing, Edge Avoidance & 0.7 m U-Turn Summary ===\n');
fprintf('Simulation Wall Time          : %.2f s (%d steps)\n', wall_time, N_t);
fprintf('Start CoM Pose [x, y, yaw]    : [%.3f m, %.3f m, %.1f deg]\n', ...
    pos_d(1,1), pos_d(1,2), yaw_deg(1));
fprintf('Peak X Reached (U-Turn Apex)  : %.3f m (target apex: %.3f m)\n', ...
    max(pos_d(:,1)), p.x_uturn_center + p.R_uturn);
fprintf('Final CoM Pose [x, y, yaw]    : [%.3f m, %.3f m, %.1f deg] (goal: [%.2f m, %.2f m, 180.0 deg])\n', ...
    pos_d(end,1), pos_d(end,2), yaw_deg(end), p.x_return_goal, 2*p.R_uturn);
fprintf('Peak Body |Roll|              : %.4f deg (RMS: %.4f deg, threshold: < 1.00 deg)\n', ...
    peak_roll_deg, rms_roll_deg);
fprintf('Peak Body |Pitch|             : %.4f deg (RMS: %.4f deg, threshold: < 1.00 deg)\n', ...
    peak_pitch_deg, rms_pitch_deg);
fprintf('Max CoM Height on 10cm Step   : %.4f m (ground nom: %.4f m, +10cm step target: %.4f m)\n', ...
    max(pos_d(:,3) - p.r_foot), p.z_nom, p.z_nom + p.h_step);
fprintf('Peak Foot World Z [FL..RR]    : [%.3f, %.3f, %.3f, %.3f] m (step=0.10m + swing=0.145m -> 0.245m)\n', ...
    max_foot_z_m(1), max_foot_z_m(2), max_foot_z_m(3), max_foot_z_m(4));
fprintf('Min Touchdown Dist to Edge    : %.4f m (safety margin requirement: >= 0.0500 m)\n', ...
    min_edge_dist_m);
fprintf('Cross-Track Path RMSE         : %.4f m (Max: %.4f m)\n', ...
    sqrt(mean(cross_track_err.^2)), max(cross_track_err));
fprintf('Mean Lateral Slip in U-Turn   : %.4f m/s\n', mean_lat_slip);

%% 4. Generate 6-Panel Stage 4 Verification Figure
[~, ~, ~, ~, ~, path_xy] = pure_pursuit_navigator(0, [0;0;p.z_nom], [0;0;0], p);

if has_gui
    fig = figure('Name', 'Stage 4: Step Climb, Foothold Edge Avoidance & 0.7m U-Turn', ...
        'Color', 'w', 'Position', [50, 50, 1340, 800]);
else
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50, 50, 1340, 800]);
end

% Panel 1: 2D Overhead Course Map & Step Platform
subplot(2, 3, 1);
hold on; grid on; box on;
patch([p.x_step_start, p.x_step_end, p.x_step_end, p.x_step_start], ...
      [p.y_step_min,   p.y_step_min, p.y_step_max, p.y_step_max], ...
      [0.93, 0.84, 0.68], 'FaceAlpha', 0.45, 'EdgeColor', [0.75, 0.50, 0.15], ...
      'LineWidth', 1.5, 'DisplayName', '10 cm Step Platform');
xline(p.x_step_start - p.d_edge_margin, 'r:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
xline(p.x_step_start + p.d_edge_margin, 'r:', 'LineWidth', 1.0, 'DisplayName', '\pm5.5 cm Edge Safety Zone');
plot(path_xy(1, :), path_xy(2, :), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference Path (R=0.7m)');
plot(pos_d(:, 1), pos_d(:, 2), 'b-', 'LineWidth', 2.2, 'DisplayName', 'Simulated CoM Trajectory');
plot(pos_d(1, 1), pos_d(1, 2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', 'Start (0, 0)');
plot(pos_d(end, 1), pos_d(end, 2), 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'End (Near Table)');
axis equal;
xlim([-0.3, 3.1]); ylim([-0.5, 1.9]);
xlabel('World X [m]'); ylabel('World Y [m]');
title('Overhead Course: 10 cm Step & R=0.70 m U-Turn');
legend('Location', 'northwest', 'FontSize', 7.5);

% Panel 2: Terrain-Aware CoM Height & Foot World Z Elevations
subplot(2, 3, 2);
hold on; grid on; box on;
plot(t_vec, (pos_d(:, 3) - p.r_foot) * 100, 'k-', 'LineWidth', 2.0, 'DisplayName', 'CoM Height z_{com}');
plot(t_vec, (nav_d(:, 4) - p.r_foot) * 100, 'r--', 'LineWidth', 1.4, 'DisplayName', 'Terrain-Aware Target z_d');
plot(t_vec, foot_pos_w(:, 3, 1) * 100, 'b-', 'LineWidth', 1.2, 'DisplayName', 'FL Foot Z');
plot(t_vec, foot_pos_w(:, 3, 2) * 100, 'Color', [0.1, 0.65, 0.25], 'LineWidth', 1.2, 'DisplayName', 'FR Foot Z');
yline(10.0, 'Color', [0.75, 0.50, 0.15], 'LineStyle', '--', 'LineWidth', 1.4, 'DisplayName', '10 cm Step Surface');
xlabel('Time [s]'); ylabel('Elevation [cm]');
title('Terrain-Aware CoM & Foot Elevations (h_{swing}=14.5 cm)');
legend('Location', 'east', 'FontSize', 7.5);

% Panel 3: Body Roll & Pitch (< 1.0 deg)
subplot(2, 3, 3);
plot(t_vec, roll_deg, 'b-', 'LineWidth', 1.7); hold on; grid on; box on;
plot(t_vec, pitch_deg, 'r-', 'LineWidth', 1.7);
yline(1.0, 'k--', 'LineWidth', 1.1); yline(-1.0, 'k--', 'LineWidth', 1.1);
ylim([-1.2, 1.2]);
xlabel('Time [s]'); ylabel('Angle [deg]');
title('Body Roll & Pitch Across Climb, U-Turn & Descent');
legend('Roll', 'Pitch', '\pm1.0^\circ Limit', 'Location', 'northeast');

% Panel 4: Yaw Heading Tracking (0 -> 180 deg U-Turn)
subplot(2, 3, 4);
plot(t_vec, rad2deg(unwrap(rpy_d(:, 3))), 'b-', 'LineWidth', 2.0); hold on; grid on; box on;
plot(t_vec, rad2deg(unwrap(nav_d(:, 3))), 'r--', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Yaw Heading \psi [deg]');
title('180^\circ Semicircular U-Turn Heading Tracking');
legend('Actual Yaw \psi', 'Pure-Pursuit \psi_{ref}', 'Location', 'southeast');

% Panel 5: Pure-Pursuit Velocity & Yaw Rate Commands
subplot(2, 3, 5);
yyaxis left;
plot(t_vec, nav_d(:, 5), 'b-', 'LineWidth', 1.8);
ylabel('Forward Speed v_{x,cmd} [m/s]');
ylim([0, 0.65]); grid on; box on;
yyaxis right;
plot(t_vec, nav_d(:, 6), 'r-', 'LineWidth', 1.6);
ylabel('Yaw Rate \omega_{z,cmd} [rad/s]');
ylim([-0.2, 0.95]);
xlabel('Time [s]');
title('Pure-Pursuit Commands (Step Transition Slowdown)');
legend('v_{x,cmd}', '\omega_{z,cmd}', 'Location', 'northeast');

% Panel 6: Foothold Distance from Step Edge (x = 1.20 m)
subplot(2, 3, 6);
hold on; grid on; box on;
for i = 1:4
    plot(t_vec, abs(nav_d(:, 6 + i) - p.x_step_start) * 100, 'LineWidth', 1.4, ...
        'DisplayName', sprintf('%s |x_{td} - x_{edge}|', p.leg_names{i}));
end
yline(5.0, 'r--', 'LineWidth', 1.8, 'DisplayName', '5.0 cm Min Safety Margin');
ylim([0, 40]);
xlabel('Time [s]'); ylabel('Distance from Step Edge [cm]');
title('Step-Edge Foothold Clearance (>= 5.0 cm)');
legend('Location', 'northeast', 'FontSize', 7.5);
drawnow;

plot_path = fullfile(pwd, 'stage4_step_uturn_navigation_plot.png');
saveas(fig, plot_path);
fprintf('Saved Stage 4 verification plot       : %s\n', plot_path);

if ~has_gui
    close(fig);
    close_system('quadruped_wbc', 0);
end

function M = fmt_mat(raw_data, n_cols, n_steps)
    M = squeeze(raw_data);
    if size(M, 1) == n_cols && size(M, 2) == n_steps && n_cols ~= n_steps
        M = M.';
    end
end
