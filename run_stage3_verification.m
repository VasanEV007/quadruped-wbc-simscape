%RUN_STAGE3_VERIFICATION Execute Stage 3 Closed-Loop Standing (1.0 s) & Diagonal Trot (2.0 s) tests.
cd('C:\Users\eluma\.gemini\antigravity\scratch\quadruped_wbc_simscape');
clear('params', 'leg_kinematics', 'vmc_stance_controller', ...
      'gait_swing_controller', 'wbc_step', 'build_quadruped_plant', 'wire_wbc_controller');
rehash;

p = params();
z_ref_world = p.z_nom + p.r_foot;

%% 1. Wire Closed-Loop WBC Controller into quadruped_wbc.slx
t_wire = tic;
mdl_path = wire_wbc_controller('quadruped_wbc', 0);
fprintf('Controller wiring time  : %.2f s\n', toc(t_wire));

%% 2. Test A: 4-Leg Standing Closed-Loop Balance (1.0 s)
fprintf('\n--- Running Test A: 4-Leg Standing Balance (1.0 s) ---\n');
set_param('quadruped_wbc/Gait_Mode_Cmd', 'Value', '0');
clear('wbc_step'); % Reset persistent cache
t_simA = tic;
simOutA = sim('quadruped_wbc', 'StopTime', '1.0');
wall_A = toc(t_simA);

[tA, posA, rpyA, qA, tauA, grfA, maskA] = extract_wbc_telemetry(simOutA);
rollA_deg  = rad2deg(rpyA(:, 1));
pitchA_deg = rad2deg(rpyA(:, 2));
yawA_deg   = rad2deg(rpyA(:, 3));
z_com_A    = posA(:, 3) - p.r_foot; % CoM height above foot contact point

peak_roll_A  = max(abs(rollA_deg));
peak_pitch_A = max(abs(pitchA_deg));
ss_roll_A    = abs(rollA_deg(end));
ss_pitch_A   = abs(pitchA_deg(end));
ss_z_err_A   = (z_com_A(end) - p.z_nom) * 1000.0; % [mm]

fprintf('Test A Wall Time        : %.2f s (%d steps)\n', wall_A, numel(tA));
fprintf('Test A Peak |Roll|      : %.4f deg\n', peak_roll_A);
fprintf('Test A Peak |Pitch|     : %.4f deg\n', peak_pitch_A);
fprintf('Test A Steady |Roll|    : %.4f deg (threshold: < 0.50 deg)\n', ss_roll_A);
fprintf('Test A Steady |Pitch|   : %.4f deg (threshold: < 0.50 deg)\n', ss_pitch_A);
fprintf('Test A Final CoM Height : %.5f m (z_nom = %.5f m, error = %.3f mm)\n', ...
    z_com_A(end), p.z_nom, ss_z_err_A);

%% 3. Test B: In-Place Diagonal Trot Closed-Loop Balance (2.0 s)
fprintf('\n--- Running Test B: In-Place Diagonal Trot (2.0 s) ---\n');
set_param('quadruped_wbc/Gait_Mode_Cmd', 'Value', '1');
clear('wbc_step'); % Reset persistent cache
t_simB = tic;
simOutB = sim('quadruped_wbc', 'StopTime', '2.0');
wall_B = toc(t_simB);

[tB, posB, rpyB, qB, tauB, grfB, maskB] = extract_wbc_telemetry(simOutB);
rollB_deg  = rad2deg(rpyB(:, 1));
pitchB_deg = rad2deg(rpyB(:, 2));
yawB_deg   = rad2deg(rpyB(:, 3));
z_com_B    = posB(:, 3) - p.r_foot;

peak_roll_B    = max(abs(rollB_deg));
peak_pitch_B   = max(abs(pitchB_deg));
max_z_err_B_mm = max(abs(z_com_B - p.z_nom)) * 1000.0;
rms_z_err_B_mm = sqrt(mean((z_com_B - p.z_nom).^2)) * 1000.0;

% Compute actual 4-foot vertical clearances above ground in Test B
foot_z_clearance_B = zeros(numel(tB), 4);
for k = 1:numel(tB)
    for i = 1:4
        q_ik = [qB(k, 2*i - 1); qB(k, 2*i)];
        [p_foot_b, ~] = leg_kinematics(q_ik, i, p);
        foot_z_clearance_B(k, i) = (posB(k, 3) - p.r_foot) + p_foot_b(3);
    end
end
max_foot_clear_mm = max(foot_z_clearance_B) * 1000.0;

fprintf('Test B Wall Time        : %.2f s (%d steps)\n', wall_B, numel(tB));
fprintf('Test B Peak |Roll|      : %.4f deg (threshold: < 1.00 deg)\n', peak_roll_B);
fprintf('Test B Peak |Pitch|     : %.4f deg (threshold: < 1.00 deg)\n', peak_pitch_B);
fprintf('Test B Max CoM Z Error  : %.3f mm (RMS: %.3f mm)\n', max_z_err_B_mm, rms_z_err_B_mm);
fprintf('Test B Peak Foot Lift   : FL=%.1f mm, FR=%.1f mm, RL=%.1f mm, RR=%.1f mm (target: 60.0 mm)\n', ...
    max_foot_clear_mm(1), max_foot_clear_mm(2), max_foot_clear_mm(3), max_foot_clear_mm(4));

% Save model with Mechanics Explorer enabled for interactive desktop viewing
set_param('quadruped_wbc', 'SimMechanicsOpenEditorOnUpdate', 'on');
save_system('quadruped_wbc', mdl_path);

%% 4. Plot 6-Panel Verification Figure
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [60, 60, 1280, 780]);

subplot(2, 3, 1);
plot(tA, rollA_deg, 'b-', 'LineWidth', 1.8); hold on; grid on; box on;
plot(tA, pitchA_deg, 'r-', 'LineWidth', 1.8);
yline(0.5, 'k--', 'LineWidth', 1.0); yline(-0.5, 'k--', 'LineWidth', 1.0);
xlabel('Time [s]'); ylabel('Angle [deg]');
title('Test A (4-Leg Stand): Roll & Pitch (< 0.5 deg)');
legend('Roll', 'Pitch', '\pm0.5^\circ Limit', 'Location', 'northeast');

subplot(2, 3, 2);
plot(tA, z_com_A * 1000, 'b-', 'LineWidth', 1.8); hold on; grid on; box on;
yline(p.z_nom * 1000, 'r--', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('CoM Height [mm]');
title('Test A (4-Leg Stand): CoM Height z_{nom}');
legend('Actual z_{com}', 'Target z_{nom} (480.8 mm)', 'Location', 'southeast');

subplot(2, 3, 3);
plot(tA, grfA(:, [3, 6, 9, 12]), 'LineWidth', 1.6); grid on; box on;
xlabel('Time [s]'); ylabel('Vertical GRF F_z [N]');
title('Test A: 4-Leg WLS Vertical Forces');
legend('FL', 'FR', 'RL', 'RR', 'Location', 'southeast');

subplot(2, 3, 4);
plot(tB, rollB_deg, 'b-', 'LineWidth', 1.8); hold on; grid on; box on;
plot(tB, pitchB_deg, 'r-', 'LineWidth', 1.8);
yline(1.0, 'k--', 'LineWidth', 1.0); yline(-1.0, 'k--', 'LineWidth', 1.0);
xlabel('Time [s]'); ylabel('Angle [deg]');
title('Test B (Diagonal Trot): Roll & Pitch (< 1.0 deg)');
legend('Roll', 'Pitch', '\pm1.0^\circ Limit', 'Location', 'northeast');

subplot(2, 3, 5);
plot(tB, foot_z_clearance_B(:, 1) * 1000, 'b-', 'LineWidth', 1.7); hold on; grid on; box on;
plot(tB, foot_z_clearance_B(:, 2) * 1000, 'r-', 'LineWidth', 1.7);
plot(tB, foot_z_clearance_B(:, 3) * 1000, 'r--', 'LineWidth', 1.4);
plot(tB, foot_z_clearance_B(:, 4) * 1000, 'b--', 'LineWidth', 1.4);
yline(p.h_swing * 1000, 'k:', 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('Foot Clearance [mm]');
title('Test B: Diagonal Trot Swing Foot Clearance');
legend('FL (Pair A)', 'FR (Pair B)', 'RL (Pair B)', 'RR (Pair A)', '60 mm Apex', 'Location', 'northeast');

subplot(2, 3, 6);
plot(tB, grfB(:, 3), 'b-', 'LineWidth', 1.5); hold on; grid on; box on;
plot(tB, grfB(:, 6), 'r-', 'LineWidth', 1.5);
plot(tB, grfB(:, 9), 'r--', 'LineWidth', 1.3);
plot(tB, grfB(:, 12), 'b--', 'LineWidth', 1.3);
xlabel('Time [s]'); ylabel('Vertical GRF F_z [N]');
title('Test B: Alternating Diagonal Stance GRFs');
legend('FL', 'FR', 'RL', 'RR', 'Location', 'northeast');

plot_path = fullfile(pwd, 'stage3_closed_loop_trot_plot.png');
saveas(fig, plot_path);
close(fig);
close_system('quadruped_wbc', 0);
fprintf('Saved Stage 3 verification plot: %s\n', plot_path);

%% ---------------- Local Helper Function ----------------
function [t_vec, pos_d, rpy_d, q_d, tau_d, grf_d, mask_d] = extract_wbc_telemetry(simOut)
    pos_ts  = simOut.get('sim_torso_pos');
    rpy_ts  = simOut.get('sim_torso_rpy');
    q_ts    = simOut.get('sim_q_legs');
    tau_ts  = simOut.get('sim_tau_legs');
    grf_ts  = simOut.get('sim_F_grf');
    mask_ts = simOut.get('sim_contact_mask');

    t_vec  = pos_ts.Time;
    pos_d  = format_ts_matrix(pos_ts.Data, 3, numel(t_vec));
    q_d    = format_ts_matrix(q_ts.Data, 8, numel(t_vec));

    % Interpolate discrete 1 kHz controller logs onto t_vec if lengths differ
    rpy_raw  = format_ts_matrix(rpy_ts.Data, 3, numel(rpy_ts.Time));
    tau_raw  = format_ts_matrix(tau_ts.Data, 8, numel(tau_ts.Time));
    grf_raw  = format_ts_matrix(grf_ts.Data, 12, numel(grf_ts.Time));
    mask_raw = format_ts_matrix(mask_ts.Data, 4, numel(mask_ts.Time));

    if numel(rpy_ts.Time) ~= numel(t_vec)
        rpy_d  = interp1(rpy_ts.Time, rpy_raw, t_vec, 'previous', 'extrap');
        tau_d  = interp1(tau_ts.Time, tau_raw, t_vec, 'previous', 'extrap');
        grf_d  = interp1(grf_ts.Time, grf_raw, t_vec, 'previous', 'extrap');
        mask_d = interp1(mask_ts.Time, mask_raw, t_vec, 'previous', 'extrap');
    else
        rpy_d  = rpy_raw;
        tau_d  = tau_raw;
        grf_d  = grf_raw;
        mask_d = mask_raw;
    end
end

function M = format_ts_matrix(raw_data, n_cols, n_steps)
    M = squeeze(raw_data);
    if size(M, 1) == n_cols && size(M, 2) == n_steps && n_cols ~= n_steps
        M = M.';
    end
end
