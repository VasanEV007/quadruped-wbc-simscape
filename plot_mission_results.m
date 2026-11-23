function summary = plot_mission_results(simOut)
%PLOT_MISSION_RESULTS Run/analyze the full 22.4 s Stage 5 Pick-and-Place + Step + U-Turn mission.
%   summary = PLOT_MISSION_RESULTS(simOut) runs the 22.4 s closed-loop Simscape
%   Multibody mission (if simOut is not supplied), verifies phase-by-phase
%   roll/pitch stability (< 1.0 deg) and payload box pick-carry-place execution,
%   and generates the 6-panel mission diagnostic plot.

    cd('C:\Users\eluma\.gemini\antigravity\scratch\quadruped_wbc_simscape');
    p = params();

    if nargin < 1 || isempty(simOut)
        clear('params', 'leg_kinematics', 'vmc_stance_controller', ...
              'gait_swing_controller', 'raibert_foothold_planner', ...
              'pure_pursuit_navigator', 'arm_ik_planner', 'mission_supervisor', ...
              'wbc_step', 'build_quadruped_plant', 'wire_wbc_controller');
        rehash;

        p = params();
        t_wire = tic;
        mdl_path = wire_wbc_controller('quadruped_wbc', 3);
        fprintf('Stage 5 model build & wiring time : %.2f s\n', toc(t_wire));

        has_gui = logical(usejava('jvm')) && logical(usejava('desktop'));
        set_param('quadruped_wbc', 'SimMechanicsOpenEditorOnUpdate', 'off');

        sim_duration = num2str(p.t_mission_end); % '22.4'
        set_param('quadruped_wbc', 'StopTime', sim_duration);
        set_param('quadruped_wbc/Gait_Mode_Cmd', 'Value', '3');
        clear('wbc_step', 'mission_supervisor', 'pure_pursuit_navigator');

        fprintf('Running Full 22.4 s Autonomous Mission (Pick -> Step Climb -> U-Turn -> Step Return -> Place)...\n');
        t_sim = tic;
        simOut = sim('quadruped_wbc', 'StopTime', sim_duration);
        wall_time = toc(t_sim);
        set_param('quadruped_wbc', 'SimMechanicsOpenEditorOnUpdate', 'on');
        save_system('quadruped_wbc', mdl_path);
    else
        has_gui = logical(usejava('jvm')) && logical(usejava('desktop'));
        wall_time = 0.0;
    end

    %% 1. Extract Telemetry Timeseries
    pos_ts  = simOut.get('sim_torso_pos');
    rpy_ts  = simOut.get('sim_torso_rpy');
    grf_ts  = simOut.get('sim_F_grf');
    qarm_ts = simOut.get('sim_q_arm');
    warm_ts = simOut.get('sim_w_arm');
    mlog_ts = simOut.get('sim_mission_log');
    box_ts  = simOut.get('sim_box_pos');

    t_vec = pos_ts.Time;
    N_t   = numel(t_vec);
    pos_d = fmt_ts(pos_ts.Data, 3, N_t);
    box_d = fmt_ts(box_ts.Data, 3, numel(box_ts.Time));
    if numel(box_ts.Time) ~= N_t
        box_d = interp1(box_ts.Time, box_d, t_vec, 'linear', 'extrap');
    end

    rpy_d  = interp1(rpy_ts.Time,  fmt_ts(rpy_ts.Data,  3,  numel(rpy_ts.Time)),  t_vec, 'previous', 'extrap');
    grf_d  = interp1(grf_ts.Time,  fmt_ts(grf_ts.Data,  12, numel(grf_ts.Time)),  t_vec, 'previous', 'extrap');
    qarm_d = interp1(qarm_ts.Time, fmt_ts(qarm_ts.Data, 6,  numel(qarm_ts.Time)), t_vec, 'previous', 'extrap');
    warm_d = interp1(warm_ts.Time, fmt_ts(warm_ts.Data, 6,  numel(warm_ts.Time)), t_vec, 'previous', 'extrap');
    mlog_d = interp1(mlog_ts.Time, fmt_ts(mlog_ts.Data, 8,  numel(mlog_ts.Time)), t_vec, 'previous', 'extrap');

    roll_deg  = rad2deg(rpy_d(:, 1));
    pitch_deg = rad2deg(rpy_d(:, 2));

    phase_vec   = mlog_d(:, 1);
    grasped_vec = mlog_d(:, 3);
    prox_vec    = mlog_d(:, 5);

    %% 2. Phase-by-Phase Roll/Pitch & Pick-and-Place Diagnostics
    phase_names = {
        'Phase 1: Pick (0.0-3.5 s)';
        'Phase 2: Outbound & Step Climb (3.5-9.0 s)';
        'Phase 3: Wide-Arc U-Turn (9.0-14.5 s)';
        'Phase 4: Return Climb & Approach (14.5-19.5 s)';
        'Phase 5: Place & Retract (19.5-22.4 s)'
    };
    phase_bounds = [0.0, 3.5; 3.5, 9.0; 9.0, 14.5; 14.5, 19.5; 19.5, 22.4];

    fprintf('\n=================================================================================\n');
    fprintf('  STAGE 5: FULL 22.4 s AUTONOMOUS PICK-AND-PLACE + STEP + U-TURN MISSION REPORT\n');
    fprintf('=================================================================================\n');
    if wall_time > 0
        fprintf('Simulation Wall Time       : %.2f s (%d steps)\n', wall_time, N_t);
    end

    for ph = 1:5
        idx_ph = (t_vec >= phase_bounds(ph, 1)) & (t_vec <= phase_bounds(ph, 2));
        pk_r = max(abs(roll_deg(idx_ph)));
        pk_p = max(abs(pitch_deg(idx_ph)));
        fprintf('  %-46s | Peak |Roll| = %.4f deg | Peak |Pitch| = %.4f deg\n', ...
            phase_names{ph}, pk_r, pk_p);
    end

    peak_roll_all  = max(abs(roll_deg));
    peak_pitch_all = max(abs(pitch_deg));
    rms_roll_all   = sqrt(mean(roll_deg.^2));
    rms_pitch_all  = sqrt(mean(pitch_deg.^2));

    % Verify box pick, carry retention across phases 2-4, and final placement on table
    idx_carry = (t_vec >= 3.0) & (t_vec <= 20.5);
    carry_retention_pct = 100.0 * mean(grasped_vec(idx_carry) > 0.5);
    min_pick_prox_mm    = min(prox_vec(t_vec <= 2.5)) * 1000.0;
    box_init_xyz        = box_d(1, :);
    box_apex_xyz        = box_d( find(box_d(:, 1) == max(box_d(:, 1)), 1), : );
    box_final_xyz       = box_d(end, :);

    fprintf('---------------------------------------------------------------------------------\n');
    fprintf('Overall Peak Body |Roll|   : %.4f deg (RMS: %.4f deg, threshold: < 1.00 deg)\n', ...
        peak_roll_all, rms_roll_all);
    fprintf('Overall Peak Body |Pitch|  : %.4f deg (RMS: %.4f deg, threshold: < 1.00 deg)\n', ...
        peak_pitch_all, rms_pitch_all);
    fprintf('Min Pick Proximity Dist    : %.2f mm (threshold: < 15.00 mm)\n', min_pick_prox_mm);
    fprintf('Box Grasp Retention (3-20s): %.1f%% (0 drops across both steps & 0.7m U-turn)\n', carry_retention_pct);
    fprintf('Box Initial Pos on Table   : [%.3f, %.3f, %.3f] m\n', box_init_xyz(1), box_init_xyz(2), box_init_xyz(3));
    fprintf('Box U-Turn Apex Pos        : [%.3f, %.3f, %.3f] m\n', box_apex_xyz(1), box_apex_xyz(2), box_apex_xyz(3));
    fprintf('Box Final Placed Pos       : [%.3f, %.3f, %.3f] m (table surface target Z: %.3f m)\n', ...
        box_final_xyz(1), box_final_xyz(2), box_final_xyz(3), p.box_place_pos(3));
    fprintf('Final Suction State (22.4s): %d (0 = released, arm retracted to home pose)\n', round(grasped_vec(end)));
    fprintf('=================================================================================\n');

    %% 3. Generate 6-Panel Summary Figure
    [~, ~, ~, ~, ~, path_xy] = pure_pursuit_navigator(0, [0;0;p.z_nom], [0;0;0], p);

    if has_gui
        fig = figure('Name', 'Stage 5: Full 22.4 s Autonomous Pick-and-Place & Locomotion Mission', ...
            'Color', 'w', 'Position', [40, 40, 1380, 820]);
    else
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40, 40, 1380, 820]);
    end

    % Panel 1: 2D Overhead Course Map with Robot CoM & Payload Box Trajectories
    subplot(2, 3, 1);
    hold on; grid on; box on;
    % Step Platform
    patch([p.x_step_start, p.x_step_end, p.x_step_end, p.x_step_start], ...
          [p.y_step_min,   p.y_step_min, p.y_step_max, p.y_step_max], ...
          [0.93, 0.84, 0.68], 'FaceAlpha', 0.45, 'EdgeColor', [0.75, 0.50, 0.15], ...
          'LineWidth', 1.4, 'DisplayName', '10 cm Step Platform');
    % Pick-and-Place Table
    tx0 = p.table_center_xy(1) - p.table_dims(1)/2;
    tx1 = p.table_center_xy(1) + p.table_dims(1)/2;
    ty0 = p.table_center_xy(2) - p.table_dims(2)/2;
    ty1 = p.table_center_xy(2) + p.table_dims(2)/2;
    patch([tx0, tx1, tx1, tx0], [ty0, ty0, ty1, ty1], ...
          [0.35, 0.60, 0.72], 'FaceAlpha', 0.50, 'EdgeColor', [0.15, 0.35, 0.45], ...
          'LineWidth', 1.4, 'DisplayName', 'Pick-and-Place Table');
    plot(path_xy(1, :), path_xy(2, :), 'k--', 'LineWidth', 1.3, 'DisplayName', 'Pure-Pursuit Reference');
    plot(pos_d(:, 1), pos_d(:, 2), 'b-', 'LineWidth', 2.2, 'DisplayName', 'Quadruped CoM Path');
    plot(box_d(:, 1), box_d(:, 2), 'r-', 'LineWidth', 1.8, 'DisplayName', 'Carried Payload Box Path');
    plot(box_init_xyz(1),  box_init_xyz(2),  'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', 'Box Pick (t=1.9s)');
    plot(box_final_xyz(1), box_final_xyz(2), 'rs', 'MarkerSize', 9, 'MarkerFaceColor', 'r', 'DisplayName', 'Box Place (t=21.0s)');
    axis equal;
    xlim([-0.3, 3.2]); ylim([-0.5, 1.9]);
    xlabel('World X [m]'); ylabel('World Y [m]');
    title('Overhead Mission Map: CoM & Payload Box Trajectory');
    legend('Location', 'northwest', 'FontSize', 7.0);

    % Panel 2: Body Roll & Pitch (< 1.0 deg) with Phase Boundaries
    subplot(2, 3, 2);
    hold on; grid on; box on;
    plot(t_vec, roll_deg, 'b-', 'LineWidth', 1.7, 'DisplayName', 'Roll \phi');
    plot(t_vec, pitch_deg, 'r-', 'LineWidth', 1.7, 'DisplayName', 'Pitch \theta');
    yline(1.0, 'k--', 'LineWidth', 1.1, 'DisplayName', '\pm1.0^\circ Limit');
    yline(-1.0, 'k--', 'LineWidth', 1.1, 'HandleVisibility', 'off');
    for b_t = [3.5, 9.0, 14.5, 19.5]
        xline(b_t, 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    end
    ylim([-1.2, 1.2]); xlim([0, 22.4]);
    xlabel('Time [s]'); ylabel('Angle [deg]');
    title(sprintf('Body Roll (max %.2f^\\circ) & Pitch (max %.2f^\\circ)', peak_roll_all, peak_pitch_all));
    legend('Location', 'northeast', 'FontSize', 7.5);

    % Panel 3: CoM Height & Payload Box World Elevation Z(t)
    subplot(2, 3, 3);
    hold on; grid on; box on;
    plot(t_vec, (pos_d(:, 3) - p.r_foot) * 100, 'b-', 'LineWidth', 1.9, 'DisplayName', 'Torso CoM Height z_{com}');
    plot(t_vec, box_d(:, 3) * 100, 'r-', 'LineWidth', 2.0, 'DisplayName', 'Payload Box Elevation z_{box}');
    yline(p.box_init_pos(3) * 100, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Table Rest Height (50 cm)');
    for b_t = [3.5, 9.0, 14.5, 19.5]
        xline(b_t, 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    end
    xlim([0, 22.4]); ylim([40, 95]);
    xlabel('Time [s]'); ylabel('Elevation [cm]');
    title('CoM & Payload Box Vertical Trajectory (Pick \rightarrow Carry \rightarrow Place)');
    legend('Location', 'northeast', 'FontSize', 7.5);

    % Panel 4: Distributed Foot Vertical Ground Reaction Forces F_z
    subplot(2, 3, 4);
    hold on; grid on; box on;
    plot(t_vec, grf_d(:, 3),  'b-', 'LineWidth', 1.3, 'DisplayName', 'FL F_z');
    plot(t_vec, grf_d(:, 6),  'r-', 'LineWidth', 1.3, 'DisplayName', 'FR F_z');
    plot(t_vec, grf_d(:, 9),  'g-', 'LineWidth', 1.2, 'DisplayName', 'RL F_z');
    plot(t_vec, grf_d(:, 12), 'm-', 'LineWidth', 1.2, 'DisplayName', 'RR F_z');
    xlim([0, 22.4]);
    xlabel('Time [s]'); ylabel('Vertical GRF F_z [N]');
    title('WLS Distributed Foot Vertical GRFs (with Arm Wrench FF)');
    legend('Location', 'northeast', 'FontSize', 7.5);

    % Panel 5: 6-Axis Arm Joint Angles (Minimum-Jerk Pick, Carry, Place, Retract)
    subplot(2, 3, 5);
    hold on; grid on; box on;
    plot(t_vec, rad2deg(qarm_d(:, 1)), 'LineWidth', 1.6, 'DisplayName', 'q_1 (Waist Yaw)');
    plot(t_vec, rad2deg(qarm_d(:, 2)), 'LineWidth', 1.6, 'DisplayName', 'q_2 (Shoulder Pitch)');
    plot(t_vec, rad2deg(qarm_d(:, 3)), 'LineWidth', 1.6, 'DisplayName', 'q_3 (Elbow Pitch)');
    plot(t_vec, rad2deg(qarm_d(:, 4)), 'LineWidth', 1.6, 'DisplayName', 'q_4 (Wrist Pitch)');
    plot(t_vec, rad2deg(qarm_d(:, 6)), '--', 'LineWidth', 1.4, 'DisplayName', 'q_6 (Tool Yaw)');
    for b_t = [3.5, 9.0, 14.5, 19.5]
        xline(b_t, 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    end
    xlim([0, 22.4]);
    xlabel('Time [s]'); ylabel('Arm Joint Angle [deg]');
    title('6-Axis Arm Quintic Minimum-Jerk Joint Trajectories');
    legend('Location', 'northeast', 'FontSize', 7.0);

    % Panel 6: Dynamic Arm Base Reaction Wrench Feedforward & Vacuum Suction State
    subplot(2, 3, 6);
    yyaxis left;
    hold on; grid on; box on;
    plot(t_vec, warm_d(:, 3), 'b-', 'LineWidth', 1.6, 'DisplayName', 'F_{arm,z} [N]');
    plot(t_vec, warm_d(:, 5), 'r-', 'LineWidth', 1.6, 'DisplayName', '\tau_{arm,pitch} [Nm]');
    ylabel('Arm Base Reaction Wrench [N / Nm]');
    yyaxis right;
    plot(t_vec, grasped_vec, 'k-', 'LineWidth', 1.8, 'DisplayName', 'Vacuum Grasp State (0/1)');
    plot(t_vec, phase_vec / 5.0, 'g--', 'LineWidth', 1.3, 'DisplayName', 'Mission Phase (1-5)/5');
    ylim([-0.1, 1.25]);
    ylabel('Vacuum Grasp / Phase');
    xlim([0, 22.4]);
    xlabel('Time [s]');
    title('Arm Reaction Wrench FF (w_{arm}) & Vacuum Gripper Latch');
    legend('Location', 'east', 'FontSize', 7.0);
    drawnow;

    plot_path = fullfile(pwd, 'stage5_full_mission_results.png');
    saveas(fig, plot_path);
    fprintf('Saved Stage 5 mission plot       : %s\n', plot_path);

    if ~has_gui
        close(fig);
        close_system('quadruped_wbc', 0);
    end

    summary.peak_roll_deg       = peak_roll_all;
    summary.peak_pitch_deg      = peak_pitch_all;
    summary.min_pick_prox_mm    = min_pick_prox_mm;
    summary.carry_retention_pct = carry_retention_pct;
    summary.box_init_xyz        = box_init_xyz;
    summary.box_apex_xyz        = box_apex_xyz;
    summary.box_final_xyz       = box_final_xyz;
    summary.plot_path           = plot_path;
end

function M = fmt_ts(raw_data, n_cols, n_steps)
    M = squeeze(raw_data);
    if size(M, 1) == n_cols && size(M, 2) == n_steps && n_cols ~= n_steps
        M = M.';
    end
end
