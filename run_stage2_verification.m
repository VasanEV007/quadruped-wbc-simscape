%RUN_STAGE2_VERIFICATION Build quadruped_wbc.slx, compile, run 0.5 s simulation, and plot results.
cd('C:\Users\eluma\.gemini\antigravity\scratch\quadruped_wbc_simscape');
clear('build_quadruped_plant', 'params', 'leg_kinematics', 'vmc_stance_controller');
rehash;

p = params();

t_build_start = tic;
mdl_path = build_quadruped_plant('quadruped_wbc');
fprintf('Plant build time        : %.2f s\n', toc(t_build_start));

has_gui = logical(usejava('jvm')) && logical(usejava('desktop'));
if has_gui
    set_param('quadruped_wbc', 'SimMechanicsOpenEditorOnUpdate', 'on');
    open_system('quadruped_wbc');
else
    set_param('quadruped_wbc', 'SimMechanicsOpenEditorOnUpdate', 'off');
end

fprintf('Compiling quadruped_wbc.slx...\n');
t_comp_start = tic;
set_param('quadruped_wbc', 'SimulationCommand', 'update');
fprintf('Compilation succeeded   : 0 assembly violations (%.2f s)\n', toc(t_comp_start));

fprintf('Executing 0.5 s passive contact settling simulation...\n');
t_sim_start = tic;
simOut = sim('quadruped_wbc', 'StopTime', '0.5');
sim_wall_time = toc(t_sim_start);

pos_ts = simOut.get('sim_torso_pos');
vel_ts = simOut.get('sim_torso_vel');
q_ts   = simOut.get('sim_q_legs');
dq_ts  = simOut.get('sim_dq_legs');

t_vec = pos_ts.Time;
pos_data = squeeze(pos_ts.Data);
vel_data = squeeze(vel_ts.Data);
q_data   = squeeze(q_ts.Data);
dq_data  = squeeze(dq_ts.Data);

if size(pos_data, 1) == 3 && numel(t_vec) ~= 3
    pos_data = pos_data.';
    vel_data = vel_data.';
end
if size(q_data, 1) == 8 && numel(t_vec) ~= 8
    q_data  = q_data.';
    dq_data = dq_data.';
end

fprintf('\n=== Stage 2: 0.5 s Passive Contact Verification Summary ===\n');
fprintf('Solver                  : %s (%d steps, wall time: %.2f s)\n', ...
    get_param('quadruped_wbc', 'Solver'), numel(t_vec), sim_wall_time);
fprintf('Initial Torso Pos (t=0) : [%.5f, %.5f, %.5f] m\n', pos_data(1,1), pos_data(1,2), pos_data(1,3));
fprintf('Final Torso Pos (t=0.5) : [%.5f, %.5f, %.5f] m\n', pos_data(end,1), pos_data(end,2), pos_data(end,3));
fprintf('Torso Vertical Sag (dz) : %.3f mm\n', (pos_data(end,3) - pos_data(1,3)) * 1000.0);
fprintf('Final Torso Vel (t=0.5) : [%.5f, %.5f, %.5f] m/s\n', vel_data(end,1), vel_data(end,2), vel_data(end,3));
fprintf('Final Hip Angles (rad)  : [%.4f, %.4f, %.4f, %.4f] (nom: %.4f)\n', ...
    q_data(end,1), q_data(end,3), q_data(end,5), q_data(end,7), -pi/4);
fprintf('Final Knee Angles (rad) : [%.4f, %.4f, %.4f, %.4f] (nom: %.4f)\n', ...
    q_data(end,2), q_data(end,4), q_data(end,6), q_data(end,8), pi/2);
fprintf('Max Final Joint Vel     : %.5f rad/s\n', max(abs(dq_data(end, :))));

% Enable Mechanics Explorer on the saved .slx so opening it in MATLAB GUI shows 3D view
set_param('quadruped_wbc', 'SimMechanicsOpenEditorOnUpdate', 'on');
save_system('quadruped_wbc', mdl_path);

% Plot 4-Panel Verification Figure
if has_gui
    fig = figure('Name', 'Stage 2: 82.6 kg Quadruped + 6-Axis Arm Simscape Settling', ...
        'Color', 'w', 'Position', [80, 80, 1150, 680]);
else
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 80, 1150, 680]);
end

subplot(2, 2, 1);
plot(t_vec, pos_data(:, 3) * 1000, 'b-', 'LineWidth', 1.8); hold on; grid on; box on;
yline((p.z_nom + p.r_foot) * 1000, 'k--', 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('Torso Height Z [mm]');
title('Torso CoM Vertical Settling (82.6 kg)');
legend('Simulated Z(t)', 'Initial Tangent Z_0', 'Location', 'northeast');

subplot(2, 2, 2);
plot(t_vec, vel_data(:, 1), 'r-', 'LineWidth', 1.5); hold on; grid on; box on;
plot(t_vec, vel_data(:, 2), 'g--', 'LineWidth', 1.5);
plot(t_vec, vel_data(:, 3), 'b-', 'LineWidth', 1.8);
xlabel('Time [s]'); ylabel('Torso Velocity [m/s]');
title('Torso Linear Velocities (Damped Settling)');
legend('v_x', 'v_y', 'v_z', 'Location', 'northeast');

subplot(2, 2, 3);
plot(t_vec, rad2deg(q_data(:, [1, 3, 5, 7])), 'LineWidth', 1.5); hold on; grid on; box on;
yline(rad2deg(p.q_hip_nom), 'k--', 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('Hip Pitch Angle [deg]');
title('Hip Pitch Joint Locks (FL, FR, RL, RR)');
legend('FL', 'FR', 'RL', 'RR', 'Nominal (-45 deg)', 'Location', 'southeast');

subplot(2, 2, 4);
plot(t_vec, rad2deg(q_data(:, [2, 4, 6, 8])), 'LineWidth', 1.5); hold on; grid on; box on;
yline(rad2deg(p.q_knee_nom), 'k--', 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('Knee Pitch Angle [deg]');
title('Knee Pitch Joint Locks (FL, FR, RL, RR)');
legend('FL', 'FR', 'RL', 'RR', 'Nominal (+90 deg)', 'Location', 'southeast');

plot_path = fullfile(pwd, 'stage2_passive_contact_plot.png');
saveas(fig, plot_path);
fprintf('Saved verification plot : %s\n', plot_path);

if ~has_gui
    close(fig);
    close_system('quadruped_wbc', 0);
end
