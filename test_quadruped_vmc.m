function tests = test_quadruped_vmc()
%TEST_QUADRUPED_VMC Unit test suite for 82.6 kg quadruped kinematics & VMC.
%   Run via: results = runtests('test_quadruped_vmc.m');
    tests = functiontests(localfunctions);
end

function test_params_mass_and_structure(testCase)
    p = params();

    % Verify total mass is exactly 82.6 kg
    verifyEqual(testCase, p.m_total, 82.6, 'AbsTol', 1.0e-12, ...
        'Total system mass (quadruped + 6-axis arm + vacuum gripper) must be 82.6 kg.');

    % Verify geometric & contact parameters
    verifyGreaterThan(testCase, p.L_thigh, 0.0);
    verifyGreaterThan(testCase, p.L_shank, 0.0);
    verifyEqual(testCase, size(p.hip_offsets), [3, 4]);
    verifyGreaterThan(testCase, p.k_contact, 0.0);
    verifyGreaterThan(testCase, p.b_contact, 0.0);
    verifyGreaterThan(testCase, p.mu_vmc, 0.0);

    % Verify VMC gain matrices are symmetric positive-definite
    verifyEqual(testCase, size(p.Kp_com), [3, 3]);
    verifyEqual(testCase, size(p.Kd_com), [3, 3]);
    verifyEqual(testCase, size(p.Kp_rot), [3, 3]);
    verifyEqual(testCase, size(p.Kd_rot), [3, 3]);
    verifyEqual(testCase, size(p.S_wrench), [6, 6]);
    verifyEqual(testCase, size(p.W_force), [12, 12]);
    verifyGreaterThan(testCase, min(eig(p.S_wrench)), 0.0);
    verifyGreaterThan(testCase, min(eig(p.W_force)), 0.0);
end

function test_leg_kinematics_nominal_and_hip_offsets(testCase)
    p = params();
    q_nom = p.q_leg_nom; % [-pi/4; pi/2]

    for leg_id = 1:4
        [p_foot_3d, J_3d, p_foot_2d, J_2d, kin_info] = leg_kinematics(q_nom, leg_id, p);

        % At q1 = -pi/4, q2 = +pi/2, foot X is directly below hip (x_rel = 0)
        verifyEqual(testCase, p_foot_2d(1), 0.0, 'AbsTol', 1.0e-12);
        verifyEqual(testCase, p_foot_2d(2), -p.z_nom, 'AbsTol', 1.0e-12);

        % 3D foot position must equal hip_offset + [x_rel; 0; z_rel]
        expected_3d = p.hip_offsets(:, leg_id) + [p_foot_2d(1); 0.0; p_foot_2d(2)];
        verifyEqual(testCase, p_foot_3d, expected_3d, 'AbsTol', 1.0e-12);

        % Lateral row of 3D Jacobian must be zero (no hip abduction)
        verifyEqual(testCase, J_3d(2, :), [0.0, 0.0], 'AbsTol', 1.0e-15);

        % Analytical determinant check
        verifyEqual(testCase, kin_info.det_J2d, det(J_2d), 'AbsTol', 1.0e-12);
    end
end

function test_leg_jacobians_against_finite_difference(testCase)
    p = params();
    rng(42); % Deterministic seed
    eps_fd = 1.0e-7;

    q_rand = [ ...
        p.q_hip_lim(1) + (p.q_hip_lim(2) - p.q_hip_lim(1)) * rand(1, 20); ...
        p.q_knee_lim(1) + (p.q_knee_lim(2) - p.q_knee_lim(1)) * rand(1, 20) ...
    ];
    q_samples = [ ...
        p.q_leg_nom, ...
        [-0.6; 1.1], ...
        [-0.3; 0.7], ...
        [-0.9; 1.8], ...
        [ 0.2; 0.9], ...
        q_rand ...
    ];

    max_err_2d = 0.0;
    max_err_3d = 0.0;

    for k = 1:size(q_samples, 2)
        q = q_samples(:, k);
        for leg_id = 1:4
            [~, J_3d, ~, J_2d] = leg_kinematics(q, leg_id, p);

            J_3d_fd = zeros(3, 2);
            J_2d_fd = zeros(2, 2);

            for j = 1:2
                dq = zeros(2, 1);
                dq(j) = eps_fd;
                [p3_plus,  ~, p2_plus,  ~] = leg_kinematics(q + dq, leg_id, p);
                [p3_minus, ~, p2_minus, ~] = leg_kinematics(q - dq, leg_id, p);

                J_3d_fd(:, j) = (p3_plus - p3_minus) / (2.0 * eps_fd);
                J_2d_fd(:, j) = (p2_plus - p2_minus) / (2.0 * eps_fd);
            end

            err_2d = norm(J_2d - J_2d_fd, 'fro');
            err_3d = norm(J_3d - J_3d_fd, 'fro');
            max_err_2d = max(max_err_2d, err_2d);
            max_err_3d = max(max_err_3d, err_3d);

            verifyLessThan(testCase, err_2d, 1.0e-8, ...
                sprintf('2D Jacobian FD mismatch at sample %d, leg %d', k, leg_id));
            verifyLessThan(testCase, err_3d, 1.0e-8, ...
                sprintf('3D Jacobian FD mismatch at sample %d, leg %d', k, leg_id));
        end
    end

    fprintf('\n  [Jacobian FD Verification] Max 2D Error: %.3e, Max 3D Error: %.3e\n', ...
        max_err_2d, max_err_3d);
end

function test_vmc_static_standing_gravity_compensation(testCase)
    p = params();
    q_legs = repmat(p.q_leg_nom, 1, 4);
    contact_mask = [1; 1; 1; 1];

    state.p_com = [0.0; 0.0; p.z_nom];
    state.v_com = [0.0; 0.0; 0.0];
    state.rpy   = [0.0; 0.0; 0.0];
    state.omega = [0.0; 0.0; 0.0];
    des_state = state;

    [tau_legs, F_grf, wrench_cmd, info] = vmc_stance_controller( ...
        state, des_state, contact_mask, q_legs, p);

    % 1. Desired vertical force must equal M_total * g = 82.6 * 9.81 = 810.306 N
    expected_Fz_total = p.m_total * p.g;
    verifyEqual(testCase, wrench_cmd, [0; 0; expected_Fz_total; 0; 0; 0], 'AbsTol', 1.0e-10);

    % 2. Symmetric 4-leg stance must distribute Fz equally across all 4 feet
    Fz_per_foot = F_grf(3, :);
    verifyEqual(testCase, Fz_per_foot, repmat(Fz_per_foot(1), 1, 4), 'AbsTol', 1.0e-7);
    verifyEqual(testCase, sum(Fz_per_foot), expected_Fz_total, 'RelTol', 1.0e-4);

    % 3. Horizontal forces Fx, Fy must be zero in symmetric static standing
    verifyEqual(testCase, F_grf(1:2, :), zeros(2, 4), 'AbsTol', 1.0e-7);

    % 4. Unconstrained WLS solution must be directly feasible in static stance
    verifyTrue(testCase, info.is_unconstrained_feasible);

    % 5. Joint torques must match -J_3d^T * F_grf and be equal across all 4 legs
    [~, J_3d] = leg_kinematics(p.q_leg_nom, 1, p);
    expected_tau_1 = -J_3d.' * F_grf(:, 1);
    for i = 1:4
        verifyEqual(testCase, tau_legs(:, i), expected_tau_1, 'AbsTol', 1.0e-7);
    end

    fprintf('  [Static Stance VMC] Total Fz = %.2f N (%.2f N/foot), Tau = [%.2f, %.2f] Nm\n', ...
        sum(Fz_per_foot), Fz_per_foot(1), tau_legs(1, 1), tau_legs(2, 1));
end

function test_vmc_friction_pyramid_and_diagonal_trot(testCase)
    p = params();
    q_legs = repmat(p.q_leg_nom, 1, 4);

    % Diagonal trot pair: FL (1) and RR (4) in stance, FR (2) and RL (3) in swing
    contact_mask = [1; 0; 0; 1];

    state.p_com = [-0.03; 0.01; p.z_nom - 0.02];
    state.v_com = [-0.20; 0.05; -0.10];
    state.rpy   = [0.04; -0.06; 0.02];
    state.omega = [0.10; -0.15; 0.05];

    des_state.p_com = [0.0; 0.0; p.z_nom];
    des_state.v_com = [0.5; 0.0; 0.0];
    des_state.rpy   = [0.0; 0.0; 0.0];
    des_state.omega = [0.0; 0.0; 0.0];
    % High forward acceleration demand to trigger friction pyramid projection
    des_state.a_com_ff = [8.5; 3.0; 0.0];

    [tau_legs, F_grf, ~, info] = vmc_stance_controller( ...
        state, des_state, contact_mask, q_legs, p);

    % 1. Swing legs (2 and 3) must have strictly zero GRF and zero stance torque
    verifyEqual(testCase, F_grf(:, [2, 3]), zeros(3, 2), 'AbsTol', 1.0e-15);
    verifyEqual(testCase, tau_legs(:, [2, 3]), zeros(2, 2), 'AbsTol', 1.0e-15);

    % 2. Stance legs (1 and 4) must strictly satisfy normal bounds and friction pyramid
    for i = [1, 4]
        fx = F_grf(1, i);
        fy = F_grf(2, i);
        fz = F_grf(3, i);

        verifyGreaterThanOrEqual(testCase, fz, p.Fz_min - 1.0e-9);
        verifyLessThanOrEqual(testCase, fz, p.Fz_max + 1.0e-9);
        verifyLessThanOrEqual(testCase, abs(fx), p.mu_vmc * fz + 1.0e-7);
        verifyLessThanOrEqual(testCase, abs(fy), p.mu_vmc * fz + 1.0e-7);
        verifyLessThanOrEqual(testCase, info.friction_ratio(i), 1.0 + 1.0e-7);
    end

    fprintf('  [Diagonal Trot & Friction Pyramid] FL F_grf = [%.1f, %.1f, %.1f] N (mu_util = %.3f)\n', ...
        F_grf(1, 1), F_grf(2, 1), F_grf(3, 1), info.friction_ratio(1));
end

function test_gait_swing_controller_trot_and_apex(testCase)
    p = params();
    q_legs  = repmat(p.q_leg_nom, 1, 4);
    dq_legs = zeros(2, 4);

    % 1. Stationary 4-leg stance mode (gait_mode = 0)
    [cm_stand, tau_sw_stand] = gait_swing_controller(0.5, q_legs, dq_legs, [0;0;0], [0;0;0], 0, p);
    verifyEqual(testCase, cm_stand, [1; 1; 1; 1]);
    verifyEqual(testCase, tau_sw_stand, zeros(2, 4), 'AbsTol', 1.0e-12);

    % 2. First half-cycle of diagonal trot: Pair A (FL, RR) stance, Pair B (FR, RL) at mid-swing apex
    t_apex_B = p.t_stand_init + 0.25 * p.T_gait;
    [cm_B, tau_sw_B, p_des_B, v_des_B] = gait_swing_controller( ...
        t_apex_B, q_legs, dq_legs, [0;0;0], [0;0;0], 1, p);
    verifyEqual(testCase, cm_B, [1; 0; 0; 1]);
    verifyEqual(testCase, tau_sw_B(:, [1, 4]), zeros(2, 2), 'AbsTol', 1.0e-12);
    verifyEqual(testCase, p_des_B(3, 2), -p.z_nom + p.h_swing, 'AbsTol', 1.0e-12);
    verifyEqual(testCase, p_des_B(3, 3), -p.z_nom + p.h_swing, 'AbsTol', 1.0e-12);
    verifyEqual(testCase, v_des_B(3, 2), 0.0, 'AbsTol', 1.0e-12);
    verifyGreaterThan(testCase, tau_sw_B(2, 2), 0.0); % Positive knee torque lifts foot upward

    % 3. Second half-cycle of diagonal trot: Pair B (FR, RL) stance, Pair A (FL, RR) at mid-swing apex
    t_apex_A = p.t_stand_init + 0.75 * p.T_gait;
    [cm_A, tau_sw_A, p_des_A] = gait_swing_controller( ...
        t_apex_A, q_legs, dq_legs, [0;0;0], [0;0;0], 1, p);
    verifyEqual(testCase, cm_A, [0; 1; 1; 0]);
    verifyEqual(testCase, tau_sw_A(:, [2, 3]), zeros(2, 2), 'AbsTol', 1.0e-12);
    verifyEqual(testCase, p_des_A(3, 1), -p.z_nom + p.h_swing, 'AbsTol', 1.0e-12);
    verifyEqual(testCase, p_des_A(3, 4), -p.z_nom + p.h_swing, 'AbsTol', 1.0e-12);
end

function test_stage4_raibert_edge_avoidance_and_pure_pursuit(testCase)
    p = params();

    % 1. Test step-edge avoidance when front hips are right near x_step_start = 1.20 m
    % Front hip offset is +0.35 m, so p_com_x = 0.84 m places front hips at x = 1.19 m
    p_com = [0.84; 0.0; p.z_nom + p.r_foot];
    v_com = [0.40; 0.0; 0.0];
    rpy   = [0.0; 0.0; 0.0];
    omega = [0.0; 0.0; 0.0];
    v_cmd = [0.40; 0.0; 0.0];

    [p_td_w, ~, z_com_des, h_swing_legs, edge_shifted] = ...
        raibert_foothold_planner(p_com, v_com, rpy, omega, v_cmd, 0.0, p);

    % Verify front feet touchdowns are shifted >= 0.05 m clear of x_step_start (1.20 m)
    for i = 1:4
        verifyGreaterThanOrEqual(testCase, abs(p_td_w(1, i) - p.x_step_start), 0.05 - 1.0e-9);
    end
    verifyTrue(testCase, any(edge_shifted(1:2)));
    verifyGreaterThanOrEqual(testCase, h_swing_legs(1), 0.14);
    verifyGreaterThan(testCase, z_com_des, p.z_nom + p.r_foot);

    % 2. Test Pure-Pursuit Navigator at U-turn entry (x = 1.95 m, y = 0.0 m, psi = 0)
    [v_cmd_pp, wz_cmd_pp, ~, ~, ~, path_xy] = ...
        pure_pursuit_navigator(4.5, [p.x_uturn_center; 0.0; 0.58], [0.0; 0.0; 0.0], p);
    verifyGreaterThan(testCase, v_cmd_pp(1), 0.0);
    verifyEqual(testCase, v_cmd_pp(2), 0.0, 'AbsTol', 1.0e-12);
    verifyGreaterThan(testCase, wz_cmd_pp, 0.1); % Positive left-turn yaw rate into R=0.70m arc
    verifyEqual(testCase, max(path_xy(2, :)), 2.0 * p.R_uturn, 'AbsTol', 1.0e-6);
end

function test_stage5_arm_ik_and_mission_supervisor(testCase)
    p = params();

    % 1. Test 6-Axis Arm FK/IK accuracy (< 1e-6 m) and dynamic reaction wrench w_arm
    [q_des, dq_des, ~, p_ee_des, ~, w_arm, arm_info] = ...
        arm_ik_planner(p.p_ee_home_b, p.p_ee_pick_b, 1.60, 1.25, 1.95, eye(3), p.m_box, p);
    verifyEqual(testCase, numel(q_des), 6);
    verifyEqual(testCase, numel(dq_des), 6);
    verifyEqual(testCase, numel(w_arm), 6);
    verifyLessThan(testCase, arm_info.ik_error_m, 1.0e-6);
    verifyEqual(testCase, arm_info.p_ee_fk_b, p_ee_des, 'AbsTol', 1.0e-6);

    % 2. Test Mission Supervisor across all 5 phases & proximity suction latch (< 0.015 m)
    p_com = [0.0; 0.0; p.z_nom + p.r_foot];
    v_com = [0.0; 0.0; 0.0];
    [ph1, gm1, q_pick, ~, ~, ~, suc1, grasp1, m_pay1, ~, ~, prox1] = ...
        mission_supervisor(1.95, p_com, v_com, eye(3), [0;0;0], q_des, zeros(6,1), p.box_init_pos, [0;0;0], p);
    verifyEqual(testCase, ph1, 1);
    verifyEqual(testCase, gm1, 0);
    verifyEqual(testCase, suc1, 1);
    verifyLessThan(testCase, prox1, 0.015);
    verifyEqual(testCase, grasp1, 1);
    verifyEqual(testCase, m_pay1, p.m_box, 'AbsTol', 1.0e-12);
    verifyEqual(testCase, numel(q_pick), 6);
end
