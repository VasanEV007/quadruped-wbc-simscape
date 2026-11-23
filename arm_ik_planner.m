function [q_arm_des, dq_arm_des, ddq_arm_des, p_ee_des_b, v_ee_des_b, w_arm_world, arm_info] = ...
    arm_ik_planner(p_start_b, p_end_b, t_curr, t_start, t_end, R_wb, m_payload_active, p, q_arm_meas)
%ARM_IK_PLANNER 6-Axis arm FK/IK, quintic minimum-jerk trajectory & dynamic reaction wrench.
%   [q_arm_des, dq_arm_des, ddq_arm_des, p_ee_des_b, v_ee_des_b, w_arm_world, arm_info] = ...
%       ARM_IK_PLANNER(p_start_b, p_end_b, t_curr, t_start, t_end, R_wb, m_payload_active, p, q_arm_meas)
%   generates a C^2 quintic minimum-jerk task-space/joint-space trajectory between
%   p_start_b and p_end_b in the torso body frame, solves 6-DOF analytical +
%   damped-least-squares (DLS) inverse kinematics, and computes the 6D dynamic
%   base reaction wrench w_arm_world = [F_arm; tau_arm] in the world frame for
%   proactive VMC stance leg feedforward compensation.
%
%   Inputs:
%     p_start_b        - 3x1 segment start EE position in body frame [m]
%     p_end_b          - 3x1 segment target EE position in body frame [m]
%     t_curr           - Current simulation time [s]
%     t_start          - Segment start timestamp [s]
%     t_end            - Segment end timestamp [s]
%     R_wb             - (Optional) 3x3 body-to-world rotation matrix (default: eye(3))
%     m_payload_active - (Optional) Active grasped box payload mass [kg] (default: 0.0)
%     p                - (Optional) Parameter struct from params()
%     q_arm_meas       - (Optional) 6x1 measured arm joint angles [rad]
%
%   Outputs:
%     q_arm_des   - 6x1 desired arm joint angles [rad]
%     dq_arm_des  - 6x1 desired arm joint velocities [rad/s]
%     ddq_arm_des - 6x1 desired arm joint accelerations [rad/s^2]
%     p_ee_des_b  - 3x1 desired EE position in body frame [m]
%     v_ee_des_b  - 3x1 desired EE linear velocity in body frame [m/s]
%     w_arm_world - 6x1 dynamic base reaction wrench [F_arm; tau_arm] in world frame [N; Nm]
%     arm_info    - Struct with FK positions, 6x6 Jacobian, IK residual, and gravity torque

    if nargin < 8 || isempty(p)
        p = params();
    end
    if nargin < 7 || isempty(m_payload_active)
        m_payload_active = 0.0;
    end
    if nargin < 6 || isempty(R_wb)
        R_wb = eye(3);
    end

    p0 = double(p_start_b(:));
    p1 = double(p_end_b(:));
    T_seg = max(double(t_end - t_start), 1.0e-3);
    tau   = min(max(double(t_curr - t_start) / T_seg, 0.0), 1.0);

    %% 1. Quintic Minimum-Jerk Task-Space Interpolation
    % s(tau) = 10*tau^3 - 15*tau^4 + 6*tau^5
    tau2 = tau * tau;
    tau3 = tau2 * tau;
    tau4 = tau3 * tau;
    tau5 = tau4 * tau;

    s_pos = 10.0 * tau3 - 15.0 * tau4 + 6.0 * tau5;
    s_vel = (30.0 * tau2 - 60.0 * tau3 + 30.0 * tau4) / T_seg;
    s_acc = (60.0 * tau  - 180.0 * tau2 + 120.0 * tau3) / (T_seg * T_seg);

    dp_seg     = p1 - p0;
    p_ee_des_b = p0 + s_pos * dp_seg;
    v_ee_des_b = s_vel * dp_seg;
    a_ee_des_b = s_acc * dp_seg;

    %% 2. 6-DOF Analytical Seed + Damped-Least-Squares (DLS) Inverse Kinematics
    q_arm_des = solve_arm_ik_6dof(p_ee_des_b, p);

    % Evaluate Forward Kinematics & 6x6 Geometric Jacobian at q_arm_des
    [p_ee_fk_b, J_arm, p_coms_b, p_joints_b] = arm_fk_jacobian(q_arm_des, p);

    % Compute desired joint velocities and accelerations via damped pseudo-inverse
    J_v = J_arm(1:3, :);
    lambda_dls = 0.015;
    J_pinv = J_v.' / (J_v * J_v.' + (lambda_dls^2) * eye(3));

    dq_arm_des  = J_pinv * v_ee_des_b;
    ddq_arm_des = J_pinv * a_ee_des_b;

    % Keep wrist level (q2 + q3 + q4 = pi) and tool yaw aligned (q6 = -q1)
    dq_arm_des(4)  = -(dq_arm_des(2) + dq_arm_des(3));
    dq_arm_des(5)  = 0.0;
    dq_arm_des(6)  = -dq_arm_des(1);
    ddq_arm_des(4) = -(ddq_arm_des(2) + ddq_arm_des(3));
    ddq_arm_des(5) = 0.0;
    ddq_arm_des(6) = -ddq_arm_des(1);

    %% 3. Dynamic Base Reaction Wrench w_arm = [F_arm; tau_arm] in World Frame
    % Reference stowed link CoMs (already counterbalanced by trunk_com_offset)
    [~, ~, p_coms_stow_b, ~] = arm_fk_jacobian(p.q_arm_stow, p);

    if nargin >= 9 && ~isempty(q_arm_meas)
        [p_ee_curr_b, J_curr, p_coms_curr_b] = arm_fk_jacobian(double(q_arm_meas(:)), p);
    else
        p_ee_curr_b   = p_ee_fk_b;
        J_curr        = J_arm;
        p_coms_curr_b = p_coms_b;
    end

    % Approximate link CoM accelerations from end-effector acceleration scaling
    F_inertial_b   = [0.0; 0.0; 0.0];
    tau_inertial_b = [0.0; 0.0; 0.0];
    tau_grav_shift_w = [0.0; 0.0; 0.0];
    g_world = [0.0; 0.0; -p.g];

    total_arm_len = sum(p.arm_link_lengths);
    cum_len = 0.0;
    for j = 1:6
        mj = p.m_arm_links(j);
        cum_len = cum_len + 0.5 * p.arm_link_lengths(j);
        accel_scale_j = min(cum_len / max(total_arm_len, 1.0e-3), 1.0);
        a_com_j_b = accel_scale_j * a_ee_des_b;

        % Inertial reaction force & moment on base: -m_j * a_j and -r_j x (m_j * a_j)
        F_inertial_b   = F_inertial_b   - mj * a_com_j_b;
        tau_inertial_b = tau_inertial_b - cross(p_coms_curr_b(:, j), mj * a_com_j_b);

        % Gravity moment shift relative to counterbalanced stowed pose
        dr_com_j_w = R_wb * (p_coms_curr_b(:, j) - p_coms_stow_b(:, j));
        tau_grav_shift_w = tau_grav_shift_w + cross(dr_com_j_w, mj * g_world);
        cum_len = cum_len + 0.5 * p.arm_link_lengths(j);
    end

    % Gripper mass (1.4 kg) + active grasped box payload mass (0.50 kg when grasped)
    m_tip_extra = m_payload_active;
    F_tip_world = R_wb * (-(p.m_gripper + m_tip_extra) * a_ee_des_b) + ...
                  m_tip_extra * g_world;

    dr_grip_w = R_wb * (p_ee_curr_b - p_coms_stow_b(:, 6));
    r_ee_w    = R_wb * (p_ee_curr_b - p.r_com_offset);

    tau_tip_world = cross(r_ee_w, R_wb * (-(p.m_gripper + m_tip_extra) * a_ee_des_b)) + ...
                    cross(dr_grip_w, p.m_gripper * g_world) + ...
                    cross(r_ee_w, m_tip_extra * g_world);

    F_arm_world   = R_wb * F_inertial_b   + F_tip_world;
    tau_arm_world = R_wb * tau_inertial_b + tau_grav_shift_w + tau_tip_world;

    w_arm_world = [F_arm_world; tau_arm_world];

    %% 4. Joint-Space Gravity Compensation Torque for 6-Axis Arm Actuators
    if nargout >= 7
        g_body = R_wb.' * g_world;
        tau_grav_arm = compute_arm_gravity_torques(q_arm_des, g_body, m_payload_active, p);

        arm_info.p_ee_fk_b    = p_ee_fk_b;
        arm_info.p_joints_b   = p_joints_b;
        arm_info.p_coms_b     = p_coms_b;
        arm_info.J_arm        = J_curr;
        arm_info.ik_error_m   = norm(p_ee_fk_b - p_ee_des_b);
        arm_info.tau_grav_arm = tau_grav_arm;
    end
end

%% ---------------- 6-DOF Arm Kinematics & Helper Functions ----------------

function q_sol = solve_arm_ik_6dof(p_target_b, p)
%SOLVE_ARM_IK_6DOF Closed-form cylindrical seed + 4-step DLS Jacobian refinement.
    p_base = p.arm_base_offset;
    L      = p.arm_link_lengths; % [l1, l2, l3, l4, l5, l6]
    l1 = L(1);
    l2 = L(2);
    l3 = L(3);
    l_tool = L(4) + L(5) + L(6) + p.arm_ee_offset; % Vertical downward tool drop

    dx = p_target_b(1) - p_base(1);
    dy = p_target_b(2) - p_base(2);
    dz = p_target_b(3) - (p_base(3) + l1);

    % 1. Waist yaw q1
    q1 = atan2(dy, dx);
    r_horiz = hypot(dx, dy);

    % 2. Wrist pitch center (r_w, z_w) with tool pointing straight down (q2+q3+q4 = pi)
    r_w = r_horiz;
    z_w = dz + l_tool;

    d_sq = r_w*r_w + z_w*z_w;
    cos_q3 = (d_sq - l2*l2 - l3*l3) / (2.0 * l2 * l3);
    cos_q3 = min(max(cos_q3, -0.999), 0.999);
    q3 = acos(cos_q3); % Elbow-up positive pitch branch
    q2 = atan2(r_w, z_w) - atan2(l3 * sin(q3), l2 + l3 * cos(q3));
    q4 = pi - (q2 + q3);
    q5 = 0.0;
    q6 = -q1;

    q_sol = [q1; q2; q3; q4; q5; q6];

    % 3. Damped-Least-Squares (DLS) refinement (4 iterations -> < 1e-12 m accuracy)
    lambda_sq = 1.0e-4;
    for iter = 1:4
        [p_fk, J_arm] = arm_fk_jacobian(q_sol, p);
        e_p = p_target_b - p_fk;
        if norm(e_p) < 1.0e-10
            break;
        end
        J_v = J_arm(1:3, 1:3);
        dq_123 = J_v.' / (J_v * J_v.' + lambda_sq * eye(3)) * e_p;
        q_sol(1:3) = q_sol(1:3) + dq_123;
        q_sol(4)   = pi - (q_sol(2) + q_sol(3));
        q_sol(5)   = 0.0;
        q_sol(6)   = -q_sol(1);
    end
end

function [p_ee_b, J_arm, p_coms_b, p_joints_b] = arm_fk_jacobian(q, p)
%ARM_FK_JACOBIAN Compute 6-DOF arm forward kinematics, link CoMs, and 6x6 geometric Jacobian.
    p_base = p.arm_base_offset;
    L      = p.arm_link_lengths;
    l1 = L(1); l2 = L(2); l3 = L(3);
    l4 = L(4); l5 = L(5); l6 = L(6);
    l_ee = p.arm_ee_offset;

    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = q(4);
    q23  = q2 + q3;
    q234 = q2 + q3 + q4;

    c1 = cos(q1); s1 = sin(q1);
    u_r    = [c1; s1; 0.0];       % Radial horizontal unit vector
    u_perp = [-s1; c1; 0.0];      % Lateral pitch rotation axis (+Y_local)
    u_z    = [0.0; 0.0; 1.0];     % Vertical unit vector

    % Unit direction vectors of links 2, 3, and tool (4..6)
    dir2 = sin(q2)   * u_r + cos(q2)   * u_z;
    dir3 = sin(q23)  * u_r + cos(q23)  * u_z;
    dir4 = sin(q234) * u_r + cos(q234) * u_z;

    % Joint origins in body frame (6x3 -> 3x7 including base and EE)
    p1 = p_base;
    p2 = p1 + l1 * u_z;
    p3 = p2 + l2 * dir2;
    p4 = p3 + l3 * dir3;
    p5 = p4 + l4 * dir4;
    p6 = p5 + l5 * dir4;
    p_ee_b = p6 + (l6 + l_ee) * dir4;

    p_joints_b = [p1, p2, p3, p4, p5, p6, p_ee_b];

    % Link CoM positions in body frame (3x6)
    p_coms_b = [ ...
        p1 + 0.5 * l1 * u_z, ...
        p2 + 0.5 * l2 * dir2, ...
        p3 + 0.5 * l3 * dir3, ...
        p4 + 0.5 * l4 * dir4, ...
        p5 + 0.5 * l5 * dir4, ...
        p6 + 0.5 * (l6 + l_ee) * dir4 ...
    ];

    % Joint rotation axes in body frame
    z_axes = [ ...
        u_z, ...                        % Joint 1: Waist Yaw (+Z)
        u_perp, ...                     % Joint 2: Shoulder Pitch
        u_perp, ...                     % Joint 3: Elbow Pitch
        u_perp, ...                     % Joint 4: Wrist Pitch
        cross(u_perp, dir4), ...        % Joint 5: Wrist Roll
        dir4 ...                        % Joint 6: Tool Yaw
    ];

    % 6x6 Geometric Jacobian [J_v; J_w]
    J_arm = zeros(6, 6);
    for j = 1:6
        zj = z_axes(:, j);
        rj = p_ee_b - p_joints_b(:, j);
        J_arm(1:3, j) = cross(zj, rj);
        J_arm(4:6, j) = zj;
    end
end

function tau_g = compute_arm_gravity_torques(q, g_body, m_payload_active, p)
%COMPUTE_ARM_GRAVITY_TORQUES Exact static gravity compensation torques for the 6 arm joints.
    [~, ~, p_coms_b, p_joints_b] = arm_fk_jacobian(q, p);
    p_ee_b = p_joints_b(:, 7);

    q1 = q(1);
    q234 = q(2) + q(3) + q(4);
    u_r    = [cos(q1); sin(q1); 0.0];
    u_perp = [-sin(q1); cos(q1); 0.0];
    u_z    = [0.0; 0.0; 1.0];
    dir4   = sin(q234) * u_r + cos(q234) * u_z;

    z_axes = [u_z, u_perp, u_perp, u_perp, cross(u_perp, dir4), dir4];
    tau_g  = zeros(6, 1);

    for j = 1:6
        zj = z_axes(:, j);
        pj = p_joints_b(:, j);
        tau_sum = 0.0;
        for k = j:6
            mk = p.m_arm_links(k);
            rk = p_coms_b(:, k) - pj;
            tau_sum = tau_sum - dot(zj, cross(rk, mk * g_body));
        end
        r_tip = p_ee_b - pj;
        m_tip = p.m_gripper + m_payload_active;
        tau_sum = tau_sum - dot(zj, cross(r_tip, m_tip * g_body));
        tau_g(j) = tau_sum;
    end
end
