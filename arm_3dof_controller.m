function [q_arm_des, dq_arm_des, ddq_arm_des, tau_arm_cmd, w_arm_world, p_ee_world_act, ee_err_world_m, arm3_info] = ...
    arm_3dof_controller(p_ee_des_w, v_ee_des_w, a_ee_des_w, p_com_w, v_com_w, R_wb, omega_w, ...
                        q_arm_meas, dq_arm_meas, m_payload_active, p)
%ARM_3DOF_CONTROLLER Decoupled 3-DOF analytical IK, base-twist velocity compensation & computed torque.
%   [q_arm_des, dq_arm_des, ddq_arm_des, tau_arm_cmd, w_arm_world, p_ee_world_act, ee_err_world_m, arm3_info] = ...
%       ARM_3DOF_CONTROLLER(p_ee_des_w, v_ee_des_w, a_ee_des_w, p_com_w, v_com_w, R_wb, omega_w, ...
%                           q_arm_meas, dq_arm_meas, m_payload_active, p)
%   implements Stage 6 decoupled 3-DOF arm control (Joint 1: base yaw, Joint 2:
%   shoulder pitch, Joint 3: elbow pitch, with wrist Joints 4-6 orientation-locked):
%     1. World-to-Mount Frame Transformation with active 6-DOF base-twist cancellation:
%          p_ee_d^B = R_wb^T * (p_ee_d^W - p_b) - p_mount^B
%          v_ee_d^B = R_wb^T * (v_ee_d^W - v_b - omega_w x (p_ee_d^W - p_b))
%     2. Closed-Form 3-DOF Analytical Elbow-Up IK & 3x3 Jacobian J_3dof with
%        smooth workspace boundary clamping and DLS velocity mapping.
%     3. Computed-Torque + Task-Space/Joint-Space PD + Gravity/Coriolis compensation
%        and 6D dynamic base reaction wrench w_arm_world = [F_arm; tau_arm].

    if nargin < 11 || isempty(p)
        p = params();
    end
    if nargin < 10 || isempty(m_payload_active)
        m_payload_active = 0.0;
    end

    p_ee_w_d = double(p_ee_des_w(:));
    v_ee_w_d = double(v_ee_des_w(:));
    a_ee_w_d = double(a_ee_des_w(:));
    p_b      = double(p_com_w(:));
    v_b      = double(v_com_w(:));
    w_b      = double(omega_w(:));
    q_meas   = double(q_arm_meas(:));
    dq_meas  = double(dq_arm_meas(:));

    p_mount_b = p.arm_base_offset(:);

    %% 1. World-to-Mount Frame Transformation & Base-Twist Cancellation
    % Position in mount frame: p_ee_d^B = R_wb^T * (p_ee_d^W - p_b) - p_mount^B
    r_ee_from_com_w = p_ee_w_d - p_b;
    p_ee_d_body     = R_wb.' * r_ee_from_com_w;
    p_ee_d_mount    = p_ee_d_body - p_mount_b;

    % Velocity in mount frame canceling base linear & angular velocity twist:
    % v_ee_d^B = R_wb^T * (v_ee_d^W - v_b - omega_w x (p_ee_d^W - p_b))
    v_rel_w      = v_ee_w_d - v_b - cross(w_b, r_ee_from_com_w);
    v_ee_d_mount = R_wb.' * v_rel_w;

    % Acceleration in mount frame canceling base centripetal & Coriolis acceleration
    a_coriolis_w = cross(w_b, cross(w_b, r_ee_from_com_w)) + 2.0 * cross(w_b, v_rel_w);
    a_ee_d_mount = R_wb.' * (a_ee_w_d - a_coriolis_w);

    %% 2. Closed-Form 3-DOF Analytical Elbow-Up IK with Smooth Workspace Clamping
    L = p.arm_link_lengths; % [l1, l2, l3, l4, l5, l6]
    l1 = L(1);
    l2 = L(2);
    l3 = L(3);
    l_tool = L(4) + L(5) + L(6) + p.arm_ee_offset; % Rigid wrist + suction cup vertical drop

    % Vector from shoulder pitch center [0; 0; l1] to wrist center [x; y; z + l_tool] in mount frame
    dp_sh = [ ...
        p_ee_d_mount(1); ...
        p_ee_d_mount(2); ...
        p_ee_d_mount(3) - l1 + l_tool ...
    ];

    % Smooth C^inf workspace boundary clamping (avoids full-stretch / folded singularities)
    r_max_hard = 0.992 * (l2 + l3);
    r_max_soft = 0.910 * (l2 + l3);
    r_min_hard = max(1.08 * abs(l2 - l3), 0.045);

    d_raw = norm(dp_sh);
    if d_raw > r_max_soft
        d_clamped = r_max_soft + (r_max_hard - r_max_soft) * ...
                    tanh((d_raw - r_max_soft) / max(r_max_hard - r_max_soft, 1.0e-4));
        dp_sh = dp_sh * (d_clamped / d_raw);
    elseif d_raw < r_min_hard && d_raw > 1.0e-6
        dp_sh = dp_sh * (r_min_hard / d_raw);
    end

    dx = dp_sh(1);
    dy = dp_sh(2);
    dz = dp_sh(3);

    % Joint 1: Base Yaw theta_1
    theta1_d = atan2(dy, dx);
    r_horiz  = hypot(dx, dy);

    % Joint 3: Elbow Pitch theta_3 (elbow-up positive branch)
    d_sq = r_horiz^2 + dz^2;
    cos_t3 = (d_sq - l2^2 - l3^2) / (2.0 * l2 * l3);
    cos_t3 = min(max(cos_t3, -0.996), 0.996);
    theta3_d = acos(cos_t3);

    % Joint 2: Shoulder Pitch theta_2
    theta2_d = atan2(r_horiz, dz) - atan2(l3 * sin(theta3_d), l2 + l3 * cos(theta3_d));

    % Wrist Joints 4-6: Orientation-compensated / locked so vacuum cup stays vertical (-Z_B)
    theta4_d = pi - (theta2_d + theta3_d);
    theta5_d = 0.0;
    theta6_d = -theta1_d;

    q_arm_des = [theta1_d; theta2_d; theta3_d; theta4_d; theta5_d; theta6_d];

    %% 3. Analytical 3x3 Jacobian J_3dof & Damped-Least-Squares Velocity/Accel IK
    [p_ee_fk_mount, J_3dof, Jdot_3dof] = fk_jacobian_3dof( ...
        q_arm_des(1:3), dq_meas(1:3), l1, l2, l3, l_tool);

    % Damped pseudo-inverse J_3dof^dagger
    lambda_dls = 0.012;
    J_dag = J_3dof.' / (J_3dof * J_3dof.' + (lambda_dls^2) * eye(3));

    dtheta_123_d  = J_dag * v_ee_d_mount;
    ddtheta_123_d = J_dag * (a_ee_d_mount - Jdot_3dof * dtheta_123_d);

    % Clamp joint velocities and accelerations to smooth physical bounds
    dtheta_123_d  = min(max(dtheta_123_d,  -6.0), 6.0);
    ddtheta_123_d = min(max(ddtheta_123_d, -45.0), 45.0);

    dq_arm_des = [ ...
        dtheta_123_d(1); ...
        dtheta_123_d(2); ...
        dtheta_123_d(3); ...
        -(dtheta_123_d(2) + dtheta_123_d(3)); ...
        0.0; ...
        -dtheta_123_d(1) ...
    ];

    ddq_arm_des = [ ...
        ddtheta_123_d(1); ...
        ddtheta_123_d(2); ...
        ddtheta_123_d(3); ...
        -(ddtheta_123_d(2) + ddtheta_123_d(3)); ...
        0.0; ...
        -ddtheta_123_d(1) ...
    ];

    %% 4. Measured Forward Kinematics & World-Frame Stabilization Error
    [p_ee_meas_mount, J_3dof_meas] = fk_jacobian_3dof( ...
        q_meas(1:3), dq_meas(1:3), l1, l2, l3, l_tool);
    v_ee_meas_mount = J_3dof_meas * dq_meas(1:3);

    p_ee_meas_body  = p_mount_b + p_ee_meas_mount;
    p_ee_world_act  = p_b + R_wb * p_ee_meas_body;
    ee_err_world_m  = norm(p_ee_world_act - p_ee_w_d);

    %% 5. Computed-Torque + Task-Space PD + Gravity/Coriolis Control (Joints 1-3 & Locked 4-6)
    g_world = [0.0; 0.0; -p.g];
    g_body  = R_wb.' * g_world;

    % Effective distal mass at wrist (Links 4-6 + Vacuum Gripper + Active Payload Box)
    m_distal = sum(p.m_arm_links(4:6)) + p.m_gripper + m_payload_active;
    m2 = p.m_arm_links(2);
    m3 = p.m_arm_links(3);

    % Diagonal approximate 3-DOF joint inertia matrix M_3dof(theta) [kg*m^2]
    t2  = q_meas(2);
    t23 = q_meas(2) + q_meas(3);
    r_arm_hor = abs(l2 * sin(t2)) + abs(l3 * sin(t23));
    I11 = 0.015 + (0.25*m2*l2^2)*sin(t2)^2 + (m3 + m_distal) * max(r_arm_hor^2, 0.02);
    I22 = 0.018 + 0.33*m2*l2^2 + m3*(l2^2 + 0.25*l3^2 + l2*l3*cos(q_meas(3))) + ...
          m_distal*(l2^2 + l3^2 + 2.0*l2*l3*cos(q_meas(3)));
    I33 = 0.010 + 0.33*m3*l3^2 + m_distal*l3^2;
    M_3dof = diag([I11, I22, I33]);

    % Centrifugal / Coriolis torque estimate C_3dof(theta, dtheta)
    h_cor = -(0.5*m3 + m_distal) * l2 * l3 * sin(q_meas(3));
    tau_cor_123 = [ ...
        2.0 * (m3 + m_distal) * r_arm_hor * (l2*cos(t2)*dq_meas(2) + l3*cos(t23)*(dq_meas(2)+dq_meas(3))) * dq_meas(1); ...
        h_cor * (2.0 * dq_meas(2) * dq_meas(3) + dq_meas(3)^2); ...
       -h_cor * (dq_meas(2)^2) ...
    ];

    % Exact 6-joint static gravity torque vector
    [tau_grav_6, p_coms_curr_b] = compute_3dof_arm_gravity(q_meas, g_body, m_payload_active, p);

    % Task-space Cartesian PD wrench mapped via J_3dof^T
    Kp_task = diag([420.0, 420.0, 520.0]);
    Kd_task = diag([38.0,  38.0,  45.0]);
    F_task_mount = Kp_task * (p_ee_d_mount - p_ee_meas_mount) + ...
                   Kd_task * (v_ee_d_mount - v_ee_meas_mount);
    tau_task_123 = J_3dof_meas.' * F_task_mount;

    % Joint-space computed-torque + PD tracking for Joints 1-3
    tau_ff_123 = M_3dof * ddtheta_123_d + tau_cor_123 + tau_grav_6(1:3);
    tau_pd_123 = p.Kp_arm(1:3) .* (q_arm_des(1:3) - q_meas(1:3)) + ...
                 p.Kd_arm(1:3) .* (dq_arm_des(1:3) - dq_meas(1:3));
    tau_123 = tau_ff_123 + tau_task_123 + tau_pd_123;

    % Stiff PD lock + gravity compensation on Wrist Joints 4-6
    tau_456 = tau_grav_6(4:6) + ...
              1.5 * p.Kp_arm(4:6) .* (q_arm_des(4:6) - q_meas(4:6)) + ...
              1.5 * p.Kd_arm(4:6) .* (dq_arm_des(4:6) - dq_meas(4:6));

    tau_arm_cmd = min(max([tau_123; tau_456], -p.tau_arm_max), p.tau_arm_max);

    %% 6. Dynamic Base Reaction Wrench w_arm_world = [F_arm; tau_arm] for VMC Feedforward
    [~, p_coms_stow_b] = compute_3dof_arm_gravity(p.q_arm_stow, [0; 0; -p.g], 0.0, p);

    a_ee_body = R_wb.' * a_ee_w_d;
    F_inertial_b     = [0.0; 0.0; 0.0];
    tau_inertial_b   = [0.0; 0.0; 0.0];
    tau_grav_shift_w = [0.0; 0.0; 0.0];

    total_len = l1 + l2 + l3 + l_tool;
    cum_len   = 0.0;
    for j = 1:6
        mj = p.m_arm_links(j);
        cum_len = cum_len + 0.5 * L(j);
        sc_j = min(cum_len / max(total_len, 1.0e-3), 1.0);
        a_cj_b = sc_j * a_ee_body;

        F_inertial_b   = F_inertial_b   - mj * a_cj_b;
        tau_inertial_b = tau_inertial_b - cross(p_coms_curr_b(:, j), mj * a_cj_b);

        dr_w = R_wb * (p_coms_curr_b(:, j) - p_coms_stow_b(:, j));
        tau_grav_shift_w = tau_grav_shift_w + cross(dr_w, mj * g_world);
        cum_len = cum_len + 0.5 * L(j);
    end

    F_tip_w = -(p.m_gripper + m_payload_active) * a_ee_w_d + m_payload_active * g_world;
    dr_grip_w = R_wb * (p_ee_meas_body - p_coms_stow_b(:, 6));
    r_ee_w    = R_wb * (p_ee_meas_body - p.r_com_offset);

    tau_tip_w = cross(r_ee_w, -(p.m_gripper + m_payload_active) * a_ee_w_d) + ...
                cross(dr_grip_w, p.m_gripper * g_world) + ...
                cross(r_ee_w, m_payload_active * g_world);

    w_arm_world = [ ...
        R_wb * F_inertial_b   + F_tip_w; ...
        R_wb * tau_inertial_b + tau_grav_shift_w + tau_tip_w ...
    ];

    if nargout >= 8
        arm3_info.p_ee_d_mount   = p_ee_d_mount;
        arm3_info.v_ee_d_mount   = v_ee_d_mount;
        arm3_info.p_ee_fk_mount  = p_ee_fk_mount;
        arm3_info.J_3dof         = J_3dof;
        arm3_info.cond_J3        = cond(J_3dof);
    end
end

%% ---------------- 3-DOF Kinematic & Dynamic Helpers ----------------

function [p_ee_mount, J_3dof, Jdot_3dof] = fk_jacobian_3dof(theta_123, dtheta_123, l1, l2, l3, l_tool)
%FK_JACOBIAN_3DOF Exact 3-DOF forward kinematics, 3x3 Jacobian, and dJ/dt in mount frame.
    t1  = theta_123(1);
    t2  = theta_123(2);
    t3  = theta_123(3);
    t23 = t2 + t3;

    dt1  = dtheta_123(1);
    dt2  = dtheta_123(2);
    dt23 = dt2 + dtheta_123(3);

    c1 = cos(t1); s1 = sin(t1);
    s2 = sin(t2); c2 = cos(t2);
    s23 = sin(t23); c23 = cos(t23);

    R_rad  = l2 * s2 + l3 * s23;
    Z_vert = l2 * c2 + l3 * c23;

    p_ee_mount = [ ...
        c1 * R_rad; ...
        s1 * R_rad; ...
        l1 + Z_vert - l_tool ...
    ];

    J_3dof = [ ...
        -s1 * R_rad,  c1 * Z_vert,  c1 * l3 * c23; ...
         c1 * R_rad,  s1 * Z_vert,  s1 * l3 * c23; ...
         0.0,        -R_rad,       -l3 * s23       ...
    ];

    if nargout >= 3
        dR_dt = l2 * c2 * dt2 + l3 * c23 * dt23;
        dZ_dt = -l2 * s2 * dt2 - l3 * s23 * dt23;

        Jdot_3dof = [ ...
            -c1*dt1*R_rad - s1*dR_dt, -s1*dt1*Z_vert + c1*dZ_dt, -s1*dt1*l3*c23 - c1*l3*s23*dt23; ...
            -s1*dt1*R_rad + c1*dR_dt,  c1*dt1*Z_vert + s1*dZ_dt,  c1*dt1*l3*c23 - s1*l3*s23*dt23; ...
             0.0,                     -dR_dt,                    -l3*c23*dt23 ...
        ];
    end
end

function [tau_g, p_coms_b] = compute_3dof_arm_gravity(q, g_body, m_payload_active, p)
%COMPUTE_3DOF_ARM_GRAVITY Compute link CoMs in body frame and 6x1 gravity torque vector.
    p_base = p.arm_base_offset(:);
    L      = p.arm_link_lengths;
    l1 = L(1); l2 = L(2); l3 = L(3);
    l4 = L(4); l5 = L(5); l6 = L(6);
    l_ee = p.arm_ee_offset;

    q1   = q(1);
    q2   = q(2);
    q23  = q(2) + q(3);
    q234 = q(2) + q(3) + q(4);

    u_r    = [cos(q1); sin(q1); 0.0];
    u_perp = [-sin(q1); cos(q1); 0.0];
    u_z    = [0.0; 0.0; 1.0];

    dir2 = sin(q2)   * u_r + cos(q2)   * u_z;
    dir3 = sin(q23)  * u_r + cos(q23)  * u_z;
    dir4 = sin(q234) * u_r + cos(q234) * u_z;

    p1 = p_base;
    p2 = p1 + l1 * u_z;
    p3 = p2 + l2 * dir2;
    p4 = p3 + l3 * dir3;
    p5 = p4 + l4 * dir4;
    p6 = p5 + l5 * dir4;
    p_ee = p6 + (l6 + l_ee) * dir4;

    p_joints_b = [p1, p2, p3, p4, p5, p6, p_ee];
    p_coms_b = [ ...
        p1 + 0.5 * l1 * u_z, ...
        p2 + 0.5 * l2 * dir2, ...
        p3 + 0.5 * l3 * dir3, ...
        p4 + 0.5 * l4 * dir4, ...
        p5 + 0.5 * l5 * dir4, ...
        p6 + 0.5 * (l6 + l_ee) * dir4 ...
    ];

    z_axes = [u_z, u_perp, u_perp, u_perp, cross(u_perp, dir4), dir4];
    tau_g  = zeros(6, 1);

    for j = 1:6
        zj = z_axes(:, j);
        pj = p_joints_b(:, j);
        s_tau = 0.0;
        for k = j:6
            rk = p_coms_b(:, k) - pj;
            s_tau = s_tau - dot(zj, cross(rk, p.m_arm_links(k) * g_body));
        end
        r_tip = p_ee - pj;
        s_tau = s_tau - dot(zj, cross(r_tip, (p.m_gripper + m_payload_active) * g_body));
        tau_g(j) = s_tau;
    end
end
