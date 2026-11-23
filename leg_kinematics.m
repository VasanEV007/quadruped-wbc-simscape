function [p_foot_3d, J_3d, p_foot_2d, J_2d, kin_info] = leg_kinematics(q_leg, leg_id, p)
%LEG_KINEMATICS Forward kinematics and analytical 2D/3D leg Jacobians.
%   [p_foot_3d, J_3d, p_foot_2d, J_2d, kin_info] = LEG_KINEMATICS(q_leg, leg_id, p)
%   computes forward kinematics and analytical linear velocity Jacobians
%   for a 2-DOF planar quadruped leg (hip pitch q1 and knee pitch q2,
%   no hip abduction/adduction).
%
%   Coordinate Convention (Right-Handed Body Frame):
%     +X : Forward (longitudinal)
%     +Y : Left (lateral, joint pitch rotation axis)
%     +Z : Up (vertical)
%     q_leg = [0; 0] corresponds to the leg stretched straight down (-Z).
%     Positive pitch rotation about +Y follows the right-hand rule:
%       p_rel = [-L1*sin(q1) - L2*sin(q1+q2);
%                0;
%               -L1*cos(q1) - L2*cos(q1+q2)]
%
%   Inputs:
%     q_leg  - 2x1 joint angle vector [q_hip_pitch; q_knee_pitch] [rad]
%     leg_id - (Optional) Leg index: 1=FL, 2=FR, 3=RL, 4=RR, or 0 for
%              hip-local frame (default: 1)
%     p      - (Optional) Parameter struct from params()
%
%   Outputs:
%     p_foot_3d - 3x1 foot position [x; y; z] in body frame [m]
%     J_3d      - 3x2 analytical 3D Jacobian d(p_foot_3d)/dq [m/rad]
%     p_foot_2d - 2x1 sagittal foot position [x_rel; z_rel] relative to hip [m]
%     J_2d      - 2x2 analytical 2D sagittal Jacobian d(p_foot_2d)/dq [m/rad]
%     kin_info  - Struct with hip/knee positions, determinant, and condition number

    if nargin < 3 || isempty(p)
        p = params();
    end
    if nargin < 2 || isempty(leg_id)
        leg_id = 1;
    end

    validateattributes(q_leg, {'numeric'}, {'numel', 2, 'real', 'finite'}, ...
        'leg_kinematics', 'q_leg');

    q1 = double(q_leg(1));
    q2 = double(q_leg(2));
    q12 = q1 + q2;

    L1 = p.L_thigh;
    L2 = p.L_shank;

    % Trigonometric terms
    s1 = sin(q1);
    c1 = cos(q1);
    s12 = sin(q12);
    c12 = cos(q12);

    % Hip offset in body frame
    if leg_id >= 1 && leg_id <= 4
        p_hip_3d = p.hip_offsets(:, leg_id);
    elseif leg_id == 0
        p_hip_3d = [0.0; 0.0; 0.0];
    else
        error('leg_kinematics:InvalidLegId', ...
            'leg_id must be 1 (FL), 2 (FR), 3 (RL), 4 (RR), or 0 (hip-local).');
    end

    % 2D Sagittal Forward Kinematics (X-Z plane relative to hip joint)
    p_knee_2d = [ ...
        -L1 * s1; ...
        -L1 * c1  ...
    ];

    p_foot_2d = [ ...
        -L1 * s1 - L2 * s12; ...
        -L1 * c1 - L2 * c12  ...
    ];

    % 3D Forward Kinematics in Body Frame
    p_knee_3d = p_hip_3d + [p_knee_2d(1); 0.0; p_knee_2d(2)];
    p_foot_3d = p_hip_3d + [p_foot_2d(1); 0.0; p_foot_2d(2)];

    % 2D Analytical Sagittal Jacobian: d([x_rel; z_rel]) / d([q1; q2])
    J_2d = [ ...
        -L1 * c1 - L2 * c12, -L2 * c12; ...
         L1 * s1 + L2 * s12,  L2 * s12  ...
    ];

    % 3D Analytical Linear Velocity Jacobian: d([x; y; z]) / d([q1; q2])
    % Lateral Y row is strictly zero because there is no hip abduction joint.
    J_3d = [ ...
        J_2d(1, 1), J_2d(1, 2); ...
        0.0,        0.0;        ...
        J_2d(2, 1), J_2d(2, 2)  ...
    ];

    if nargout >= 5
        kin_info.p_hip_3d = p_hip_3d;
        kin_info.p_knee_3d = p_knee_3d;
        kin_info.p_knee_2d = p_knee_2d;
        kin_info.det_J2d = -L1 * L2 * sin(q2);
        kin_info.cond_J2d = cond(J_2d);
    end
end
