function [C, dC_dtheta, S] = orthotropic_constitutive_matrix(theta, E_L, E_T, nu_LT, G_LT)
% 计算正交各向异性层合单元在角度theta下的面内本构矩阵及其角度导数

    nu_TL = nu_LT * E_T / E_L;
    denom = 1 - nu_LT * nu_TL;
    Q11 = E_L / denom;
    Q22 = E_T / denom;
    Q12 = nu_LT * E_T / denom;
    Q66 = G_LT;

    Q = [Q11, Q12, 0;
         Q12, Q22, 0;
         0,   0,   Q66];

    c = cos(theta);
    s = sin(theta);

    T = [c^2,  s^2,   2*s*c;
         s^2,  c^2,  -2*s*c;
         -s*c, s*c,  c^2-s^2];

    dT_dtheta = [-2*c*s,    2*c*s,    2*(c^2-s^2);
                  2*c*s,   -2*c*s,   -2*(c^2-s^2);
                 (s^2-c^2), (c^2-s^2), -4*c*s];

    C = T' * Q * T;
    dC_dtheta = dT_dtheta' * Q * T + T' * Q * dT_dtheta;

    rc = rcond(C);
    if ~isfinite(rc) || rc < 1e-12
        C = C + 1e-9 * eye(3);
    end

    if nargout >= 3
        S = C \ eye(3);
    end
end
