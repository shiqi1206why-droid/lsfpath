function Ke = element_stiffness(theta, E_L, E_T, nu_LT, G_LT, t, dx, dy)
    % 计算正交各向异性材料在角度 theta 下的单元刚度矩阵

    nu_TL = nu_LT * E_T / E_L;
    Q11 = E_L / (1 - nu_LT * nu_TL);
    Q22 = E_T / (1 - nu_LT * nu_TL);
    Q12 = nu_LT * E_T / (1 - nu_LT * nu_TL);
    Q66 = G_LT;

    Q = [Q11, Q12, 0;
         Q12, Q22, 0;
         0,   0,   Q66];

    c = cos(theta);
    s = sin(theta);
    T = [c^2,  s^2,   2*s*c;
         s^2,  c^2,  -2*s*c;
         -s*c, s*c,  c^2-s^2];

    C = T' * Q * T;
 
    gauss_points = [-1/sqrt(3), 1/sqrt(3)];
    gauss_weights = [1, 1];
    Ke = zeros(8, 8);

    for i = 1:2
        for j = 1:2
            xi = gauss_points(i);
            eta = gauss_points(j);
            w = gauss_weights(i) * gauss_weights(j);

            dN_dxi = 0.25 * [-(1-eta), (1-eta), (1+eta), -(1+eta)];
            dN_deta = 0.25 * [-(1-xi), -(1+xi), (1+xi), (1-xi)];

            J = [dx/2, 0; 0, dy/2];
            detJ = det(J);
            dN_dxy = J \ [dN_dxi; dN_deta];

            B = zeros(3, 8);
            for k = 1:4
                B(1, 2*k-1) = dN_dxy(1, k);
                B(2, 2*k) = dN_dxy(2, k);
                B(3, 2*k-1) = dN_dxy(2, k);
                B(3, 2*k) = dN_dxy(1, k);
            end

            Ke = Ke + B' * C * B * detJ * w * t;
        end
    end
end

