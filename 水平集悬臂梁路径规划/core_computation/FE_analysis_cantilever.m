function [U, K, F] = FE_analysis_cantilever(nelx, nely, theta_e, E_L, E_T, nu_LT, G_LT, t, F_mag, dx, dy)
    % 悬臂梁边界条件下的有限元分析
    % 返回值：U（位移），K（刚度矩阵），F（载荷向量）

    ndof = 2*(nelx+1)*(nely+1);
    K = sparse(ndof, ndof);
    F = sparse(ndof, 1);

    for elx = 1:nelx
        for ely = 1:nely
            n1 = (nely+1)*(elx-1) + ely;
            n2 = (nely+1)*elx + ely;
            n3 = n2 + 1;
            n4 = n1 + 1;
            edof = [2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n3-1, 2*n3, 2*n4-1, 2*n4];
            Ke = element_stiffness(theta_e(ely, elx), E_L, E_T, nu_LT, G_LT, t, dx, dy);
            K(edof, edof) = K(edof, edof) + Ke;
        end
    end

    F(2*(nely+1)*nelx+nely+2,1) = F_mag;
    fixeddofs=1:2*(nely+1);
    

    alldofs = 1:ndof;
    freedofs = setdiff(alldofs, fixeddofs);
    K_ff = K(freedofs, freedofs) + 1e-10*speye(length(freedofs));
    U = zeros(ndof, 1);
    U(freedofs) = K_ff \ F(freedofs);
end

