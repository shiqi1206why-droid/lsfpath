function [U, K, F] = FE_analysis_cantilever(nelx, nely, theta_e, E_L, E_T, nu_LT, G_LT, t, F_mag, dx, dy, material_mask)
    % 悬臂梁边界条件下的有限元分析
    % 返回值：U（位移），K（刚度矩阵），F（载荷向量）

    if nargin < 12 || isempty(material_mask)
        material_mask = true(nely, nelx);
    end
    if ~isequal(size(material_mask), [nely, nelx])
        error('material_mask尺寸错误：期望[%d,%d]，实际[%d,%d]。', ...
            nely, nelx, size(material_mask, 1), size(material_mask, 2));
    end
    material_mask = logical(material_mask);

    ndof = 2*(nelx+1)*(nely+1);
    K = sparse(ndof, ndof);
    F = sparse(ndof, 1);

    for elx = 1:nelx
        for ely = 1:nely
            if ~material_mask(ely, elx)
                continue;
            end
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
    K_ff = K(freedofs, freedofs);
    diag_vals = abs(diag(K_ff));
    diag_vals = diag_vals(isfinite(diag_vals) & diag_vals > 0);
    if isempty(diag_vals)
        diag_scale = 1.0;
    else
        diag_scale = median(diag_vals);
    end
    reg = max(1e-10, 1e-8 * diag_scale);
    K_ff = K_ff + reg * speye(length(freedofs));
    K(freedofs, freedofs) = K_ff;  % 保持后续能量一致性检查与求解矩阵一致

    U = zeros(ndof, 1);
    
    warn_near = warning('query', 'MATLAB:nearlySingularMatrix');
    warn_sing = warning('query', 'MATLAB:singularMatrix');
    warning('off', 'MATLAB:nearlySingularMatrix');
    warning('off', 'MATLAB:singularMatrix');
    try
        U(freedofs) = K_ff \ F(freedofs);
    catch ME
        warning(warn_near.state, 'MATLAB:nearlySingularMatrix');
        warning(warn_sing.state, 'MATLAB:singularMatrix');
        rethrow(ME);
    end
    warning(warn_near.state, 'MATLAB:nearlySingularMatrix');
    warning(warn_sing.state, 'MATLAB:singularMatrix');
end

