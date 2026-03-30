function strain_energy = compute_strain_energy(nelx, nely, U, theta_e, E_L, E_T, nu_LT, G_LT, t, dx, dy, material_mask)
    % 计算单元应变能密度

    if nargin < 12 || isempty(material_mask)
        material_mask = true(nely, nelx);
    end
    if ~isequal(size(material_mask), [nely, nelx])
        error('material_mask尺寸错误：期望[%d,%d]，实际[%d,%d]。', ...
            nely, nelx, size(material_mask, 1), size(material_mask, 2));
    end
    material_mask = logical(material_mask);

    strain_energy = zeros(nely, nelx);
    for ely = 1:nely
        for elx = 1:nelx
            if ~material_mask(ely, elx)
                strain_energy(ely, elx) = 0;
                continue;
            end
            n1 = (nely+1)*(elx-1) + ely;
            n2 = (nely+1)*elx + ely;
            n3 = n2 + 1;
            n4 = n1 + 1;
            nodes = [n1, n2, n3, n4];

            edof = [];
            for n = nodes
                edof = [edof, 2*n-1, 2*n]; %#ok<AGROW>
            end
            Ue = U(edof);

            Ke = element_stiffness(theta_e(ely, elx), E_L, E_T, nu_LT, G_LT, t, dx, dy);
            strain_energy(ely, elx) = 0.5 * Ue' * Ke * Ue / (dx * dy * t);
        end
    end
end

