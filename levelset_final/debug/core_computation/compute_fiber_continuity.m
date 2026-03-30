function FCS = compute_fiber_continuity(theta_e, material_mask_core)
    % 根据相邻单元角度差计算纤维连续性指标

    [nely, nelx] = size(theta_e);
    continuous = 0;
    total = 0;
    angle_threshold = 10 * pi/180;

    if nargin < 2 || isempty(material_mask_core)
        material_mask_core = true(nely, nelx);
    elseif ~isequal(size(material_mask_core), [nely, nelx])
        error('material_mask_core尺寸错误：期望[%d,%d]，实际[%d,%d]。', ...
            nely, nelx, size(material_mask_core, 1), size(material_mask_core, 2));
    else
        material_mask_core = logical(material_mask_core);
    end

    for i = 1:nely
        for j = 1:nelx
            if j < nelx
                if ~(material_mask_core(i, j) && material_mask_core(i, j+1))
                    continue;
                end
                if ~(isfinite(theta_e(i, j)) && isfinite(theta_e(i, j+1)))
                    continue;
                end
                diff = abs(theta_e(i, j) - theta_e(i, j+1));
                diff = min(diff, pi - diff);
                if diff <= angle_threshold
                    continuous = continuous + 1;
                end
                total = total + 1;
            end
            if i < nely
                if ~(material_mask_core(i, j) && material_mask_core(i+1, j))
                    continue;
                end
                if ~(isfinite(theta_e(i, j)) && isfinite(theta_e(i+1, j)))
                    continue;
                end
                diff = abs(theta_e(i, j) - theta_e(i+1, j));
                diff = min(diff, pi - diff);
                if diff <= angle_threshold
                    continuous = continuous + 1;
                end
                total = total + 1;
            end
        end
    end

    if total > 0
        FCS = continuous / total;
    else
        FCS = 0;
    end
end

