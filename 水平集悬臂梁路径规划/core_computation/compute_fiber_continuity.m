function FCS = compute_fiber_continuity(theta_e)
    % 根据相邻单元角度差计算纤维连续性指标

    [nely, nelx] = size(theta_e);
    continuous = 0;
    total = 0;
    angle_threshold = 10 * pi/180;

    for i = 1:nely
        for j = 1:nelx
            if j < nelx
                diff = abs(theta_e(i, j) - theta_e(i, j+1));
                diff = min(diff, pi - diff);
                if diff <= angle_threshold
                    continuous = continuous + 1;
                end
                total = total + 1;
            end
            if i < nely
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

