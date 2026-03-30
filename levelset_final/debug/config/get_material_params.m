function mat = get_material_params(material_type)
    % 材料参数库
    % 
    % 输入：
    %   material_type - 材料类型名称
    % 
    % 输出：
    %   mat - 材料参数结构体
    %
    % 支持的材料：
    %   'carbon_fiber' - 碳纤维 T300/5208（默认）
    %   'glass_fiber'  - 玻璃纤维 E-Glass/Epoxy
    %   'custom'       - 自定义材料（返回模板）
    
    if nargin < 1
        material_type = 'carbon_fiber';
    end
    
    switch lower(material_type)
        case 'carbon_fiber'
            mat.name = 'Carbon Fiber T300/5208';
            mat.E_L = 137.9e9;      % 纵向杨氏模量 (Pa)
            mat.E_T = 10.34e9;      % 横向杨氏模量 (Pa)
            mat.nu_LT = 0.29;       % 纵向泊松比
            mat.G_LT = 6.89e9;      % 面内剪切模量 (Pa)
            mat.G_LW = 6.89e9;      % 面外剪切模量 (Pa)
            mat.G_TW = 3.7e9;       % 横向剪切模量 (Pa)
            mat.thickness = 0.001;  % 板厚 (m)
            
        case 'glass_fiber'
            mat.name = 'Glass Fiber E-Glass/Epoxy';
            mat.E_L = 45e9;         % 纵向杨氏模量 (Pa)
            mat.E_T = 12e9;         % 横向杨氏模量 (Pa)
            mat.nu_LT = 0.28;       % 纵向泊松比
            mat.G_LT = 5.5e9;       % 面内剪切模量 (Pa)
            mat.G_LW = 5.5e9;       % 面外剪切模量 (Pa)
            mat.G_TW = 3.5e9;       % 横向剪切模量 (Pa)
            mat.thickness = 0.001;  % 板厚 (m)
            
        case 'custom'
            % 返回模板供用户填写
            mat.name = 'Custom Material';
            mat.E_L = 0;
            mat.E_T = 0;
            mat.nu_LT = 0;
            mat.G_LT = 0;
            mat.G_LW = 0;
            mat.G_TW = 0;
            mat.thickness = 0;
            warning('返回自定义材料模板，请填写所有参数');
            
        otherwise
            error('未知材料类型: %s\n支持的类型: carbon_fiber, glass_fiber, custom', material_type);
    end
    
    % 计算派生参数
    mat.nu_TL = mat.nu_LT * mat.E_T / mat.E_L;
end

