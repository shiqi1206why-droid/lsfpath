%% TOPOLOGY OPTIMIZATION USING THE LEVEL-SET METHOD, VIVIEN J. CHALLIS 2009
function [struc] = shuipingji2(nelx,nely,volReq,stepLength,numReinit,topWeight)
%Initialization
nelx = 80;
nely = 50;
volReq = 0.5;
stepLength = 3;
numReinit = 2;
topWeight = 4;
struc = ones(nely,nelx);
% 初始化水平集函数，使lsf满足水平集函数的定义
[lsf] = reinit(struc);
shapeSens = zeros(nely,nelx); 
topSens = zeros(nely,nelx);
[KE,KTr,lambda,mu] = materialInfo();
% Main loop:
 for iterNum = 1:200
     % FE-analysis, calculate sensitivities
     [U] = FE(struc,KE);
     % 计算形状灵敏度和拓扑灵敏度
     for ely = 1:nely
         for elx = 1:nelx
              n1 = (nely+1)*(elx-1)+ely;
              n2 = (nely+1)* elx +ely;
              Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
               shapeSens(ely,elx) = -max(struc(ely,elx),0.0001)*Ue'*KE*Ue;
               topSens(ely,elx) = struc(ely,elx)*pi/2*(lambda+2*mu)/mu/(lambda+mu)*(4*mu*Ue'*KE*Ue+(lambda-mu)*Ue'*KTr*Ue);
         end
     end
     % Store data, print & plot information
     % 将在main loop中计算的灵敏度储存在objective中
     objective(iterNum) = -sum(shapeSens(:)); 
     % 计算总体体积分数
     volCurr = sum(struc(:))/(nelx*nely);
     % 绘制结构
     disp([' It.: ' num2str(iterNum) 'Compl.: ' sprintf('%10.4f',objective(iterNum)) 'Vol.:' sprintf('%6.3f ',volCurr)])
     colormap(gray); 
     imagesc(-struc,[-1,0]); axis equal; axis tight; axis off; drawnow;
     % 检查是否收敛，必须同时满足三个条件
     if iterNum > 5 && ( abs(volCurr-volReq) < 0.005 ) && all( abs(objective(end)-objective(end-5:end-1) ) < 0.01*abs(objective(end)) )
      return;
     end
     % 更新拉格朗日函数
     if iterNum == 1
        la = -0.01; La = 1000; alpha = 0.9;
     else
         la = la - 1/La * (volCurr - volReq);
         La = alpha * La;
     end
     % 添加含拉格朗日的灵敏度，根据公式6,8,9,13,14，两项均少个负号是由于作者定义
     % la为负数，如果把la改成正数，符号就能完全与公式相对应
     shapeSens = shapeSens - la + 1/La*(volCurr - volReq);
     topSens = topSens + pi*(la - 1/La*(volCurr - volReq));
     % 调用update函数更新结构和lsf函数矩阵
     [struc,lsf] = updateStep(lsf,shapeSens,topSens,stepLength,topWeight);
     % 调用reinit函数重新初始化lsf函数使其不偏离水平集函数定义
      if ~mod(iterNum,numReinit)
          [lsf] = reinit(struc);
      end
      fprintf('主循环已结束，准备保存文件...\n');
      %% ========= 在主循环结束后添加 =========
      % 显示最终结果
      figure(2);
      imagesc(1-struc);
      colormap(gray);
      axis equal tight;
      title(sprintf('拓扑优化最终结果 - 第%d次迭代', iterNum));

      % 保存结果到文件
      save('topo_result.mat', 'struc', 'lsf', 'nelx', 'nely');
      fprintf('\n=================================\n');
      fprintf('拓扑优化完成！\n');
      fprintf('结果已保存到 topo_result.mat\n');
      fprintf('结构中材料体积分数：%.2f%%\n', sum(struc(:))/(nelx*nely)*100);
      fprintf('=================================\n');
 end

 %%---- 初始化水平集函数  ----
 function [lsf] = reinit(struc)
 % 创建一个新数组strucfull，使其比原数组struc周围多一圈0
 strucFull = zeros(size(struc)+2); 
 strucFull(2:end-1,2:end-1) = struc;
 % Use "bwdist" (Image Processing Toolbox)
 lsf = (~strucFull).*(bwdist(strucFull)-0.5) - strucFull.*(bwdist(strucFull-1)-0.5);
 %%-设计更新函数
 function [struc,lsf] = updateStep(lsf,shapeSens,topSens,stepLength,topWeight)
  % Smooth the sensitivities 
  [shapeSens] = conv2(padarray(shapeSens,[1,1],'replicate'),1/6*[ 0 1 0 ; 1 2 1 ;0 1 0],'valid');
 [topSens] = conv2(padarray(topSens,[1,1],'replicate'),1/6*[ 0 1 0 ; 1 2 1 ; 0 1 0 ] ,'valid');
 % 将试件底部两端和中间的灵敏度置零，使其不被优化
 % shapeSens(end,[1,round(end/2):round(end/2+1),end]) = 0;
 %  topSens(end,[1,round(end/2):round(end/2+1),end]) = 0;
  % Design update via evolution
  [struc,lsf] = evolve(-shapeSens,topSens.*(lsf(2:end-1,2:end-1)<0),lsf,stepLength,topWeight);
  %%---- 水平集函数演化----
   function [struc,lsf] = evolve(v,g,lsf,stepLength,w)
  % Extend sensitivites using a zero border
  vFull = zeros(size(v)+2); 
  vFull(2:end-1,2:end-1) = v;
  gFull = zeros(size(g)+2);
  gFull(2:end -1,2:end-1) = g;
  % 公式5的具体体现计算时间步长
  dt = 0.1/max(abs(v(:)));
 % Evolve for total time stepLength * CFL value:
 for i = 1:(10*stepLength)
   % 使用迎风有限差分来计算公式3中的模，先计算xy的正负向差分，p左或上，m右或下
   dpx = circshift(lsf,[0,-1])-lsf;
   dmx = lsf - circshift(lsf,[0,1]);
   dpy = circshift(lsf,[-1,0]) - lsf;
   dmy = lsf - circshift(lsf,[1,0]);
   % 使用迎风格式更新水平集函数，公式3的具体体现。v>0即max(vfull,0)向外扩张，边界为0，
   % 边界外大于0，边界内小于0，取向右差分，v<0差分方向相反
   % 向左向右是取决v的具体方向
   lsf = lsf - dt * min(vFull,0).* sqrt( min(dmx,0).^2+max(dpx,0).^2+min(dmy,0).^2+max(dpy,0).^2 ) - dt * max(vFull,0) .*  sqrt( max(dmx,0).^2+min(dpx,0).^2+max(dmy,0).^2+min(dpy,0).^2 )- w*dt*gFull;
 end
  % New structure obtained from lsf
  strucFull = (lsf<0);
  struc = strucFull(2:end-1,2:end-1);
  %%有限元分析
  function [U] = FE(struc,KE)
 [nely,nelx] = size(struc);
  K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
  F = sparse(2*(nely+1)*(nelx+1),1); U = zeros(2*(nely+1)*(nelx+1),1);
  for elx = 1:nelx
      for ely = 1:nely
          n1 = (nely+1)*(elx-1)+ely;
          n2 = (nely+1)* elx +ely;
          edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
          K(edof,edof) = K(edof,edof) + max(struc(ely,elx),0.0001)*KE;
      end
  end
  %定义负载和支撑
  F(2*(nely+1)*nelx+nely+2,1) = -1;
  fixeddofs=1:2*(nely+1);
  % Solving
  alldofs = 1:2*(nely+1)*(nelx+1);%所有节点
  freedofs = setdiff(alldofs,fixeddofs);%将所有节点与固定节点进行对比以得出自由节点
  U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
  %%---- 定义材料属性 ----
  function [KE,KTr,lambda,mu] = materialInfo()
  % Set material parameters, find Lame values
  E = 1.; 
  nu = 0.3;
  lambda = E*nu/((1+nu)*(1-nu)); mu = E/(2*(1+nu));
  % 使用stiffness函数生成单元刚度矩阵
  k =[ 1/2-nu/6 1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 -1/4+nu/12 -1/8-nu/8 nu/6 1/8-3*nu/8];
  KE = E/(1-nu^2)*stiffnessMatrix(k);   
  % Find "trace" matrix "KTr"
  k=[1/3 1/4 -1/3 1/4 -1/6 -1/4 1/6 -1/4];
  KTr = E/(1-nu)*stiffnessMatrix(k);
  %%---- ELEMENT STIFFNESS MATRIX ----
  function [K] = stiffnessMatrix(k)
  % Forms stiffness matrix from first row   
   K=[k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
    k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
    k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
    k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
    k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
    k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
    k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
    k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];   
% Gives Euclidean distance to the closestnon-zero entry in "image"
dist=zeros(size(image));
[yinds,xinds]=find(image~=0);
for i=1:size(image,2)
    for j=1:size(image,1)
         dist(j,i)=sqrt(min((yinds-j).^2+(xinds-i).^2));
    end
end

