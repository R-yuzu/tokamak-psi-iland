clear;clc;format longE;
cd 'C:\Users\wudou\zhz\zhz\clt\'
mx=256;
mz=256;
my=64;
aa=1;
r0=1;
nst_start=0;
nst_step=1;
nst_end=1;
phi_slide = 1;
room=0.0545069847;
Raxis = 1.84699652;  
Zaxis = -0.0321515666;  
tkdir = './data/exam8-8-8-5-15-5/tk000000';
bsdir = './data/exam8-8-8-5-15-5/eq_pgfile_rzbs.dat';
psi0dir='./data/exam8-8-8-5-15-5/eq_pgfile_2d.dat';
function [RBS, ZBS] = read_rzbs(dir)
    fid = fopen(dir, 'rt');
    if fid == -1
        error('无法打开文件: %s', dir);
    end
    data = textscan(fid, '%f %f', 'CollectOutput', true);
    fclose(fid);
    if isempty(data) || isempty(data{1})
        error('文件为空: %s', dir);
    end
    RBS = data{1}(:, 1);
    ZBS = data{1}(:, 2);
end
function [ncase,nstep,time,xtotal]=load_tkfile(dir)
% ref: https://www.cnblogs.com/jun-phy/p/12771979.html
% xtotal(rho,pressure,v_R,v_phi,v_Z,B_R,B_phi,B_Z)
    fid=fopen(dir,'rb');
    fseek(fid,4,'cof');
    ncase=fread(fid,1,'int32');
    nstep=fread(fid,1,'int32');
    time=fread(fid,1,'float64');
    fseek(fid,4,'cof');
    fseek(fid,4,'cof');
    xtotal=fread(fid,'float64');
    fclose(fid);
end
function [psi, R, Z, br_gfile, by_gfile, bz_gfile] = read_2d(dir)
    %RR, ZZ, THT, PSI, Q, NE, NI, PRES, PE, PI, VR, VT, VZ, BR, BT, BZ, JR, JT, JZ, 
    % [d/dR and d/dZ of above variables(NE-JZ)], PSI_DRRZ,PSI_DZRZ, ERR_R, ERR_T, ERR_Z, DIVB
    fid = fopen(dir, 'r');
    if fid == -1
        error('无法打开文件: %s', dir);
    end
    header = fgetl(fid);
    data = textscan(fid, '%f %f %*f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %*[^\n]', ...
        'Delimiter', ' ', 'MultipleDelimsAsOne', true);
    fclose(fid);
    Rcol = data{1};
    Zcol = data{2};
    PSIcol = data{3};
    BRcol= data{4};
    BTcol= data{5};
    BZcol= data{6};
    R_unique = unique(Rcol);
    Z_unique = unique(Zcol);
    R_unique = sort(R_unique);
    Z_unique = sort(Z_unique);
    [~, R_idx] = ismember(Rcol, R_unique);
    [~, Z_idx] = ismember(Zcol, Z_unique);
    nR = length(R_unique);
    nZ = length(Z_unique);
    psi = NaN(nZ, nR);
    br_gfile= NaN(nZ, nR);
    by_gfile= NaN(nZ, nR);
    bz_gfile= NaN(nZ, nR);
    for k = 1:length(PSIcol)
        i = Z_idx(k);
        j = R_idx(k);
        psi(i, j) = PSIcol(k);
        br_gfile(i, j)=BRcol(k);
        by_gfile(i, j)=BTcol(k);
        bz_gfile(i, j)=BZcol(k);
    end
    R = R_unique;
    Z = Z_unique;
end
function [br,bz]=psi_to_br_and_bz(mx,mz,psi,dz,dr,r_inv)
    br=zeros(mz,mx);
    bz=zeros(mz,mx);
    for i = 3:mz-2
    for j = 3:mx-2
        dpsidz = ( -psi(i+2, j) + 8*psi(i+1, j) - 8*psi(i-1, j) + psi(i-2, j) ) / (12 * dz);
        dpsidr = ( -psi(i, j+2) + 8*psi(i, j+1) - 8*psi(i, j-1) + psi(i, j-2) ) / (12 * dr);
        br(i, j) = - r_inv(j) * dpsidz;
        bz(i, j) =   r_inv(j) * dpsidr;
    end
    end
    scale=0.0545069847;%来源于归一化磁通
    br=br(3:end-2,3:end-2).*(scale);
    bz=bz(3:end-2,3:end-2).*(scale);
end
function psi1 = generate_psi_v1(r, br, bz, dr, dz, mx, mz, psi, fixed)
mu0 = 1;
omega = 1.6;     % SOR 松弛因子
scale = 1;
maxIter = 100000;
tol = 1e-8;
psi1=zeros(mz+4,mx+4);%zeros
psi0=psi_initial(psi);%initial
psi1(3:end-2,3:end-2) = psi0./(scale);   % mz+4, mx+4
for iter = 1:maxIter    
    err = 0;    
    for j = 3:mz+2
        for i = 3:mx+2  
            if fixed(j-2, i-2)
                continue;   % 跳过固定点
            end
            ri = r(i-2);   % 对应物理r坐标            
            dBr_dz = (br(j+1,i)-br(j-1,i))/(2*dz);
            dBz_dr = (bz(j,i+1)-bz(j,i-1))/(2*dr);
            jphi = (dBr_dz - dBz_dr)/mu0;            
            A = 2/dr^2 + 2/dz^2;
            rhs = (psi1(j,i+1)+psi1(j,i-1))/dr^2 ...
                + (psi1(j+1,i)+psi1(j-1,i))/dz^2 ...
                + (1/ri)*(psi1(j,i+1)-psi1(j,i-1))/(2*dr) ...
                - mu0 * ri * jphi;            
            psi_new = rhs / A;            
            psi1(j,i) = (1-omega)*psi1(j,i) + omega*psi_new;            
            err = max(err, abs(psi_new - psi1(j,i)));
        end
    end     
    if mod(iter,100)==0
        fprintf('Iter=%d, err=%e\n',iter,err);
    end    
    if err < tol
        fprintf('Converged at %d iterations\n',iter);
        break
    end
end
psi1=psi1./(0.0545069847);
psi1=psi1.*scale;
end
function [i,j] = index(r,z,rpoint,zpoint)
    [Rgrid, Zgrid] = meshgrid(r, z); 
    dist = (Rgrid - rpoint).^2 + (Zgrid - zpoint).^2;
    [~, idx] = min(dist(:));
    j = mod(idx-1, size(Rgrid,1)) + 1; %z
    i = floor((idx-1)/size(Rgrid,1)) + 1; %r
end
function boundary_mask=boundry_curve(rbs,zbs,r,z,mx,mz)
% ========================
% 将离散边界点拟合为闭合曲线，并获取曲线上的网格点索引
% ========================
% 参数化：按累积弧长归一化
tau = [0; cumsum(sqrt(diff(rbs).^2 + diff(zbs).^2))];
tau = tau / tau(end);  % 归一化到 [0,1]
% 使用周期样条插值（确保闭合曲线的光滑性）
pp_r = csape(tau, rbs, 'periodic');
pp_z = csape(tau, zbs, 'periodic');
% 生成密集插值点（数量根据网格密度调整，建议 10×网格点数）
nDense = 1000;  % 可调
t_dense = linspace(0, 1, nDense)';
r_dense = fnval(pp_r, t_dense);
z_dense = fnval(pp_z, t_dense);
[Rgrid, Zgrid] = meshgrid(r, z); 
% 获取完整网格点坐标（含 ghost cells）
grid_points = [Rgrid(:), Zgrid(:)];   % 所有网格点展平为 N×2 矩阵
% 找到每个密集点在完整网格上的最近邻（欧氏距离）
idx_nearest = dsearchn(grid_points, [r_dense, z_dense]);
% 去重，得到唯一的线性索引
idx_unique = unique(idx_nearest);
% 转换为二维下标（完整网格的行、列）
[i_full, j_full] = ind2sub(size(Rgrid), idx_unique);
% 转换为内部网格索引（去掉 ghost cells，用于 psi 数组）
i_internal = i_full ;
j_internal = j_full ;
% 剔除超出内部网格范围的点（ghost cells 边界可能被误选，但通常不会）
valid = (i_internal >= 1) & (i_internal <= mz) & ...
        (j_internal >= 1) & (j_internal <= mx);
i_boundary = i_internal(valid);
j_boundary = j_internal(valid);
% 创建边界掩码（内部网格尺寸 mz×mx）
boundary_mask = false(mz, mx);
boundary_mask(sub2ind([mz, mx], i_boundary, j_boundary)) = true;
end
function fixed=fix(boundary_mask,i0,j0)
    fixed = false(size(boundary_mask));
    %fixed = boundary_mask;
    fixed(1:end,1)=true;fixed(1:end,end)=true;
    fixed(1,1:end)=true;fixed(end,1:end)=true;
    fixed(2:end-1,2)=true;fixed(2:end-1,end-1)=true;
    fixed(2,2:end-1)=true;fixed(end-1,2:end-1)=true;
    fixed(3:end-2,3)=true;fixed(3:end-2,end-2)=true;
    fixed(3,3:end-2)=true;fixed(end-2,3:end-2)=true;
    fixed(j0,i0)=true;%z,r-mz,mx
end
function psi0=psi_initial(psi)
    psi0=zeros(size(psi));
    psi0(1:end,1)=psi(1:end,1);psi0(1:end,end)=psi(1:end,end);
    psi0(1,1:end)=psi(1,1:end);psi0(end,1:end)=psi(end,1:end);
    psi0(2:end-1,2)=psi(2:end-1,2);psi0(2:end-1,end-1)=psi(2:end-1,end-1);
    psi0(2,2:end-1)=psi(2,2:end-1);psi0(end-1,2:end-1)=psi(end-1,2:end-1);
    psi0(3:end-2,3)=psi(3:end-2,3);psi0(3:end-2,end-2)=psi(3:end-2,end-2);
    psi0(3,3:end-2)=psi(3,3:end-2);psi0(end-2,3:end-2)=psi(end-2,3:end-2);
    psi0=psi0.*(0.0545069847);
    %psi0=psi.*(0.0545069847);
end
function psi1 = generate_psi_v2(r, br, bz, dr, dz, mx, mz, psi, fixed)
omega = 1.7;     
maxIter = 100000;
tol = 1e-15;
psi1 = zeros(mz+4, mx+4);
psi0 = psi_initial(psi);
psi1(3:end-2,3:end-2) = psi0;
for iter = 1:maxIter    
    err = 0;    
    for j = 3:mz+2
        for i = 3:mx+2            
            if fixed(j-2, i-2)
                continue;
            end            
            ri = r(i-2);
            d_rBz_dr = ( r(i-1)*bz(j,i+1) - r(i-3)*bz(j,i-1) ) / (2*dr);
            d_rBr_dz = ( ri*br(j+1,i) - ri*br(j-1,i) ) / (2*dz);            
            rhs = d_rBz_dr - d_rBr_dz;
            A = 2/dr^2 + 2/dz^2;            
            psi_new = ( ...
                (psi1(j,i+1)+psi1(j,i-1))/dr^2 ...
              + (psi1(j+1,i)+psi1(j-1,i))/dz^2 ...
              - rhs ) / A;            
            psi1(j,i) = (1-omega)*psi1(j,i) + omega*psi_new;            
            err = max(err, abs(psi_new - psi1(j,i)));
        end
    end    
    if mod(iter,200)==0
        fprintf('Iter=%d, err=%e\n',iter,err);
    end    
    if err < tol
        fprintf('Converged at %d iterations\n',iter);
        break
    end
end
psi1=psi1./(0.0545069847);
end
function theta=generate_theta(r,z,Raxis,Zaxis,mx,mz)
    theta=zeros(mz,mx);
    zl_divide_rl=zeros(mz,mx);
    for i=1:mz
        for j=1:mx
            rlength=r(j)-Raxis;
            zlength=z(i)-Zaxis;
                if (-0.000000000000003<=rlength)&&(rlength<=0.000000000000003)
                    if zlength>0    
                        theta(i,j)=0.5*pi;
                    end
                    if zlength<0
                        theta(i,j)=1.5*pi;
                    end
                    continue;
                end
            zl_divide_rl(i,j)=zlength/rlength;
            theta(i,j)=atan(zl_divide_rl(i,j));
            if (rlength>0)&&(zlength<0)
                theta(i,j)=theta(i,j)+2*pi;
            end
            if rlength<0
                theta(i,j)=theta(i,j)+pi;
            end
        end
    end
end
function psi_theta=psitheta(psi10,theta,mx,mz)
    psi_theta=zeros(mz,mx,2);
    for i=1:mz
        for j=1:mx
            psi_theta(i,j,1)=psi10(i,j);
            psi_theta(i,j,2)=theta(i,j);
        end
    end
end
function brzmod_psi_theta=mod_psi_theta(mx,mz,br,bz,psi_theta)
    br_in=br(3:mz+2,3:mx+2);
    bz_in=bz(3:mz+2,3:mx+2);
    br_square=br_in.^2;
    bz_square=bz_in.^2;
    brzmod_square=bz_square+br_square;
    brzmod=sqrt(brzmod_square);
    brzmod_psi_theta=zeros(mz,mx,3);
    for i=1:mz
        for j=1:mx
            brzmod_psi_theta(i,j,1)=brzmod(i,j);
            brzmod_psi_theta(i,j,2)=psi_theta(i,j,1);
            brzmod_psi_theta(i,j,3)=psi_theta(i,j,2);
        end
    end
end
function brzmod_sqrtpsi_theta=mod_sqrtpsi_theta(brzmod_psi_theta,mz,mx)
    brzmod_sqrtpsi_theta=zeros(mz,mx,3);
    brzmod_sqrtpsi_theta(:,:,1)=brzmod_psi_theta(:,:,1);
    brzmod_sqrtpsi_theta(:,:,3)=brzmod_psi_theta(:,:,3);
    brzmod_sqrtpsi_theta(:,:,2)=sqrt(brzmod_psi_theta(:,:,2));
end
function figure1(r,z,psi,Raxis,Zaxis,rbs,zbs)
    figure;
    contour(r, z, psi, 2000); %mz,mx
    xlabel('R (m)'); ylabel('Z (m)');
    title('Magnetic Surfaces (\psi = const)');
    axis equal;
    colorbar;
    hold on;
    scatter(Raxis,Zaxis,5,'d');
    hold on;
    scatter(rbs,zbs,5,'o');
    hold off;
end
function figure2(r,z,br,mx,mz)
    figure;
    contour(r, z, br(3:mz+2,3:mx+2), 1000);
    xlabel('R (m)'); ylabel('Z (m)');
    title('br');
    axis equal;
    colorbar;
end
function figure3(r,z,bz,mx,mz)
    figure;
    contour(r, z, bz(3:mz+2,3:mx+2), 1000);
    xlabel('R (m)'); ylabel('Z (m)');
    title('bz');
    axis equal;
    colorbar;
end
function figure4(mx,mz,br_psi,bz_psi,r,z)
    figure;
    contour(r(3:mx-2), z(3:mz-2), br_psi, 1000);
    xlabel('R (m)'); ylabel('Z (m)');
    title('br-psi');
    axis equal;
    colorbar;
    figure;
    contour(r(3:mx-2), z(3:mz-2), bz_psi, 1000);
    xlabel('R (m)'); ylabel('Z (m)');
    title('bz-psi');
    axis equal;
    colorbar;
end
function figure5(r,z,br_gfile,bz_gfile)
    figure;
    contour(r, z, br_gfile, 1000);
    xlabel('R (m)'); ylabel('Z (m)');
    title('br-gfile');
    axis equal;
    colorbar;
    figure;
    contour(r, z, bz_gfile, 1000);
    xlabel('R (m)'); ylabel('Z (m)');
    title('bz-gfile');
    axis equal;
    colorbar;
end
function figure6(br,bz,br_psi,bz_psi,br_gfile,bz_gfile,mx,mz,r,z)
    figure;
    subplot(2,3,1);
    contour(r(3:mx-2), z(3:mz-2), br(5:mz,5:mx)-br_psi, 1000);
    xlabel('R (m)'); ylabel('Z (m)');
    title('br between br-psi');
    axis equal;
    colorbar;
    subplot(2,3,2);
    contour(r(3:mx-2), z(3:mz-2), br(5:mz,5:mx)-br_gfile(3:mz-2,3:mx-2), 1000);
    xlabel('R (m)'); ylabel('Z (m)');
    title('br between br-gfile');
    axis equal;
    colorbar;
    subplot(2,3,3);
    contour(r(3:mx-2), z(3:mz-2), br_gfile(3:mz-2,3:mx-2)-br_psi, 1000);
    xlabel('R (m)'); ylabel('Z (m)');
    title('br-gfile between br-psi');
    axis equal;
    colorbar;
    subplot(2,3,4);
    contour(r(3:mx-2), z(3:mz-2), bz(5:mz,5:mx)-bz_psi, 1000);
    xlabel('R (m)'); ylabel('Z (m)');
    title('bz between bz-psi');
    axis equal;
    colorbar;
    subplot(2,3,5);
    contour(r(3:mx-2), z(3:mz-2), bz(5:mz,5:mx)-bz_gfile(3:mz-2,3:mx-2), 1000);
    xlabel('R (m)'); ylabel('Z (m)');
    title('bz between bz-gfile');
    axis equal;
    colorbar;
    subplot(2,3,6);
    contour(r(3:mx-2), z(3:mz-2), bz_gfile(3:mz-2,3:mx-2)-bz_psi, 1000);
    xlabel('R (m)'); ylabel('Z (m)');
    title('bz-gfile between bz-psi');
    axis equal;
    colorbar;
end
function figure7(r,z,varpsi,Raxis,Zaxis,rbs,zbs)
    figure;
    contour(r, z, varpsi(3:end-2,3:end-2), 2000); %mz,mx
    xlabel('R (m)'); ylabel('Z (m)');
    title('Magnetic Surfaces (\psi = const)');
    axis equal;
    colorbar;
    hold on;
    scatter(Raxis,Zaxis,5,'d');
    hold on;
    scatter(rbs,zbs,5,'o');
    hold off;
end
function figure8(r,z,psi,i0,j0,rbs,zbs,boundary_mask)
    figure;
    contour(r, z, psi, 500); %mz,mx
    xlabel('R (m)'); ylabel('Z (m)');
    title('Magnetic Surfaces (\psi = const)');
    axis equal;
    colorbar;
    hold on;
    scatter(rbs,zbs,18,'o');
    hold on;
    scatter(r(i0),z(j0),5,'d');
    hold on;
    [row, col] = find(boundary_mask);
    z_coords = z(row);   % 行下标对应z坐标
    r_coords = r(col);   % 列下标对应x坐标
    scatter(r_coords, z_coords,5, '*');
end
function figure9(fixed,r,z,psi,rbs,zbs)
    figure;
    contour(r, z, psi, 500); %mz,mx
    xlabel('R (m)'); ylabel('Z (m)');
    title('Magnetic Surfaces (\psi = const)');
    axis equal;
    colorbar;
    hold on;
    scatter(rbs,zbs,18,'o');
    hold on;
    [row, col] = find(fixed);
    z_coords = z(row);   % 行下标对应z坐标
    r_coords = r(col);   % 列下标对应x坐标
    scatter(r_coords, z_coords,5, '*');
end
function figure10(r,z,varpsi_gfile,Raxis,Zaxis,rbs,zbs)
    figure;
    contour(r, z, varpsi_gfile(3:end-2,3:end-2), 2000); %mz,mx
    xlabel('R (m)'); ylabel('Z (m)');
    title('Magnetic Surfaces (\psi = const)');
    axis equal;
    colorbar;
    hold on;
    scatter(Raxis,Zaxis,5,'d');
    hold on;
    scatter(rbs,zbs,5,'o');
    hold off;
end
function figure11(r,z,rbs,zbs,Raxis,Zaxis,deltapsi)
    figure;
    contour(r, z, deltapsi, 2000); %mz,mx
    xlabel('R (m)'); ylabel('Z (m)');
    title('Magnetic Surfaces (\psi = const)');
    axis equal;
    colorbar;
    hold on;
    scatter(Raxis,Zaxis,5,'d');
    hold on;
    scatter(rbs,zbs,5,'o');
    hold off;
end
function figure12(r,z,rbs,zbs,Raxis,Zaxis,psi_theta)
    figure;
    subplot(1,2,1);
    contour(r, z, psi_theta(:,:,1), 500); %mz,mx
    xlabel('R (m)'); ylabel('Z (m)');
    title('Magnetic Surfaces (\psi = const)');
    axis equal;
    colorbar;
    hold on;
    scatter(Raxis,Zaxis,5,'d');
    hold on;
    scatter(rbs,zbs,5,'o');
    hold off;
    subplot(1,2,2);
    contour(r,z,psi_theta(:,:,2),500);
    xlabel('R (m)'); ylabel('Z (m)');
    title('Theta (\theta)');
    axis equal;
    colorbar;
end
function figure13(data)
data1 = data(:,:,2);
data2 = data(:,:,3);
data3 = data(:,:,1);
mask = data1 <= 1.01;
data1   = data1(mask);
data2   = data2(mask);
data3   = data3(mask);
scatter(data1(:), data2(:), 15, data3(:), '.');
axis([-0.01 1.01 -0.01 2*pi+0.1])
colorbar;               % 显示颜色条
xlabel('\psi'); ylabel('\theta');
title('B_{rz,mod}');
end
function figure14(data)
data1 = data(:,:,2);
data2 = data(:,:,3);
data3 = data(:,:,1);
mask = data1 <= 1.01;
data1   = data1(mask);
data2   = data2(mask);
data3   = data3(mask);
scatter(data1(:), data2(:), 200, data3(:), '.');
axis([0.20 1.01 -0.01 2*pi+0.1])
colorbar;               % 显示颜色条
xlabel('\psi^{0.5}'); ylabel('\theta');
title('B_{rz,mod}');
end
[ncase,nstep,time,xtotal]=load_tkfile(tkdir);%tkfile-br,bz
[rbs,zbs] =read_rzbs(bsdir);%boundary
[psi,r,z,br_gfile,by_gfile,bz_gfile] =read_2d(psi0dir);%gfile
dr=(r(end)-r(1))/255.0;
dz=(z(end)-z(1))/255.0;
xtk=reshape(xtotal,mx+4,mz+4,my+4,8);%x-parametres
xslide=xtk(:,:,34,:);
xphi=squeeze(xslide);%slide
rho=xphi(:,:,1)';pre=xphi(:,:,2)';vr=xphi(:,:,3)';vy=xphi(:,:,4)';vz=xphi(:,:,5)';br=xphi(:,:,6)';by=xphi(:,:,7)';bz=xphi(:,:,8)';
r_inv=1./r';%r-inverse
[br_psi,bz_psi]=psi_to_br_and_bz(mx,mz,psi,dz,dr,r_inv);%psi
boundary_mask=boundry_curve(rbs,zbs,r,z,mx,mz);%boundary
[i0,j0]=index(r,z,Raxis,Zaxis);%mag-axis
fixed=fix(boundary_mask,i0,j0);%fixed-point
%varpsi=generate_psi_v1(r, br, bz, dr, dz, mx, mz, psi, fixed);%psi1
br_gfile_full=zeros(mz+4,mx+4);
bz_gfile_full=zeros(mz+4,mx+4);
br_gfile_full(3:end-2,3:end-2)=br_gfile;
bz_gfile_full(3:end-2,3:end-2)=bz_gfile;
%varpsi_gfile=generate_psi_v1(r, br_gfile_full, bz_gfile_full, dr, dz, mx, mz, psi, fixed);%psi1-gfile
varpsi_v2=generate_psi_v2(r, br, bz, dr, dz, mx, mz, psi, fixed);%psi1
%varpsi_gfile_v2=generate_psi_v2(r, br_gfile_full, bz_gfile_full, dr, dz, mx, mz, psi, fixed);%psi1-gfile
deltapsi=psi-varpsi_v2(3:end-2,3:end-2);%error-psi
%deltapsi_gfile=psi-varpsi_gfile_v2(3:end-2,3:end-2);%error-psi-gfile
psi10=varpsi_v2(3:end-2,3:end-2);
theta=generate_theta(r,z,Raxis,Zaxis,mx,mz);
psi_theta=psitheta(psi10,theta,mx,mz);
brzmod_psi_theta=mod_psi_theta(mx,mz,br,bz,psi_theta);%mod-psi-theta
brzmod_sqrtpsi_theta=mod_sqrtpsi_theta(brzmod_psi_theta,mz,mx);%mod-sqrt-theta
% ========================%
% 绘图 %
% ========================%
%figure1(r,z,psi,Raxis,Zaxis,rbs,zbs);
%figure2(r,z,br,mx,mz);
%figure3(r,z,bz,mx,mz);
%figure4(mx,mz,br_psi,bz_psi,r,z);
%figure5(r,z,br_gfile,bz_gfile);
%figure6(br,bz,br_psi,bz_psi,br_gfile,bz_gfile,mx,mz,r,z);
%figure7(r,z,varpsi_v2,Raxis,Zaxis,rbs,zbs);
%figure8(r,z,psi,i0,j0,rbs,zbs,boundary_mask);
%figure9(fixed,r,z,psi,rbs,zbs);
%figure10(r,z,varpsi_gfile,Raxis,Zaxis,rbs,zbs);
%figure11(r,z,rbs,zbs,Raxis,Zaxis,deltapsi);
%figure12(r,z,rbs,zbs,Raxis,Zaxis,psi_theta);
%figure13(brzmod_psi_theta);
disp('complete');figure14(brzmod_sqrtpsi_theta);
%B_r(r,z) 
%B_z(r,z)
%\psi(r,z)
%B_p=\frac{1}{r}\nabla\psi\times\hat(e)_{\phi}
%\nabla\psi=(\frac{\partial\psi}{\partial r},0,\frac{\partial\psi}{\partial z})
%\nabla\psi\times\hat{e}_{\phi}
%B_r=-\frac{1}{r}\frac{\partial\psi}{\partial z}
%B_z=\frac{1}{r}\frac{\partial\psi}{\partial r}
%\partial_r\psi=rB_z
%\partial_z\psi=-rB_r
%F=(rB_z,-rB_r)
%\nabla\psi=\vec{F}
%J[\psi]=\int_{\Omega}[(\partial_r\psi-rB_z)^2+(\partial_z\psi+rB_r)^2]drdz
%\psi=\psi+\varepsilon\eta
%\delta J=0
%\delta J=2\int(\partial_r\psi-rB_z)\partial_r\eta+(\partial_z\psi+rB_r)\partial_z\eta drdz
%\partial_r(\partial_r\psi-rB_z)+\partial_z(\partial_z\psi+rB_r)=0
%\partial^2_r\psi+\partial^2_z\psi=\partial_r(rB_z)-\partial_z(rB_r)
%\Delta\psi=\partial_r(rB_z)-\partial_z(rB_r)