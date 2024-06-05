%% 16个格点的Kitaev模型
clear; clc;
N=4*4;  %总自旋数目

%% 生成哈密顿量矩阵
[H_x,H_y,H_z,H_h]=createHamiltonian();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jx=1;
Jy=1;
Jz=1;
h = 0.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=Jx*H_x+Jy*H_y+Jz*H_z+h*H_h;  

tic
[V,D]=eigs(H,1,'sa');   % 稀疏矩阵对角化
toc

sigma_z = sparse([1, 0; 0, -1]);
sigma_x=sparse([0 1; 1 0]); 

sigma_z_prod = sigma_z;
for i = 1:15
    sigma_z_prod = kron(sigma_z_prod, sigma_z);
end

sigma_x_prod = sigma_x;
for i = 1:15
    sigma_x_prod = kron(sigma_x_prod, sigma_x);
end

% V' * sigma_z_prod * V


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H_x,H_y,H_z,H_h] = createHamiltonian(~)
%% 16个格点的Kitaev模型
sigma_x=sparse([0 1; 1 0]);  
sigma_y=sparse([0 -1i;1i 0]);
sigma_z=sparse([1 0; 0 -1]); 


%% x-bond
   H_x=+kron(speye(2^0), kron(kron(kron(sigma_x, speye(2^3)), sigma_x), speye(2^11)))+... %1-5 site
       +kron(speye(2^2), kron(kron(kron(sigma_x, speye(2^3)), sigma_x), speye(2^9)))+...  %3-7 site
       +kron(speye(2^5), kron(kron(kron(sigma_x, speye(2^3)), sigma_x), speye(2^6)))+...  %6-10 site
       +kron(speye(2^7), kron(kron(kron(sigma_x, speye(2^3)), sigma_x), speye(2^4)))+...  %8-12 site
       +kron(speye(2^8), kron(kron(kron(sigma_x, speye(2^3)), sigma_x), speye(2^3)))+...  %9-13 site
       +kron(speye(2^10), kron(kron(kron(sigma_x, speye(2^3)), sigma_x), speye(2^1)));    %11-15 site
       
%% y-bond
   H_y=+kron(speye(2^1), kron(kron(sigma_y, sigma_y), speye(2^13)))+... %2-3 site
       +kron(speye(2^4), kron(kron(sigma_y, sigma_y), speye(2^10)))+... %5-6 site
       +kron(speye(2^6), kron(kron(sigma_y, sigma_y), speye(2^8)))+...  %7-8 site
       +kron(speye(2^9), kron(kron(sigma_y, sigma_y), speye(2^5)))+...  %10-11 site
       +kron(speye(2^12), kron(kron(sigma_y, sigma_y), speye(2^2)))+... %13-14 site
       +kron(speye(2^14), kron(kron(sigma_y, sigma_y), speye(2^0)));    %15-16 site
       
%% z-bond
    H_z=+kron(speye(2^0), kron(kron(sigma_z, sigma_z), speye(2^14)))+... %1-2 site
        +kron(speye(2^2), kron(kron(sigma_z, sigma_z), speye(2^12)))+... %3-4 site
        +kron(speye(2^5), kron(kron(sigma_z, sigma_z), speye(2^9)))+...  %6-7 site
        +kron(speye(2^8), kron(kron(sigma_z, sigma_z), speye(2^6)))+...  %9-10 site
        +kron(speye(2^10), kron(kron(sigma_z, sigma_z), speye(2^4)))+... %11-12 site
        +kron(speye(2^13), kron(kron(sigma_z, sigma_z), speye(2^1)));    %14-15 site
        
%% magnatic field
H_h = sparse(0);
    for i = 1:16
        H_h = H_h + (kron(kron(speye(2^(i-1)), sigma_x), speye(2^(16-i))) + ...
            kron(kron(speye(2^(i-1)), sigma_y), speye(2^(16-i))) + ...
            kron(kron(speye(2^(i-1)), sigma_z), speye(2^(16-i))));
        
    end
end