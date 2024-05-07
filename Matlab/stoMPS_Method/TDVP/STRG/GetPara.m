function [ Para ] = GetPara( beta )

Para.beta = beta;
fprintf('beta = %.3f \n ', Para.beta);
% // ==================== Model ====================

Para.ModelName = 'XY_chain';
Para.d = 2; % local spin; 2 -> spin-1/2, 3 -> spin-1

% Geometry
Para.L = 10;
Para.BCX = 'OBC';

% Interaction
Para.Model.J = 1;


% // ==================== Field ====================
Para.Field.h = [0,0,0];


% // ================= Statistics ==================
Para.Step_max = 5000; % Number of samples


% // ===================== RG ======================
Para.MCrit = 40; % Bond dimension
Para.It_max = 10000; % Highest order of expansion
Para.VariProd_step_max = 10000; 
Para.VariSum_step_max = 10000;

addpath(['./Model/', Para.ModelName])

Para.SamRslt_path = ['SamRslt/', Para.ModelName, '_beta-', num2str(Para.beta), '_D-', num2str(Para.MCrit), '_L-', num2str(Para.L)];
    
if ~exist(Para.SamRslt_path, 'dir')
    mkdir(Para.SamRslt_path);
end
% // ============ Entanglement entropy==============

Para.Save_EE = 1;

if ~exist(['anaRslt/EE_Rslt/', Para.ModelName], 'dir')
    mkdir(['anaRslt/EE_Rslt/', Para.ModelName]);
end

Para.EE_name = [Para.ModelName, '/beta-',num2str(Para.beta),'_D-',num2str(Para.MCrit),'_L-',num2str(Para.L),'_i-'];
end

