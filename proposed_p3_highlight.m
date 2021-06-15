% In this version, both S and Si_s are constrained with inequality constraint

clear
clc
warning off;

para = xlsread('chooseBest_p3.xlsx');

DataName = cell(13 , 1);
DataName{1} = 'bbcsport';
DataName{2} = 'bbcsport2view';
DataName{3} = 'proteinFold';
DataName{4} = 'caltech101_mit';
DataName{5} = 'CCV';
DataName{6} = 'flower17_DL_fea';
DataName{7} = 'mfeat';
DataName{8} = 'plant';
DataName{9} = 'psortPos';
DataName{10} = 'UCI_DIGIT';
DataName{11} = 'flower102';
DataName{12} = 'nonpl';
DataName{13} = 'flower17';

% DataName{1} = 'caltech101_mit';
% DataName{2} = 'CCV';
% DataName{3} = 'flower17_DL_fea';
% DataName{4} = 'flower102';
% DataName{5} = 'plant';
% DataName{6} = 'flower17';

for ICount = 1 : 13
    
    % Set path
    path = '../0_KernelProcessingandClustering';
    %pathdata = '/media/sihang/8ffb0279-b7eb-473a-b956-305a3b433e0c/Project/Kernel/datasets';
    pathdata='/media/maggie/Data/papers/sh/datasets';
    addpath(genpath(path));
    dataName = DataName{ICount};
    load([pathdata,'/',dataName,'_Kmatrix'],'KH','Y');
    
    
    % Preprocessing
    CluNum = length(unique(Y));  % num of classes
    ker_num = size(KH,3); %num of views
    sample_num = size(KH,1);
    KH = kcenter(KH);
    KH = knorm(KH);
    f_num = CluNum * 4;
    H = zeros(f_num, sample_num, ker_num);
    
    opt.disp = 0;
    for p=1:ker_num % m - kernels
        KH(:,:,p) = (KH(:,:,p)+KH(:,:,p)')/2;
        [Hp, ~] = eigs(KH(:,:,p), f_num, 'la', opt);
        H(:,:,p) = Hp';
    end
    H_total = mean(H,3);
    
    
    %% H'H
    %%%%%%%%%%%%%%%%%%%  Proposed Method  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H2 = zeros(size(KH));
    for i = 1 : ker_num
        H2(:,:,i) = H(:,:,i)' * H(:,:,i);
    end
    
    %% Try different hyper-parameters
%     for alpha = 2.^[-15 : 2: 15]
    idx_acc_nmi_prt = 0;
    for paraCount = 1:2:6
        idx_acc_nmi_prt = idx_acc_nmi_prt + 1;
        ai = para(ICount, paraCount + 1);
    parameters_alpha = 2.^[-15 : 2: 15];
        for alpha = parameters_alpha(ai)
%         tic;
        %% Initialization
        % Initialization of P
        P = eye(sample_num);
        
        % Initialization of S
        S = zeros(sample_num, sample_num);
        S_record = zeros(sample_num, sample_num, ker_num);
        
        % Initialization of Si
        for i = 1:ker_num
            A = H2(:,:,i) * P;
            [Si] = qp(S, A, alpha);
            S_record(:,:,i) = Si;
        end
        
        
        %% Optimization
        obj_record = [];
        it_count = 1;
        flag = 1;
        while it_count < 50 && flag == 1
            
            fprintf('DataName : %s , it_count  %f \n', dataName,  it_count);
            
            %% Update Si, namely S_record
            for i = 1:ker_num
                A = H2(:,:,i) * P;
                [Si] = qp(S, A, alpha);
                S_record(:,:,i) = Si;
            end
            %                 Si_Value = cal_obj(H2, P, S, S_record, alpha)
            
            %% Update S
            S1 = mean(S_record, 3);
            S1(S1 < 0) = 0;
            S1(S1 > 1) = 1;
            %                 S_Value = cal_obj(H2, P, S, S_record, alpha)
            S_diff = norm(S-S1, 'fro');
            S = S1;
            
            % objective value calculation
            it_count = it_count+1;
%             cur_obj = cal_obj(H2, P, S, S_record, alpha);
%             obj_record = [obj_record, cur_obj];
            
%             if it_count > 2 && abs((obj_record(it_count-1)-obj_record(it_count-2))/ obj_record(it_count-2)) < 10^(-4)
%                 flag = 0;
%             end
            if it_count > 2 && S_diff < 10^(-3)
                flag = 0;
            end
        end
        
%         [res, idx]= myNMIACC(S,Y,CluNum);
%         m_i = 0;
        
%         res_time = toc
            
            res_tsne = tsne(S,Y,2);
            print(gcf,'-dpng',['./highlightPNG_3/',dataName,'highlight',num2str(idx_acc_nmi_prt),'.png'])

        
        end
    end
end


%% The optimization process of Si
% min \alpha * ||S||_F^2  - 2 Tr((A+alpha * B)*S^T), s.t. 0<= S <=1
function [Si] = qp(A, B, alpha)
Si = (2 * alpha * A + B) / (2*alpha);
Si(Si > 1) = 1;
Si(Si < 0) = 0;
end


%% Calculate the objective function
function obj = cal_obj(H2, P, S, Si_s, alpha)
obj = 0;
ker_num = size(H2, 3);

%part 1
for i = 1 : ker_num
    obj = obj + trace(Si_s(:,:,i) * P' * H2(:,:,i));
end

obj = -obj;
% part 2
for i = 1 : ker_num
    obj = obj + alpha * norm(S-Si_s(:,:,i), 'fro')^2;
end

end
