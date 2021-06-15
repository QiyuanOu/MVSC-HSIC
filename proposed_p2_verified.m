% In this version, both S and Si_s are constrained with inequality constraint

clear
clc
warning off;


DataName = cell(13 , 1);
% DataName{1} = 'bbcsport';
% DataName{2} = 'bbcsport2view';
% DataName{3} = 'proteinFold';
% DataName{4} = 'caltech101_mit';
% DataName{5} = 'CCV';
% DataName{6} = 'flower17_DL_fea';
% DataName{7} = 'mfeat';
% DataName{8} = 'plant';
% DataName{9} = 'psortPos';
% DataName{10} = 'UCI_DIGIT';
% DataName{11} = 'flower102';
% DataName{12} = 'nonpl';
% DataName{13} = 'flower17';

DataName{1} = 'caltech101_mit';
DataName{2} = 'CCV';
DataName{3} = 'flower17_DL_fea';
DataName{4} = 'flower102';
DataName{5} = 'plant';
DataName{6} = 'flower17';


for ICount =  1 : 6
    
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
    sample_num = size(KH,1); %num of data points
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
    
    
    %%
    %%%%%%%%%%%%%%%%%%%  Proposed Method  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H2 = zeros(size(KH));
    for i = 1 : ker_num
        H2(:,:,i) = H(:,:,i)' * H(:,:,i);
    end
%     H2 = single(H2);          % gpu
%     H2 = gpuArray(H2);  
    
    p_num1 = max(50, CluNum*2);
    p_num2 = max(100, CluNum*4);
    parameters = [CluNum, p_num1, p_num2];
    
    %% Try different hyper-parameters
    for m_i = 1 : 3
        for alpha = 2.^[-15 : 2: 15]
            tic;
            % parameters
            % m -- number of anchor points
            m = parameters(m_i);
            
            %% Initialization
            % Initialization of P, strategy 1, for large scale dataset
            P = randn(sample_num, m, ker_num);
            for i = 1:ker_num
                P(:,:,i) = orth(P(:,:,i));
            end
            % Initialization of P, strategy 2
            %             [U, D, V] = svd((mean(H2, 3) + eye(size(H2, 1))));
            %             D = D(:, 1:m);
            %             P = U * (D~=0);
            
            % Initialization of S
            S = zeros(sample_num, m);
%             S = gpuArray(S);          % gpu
            S_record = zeros(sample_num, m, ker_num);
%             S_record = gpuArray(S_record);            % gpu
            
            % Initialization of Si
            for i = 1:ker_num
                A = H2(:,:,i) * P(:,:,i);
                [Si] = qp(S, A, alpha);
                S_record(:,:,i) = Si;
            end
            
            %% Optimization
            obj_record = [];
            it_count = 1;
            flag = 1;
            P1 = zeros(size(P));
%             P1 = gpuArray(P1);            % gpu
            while it_count < 50 && flag == 1
                
                fprintf('DataName : %s , it_count  %f \n', dataName,  it_count);
                
                %% Update P
                B = zeros(sample_num, m);
                p_diff = 0;
                for i = 1:ker_num
                    B = H2(:,:,i) * S_record(:,:,i);
                    [U, D, V] = svd(B, 'econ');
                    C = zeros(m);
%                     C = gpuArray(C);          % gpu
                    for j = 1 : m
                        C(j,j) = sign(D(j,j));
                    end
                    P1(:,:,i) = U * C' * V';
                    p_diff = p_diff+norm(P(:,:,i)-P1(:,:,i),'fro');
                end
                P = P1;
%                 
%                 [U1,D1,V1] = svd(B);
%                 C1 = zeros(size(D1));
%                 for i = 1 : m
%                     C1(i,i) = sign(D1(i,i));
%                 end
%                 P1 = U1 * C1 * V1'; 
%                 haha = norm(P-P1,'fro')
                
                %% Update Si, namely S_record
                for i = 1:ker_num
                    A = H2(:,:,i) * P(:,:,i);
                    [Si] = qp(S, A, alpha);
                    S_record(:,:,i) = Si;
                end
                %                 Si_Value = cal_obj(H2, P, S, S_record, alpha)
                
                
                %% Update S
                S = mean(S_record, 3);
                S(S < 0) = 0;
                S(S > 1) = 1;
                %                 S_Value = cal_obj(H2, P, S, S_record, alpha)
                
                % objective value calculation
                it_count = it_count+1;
%                 cur_obj = cal_obj(H2, P, S, S_record, alpha);
%                 obj_record = [obj_record, cur_obj];
%                 
%                 if it_count > 2 && abs((obj_record(it_count-1)-obj_record(it_count-2))/ obj_record(it_count-2)) < 10^(-4)
%                     flag = 0;
%                 end
                if it_count > 2 && p_diff < 1
                     flag = 0;
                end
            end
            
            % Result using S for the best
            %             labels=litekmeans(S, CluNum, 'MaxIter', 100,'Replicates',10); %kmeans(U, c, 'emptyaction', 'singleton', 'replicates', 100, 'display', 'off');
            %             result=ClusteringMeasure(labels,Y);
            %             fprintf('result  %f\n', result);
            
%             S = gather(S);  % gpu
            res = myNMIACC(S,Y,CluNum);
            res_time = toc
            save(['./result_record_v2time/', dataName, 'result', num2str(m_i),'_', num2str(log2(alpha)),'.mat'], ...
                 'res', 'res_time')
        end
    end
end


%% The optimization process of Si
% min \alpha * ||S||_F^2  - 2 Tr((A+alpha * B)*S^T), s.t. 0<= S <=1
function [Si] = qp(A, B, alpha)   %A:S, B:A
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
    obj = obj + trace(Si_s(:,:,i) * P(:,:,i)' * H2(:,:,i));
end

obj = -obj;
% part 2
for i = 1 : ker_num
    obj = obj + alpha * norm(S-Si_s(:,:,i), 'fro')^2;
end

end
