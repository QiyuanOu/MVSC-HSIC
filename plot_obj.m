% Note that: the repo './vresult_record' is the measurements of 13
% datasetes using ablation2; and that './result_record' is the same for ab3.

clc
clear


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


total_result = zeros(12, 3);
% parameters1 = 1 : 3;
parameters1 = [0,0,0];
parameters2 = -15 : 2: 15;
result_record = zeros(length(parameters1), length(parameters2), 3);

for i_count = 1 : 13
    dataName = DataName{i_count};
    for m_i = 1 : length(parameters1)
        m = parameters1(m_i);
        for alpha_i = 1 : length(parameters2)
            alpha = parameters2(alpha_i);
            load(['./result_record/', dataName, 'result', num2str(m),...
                '_', num2str(alpha),'.mat'], 'res')
%             plot(obj_record)
            result_record(m_i, alpha_i, 1) = res(1);
            result_record(m_i, alpha_i, 2) = res(2);
            result_record(m_i, alpha_i, 3) = res(3);
        end
    end
%     save([dataName, '_TotalResult.mat'], 'result_record')
    best1 = max(max(result_record(:,:,1)));
    best2 = max(max(result_record(:,:,2)));
    best3 = max(max(result_record(:,:,3)));
    [posx1,posy1] = find(result_record(:,:,1) == best1)
    [posx2,posy2] = find(result_record(:,:,2) == best2)
    [posx3,posy3] = find(result_record(:,:,3) == best3)
    
    total_result(i_count, 1) = best1;
    total_result(i_count, 2) = best2;
    total_result(i_count, 3) = best3;
%     save([dataName, '_TotalResult.mat'], 'result_record')
%     save('Total_result.mat', 'total_result')
end