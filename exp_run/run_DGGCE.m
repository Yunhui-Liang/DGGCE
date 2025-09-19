%
%
%

clear;
clc;

data_path = fullfile(pwd, '..',  filesep, "data_BPs", filesep);
addpath(data_path);
lib_path = fullfile(pwd, '..',  filesep, "lib", filesep);
addpath(lib_path);
code_path = genpath(fullfile(pwd, '..',  filesep, 'DGGCE_ICASSP_2026', filesep));
addpath(code_path);


dirop = dir(fullfile(data_path, '*.mat'));
datasetCandi = {dirop.name};


exp_n = 'DGGCE_RES';
for i1 =1 : length(datasetCandi)
    data_name = datasetCandi{i1}(1:end-4);
    dir_name = [pwd, filesep, exp_n, filesep, data_name];
    create_dir(dir_name);
    clear BPs Y;
    load(strcat(data_path, datasetCandi{i1}));
    assert(size(BPs, 1) == size(Y, 1));
    nSmp = size(BPs, 1);
    nCluster = length(unique(Y));
    
    nBase = 20;
    
    %*********************************************************************
    % DGGCE_RES
    %*********************************************************************
    fname2 = fullfile(dir_name, [data_name, '_', exp_n, '.mat']);
    if ~exist(fname2, 'file')
        nRepeat = 10;
        
        seed = 2024;
        rng(seed, 'twister')
        
        % Generate 50 random seeds
        random_seeds = randi([0, 1000000], 1, nRepeat * nRepeat );
        
        % Store the original state of the random number generator
        original_rng_state = rng;
        
        nMeasure = 15;
        
        DGGCE_RES_result = zeros(nRepeat, nMeasure);
        DGGCE_RES_time = zeros(nRepeat, 1);
        
        for iRepeat = 1:nRepeat
            idx = (iRepeat - 1) * nBase + 1 : iRepeat * nBase;
            BPi = BPs(:, idx);
                        
            t1_s = tic;
           
            nSub = 6;
            nSubBase = 5;
            seed = 2026;
            label = DGGCE(BPi, nSub, 5, seed, nCluster);
                        
            t1 = toc(t1_s);
            result_10 = my_eval_y(label, Y);
            DGGCE_RES_result(iRepeat, :) = [result_10', t1];
        end
        DGGCE_RES_result_summary = mean(DGGCE_RES_result, 1);
        DGGCE_RES_result_summary_std = std(DGGCE_RES_result);
        save(fname2, 'DGGCE_RES_result', 'DGGCE_RES_result_summary', 'DGGCE_RES_result_summary_std');
        
        disp([data_name, ' has been completed!']);
    end
end
rmpath(data_path);
rmpath(lib_path);
rmpath(code_path);