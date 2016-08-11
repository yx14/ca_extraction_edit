function run_cnmf_on_cwd()
%%
data_path = '\\nerffs01\mouselab\data\2photon\reg\140808_KS092_2P_KS\run02_ori_ds_V1'; 
data2 = readtiff(data_path, 1:32);

%% crop data
rect = [128, 140, 140, 140];
 
temp = imcrop(data2(:, :, 1), rect);
data = zeros(size(temp, 1), size(temp, 2), size(data2, 3)); 
for i = 1: size(data2, 3)
    temp = imcrop(data2(:, :, i), rect);
    data(1:size(temp, 1), 1:size(temp, 2), i) = temp;
end

clear data2;

%% Set parameters
sizY = size(data);                  % size of data matrix
patch_size = [70, 70];  %[50]                 % size of each patch along each dimension (optional, default: [32,32])
overlap = [18,18];  %[15] [10, 10]                      % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1: end - 1),patch_size,overlap);

d1 = sizY(1);
d2 = sizY(2);
K = 8;                                           % number of components to be found
tau = 2;                                          % std of gaussian kernel (size of neuron) 
p = 2;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.75;                                  % merging threshold

options = CNMFSetParms(...                      
    'd1',d1,'d2',d2,'ssub', 1, 'tsub', 8, ...                         % dimensions of datasets
    'init_method', 'HALS', ... 
    'eta', 0.7, ... %for sparse NMF initialization
    'beta', 0.5, ... %for sparse NMF initialization 
    'search_method','dilate', 'se', strel('disk', 2),... %'ellipse','dist',3,...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'fudge_factor',0.98,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr, 'sx', max(d1, d2)/2,...                    % merging threshold
    'thr2', 0.9, ... 
    'medw', [3, 3], ...
    'clos_op', strel('square', 3),...
    'gSig',tau, ...
    'px_min', 100, ... %minimum number of pixels for an ordered component
    'px_max', 10000, ... %maximum number of pixels for an ordered component
    'max_df_f', 0.75 ... %highest Df/F peak must exceed max_df_f for an ordered component
    );
disp('params updated'); 

%% Run on patches

[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,p,options); 

% temporal merge (identify multiple components of the same axon based on time alone)
% create Y with size dxT
data_res = reshape(data, [options.d1*options.d2, size(data, 3)]); 
[A_comb, C_comb, S_comb, f, P, YrA] = temporal_merge(data_res, A, b, C, f, P, options); 
disp('done'); 
%% order components

[A_or,C_or,S_or,~] = order_ROIs(data_res, A_comb,C_comb,b, f, S_comb,P,options);  

%% save results
save('results', 'options', 'rect', 'A_comb', 'C_comb', 'S_comb', 'f', 'P', 'YrA', 'A_or', 'C_or', 'S_or', '-v7.3');

