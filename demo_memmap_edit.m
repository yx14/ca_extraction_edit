clear;
%% load file

%path_to_package = 'C:\Users\julia\Downloads\ca_source_extraction-master_mem\ca_source_extraction-master';   % path to the folder that contains the package
%addpath(genpath(path_to_package));
             
%filename = '';      % path to stack tiff file
%foldername = 'Q:\data\2photon\reg\150911_KS145_2P_KS';    % path to folder that contains a sequence of tiff files

%if exist([filename(1:end-3),'mat'],'file')
%    data = matfile([filename(1:end-3),'mat'],'Writable',true);
%else
%    sframe=1;						% user input: first frame to read (optional, default 1)
%    num2read=[];					% user input: how many frames to read   (optional, default until the end)
%    chunksize=1000;                 % user input: read and map input in chunks (optional, default read all at once)
%    %data = memmap_file(filename,sframe,num2read,chunksize);
%    data = memmap_file_sequence(foldername);
%end

%data2 = readtiff('Q:\data\2photon\reg\150911_KS145_2P_KS');
%data2 = readtiff('Q:\data\2photon\reg\160607_KS166_2P_KS\run03_ori12_V1', 1:96);
data2 = readtiff('\\nerffs01\mouselab\data\2photon\reg\140808_KS092_2P_KS\run02_ori_ds_V1', 1:32);
%%

%data is the cropped set of images 
%temp = imcrop(data2(:, :, 1), [253, 0.5, 261, 414]); first dataset
%temp = imcrop(data2(:, :, 1), [204,  240,  329,  141] [9,  169,  223,  185]);
temp = imcrop(data2(:, :, 1), rect);
data = zeros(size(temp, 1), size(temp, 2), size(data2, 3)); 
for i = 1: size(data2, 3)
    temp = imcrop(data2(:, :, i), rect);
    data(1:size(temp, 1), 1:size(temp, 2), i) = temp;
end

clear data2;
%%
%data = stackGroupProject(data(:, :, 1:8000), 2, 'sum');  %time downsampling 
data = data(:, :, 1:4000); 
disp('stacked');
%%
data = data(1:125, 1:125, 1:4096);

%% Set parameters
sizY = size(data);                  % size of data matrix
patch_size = [70, 70];  %[50]                 % size of each patch along each dimension (optional, default: [32,32])
overlap = [15,15];  %[15] [10, 10]                      % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1: end - 1),patch_size,overlap);

d1 = sizY(1);
d2 = sizY(2);
K = 15;                                           % number of components to be found
tau = 2;                                          % std of gaussian kernel (size of neuron) 
p = 2;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.75;                                  % merging threshold

options = CNMFSetParms(...                      
    'd1',d1,'d2',d2,'ssub', 1, 'tsub', 8, ...                         % dimensions of datasets
    'init_method', 'HALS', ... 
    'eta', 0.7, ...
    'beta', 0.5, ...
    'search_method','dilate', 'se', strel('disk', 2),... %'ellipse','dist',3,...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'fudge_factor',0.98,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr, 'sx', max(d1, d2)/2,...                    % merging threshold
    'thr2', 0.9, ... 
    'medw', [3, 3], ...
    'clos_op', strel('square', 3),...
    'gSig',tau...
    );
disp('params updated'); 
%% Run on patches

[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,p,options);

% temporal merge (identify multiple components of the same axon based on time alone)
%create Y with size dxT
data_res = reshape(data, [options.d1*options.d2, size(data, 3)]); 
[A_comb, C_comb, S_comb, f, P, YrA] = temporal_merge(data_res, A, b, C, f, P, options); 
disp('done'); 

% order and plot

[A_or,C_or,S_or,P] = order_ROIs(A_comb,C_comb,S_comb,P,options); % order components by linearity 

%% contour 
contour_threshold = 0.95;                 % amount of energy used for each component to construct contour plot
figure;
[Coor,json_file] = plot_contours(A_or,reshape(P.sn,sizY(1),sizY(2)),contour_threshold,1); % contour plot of spatial footprints
%savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)

%% display components
%(A'*Y- (A'*A)*C - (A'*full(b))*f) + C;

Yr2 = double(reshape(data, d1*d2,size(data, 3)));
Cn2 =  reshape(P.sn,d1,d2);

C_or2 = C_or(:, 1:4096); 
C_comb2 = C_comb(:, 1: 4096);
f2 = f(1:4096); 
Yr3 = Yr2(:, 1:4096);
plot_components_GUI(Yr2,A_or,C_or,b,f,Cn2,options)
%% cropping coordinates 

test = data2(:, :, 1);

figure;
[i2, rect] = imcrop(test);
 