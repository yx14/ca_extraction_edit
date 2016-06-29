% julia demo

%% read the files
folder = 'Q:\data\2photon\reg\150911_KS145_2P_KS';
%addpath(genpath('/Users/epnevmatikakis/Documents/MATLAB/github/ca_source_extraction'));
filenames = dir([folder,'/*.tif']);

max_num_files = 35; % maximum number of files to read
num_files = min(max_num_files,length(filenames));

Tmax = 10000;
temp_filename = [folder,'/',filenames(1).name];
Y_temp = bigread2(temp_filename);
rm_pix = 50;  % remove first 50 rows
Y_temp = Y_temp(rm_pix+1:end,:,:);
Y_temp = imresize(Y_temp,0.5);
[d1,d2,~] = size(Y_temp);
d = d1*d2;
Y = zeros(d1,d2,Tmax);
cnt = 0;
for i = 1:num_files;
    if i > 1
        temp_filename = [folder,'/',filenames(i).name];
        Y_temp = bigread2(temp_filename);
        Y_temp = Y_temp(rm_pix+1:end,:,:);
        Y_temp = imresize(Y_temp,0.5); %dsData(Y_temp,options);
    end        
    Y(:,:,cnt+(1:size(Y_temp,3))) = Y_temp;
    cnt = cnt + size(Y_temp,3);
    if cnt + size(Y_temp,3) > Tmax
        break;
    end
end
Y(:,:,cnt+1:end) = [];
T = cnt;

%% set up options
K = 100;                                          % number of components to be found
tau = [];                                         % std of gaussian kernel (size of neuron) 
p = 1;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                                  % merging threshold

options = CNMFSetParms(...                      
    'd1',d1,'d2',d2,...                         % dimensions of datasets
    'search_method','dilate',...                % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'fudge_factor',0.98,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr,...                    % merging threshold
    'gSig',tau,...
    'init_method','HALS',...
    'tsub',2,'ssub',2,...
    'temporal_parallel',0 ...
    );


%% pre-process and initialize
[P,Y] = preprocess_data(Y,p);
[Ain, Cin, bin, fin] = initialize_components(Y, K, tau, options);

%% clean up spatial components
th = 1e-2;  % for each component set pixels below th*max_value to 0
for i = 1:K
    Ain(Ain(:,i)<th*max(Ain(:,i)),i) = 0;
end
Ain = sparse(Ain);

%% deconvolution
Cn =  reshape(P.sn,d1,d2);
Yr = reshape(Y,d,T);
%clear Y;
 
[C,f,P,S,YrA] = update_temporal_components(Yr,Ain,bin,Cin,fin,P,options);
[Am,Cm,K_m,merged_ROIs,P,Sm] = merge_components(Yr,Ain,bin,C,f,P,S,options);
 
[C2,f2,P,S2] = update_temporal_components(Yr,Am,bin,Cm,f,P,options);
%% do some plotting

[A_or,C_or,S_or,P] = order_ROIs(Am,C2,S2,P); % order components
figure;
plot_components_GUI(Yr,A_or,C_or,bin,f2,Cn,options)

%%
figure;
for i = 8:8
    %subplot(131);imagesc(reshape(Ain(:,i),d1,d2)); title(num2str(i)); colorbar;
    if i <= K
        subplot(3,1,[1,2]);imagesc(reshape(A_or(:,i),d1,d2)); title(num2str(i)); axis equal; axis tight;
        subplot(3,1,3); plot(C_or(i,:)+YrA(i,:),'r'); hold on; plot(C(i,:),'b','linewidth',2); hold off;
        %hold all;plot(C2(i,:));
    else
        subplot(3,1,[1,2]);imagesc(reshape(b,d1,d2)); title('background'); axis equal; axis tight;
        subplot(3,1,3);  plot(f,'b','linewidth',2); hold off;
    end
    hold off;
    pause;
end