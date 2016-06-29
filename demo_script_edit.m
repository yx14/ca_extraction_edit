clear;
%% load file

addpath(genpath('utilities'));
             
%nam = 'demoMovie.tif';          % insert path to tiff stack here
%sframe=1;						% user input: first frame to read (optional, default 1)
%num2read=2000;					% user input: how many frames to read   (optional, default until the end)

%Y = bigread2(nam,sframe,num2read);

%Y = readtiff('Q:\data\2photon\reg\150911_KS145_2P_KS', 1:100);
Y = readtiff('Q:\data\2photon\reg\160607_KS166_2P_KS\run03_ori12_V1', 1:32);
%Y = stackGroupProject(Y, 8, 'sum');
Y = Y - min(Y(:)); 
if ~isa(Y,'double');    Y = double(Y);  end         % convert to double

[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;                                          % total number of pixels
%% Set parameters

K = 100;                                           % number of components to be found
tau = 2;                                          % std of gaussian kernel (size of neuron) 
p = 2;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                                  % merging threshold

options = CNMFSetParms(...                      
    'd1',d1,'d2',d2, 'ssub', 2, 'tsub', 2, ...                         % dimensions of datasets
    'init_method', 'HALS', ... 
    'eta', 1, ...
    'beta', 0.5, ...
    'search_method','dilate', 'se', strel('disk', 3, 0),... %'ellipse','dist',3,...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'fudge_factor',0.98,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr, 'sx', max(d1, d2)/2,...                    % merging threshold
    'medw', [3, 3], ...
    'clos_op', strel('square', 3),...
    'gSig',tau...
    );
%% Data pre-processing
Y = double(Y);
Y = reshape(data, [options.d1, options.d2, size(data, 3)]); 
[P,Y] = preprocess_data(Y,p);
disp('preprocessing complete');
%%
Y = double(Y);
%% fast initialization of spatial components using greedyROI and HALS

[Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options);  % initialize

% display centers of found components
Cn =  reshape(P.sn,d1,d2); %correlation_image(Y); %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
%%
figure;imagesc(Cn);
    axis equal; axis tight; hold all;
    scatter(center(:,2),center(:,1),'mo');
    title('Center of ROIs found from initialization algorithm');
    drawnow;

%% manually refine components (optional)
refine_components = false;  % flag for manual refinement
if refine_components
    [Ain,Cin,center] = manually_refine_components(Y,Ain,Cin,center,Cn,tau,options);
end
    
%% update spatial components
 
Yr = reshape(Y, [size(Y, 1)*size(Y, 2), size(Y, 3)]); 
%clear Y;
%[A,b,Cin] = update_spatial_components(Yr,Cin,fin,Ain,P,options);

%% clean up spatial components
th = 1e-2;  % for each component set pixels below th*max_value to 0
for i = 1:K
    Ain(Ain(:,i)<th*max(Ain(:,i)),i) = 0;
end
Ain = sparse(Ain);

%% update temporal components
P.p = 0;    % set AR temporarily to zero for speed
[C,f,P,S] = update_temporal_components(Yr,Ain,bin,Cin,fin,P,options);

%% merge found components
[Am,Cm,K_m,merged_ROIs,P,Sm] = merge_components(Yr,Ain,bin,C,f,P,S,options);

%%
display_merging = 1; % flag for displaying merging example
if and(display_merging, ~isempty(merged_ROIs))
    i = 1; %randi(length(merged_ROIs));
    ln = length(merged_ROIs{i});
    figure;
        set(gcf,'Position',[300,300,(ln+2)*300,300]);
        for j = 1:ln
            subplot(1,ln+2,j); imagesc(reshape(A(:,merged_ROIs{i}(j)),d1,d2)); 
                title(sprintf('Component %i',j),'fontsize',16,'fontweight','bold'); axis equal; axis tight;
        end
        subplot(1,ln+2,ln+1); imagesc(reshape(Am(:,K_m-length(merged_ROIs)+i),d1,d2));
                title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight; 
        subplot(1,ln+2,ln+2);
            plot(1:T,(diag(max(C(merged_ROIs{i},:),[],2))\C(merged_ROIs{i},:))'); 
            hold all; plot(1:T,Cm(K_m-length(merged_ROIs)+i,:)/max(Cm(K_m-length(merged_ROIs)+i,:)),'--k')
            title('Temporal Components','fontsize',16,'fontweight','bold')
        drawnow;
end

%% repeat
P.p = p;    % restore AR value
%[A2,b2,Cm] = update_spatial_components(Yr,Cm,f,Am,P,options);
[C2,f2,P,S2] = update_temporal_components(Yr,Am,bin,Cm,f,P,options);

%% do some plotting

[A_or,C_or,S_or,P] = order_ROIs(Am,Cm,S,P, options); % order components
%%
K_m = size(C_or,1);
[C_df,~] = extract_DF_F(Yr,[A_or,b2],[C_or;f2],K_m+1); % extract DF/F values (optional)

contour_threshold = 0.95;                       % amount of energy used for each component to construct contour plot
figure;
[Coor,json_file] = plot_contours(A_or,reshape(P.sn,d1,d2),contour_threshold,1); % contour plot of spatial footprints
%savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)
%% display components


plot_components_GUI(Yr,A_or,C_or,bin,f2,Cn,options);
%% make movie

make_patch_video(A_or,C_or,bin,f,Yr,Coor,options)