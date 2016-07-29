function [A_comb2,C_comb,S_comb,f, P, YrA] = temporal_merge(Y, A,b, C,f, P, options)

%merges components if they have strong temporal correlation regardless of
%spatial overlap 

%outputs: 
%A_comb2 is A with components merged and others removed
%C and S just keep the components from merge_x, NEED TO FIX: track components from
%merge_y 

thr = options.thr2; 
nr = size(A, 2); 
%produce a correlation matrix where each element is the temporal correlation between
%two components
corr_m = corr(full(C(1:nr, :)'));
corr_m = triu(corr_m) >= thr; 
%remove diagonal
corr_m = corr_m - diag(diag(corr_m));
[merge_x, merge_y, ~] = find(corr_m); 
merge_t = unique(merge_x);
merge_u = size(merge_t, 1); 
A_comb = zeros(size(A, 1), size(A, 2));
C_comb = zeros(size(A, 2) - merge_u, size(C, 2)); 
 
  
%A_comb has components merged within merge_x elements
for i = 1: size(A, 2) 
    if any(i == merge_x)
        [ind, ~] = find(merge_x == i); 
        for j = ind'
            A_comb(:, i) = A(:, i) + A(:, merge_y(j));
        end
    else
        A_comb(:, i) = A(:, i); 
    end
end
merge_ct = 1;
A_comb2 = zeros(size(A, 1), size(A, 2) - merge_u);
%A_comb2 drops merge_y components 
for i = 1: size(A, 2) 
    if ~any(i == merge_y)
        A_comb2(:, merge_ct) = A_comb(:, i);
        %initialize with the first component's temporal footprint and update
        %afterward 
        C_comb(merge_ct, :) = C(i, :);
        merge_ct = merge_ct + 1; 
    end
end
%update temporal components
 
options.temporal_iter = 1;
[C_comb,f,P,S_comb,YrA] = update_temporal_components(Y,A_comb2,b,C_comb,f,P,options);

disp(merge_u);
end