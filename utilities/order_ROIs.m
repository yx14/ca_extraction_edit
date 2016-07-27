function [A_or,C_or,S_or,P_or,srt] = order_ROIs(A,C,S,P,options, srt)

% ordering of the found components based on their maximum temporal
% activation and their size (through their l_inf norm)
% you can also pre-specify the ordering sequence


 
score_arr = zeros(5, size(A, 2) );
%remove if round, based on distribution of x and y coordinates 
%remove if most of the pixels [x, y] in [range(x), range(y)] are filled 
for i = 1: size(A, 2)
    res = reshape(A(:, i), [options.d1, options.d2]); 
    [xpos, ypos, ~] = find(res); 
     
    x_val = std(xpos)/range(xpos);
    y_val = std(ypos)/range(ypos);
    if x_val/y_val > 1.1 || x_val/y_val < 0.9
        score_arr(1, i) = 1;
    end
   
    filled = length(find(res(min(xpos):max(xpos) , min(ypos):max(ypos))));
    if filled/(range(xpos)*range(ypos)) < 0.5;
        score_arr(3, i) = 1;
    end
    
                
end


%remove if pixels out of range (px_min, px_max)
px_min = 100; 
px_max = 10000; 
for i = 1: size(A, 2) 
    [x, ~, ~] = find(A(:, i));
    if length(x) < px_max && length(x) > px_min
        score_arr(2, i) =  1;
    end
end

%remove components, assuming there's at least one 'real' axon 
for i = size(A, 2): -1: 1
    temp = sum(score_arr); 
    if temp(i) ~= 3
        A(:, i) = []; 
        C(i, :) = [];
        S(i, :) = []; 
    end
end

%reshape A to d1*d2*T;  
A_or3 = reshape(full(A), [options.d1, options.d2, size(A, 2)]);
%calculate the R^2 values of x and y
r_vals = zeros(size(A, 2), 1); 
for i = 1: size(A, 2);
    A_vec = A_or3(:, :, i);
    %get x,y coords and values of nonzero indices 
    [x_coords, y_coords, ~] = find(A_vec);
    x_coords(:, 2) = ones(size(y_coords, 1), 1);
    [~, ~, ~, ~, stats1] = regress(y_coords, x_coords); 
    r_val1 = stats1(1); 
    %rotate it 30 degrees and calculate the correlation again (in case it's
    %vertical/horizontal the first time). appends maximum of the two
    %correlations to r_vals. 
    x_coords2 = cos(30)*x_coords(:, 1) - sin(30)*y_coords;
    y_coords2 = sin(30)*x_coords(:, 1) + cos(30)*y_coords;
    x_coords2(:, 2) = ones(size(y_coords, 1), 1); 
    [~, ~, ~, ~, stats2] = regress(y_coords2, x_coords2);
    r_val2 = stats2(1); 
    r_vals(i, 1) = max(r_val1, r_val2);  
end

nA = sqrt(sum(A.^2));
nr = length(nA);
%A = A/spdiags(nA(:),0,nr,nr);
%C = spdiags(nA(:),0,nr,nr)*C;
mA = sum(A.^4).^(1/4);
%sA = sum(A); initially blocked out 
mC = max(C,[],2);

if ~exist('srt', 'var')||isempty(srt)
    %[~, srt] = sort(r_vals.*pixl_dif, 'descend');
    [~,srt] = sort(r_vals.*mC.*mA','descend');
end

A_or = A(:,srt);
C_or = C(srt,:);
if nargin < 4
    P_or = [];
else
    P_or = P;
    %JUST FOR NOW 
    %{
    if isfield(P,'gn'); P_or.gn=P.gn(srt); end
    if isfield(P,'b'); P_or.b=P.b(srt); end
    if isfield(P,'c1'); P_or.c1=P.c1(srt); end
    if isfield(P,'neuron_sn'); P_or.neuron_sn=P.neuron_sn(srt); end
    %}
end

if nargin < 3 || isempty(S)
    S_or = [];
else
    S_or = S(srt,:);
end