% Currency Exchange Example
%   Reference: Xiaoye Jiang, Lek-Heng Lim, Yuan Yao, Yinyu Ye, Statistical
%   Ranking and combinatorial Hodge Theory, Mathematical Programming, 2010,
%   Section 7.2
%

%   Yuan Yao
%   Peking University
%   Beijing, 100871
%   China
%   yyao04@gmail.com
%
%   Last updated on 2010-10-27
%   From early versions at Stanford University, 2006-2009

% Yahoo!Finance currency exchange data collected on November 6, 2007
%   - A: is a currency exchange matrix such that A(i,j) is the amount of
%   money j in exchange of 1 unit of money i.
load currency.mat

% Pairwise comparison matrix (skew-symmetric)
W = log(A);

% Pairwise comparison graph (complete graph here)
G = (abs(W)>0);

% 0-order coboundary operator
d0 = [];
for i=1:size(G,1),
    for j=i+1:size(G,1),
        d0_row = zeros(1,length(G));
        d0_row(i) = 1;
        d0_row(j) = -1;
        d0 = [d0; d0_row];
    end
end

% vectorize the pairwise comparison edge flow
w = squareform(W)';

% 0-Laplacian or Graph Laplacian
L0 = d0'*d0;

% Divergence
div = d0'*w;

% Global ranking score
x_global = lsqr(L0,-div);

disp('*************************')
disp('An Universal Equivalent for 7 Currencies')
exp(x_global)   % Table 2, Section 7.2 in reference.