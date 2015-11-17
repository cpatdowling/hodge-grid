% UK University Ranking Example
%   Reference: Xiaoye Jiang, Lek-Heng Lim, Yuan Yao, Yinyu Ye, Statistical
%   Ranking and combinatorial Hodge Theory, Mathematical Programming, 2010,
%   Section 7.3
%

%   Yuan Yao
%   Peking University
%   Beijing, 100871
%   China
%   yyao04@gmail.com
%
%   Last updated on 2010-10-27
%   From early versions at Stanford University, 2006-2009

%   UK Universities Web Link Structure dataset
%   http://www.cybermetrics.wlv.ac.uk/database/stats/data
%   data_uk.mat
%       - W: Web link structure
%       - webpage: name of university webpages
%   RAE_rank.mat: Research Assessment Exercise (http://www.rae.ac.uk) 
%       - RAE2001: RAE Score in 2001
%       - RAE1996: RAE Score in 1996
load data_uk W v webpage
load RAE_rank RAE2001 RAE1996

W_uk=W;
v_uk=v;

D = sum(W,2);
id = find(D>0 & RAE2001>0);
n = length(id);

D = D(id);
W = W(id,id);
T = zeros(n,n);
T = diag(1./D) * W;

alpha = .99;

% T1: PageRank with alpha=0.99 as random surfer
T1 = alpha * T + (1-alpha)*ones(n,1)*ones(1,n)/n;
% Diffused PageRank
T2 = T1^2;
T4 = T2^2;

% Load precomputed edges and triangles
load simplices_uk edges triangs        
n0 = length(T1);        
n1 = max(size(edges));        
d0 = zeros(n1,n0);        
w0 = zeros(n1,1);  
w1 = w0;     
w4 = w0;
for i=1:n1,                
    d0(i,edges(1,i))=-1;                
    d0(i,edges(2,i))=1;                
    w0(i) = log(T1(edges(1,i),edges(2,i)))-log(T1(edges(2,i),edges(1,i)));  
    w1(i) = log(T2(edges(1,i),edges(2,i)))-log(T2(edges(2,i),edges(1,i)));
    w4(i) = log(T4(edges(1,i),edges(2,i)))-log(T4(edges(2,i),edges(1,i)));
end        
L0 = d0'*d0;        
x0 = lsqr(L0,d0'*w0);
x1 = lsqr(L0,d0'*w1);
x4 = lsqr(L0,d0'*w4);

[evec,d]=eigs(T1',1);  % PageRank score

[U,S,V]=svds(W,1);  % SVD for HITS

score_page = abs(evec);         % PageRank
score_HITSauth= abs(V(:,1));    % HITS authority
score_HITShub= abs(U(:,1));     % HITS hub
score_out = D;                  % out-degree score
score_in = sum(W,1)';           % in-degree score
score_research01 = RAE2001(id); % research 2001 score
score_research96 = RAE1996(id); % research 1996 score
score_avg = (score_in + score_out)/2;   % total link
score_diff= score_in - score_out;       % in-out difference
score_hodge0= x0;                       % hodge rank score
score_hodge1= x1;                       % hodge rank with 2-iteration T2 
score_hodge4= x4;                       % hodge rank with 4-iteration T4

% Compute Kendall's tau-distance among 5 rankings
%s=[score_research01,score_research96,score_page,score_in,score_out,score_avg,score_diff,score_HITSauth,score_HITShub, score_hodge0,score_hodge1];
s=[score_research01,score_in, score_out, score_HITSauth, score_HITShub, score_page, score_hodge0, score_hodge1, score_hodge4];
N=size(s,2);
S=zeros(N,N);
for i=1:N,
    for j=i+1:N, 
        S(i,j)=kendall(s(:,i),s(:,j),'tau');
    end
end
S = S+S';
%disp('RAE2001  RAE1996  pagerank  in-degree  out-degree avg_link in-out_diff  HITS_auth HITS_hub hodge0 hodge1 hodge4');
disp('RAE2001 in-degree  out-degree HITS_auth HITS_hub  pagerank     hodge0  hodge1   hodge4');
disp(S)

% This shows the following Table 3, Section 7.3 in Reference
% Kendall's tau distance of 4 rankings
%  research  pagerank  in-degree  out-degree score_link in-out-difference hodgerank
%         0    0.4029    0.0986    0.1145    0.1016    0.2314    0.2299
%    0.4029         0    0.4939    0.4382    0.4693    0.2708    0.2722
%    0.0986    0.4939         0    0.0622    0.0310    0.2283    0.2274
%    0.1145    0.4382    0.0622         0    0.0313    0.2906    0.2894
%    0.1016    0.4693    0.0310    0.0313         0    0.2593    0.2581
%    0.2314    0.2708    0.2283    0.2906    0.2593         0    0.0094
%    0.2299    0.2722    0.2274    0.2894    0.2581    0.0094         0

