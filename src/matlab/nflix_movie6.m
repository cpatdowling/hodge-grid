% Netflix Example
%   Reference: Xiaoye Jiang, Lek-Heng Lim, Yuan Yao, Yinyu Ye, Statistical
%   Ranking and combinatorial Hodge Theory, Mathematical Programming, 2010,
%   Section 7.1
%

%   Yuan Yao
%   Peking University
%   Beijing, 100871
%   China
%   yyao04@gmail.com
%
%   Last updated on 2010-10-27
%   From early versions at Stanford University, 2006-2009

% Load Netflix data for 6 movies ratings over 74 months, 1999-11 ~ 2005-12
%   - movieTitle contains:
%       'Shakespeare in Love'
%       'Witness'
%       'October Sky'
%       'The Waterboy'
%       'Interview with the Vampire'
%       'Dune'
%   - movieID contains movie ID in netflix
%   - movieYear is the DVD year 
%   - Y: ratings for each of 74 months
%   - numRates: number of ratings for each movie over 74 months
%   - scoreMean: mean score of each movie over 74 months
%   - mrqe=[85 77 76 44 65 66]: movie review query engine scores in July
%   2009, available at http://www.mrqe.com
%
load nflix_movie6.mat

% Compute the mean score for each of 74 months
mn=[];
for i=1:74, mn(i,:)=sum(Y{i},1)./sum(Y{i}>0,1);end
num=[];
scoreSum=[];
for i=1:74, scoreSum(i,:)=sum(Y{i},1); num(i,:)=sum(Y{i}>0,1);end 
% Compute the global mean score
mm0 = sum(scoreSum,1)./sum(num,1);
 
% Draw picture to show the temporal drifts of mean scores over 74 months
% Figure 4, Section 7.1, in JLYY 2010. 
for k=1:6,
	subplot(2,3,k),
	b(:,k) = regress(mn(:,k),[1:74; ones(1,74)]');
	plot(1:74,mn(:,k),'b*',1:74,[1:74; ones(1,74)]'*b(:,k),'r-.'),
	axis([1 75 2 5]);
end
for k=1:6,subplot(2,3,k),title(sprintf('%s',char(movieTitle{k})));end                                      
subplot(2,3,5),xlabel('Month Number')                                      
subplot(2,3,4),ylabel('Average Score')
subplot(2,3,1),ylabel('Average Score')
saveas(gcf,'Drift_movie6.jpg','jpg');

eps=1e-6;
G=zeros(6,6,74);    % a graph for each of 74 months
GG=zeros(6,6);      % average graph
WW=zeros(6,6);      % pairwise comparison skew-symmetric matrix
WW1 = WW;
WW2 = WW;
WW3 = WW;
for l=1:74,
	edges{l}=[];
	W{l}=[];
	d0{l}=[];
	for i=1:6,
		for j=i+1:6,
			k=find(Y{l}(:,i)>0 & Y{l}(:,j)>0); 
			if k>0,
				edges{l}=[edges{l},[i;j]];
				G(i,j,l)=1;
				G(j,i,l)=1;
				GG(i,j) = GG(i,j) +length(k);
				GG(j,i) = GG(j,i) +length(k);
                
                % Difference model
				WW(i,j) = WW(i,j) + full(sum(Y{l}(k,i)-Y{l}(k,j),1));   
				WW(j,i) = WW(j,i) + full(sum(Y{l}(k,j)-Y{l}(k,i),1));
				
                % Log ratio
				WW1(i,j) = WW1(i,j) + full(sum(log(Y{l}(k,i)) - log(Y{l}(k,j)),1)); 
				WW1(j,i) = WW1(j,i) + full(sum(log(Y{l}(k,j)) - log(Y{l}(k,i)),1));
				
                % Binary Model (Uniform Probability Model)
				WW2(i,j) = WW2(i,j) + full(sum(Y{l}(k,i)>Y{l}(k,j),1) - sum(Y{l}(k,i)<Y{l}(k,j),1));    
				WW2(j,i) = - WW2(i,j);
				
                % Bradley-Terry model or Log Odds. (This model is not shown
                % in the paper due to its poor performance. It can be
                % updated in the future)
				WW3(i,j) = WW3(i,j) + full(log(eps+sum(Y{l}(k,i)>Y{l}(k,j),1))-log(eps+sum(Y{l}(k,j)>Y{l}(k,i),1)));   
				WW3(j,i) = -WW3(i,j);
			end
		end
	end
end

% Compute the 0-order coboundary operator dd0
dd0 = [];

% Compute the 1-dimensional simplices, i.e. edges
ss1 = [];
for i=1:length(GG),
    for j=(i+1):length(GG),
        if G(i,j)>0, 
            ss1 = [ss1; [i,j]];
            dd_row = zeros(1,length(GG));
            dd_row(i) = 1;
            dd_row(j) = -1;
            dd0 = [dd0; dd_row];
        end
    end
end
ss1 = ss1'; % edge matrix of 2-by-(# edges), each column is a node pair  

% Compute the 1st-order coboundary operator dd1 and 2-dimensional 
% simplices, i.e triangles.
dd1 = [];
ss2 = [];
for i=1:length(GG),
    for j=(i+1):length(GG),
        for k=(j+1):length(GG),
            if ((GG(i,j)>0)&(GG(j,k)>0)&(GG(k,i)>0)),
                ss2 = [ss2; [i,j,k]];
                dd_row = zeros(1,size(dd0,1));
                dd_row(find(ss1(1,:)==i&ss1(2,:)==j))=1;
                dd_row(find(ss1(1,:)==j&ss1(2,:)==k))=1;
                dd_row(find(ss1(1,:)==i&ss1(2,:)==k))=-1;
                dd1 = [dd1; dd_row];
            end
        end
    end
end
% triangle matrix of 3-by-(# triangles), each column is a 
%   node triple for triangle
ss2 = ss2'; 

% Compute graph laplacian 'LL?', divergence 'div?', and global ranking 'xx?'
% Cyclicity Ratio as residues are computed, 'res?' 
LL0 = dd0'*diag(squareform(GG))*dd0;
ww0 = squareform(WW')';
div0 = dd0'*ww0;
xx0 = lsqr(LL0,div0);
res0 = (diag(1./squareform(GG))*ww0-dd0*xx0)'*diag(squareform(GG))*(diag(1./squareform(GG))*ww0-dd0*xx0)/(ww0'*diag(1./squareform(GG))*ww0);

ww1 = squareform(WW1')';
div1 = dd0'*ww1;
xx1 = lsqr(LL0,div1);
res1 = (diag(1./squareform(GG))*ww1-dd0*xx1)'*diag(squareform(GG))*(diag(1./squareform(GG))*ww1-dd0*xx1)/(ww1'*diag(1./squareform(GG))*ww1);

ww2 = squareform(WW2')';
div2 = dd0'*ww2;
xx2 = lsqr(LL0,div2);
res2 = (diag(1./squareform(GG))*ww2-dd0*xx2)'*diag(squareform(GG))*(diag(1./squareform(GG))*ww2-dd0*xx2)/(ww2'*diag(1./squareform(GG))*ww2);

LL3 = dd0'*dd0;
ww3 = squareform(WW3')';
div3 = dd0'*ww3;
xx3 = lsqr(LL3,div3);
res3 = norm(ww3- dd0*xx3)^2/norm(ww3)^2;

disp('Cyclicity Ratio (%)')
disp(sprintf('Hodge-difference: %1.2f%%; Hodge-ratio: %1.2f%%; Hodge-binary: %1.2f%%',res0*100,res1*100,res2*100))
disp('---- Hodge-binary has the smallest cyclicity ratio! (Table 1, Section 7.1 in Reference)')

% Find the ranking
[ig,rkm]=sort(mm0,'descend');
[ig,rk0]=sort(xx0,'descend');
[ig,rk1]=sort(xx1,'descend');
[ig,rk2]=sort(xx2,'descend');
[ig,rk3]=sort(xx3,'descend');

% Compute trianglular curls
curl0=dd1*diag(1./squareform(GG))*ww0;
curl1=dd1*diag(1./squareform(GG))*ww1;
curl2=dd1*diag(1./squareform(GG))*ww2;
curl3=dd1*ww3;
% relative curl on each edge
for i=1:size(ss2,2),
    curl_rel0(:,i)=curl0(i)/3./[xx0(ss2(1,i))-xx0(ss2(2,i)), xx0(ss2(2,i))-xx0(ss2(3,i)), xx0(ss2(1,i))-xx0(ss2(3,i))]';
	curl_rel1(:,i)=curl1(i)/3./[xx1(ss2(1,i))-xx1(ss2(2,i)), xx1(ss2(2,i))-xx1(ss2(3,i)), xx1(ss2(1,i))-xx1(ss2(3,i))]';
    curl_rel2(:,i)=curl2(i)/3./[xx2(ss2(1,i))-xx2(ss2(2,i)), xx2(ss2(2,i))-xx2(ss2(3,i)), xx2(ss2(1,i))-xx2(ss2(3,i))]';
    curl_rel3(:,i)=curl3(i)/3./[xx3(ss2(1,i))-xx3(ss2(2,i)), xx3(ss2(2,i))-xx3(ss2(3,i)), xx3(ss2(1,i))-xx3(ss2(3,i))]';
end

disp('Edges with large relative curls (0: Hodge-difference):')
[ii,jj]=find(abs(curl_rel0)>1);
for i=1:length(ii),
    node_a = ss2(ii(i),jj(i));
    if ii(i)>= 3, 
        node_b = ss2(1, jj(i));
    else
        node_b = ss2(ii(i)+1,jj(i));
    end
    
    disp([node_a,node_b, curl_rel0(ii(i),jj(i))]);
end
% Edges with large relative curls (0):
%     1.0000    3.0000   -1.0709
% 
%     2.0000    3.0000   -4.6041
% 
%     2.0000    3.0000   -4.7994
% 
%     2.0000    3.0000    3.9868
% 

disp('Edges with large relative curls (1: Hodge-ratio):')
[ii,jj]=find(abs(curl_rel1)>1);
for i=1:length(ii),
    node_a = ss2(ii(i),jj(i));
    if ii(i)>= 3, 
        node_b = ss2(1, jj(i));
    else
        node_b = ss2(ii(i)+1,jj(i));
    end
    
    disp([node_a,node_b, curl_rel1(ii(i),jj(i))]);
end
% Edges with large relative curls (1):
%     1.0000    3.0000  -27.4334
% 
%     1.0000    3.0000   -1.5573
% 
%     1.0000    3.0000    7.1993
% 


disp('Edges with large relative curls (2: Hodge-binary):')
[ii,jj]=find(abs(curl_rel2)>1);
for i=1:length(ii),
    node_a = ss2(ii(i),jj(i));
    if ii(i)>= 3, 
        node_b = ss2(1, jj(i));
    else
        node_b = ss2(ii(i)+1,jj(i));
    end
    
    disp([node_a,node_b, curl_rel2(ii(i),jj(i))]);
end
% Edges with large relative curls (2):
%     2.0000    3.0000    3.6039
% 
%     2.0000    3.0000    4.1338

disp('Edges with large relative curls (3: Hodge-Bradley-Terry):')
[ii,jj]=find(abs(curl_rel3)>1);
for i=1:length(ii),
    node_a = ss2(ii(i),jj(i));
    if ii(i)>= 3, 
        node_b = ss2(1, jj(i));
    else
        node_b = ss2(ii(i)+1,jj(i));
    end
    
    disp([node_a,node_b, curl_rel3(ii(i),jj(i))]);
end
% Edges with large relative curls (3):
%     2.0000    3.0000    3.2007
% 
%     1.0000    2.0000   -6.2463
% 
%     1.0000    3.0000   -2.8604
% 
%     1.0000    3.0000   -3.8011
% 
%     1.0000    3.0000   -4.5732
% 
%     5.0000    6.0000   -1.1601
% 
%     2.0000    3.0000  -14.3464
% 
%     2.0000    3.0000   16.7838
% 
%     2.0000    3.0000   17.1338
% 
%     5.0000    6.0000   -1.3125
% 
%     5.0000    6.0000   -1.3365

