function dist=kendall(x,y,opt,D)
% Compute Kendall's Rank Correlation between two global rankings x,y
%	function dist = kendall(x,y,opt,D)
%		opt = ['correlation']|'tau'
%

if nargin<3, 
	opt = 'correlation';
end

n=size(x,1);
if n==1, x=x'; n=length(x); end

n=size(y,1);
if n==1, y=y'; n=length(y); end

e=ones(n,1);
if strncmp(opt,'correlation',3),
	xx = x*e'-e*x';
	yy = y*e'-e*y';
	xx = squareform(xx)';
	yy = squareform(yy)';
	if nargin < 4,
		dist = (xx'*yy)/sqrt(xx'*xx)/sqrt(yy'*yy);
	else
		dist = (xx'*(D.*yy))/sqrt(xx'*(D.*xx))/sqrt(yy'*(D.*yy));
	end
else
	xx=squareform((x*e'-e*x')>0);
	yy=squareform((y*e'-e*y')>0);
	dist=sum(abs(xx-yy))/(n*(n-1));
end


% disp('Cheers!');
