x = 0:10;  y = sin(x);
xx = 0:.25:10;
yy = spline(x,y,xx);
plot(x,y,'o',xx,yy)

%%

% make sure X is a vector:
if length(find(size(x)>1))>1 
  error('MATLAB:chckxy:XNotVector','X must be a vector.') 
end

% ensure X is real
if any(~isreal(x)) 
  error('MATLAB:chckxy:XComplex','The X vector should have real elements.') 
end

% deal with NaN's among the sites:
nanx = find(isnan(x));
if ~isempty(nanx)
   x(nanx) = [];
   warning('MATLAB:chckxy:nan','All data points with NaN as their site will be ignored.')
end

n=length(x);
if n<2 
  error('MATLAB:chckxy:NotEnoughPts','There should be at least two data points.') 
end

% re-sort, if needed, to ensure strictly increasing site sequence:
x=x(:).'; 
dx = diff(x);

if any(dx<0), [x,ind] = sort(x); dx = diff(x); else ind=1:n; end

if ~all(dx), error('MATLAB:chckxy:RepeatedSites','The data sites should be distinct.'), end

% if Y is ND, reshape it to a matrix by combining all dimensions but the last:
sizey = size(y);


while length(sizey)>2&&sizey(end)==1, sizey(end) = []; end


yn = sizey(end); 
sizey(end)=[]; 
yd = prod(sizey);

if length(sizey)>1
   y = reshape(y,yd,yn);
else
   % if Y happens to be a column matrix, change it to the expected row matrix.
   if yn==1
       yn = yd;
       y = reshape(y,1,yn); 
       yd = 1; 
       sizey = yd;
   end
end

% determine whether not-a-knot or clamped end conditions are to be used:
nstart = n+length(nanx);
if yn==nstart
   endslopes = [];
elseif nargout==4&&yn==nstart+2
   endslopes = y(:,[1 n+2]); y(:,[1 n+2])=[];
   if any(isnan(endslopes))
      error('MATLAB:chckxy:EndslopeNaN','The endslopes cannot be NaN.')
   end
   if any(isinf(endslopes))
       error('MATLAB:chckxy:EndslopeInf','The endslopes cannot be Inf.')
   end
else
   error('MATLAB:chckxy:NumSitesMismatchValues',...
        ['The number of sites, ' int2str(nstart), ...
        ', is incompatible with the number of values, ' int2str(yn) '.'])
end

% deal with NaN's among the values:
if ~isempty(nanx)
    y(:,nanx) = [];
end

y=y(:,ind);
nany = find(sum(isnan(y),1));
if ~isempty(nany)
   y(:,nany) = []; x(nany) = [];
   warning('MATLAB:chckxy:IgnoreNaN','All data points with NaN in their value will be ignored.')
   n = length(x);
   if n<2 
     error('MATLAB:chckxy:NotEnoughPts', 'There should be at least two data points.') 
   end
end


%%
n = length(x); yd = prod(sizey);

% Generate the cubic spline interpolant in ppform

dd = ones(yd,1); dx = diff(x); divdif = diff(y,[],2)./dx(dd,:); 
if n==2
   if isempty(endslopes) % the interpolant is a straight line
      pp=mkpp(x,[divdif y(:,1)],sizey);
   else         % the interpolant is the cubic Hermite polynomial
      pp = pwch(x,y,endslopes,dx,divdif); pp.dim = sizey;
   end
elseif n==3&&isempty(endslopes) % the interpolant is a parabola
   y(:,2:3)=divdif;
   y(:,3)=diff(divdif')'/(x(3)-x(1));
   y(:,2)=y(:,2)-y(:,3)*dx(1); 
   pp = mkpp(x([1,3]),y(:,[3 2 1]),sizey);
else % set up the sparse, tridiagonal, linear system b = ?*c for the slopes
   b=zeros(yd,n);
   b(:,2:n-1)=3*(dx(dd,2:n-1).*divdif(:,1:n-2)+dx(dd,1:n-2).*divdif(:,2:n-1));
   if isempty(endslopes)
      x31=x(3)-x(1);xn=x(n)-x(n-2);
      b(:,1)=((dx(1)+2*x31)*dx(2)*divdif(:,1)+dx(1)^2*divdif(:,2))/x31;
      b(:,n)=...
      (dx(n-1)^2*divdif(:,n-2)+(2*xn+dx(n-1))*dx(n-2)*divdif(:,n-1))/xn;
   else
      x31 = 0; xn = 0; b(:,[1 n]) = dx(dd,[2 n-2]).*endslopes;
   end
   dxt = dx(:);
   c = spdiags([ [x31;dxt(1:n-2);0] ...
        [dxt(2);2*[dxt(2:n-1)+dxt(1:n-2)];dxt(n-2)] ...
        [0;dxt(2:n-1);xn] ],[-1 0 1],n,n);

   % sparse linear equation solution for the slopes
   mmdflag = spparms('autommd');
   spparms('autommd',0);
   s=b/c;
   spparms('autommd',mmdflag);
   
end


%% 

d=size(y,1);
dxd=repmat(dx,d,1);

n=numel(x);

dzzdx=(divdif-s(:,1:n-1))./dxd;
dzdxdx=(s(:,2:n)-divdif)./dxd;

dnm1=d*(n-1)

%%
%reshape((dzdxdx-dzzdx)./dxd,dnm1,1)
[reshape((dzdxdx-dzzdx)./dxd,dnm1,1)...
 reshape(2*dzzdx-dzdxdx,dnm1,1)...
 reshape(s(:,1:n-1),dnm1,1)...
 reshape(y(:,1:n-1),dnm1,1)]