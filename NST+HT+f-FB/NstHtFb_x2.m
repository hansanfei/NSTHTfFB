function [mu_new,NbIter] =NstHtFb_x2(y,y0,A,phi,MaxNbIter,epsilon)
% AdptNST+HT+f-FB with f(x)=x^2;
% Input:
%  observed data: y,   pre-computed vector: y0=A'*((A*A')\y);
%  measurement matrix: A,  pre-computed matrix: phi=A'*((A*A') \A);
%  maximum number of iteration: MaxNbIter, tolerance on the relative error: epsilon.
% Output:
% recovered sparse signal: mu_new
N=size(A,2);

%% initialization
% mu=ones(N,1);
mu_new=zeros(N,1);
NbIter=0;

%% main loop
while ( NbIter < MaxNbIter  && (norm(y-A*mu_new) > epsilon) )
  xnew=mu_new+y0-phi*mu_new;
  absv=abs(xnew);
  %     zero_idx=find(absv<Eps*max(absv));
  %     absv(zero_idx)=zeros(size(zero_idx));
  [~,sorted_idx]=sort(absv,'descend');
  Indice=min((NbIter+1)^2,N-1);
  Snew=sort(sorted_idx(1:Indice));
  Snew_c=sort(sorted_idx((Indice+1):N));
  mu_new=zeros(N,1);
  mu_new(Snew)=xnew(Snew)+(A(:,Snew)'*A(:,Snew))\A(:,Snew)'*(A(:,Snew_c)*xnew(Snew_c));
  NbIter=NbIter+1;
end
% x=mu_new;
