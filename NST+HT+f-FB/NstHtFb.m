function [x,NbIter] =NstHtFb(y,y0,A,phi,s,MaxNbIter,TolRes)
% AdptNST+HT+f-FB with f(x)=s; i.e, NST+HT+FB
% Input: 
%  observed data: y,   pre-computed vector: y0=A'*((A*A')\y);
%  measurement matrix: A,  pre-computed matrix: phi=A'*((A*A') \A)£»
%  spareity: s, maximum number of iteration: MaxNbIter, tolerance on the relative error: epsilon.
% Output:
% recovered sparse signal: x
N=size(A,2);
%% Initialization
% mu=ones(N,1);
mu_new=zeros(N,1);
NbIter=0;
[~, idx] = sort(abs(mu_new));
S = idx(1:s);
Snew = idx(end-s+1:end); % Assign something completely different to make sure we get at least once in the for loop
%% main loop
while norm( y-A*mu_new) > TolRes && (sum(S==Snew)<s) &&  NbIter < MaxNbIter
    xnew=mu_new+y0-phi*mu_new;
    absv=abs(xnew);
%     zero_idx=find(absv<Eps*max(absv));
%     absv(zero_idx)=zeros(size(zero_idx));
    [~,sorted_idx]=sort(absv,'descend');
    Snew=sort(sorted_idx(1:s)); 
    Snew_c=sort(sorted_idx(s+1:N)); 
    mu_new=zeros(N,1);
    mu_new(Snew)=xnew(Snew)+(A(:,Snew)'*A(:,Snew))\A(:,Snew)'*(A(:,Snew_c)*xnew(Snew_c));
	NbIter=NbIter+1;
    aux = S;
    S = Snew;
    Snew = aux;
end
x=mu_new;
end
