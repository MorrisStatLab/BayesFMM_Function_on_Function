function H1A=get_H1A(n,h,extend)

%%%% Get H1A filter
%%%% Input: n=length of signal
%%%%        2*N=length of filter
%%%%        h= filter
%%%%        extend = if 1, keep all coeffs, if 0, truncate
%%%%

if nargin==2
    extend=1;
end;

N=length(h)/2;

ictr=1;
H1A=repmat(0,floor((n-1)/2)+N,n);
for (i=1:(N-1))
    H1A(ictr,(1:2*i))=h((2*(N-i)+1):(2*N));
    ictr=ictr+1;
end;
for (i=1:(floor(n/2)-N+1))
    H1A(ictr,(2*(i-1)+1):(2*(i-1)+2*N))=h;
    ictr=ictr+1;
end;
for (i=1:(N-1+mod(n-2*N,2)))
    H1A(ictr,(2*(floor(n/2)-N+i)+1):n)=h(1:2*(N-i)+mod(n-2*N,2));
    ictr=ictr+1;
end;

%%%% 
%%%%  Some arbitrary decisions must be made regarding which of the
%%%%  coefficients to throw away:
%%%%
%%%% 1. Throw away first floor((N-1)/2) rows and last ceil((N-1)/2) rows.
%%%%    Note: When N is even, then odd number of rows must be thrown away.
%%%%    Here we throw away an extra row at the end.
%%%% 2. If scaling coefficent (abs(sum(h))>1e-6) and N even, then keep extra
%%%%     row at the end.
%%%%
eps=1e-6;
if extend==0
    H1A=H1A((floor((N-1)/2)+1):(floor((N-1)/2)+floor(n/2)+mod(n,2)*(abs(sum(h)>eps))),:);
end;
