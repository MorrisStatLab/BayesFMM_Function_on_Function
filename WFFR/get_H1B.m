function H1B=get_H1B(n,h,boundary,extend)

%%%% Get H1B filter
%%%% Input: n=length of signal
%%%%        2*N=length of filter
%%%%        h= filter
%%%%        boundary = type of boundary condition
%%%%        extend = 1 means keep all floor((n-1)/2+N) coefficients, like
%%%%                 Matlab
%%%% Assumes reflection boundary conditions -- can be adapted to include
%%%% others.

if nargin==2
    boundary='reflection';
end;

if nargin==3
    extend=1;
end;

if boundary=='zeros     '
    H1B=0;
    G1B=0;
else
N=length(h)/2;

ictr=1;
H1B=repmat(0,floor((n-1)/2)+N,n);
for (i=1:(N-1))
    if boundary=='reflection'
        H1B(ictr,(1:(2*(N-i))))=h((2*(N-i)):-1:1);
    else
        H1B(ictr,(n-2*(N-i)+1):n)=h(1:2*(N-i));
    end;
    ictr=ictr+1;
end;
for (i=1:(floor(n/2)-N+1))
    ictr=ictr+1;
end;
for (i=1:(N-1+mod(n-2*N,2)))
    if boundary=='reflection'
        H1B(ictr,(((n-2*(i-1)-1)+mod(n-2*N,2)):(n)))=h((2*N):-1:(2*(N-i)+1+mod(n-2*N,2)));
    else
        H1B(ictr,(1:(2*i-mod(n-2*N,2))))=h((2*(N-i)+1+mod(n-2*N,2)):end);
    end;
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
    H1B=H1B((floor((N-1)/2)+1):(floor((N-1)/2)+floor(n/2)+mod(n,2)*(abs(sum(h)>eps))),:);
end;
end;
