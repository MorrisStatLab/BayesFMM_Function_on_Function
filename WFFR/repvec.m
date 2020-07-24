function y=repvec(a,n)

%%% repvec(a,n) creates a vector consisting of repeating each element of a
%%% (a_i) n_i times.  It behaves like Splus function rep(a,n) when a and n 
%%% are both vectors of the same length.
%%%
%%% Input: row vectors a and n, each of length k
%%% Output: row vector y, of length sum(n).

k=max(size(a));
y=repmat(a(1),1,n(1));
if k>1
    for i=2:k
        y=[y,repmat(a(i),1,n(i))];
    end % for loop
end % if
    
    
