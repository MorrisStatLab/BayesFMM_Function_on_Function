function y=reverse(x)

%% y=reverse(x): reverses order of x, and multiply by (-1 1 -1 1 ..)
%%

n=length(x);
y=x;
for (i=1:floor(n/2))
    temp=x(i);
    y(i)=y(n-i+1);
    y(n-i+1)=temp;
end;

y=y.*repmat([1,-1],1,n/2);

