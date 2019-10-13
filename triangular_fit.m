%fit triangular pdf function
function p = triangular_fit(X,lam)
if X(1)>X(2)
    p = 1000000;
elseif X(1)>X(3)
    p = 1000000;
elseif X(3)>X(2)
    p = 1000000;
elseif X(1)<0
    p = 1000000;
else
pd3 = makedist('Triangular','a',X(1),'b',X(3),'c',X(2));
cdf3 = cdf(pd3,lam);
p = cdf3;
end
end