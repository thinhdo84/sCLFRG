function check = isRedundant(A,b,ai,bi)
x_max = linprog(-ai',A,b,...
    [],[],[],[],...
    optimoptions('linprog','Algorithm','interior-point','Display','off',...
    'MaxIterations',1000))';

if  numel(x_max)==0
    check = false;
else
    if ai*x_max' <= bi
        check = true;
    else
        check = false;
    end
end
end