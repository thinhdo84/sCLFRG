function u= explicit_mP(sol, x)
[~,j] = isinside(sol{1}.Pn,x);
u=sol{1}.Fi{j}*x + sol{1}.Gi{j};
end