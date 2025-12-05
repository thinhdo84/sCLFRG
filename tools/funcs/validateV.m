function dV = validateV(V,f,g,omega,ell,x)
u = omega(x);
dV = V.grad(x)*(f(x) + g(x)*u) + ell(x,u);
end