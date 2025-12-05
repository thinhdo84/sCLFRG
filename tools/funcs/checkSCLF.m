function bool = checkSCLF(V,x,u,g,sys)
bool = V.grad(x) * (sys.A*x + sys.B*u) +g <= 1e-8;
end