function kappa =findLevelSet6D(Vx,beta,pmax,vmax,N,sys)

xGrid = linspace(-pmax, pmax, N);
yGrid=xGrid;
zGrid=xGrid;

vxGrid  = linspace(-vmax, vmax, N);
vyGrid  = vxGrid ;
vzGrid = vxGrid ;

kappMin = +10e6;
for x = xGrid
    for y = yGrid
        for z = zGrid
            for vx = vxGrid
                for vy = vyGrid
                    for vz = vzGrid
                        Xtmp=[x,y,z,vx,vy,vz]';
                        if (~checkCLF(Vx,sys.A,sys.B,beta,sys.U,Xtmp))&&(Vx.eval(Xtmp)<=kappMin)
%                             Vxtmp =Vx.eval(Xtmp);
%                             if Vxtmp <=  kappMin
                                kappMin = Vx.eval(Xtmp);
%                             end
                        end
                    end
                end
            end
        end
    end
end

kappa = kappMin;





end