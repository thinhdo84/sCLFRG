function y=Ellipsoid_plot2D(P,e,c,varargin)
%
%  Draw a patch region inside the following ellipsoid in 2D:
%  E =  {x : (x-c)'P(x-c)<=e  }
%
%  DO Huu-Thinh
%
th = linspace(0,2*pi,100);
xunit = cos(th);
yunit = sin(th);
options = optimoptions('fsolve','Display','off');
B = fsolve(@(X) X^-2-P,eye(2)*0.01,options);
a=sqrt(e);
X_unit=[xunit;yunit];
X=a*B*X_unit+c;

if nargin ==3
    y=plot(X(1,:), X(2,:),'k','linewidth',1.2);
else
    if nargin >3
        y=plot(X(1,:), X(2,:),varargin{:});
    end
end

end