function drawHyperPlane2D(a,b,X, varargin)

if nargin ==3
    if a(2)==0
        if a(1)==0
            error('Cannot draw the hyperplane')
        else
            xline(b/a(1),'color','k','linestyle','-.','linewidth',2);
        end
    else
        plot(X, [(b-a(1)*X(1))/a(2), (b-a(1)*X(2))/a(2)], 'k-.','linewidth',2);
    end
else
    
    if a(2)==0
        if a(1)==0
            error('Cannot draw the hyperplane')
        else
            xline(b/a(1),varargin{:});
        end
    else
        plot(X, [(b-a(1)*X(1))/a(2), (b-a(1)*X(2))/a(2)], varargin{:});
    end
    
    
    
    
end