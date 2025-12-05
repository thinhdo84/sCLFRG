function varargout = checkInvariance(V,center, levelset, fx, gx, U, NSample)
% Check the invariance condition inside the level set:
% {x: V(x) <= alpha}
% with the vector field dx = f(x) + g(u)u;
% input constraint set U (polytopic)

% We first define the outerbox of the set, then sample randomly inside that
% box with Nsamples points
tic;
disp('Computing the bounding box of the level set...')
[ub, lb] = bounding_box_levelset(V, center, levelset);
SampleBox = Polyhedron('ub', ub, 'lb', lb);
Samples = cprnd(NSample,SampleBox.A,SampleBox.b)'; % NSample vector stacked side-by-side
fprintf('%d random points sampled...\n', NSample)

falseCnt = 0;
interiorCnt = 0;
for i = 1:NSample
    X = Samples(:,i);
    if V.eval(X)<=levelset
        interiorCnt = interiorCnt +1;
        checkResult = InvariancePoint(V,fx(X),gx(X),U,X);
        if ~checkResult
            falseCnt = falseCnt +1;
        end
    end
end
fprintf('%d points inside the level set...\n', interiorCnt)

if falseCnt > 0
warning('Inadmissible point detected')
checkResult = false;
end
successRate = falseCnt/NSample;
if nargout == 1
    varargout{1} = checkResult;
elseif nargout == 2
    varargout{1} = checkResult;
    varargout{2} = successRate;
else
    error('Function expects 1 or 2 outputs.');
end
fprintf('Computation time %.3f(s)\n', toc)
end