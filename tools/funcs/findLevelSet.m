function contour_level = findLevelSet(Vx,surface_data,beta,sys)
ws_pts=[surface_data.X;surface_data.Y];
V_grid = [];
for x1 =surface_data.X
    for x2 = surface_data.Y
         V_grid = [V_grid Vx.eval([x1;x2])];
    end
    
end
V_grid = reshape(V_grid,[numel(surface_data.X), numel(surface_data.Y)])';
c = 0.001;
% c=0.001;
c_pre = c;
% delta_c = 100e-4;
delta_c=0.02;
in_ws_flag = true;
contour_pts = contourc(surface_data.X,surface_data.Y,V_grid,[1 1]*c);
% contour_pts
contour_pts(:,1) = [];
for j=1:size(contour_pts,2)
    if any(surface_data.ws.A*contour_pts(:,j)>surface_data.ws.b)
        in_ws_flag = false;
        error('initial level set is too large')
    end
end
finish_flag = false;
while in_ws_flag && ~finish_flag
    level_pts = findPoints_inlevel(c, Vx, ws_pts);
    contour_pts = contourc(surface_data.X,surface_data.Y,V_grid,[1 1]*c);
    contour_pts(:,1) = [];
    for j=1:size(contour_pts,2)
        if any(surface_data.ws.A*contour_pts(:,j)>surface_data.ws.b)
            in_ws_flag = false;
        end
    end
    
    for i_cnt =1:size(level_pts,2)
        if ~checkCLF(Vx,sys.A,sys.B,beta,sys.U, level_pts(:,i_cnt) )
            finish_flag = true;
        end
    end
    c_pre=c;
    c=delta_c+c;
    contour_level = c_pre-delta_c;
end
end


function pts = findPoints_inlevel(level, V, check_pts)
n_pts = size(check_pts,2);
pts=[];
for i =1:n_pts
    if V.eval(check_pts(:,i)) <= level
        pts = [pts, check_pts(:,i)] ;
    end
end
end