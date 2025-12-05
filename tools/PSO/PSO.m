%Class of particles
classdef PSO < handle
    properties (SetAccess = private)
        N,np,bounds,g_best_position, g_best_objective,...
            acc_constant,particle_list,terminate_check,iter,iter_max ,init,list_pos_init;
    end
    methods (Access = public)
        function obj = PSO(init)
            obj.iter = 0;
            obj.iter_max = init.K_max;
            obj.N = init.N; % number of particle
            obj.np = size(init.bounds.A,2); % dimension of the particle (size)
            obj.acc_constant = init.acc_constant;
            obj.particle_list = cell(1,obj.N);
            obj.bounds = init.bounds;
            obj.init = init;
            
            
            if isfield(obj.init,'box')
                if (obj.init.box)&&(size(obj.bounds.A,2)==2)
                    disp('dim=2 and box bound');
                    ns = floor(sqrt(obj.N))+1;
                    rand_init_pts=zeros(2,ns^2);
                    i=0;
                    max_x = max(obj.bounds.V(:,1));
                    min_x = min(obj.bounds.V(:,1));
                    min_y = min(obj.bounds.V(:,2));
                    max_y = max(obj.bounds.V(:,2));
                    for xtmp = linspace(min_x,max_x,ns)
                        for ytmp = linspace(min_y,max_y,ns)
                            i=i+1;
                            rand_init_pts(:,i) =[xtmp; ytmp];
                        end
                    end
                else
                    rand_init_pts = cprnd(obj.N, obj.bounds.A,obj.bounds.b)';
                end
            else
                rand_init_pts = cprnd(obj.N, obj.bounds.A,obj.bounds.b)';
            end
            obj.list_pos_init = rand_init_pts;
%             size(rand_init_pts(:,i))
%             rand_init_pts
            obj.g_best_position = [];
            obj.g_best_objective = [];
            
            for i = 1:obj.N
                obj.particle_list{i} = particle(init.bounds, obj.np,...
                    rand_init_pts(:,i),obj.acc_constant,init.K_max,...
                    init.RNN,init.TrueFunc,init.epsilon,init);
                if i ==1
                    obj.g_best_position = obj.particle_list{i}.p_best_position;
                    obj.g_best_objective = obj.particle_list{i}.p_best_objective;
                else
                    if  obj.g_best_objective <= obj.particle_list{i}.p_best_objective
                        obj.g_best_position = obj.particle_list{i}.p_best_position;
                        obj.g_best_objective = obj.particle_list{i}.p_best_objective;
                    end
                    
                end
                
            end
        end
        
        function runPSO(obj)
            obj.iter = obj.iter+1;
            for i = 1:obj.N
                obj.particle_list{i}.update_velocity(obj.g_best_position);
                obj.particle_list{i}.update_position();
                obj.particle_list{i}.eval();
                if  obj.g_best_objective <= obj.particle_list{i}.p_best_objective
                    obj.g_best_position = obj.particle_list{i}.p_best_position;
                    obj.g_best_objective = obj.particle_list{i}.p_best_objective;
                end
            end
            
        end
    end
end

