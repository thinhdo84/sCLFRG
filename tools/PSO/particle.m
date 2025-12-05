% Class of particles
% Objective: maximize the fitness/objective
% A Generic Particle Swarm Optimization Matlab Function
% simplified version to speedup


classdef particle< handle
    properties (SetAccess = private)
        bounds, np, position, velocity, p_best_position, p_best_objective, ...
            objective, k;
        %         acc_constant;distance_to_global_best,
        phi_a, phi_b , K_max;
        RNN,TrueFunction;
        reInitFlag, reInitRadius, arrive2Global;
        para,init;
        currObj,currPosition,currVelo; %step k
    end
    methods (Access = public)
        function obj = particle(bounds, np, position_init,...
                acc_constant,K_max, TrueFunction,epsilon,init) % initialization
            obj.init = init;
            % constant
            obj.para.g =9.81;
            obj.para.m =0.25;
            obj.para.l =0.5;
            %
            obj.np = np; %the vector size
            obj.bounds = bounds; % linear wall for particles
            % particle's props
            % functions
            obj.TrueFunction = TrueFunction;
            obj.position = position_init;
            obj.velocity = position_init*0;
            obj.p_best_position = position_init;
            obj.objective =  obj.fitness(obj.position(:,end));
            obj.p_best_objective = obj.fitness(obj.p_best_position);
            obj.k = 0;
            
            % simplified version
            obj.currPosition = position_init;
            obj.currVelo = cprnd(1,bounds.A,bounds.b)'*0;
            obj.currObj = obj.fitness(obj.currPosition);
            
            % for reinitialization
            obj.arrive2Global=false;
            if isfield(init, 'reInitFlag')
                if init.reInitFlag
                    obj.reInitFlag=true;
                    obj.reInitRadius=init.reInitRadius;
                else
                    obj.reInitFlag=false;
                end
            end
            
            % dynamics
            obj.phi_a = 1-epsilon;
            obj.phi_b = mean(acc_constant)-1+epsilon;
            obj.K_max = K_max;
        end
        
        function eval(obj)
            if obj.p_best_objective < obj.fitness(obj.currPosition)
                obj.p_best_position = obj.currPosition;
                obj.p_best_objective = obj.fitness(obj.currPosition);
            end
        end
        
        function update_velocity(obj, global_best_position)
            obj.k = obj.k+1;
            gamma1 = rand(1,obj.np);
            gamma2 = rand(1,obj.np);
            % random for each dimenstion >>> diag
            v_cognitive = obj.acc_constant(1)*diag(gamma1)*(obj.p_best_position - obj.currPosition);
            v_social = obj.acc_constant(2)*diag(gamma2)*(global_best_position-obj.currPosition);
            phi_k = (obj.phi_b-obj.phi_a)*(obj.k-1)/obj.K_max + obj.phi_a;
            vel = phi_k * obj.currVelo + v_cognitive + v_social;
            
            % reinitialize the position and velocity
            if obj.reInitFlag
                if  norm(obj.currPosition-global_best_position)<=obj.reInitRadius
                    obj.arrive2Global = true;
                end
            end
            % update the velocity
            if obj.arrive2Global
                obj.currVelo = obj.currVelo*0.5;
            else
                obj.currVelo = vel;
            end
        end
        
        function update_position(obj)
            xtmp = obj.currPosition + obj.currVelo;
            if ~obj.bounds.contains(xtmp)
                xtmp = find_close_pts(obj.bounds.A,obj.bounds.b,...
                    obj.bounds.chebyCenter.x,xtmp);
            end
            if ~isnan(obj.fitness(xtmp))
                if obj.arrive2Global
                    r = rand;
                    ri=randi([1 size(obj.bounds.V,1)],1,1);
                    obj.currPosition=rand*(r*obj.position(:, 1) + (1-r)*obj.currPosition-obj.bounds.V(ri,:)')+...
                        obj.bounds.V(ri,:)';
                else
                    obj.currPosition = xtmp;
                end
            end
        end
        function y = fitness(obj, X)
            y = obj.TrueFunction(X);
        end
    end
end

function xsat= find_close_pts(A,b,c,x)
co= [A*(x-c);1;-1];
bo=[b;1;0];
max_lam = 100;
for i=1:numel(co)
    if co(i)>0
        max_lam =  min(bo(i)/co(i),max_lam);
    end
end
if max_lam>1+1e-5 || max_lam<0
    disp(max_lam)
    error( 'wrong scaling')
end
xsat=x*max_lam;
end