classdef PMSM < handle
    % Properties of the simulation model
    properties (SetAccess = private)
        % fidelity: this string may take three values 'low', 'med', 'high",
        % indicating the level of fidelity of the model.
        fidelity;  
        parameters; % System parameters
        SamlingTimes;
        SwitchingPeriod;
    end   
    
    
    methods (Access = public)
        function simulate(obj,
        
    end
end