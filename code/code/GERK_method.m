
% GERK method (Algorithm 1) and different instances

classdef GERK_method
    
    properties
        name;
        grad_fstar;
        L_fstar;
        grad_gstar;
        L_gstar;
    end
    
    methods
        function obj = GERK_method(name, grad_fstar, L_fstar, grad_gstar, L_gstar)
            if nargin == 5
                obj.name = name;
                obj.grad_fstar = grad_fstar;
                obj.L_fstar = L_fstar;
                obj.grad_gstar = grad_gstar;
                obj.L_gstar = L_gstar;
            elseif nargin == 3
                obj.name = name;
                obj.grad_fstar = grad_fstar;
                obj.L_fstar = L_fstar;
            else
                disp('Error in creating GERK method object: Need 5 parameters!')
            end
        end
    end
   
    
end






