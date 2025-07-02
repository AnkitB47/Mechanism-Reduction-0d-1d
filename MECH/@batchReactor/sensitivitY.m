function [out1,out2] = sensitivitY(obj,y)
% dfdy = obj.sensitivitY(y) computes the Jacobian of the rhs at y for the 
% RHS of the batchRector class with respect to Y.
% Usuage: Jac = obj.sensitivitY(Z])  Z must befit the requrenments of
% obj.reactor , eg. [E,Vol,Y] for common reactor, [T,Y] for all other
% reactors.
% It uses Matlab's numJacODE function, and the set of arguments used by
% conhp, conuv etc.
% It checks the reactor propertie of the batchReactor object to deside wich
% RHS is valid. Possible Values are 'constPressReactor',  'constPressNetReactor',
% 'constVolumeReactor', and  'constVolumeNetReactor'.
% 'commonReactor' and 'tibox' are NOT (jet) valid.
%
% (c) 2011 U.P. for VIRTUHCON

options = odeset('RelTol',1.e-5,'AbsTol',1.e-12,'Stats','off');
switch obj.reactor
    case {'constPressReactor'  'constPressNetReactor'}
        dfdy = numJacODE(@obj.conhp,0,y,options);
        f = obj.conhp(0,y);
        dfdy(:,1) = [];  dfdy(1,:) = [];
        f(1) = [];
        
    case {'constVolumeReactor'  'constVolumeNetReactor' }
         dfdy = numJacODE(@obj.conuv,0,y,options);
         f = obj.conuv(0,y);
        dfdy(:,1) = [];  dfdy(1,:) = [];
        f(1) = [];
    case 'constTempReactor'
        % remove temperature from solution vector, this is for
        % compatibility reasons...
        y = y(2:end);
        dfdy = numJacODE(@obj.conht,0,y,options);
        f = obj.conht(0,y);    
    case 'commonReactor'
        % here y must be [E,Vol,Y] !!!!
        setMassFractions(obj.gas,y(2:end));
        set(obj.gas,'T',y(1));
        y(1) = intEnergy_mass(obj.gas);
        y(3:end+1) = y(2:end);
        y(2) = 1./density(obj.gas);       
        
        dfdy = numJacODE(@obj.creactor,0,y,options);
        f = obj.creactor(0,y);     
        dfdy(:,1:2) = [];  dfdy(1:2,:) = [];
        f(1:2) = [];
    otherwise
        error('batchReactor:jacobi:unknownReactor',...
            'The type of the reactor is unknown.')
end
switch nargout
    case  0
      % do nothing
    case 1
        out1 = dfdy;
    case 2
         out1 = dfdy;
         out2 =  f; 
    otherwise
        % matlab error?

end


end
