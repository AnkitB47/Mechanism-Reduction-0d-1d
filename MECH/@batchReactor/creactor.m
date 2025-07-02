function dydt = creactor(obj,t,y)
% REACTOR ODE system for a generic zero-dimensional reactor.
%
%    Function REACTOR evaluates the system of ordinary differential
%    equations for a zero-dimensional reactor with arbitrary heat
%    transfer and  volume change. 
%
% Solution vector components:
%    y(1)   Total internal energy U
%    y(2)   Volume V
%    y(3)   Mass of species 1
%    ....
%    y(2+nsp) Mass of last species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
%                  *  IMPORTANT NOTE *
%          Here we claim that Y is given already in
%                 (NORMED) MASSFRACTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
%  y = y(:,j);
  int_energy =y(1);
  vol =y(2);
  masses =y(3:end);
  
  % set the state of the obj.gas by specifying (u,v,{Y_k})
  setMassFractions(obj.gas,masses);
  setState_UV(obj.gas, [int_energy vol]);
  p = pressure(obj.gas);

  % volume equation
  vdt = feval(@obj.vdot, t, vol, obj.gas);

  % energy equation
  a = feval(@obj.area, t, vol);
  q = feval(@obj.heatflux, t, obj.gas);  
  udt = -p * vdt + a * q;

  % species equations
  ydt =  ydot(obj.gas);

  % set up column vector for dydt
  dydt = [ udt
		vdt
		ydt ];
end
