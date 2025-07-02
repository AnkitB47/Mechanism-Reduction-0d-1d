function obj = burnerFlame(obj)


% parameter values


if isempty(obj.space)
    initial_grid = [0.0 0.02 0.04 0.06 0.08 0.1 ...
		0.15 0.2 0.4 0.49 0.5];  % m
else
    initial_grid = obj.space;
end

tol_ss    = [1.0e-5 1.0e-13];       % [rtol atol] for steady-state
                                    % problem
tol_ts    = [1.0e-4 1.0e-9];        % [rtol atol] for time stepping

loglevel  = 1;                      % amount of diagnostic output (0
                                    % to 5)
				    
refine_grid = 1;                    % 1 to enable refinement, 0 to
                                    % disable 				   
max_jacobian_age = [5, 10];        
				   
%%%%%%%%%%%%%%%% create the gas object %%%%%%%%%%%%%%%%%%%%%%%%
%
% This object will be used to evaluate all thermodynamic, kinetic,
% and transport properties
%
gas = obj.gas;

% set its state to that of the unburned gas at the burner
set(gas,'T', obj.initialTemperature, 'P',obj.initialPressure, 'X',obj.compLeft);



%%%%%%%%%%%%%%%% create the flow object %%%%%%%%%%%%%%%%%%%%%%%

f = AxisymmetricFlow(gas,'flow');
set(f, 'P',obj.initialPressure, 'grid', initial_grid);
set(f, 'tol', tol_ss, 'tol-time', tol_ts);



%%%%%%%%%%%%%%% create the burner %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  The burner is an Inlet object. The temperature, mass flux, 
%  and composition (relative molar) may be specified.
%
burner = Inlet('burner');
set(burner, 'T',obj.initialTemperature, 'MassFlux',obj.mdotLeft, 'X', obj.compLeft);



%%%%%%%%%%%%%% create the outlet %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  The type of flame is determined by the object that terminates
%  the domain. An Outlet object imposes zero gradient boundary
%  conditions for the temperature and mass fractions, and zero
%  radial velocity and radial pressure gradient.
%
s = Outlet('out');


%%%%%%%%%%%%% create the flame object  %%%%%%%%%%%%
%
% Once the component parts have been created, they can be assembled
% to create the flame object.
%
fl = flame(gas, burner, f, s);
setMaxJacAge(fl, max_jacobian_age(1),  max_jacobian_age(2));

% if the starting solution is to be read from a previously-saved
% solution, uncomment this line and edit the file name and solution id.

solve(fl, loglevel, refine_grid);

%%%%%%%%%%%% enable the energy equation %%%%%%%%%%%%%%%%%%%%%
%
%  The energy equation will now be solved to compute the
%  temperature profile. We also tighten the grid refinement
%  criteria to get an accurate final solution.
%

enableEnergy(f);
setRefineCriteria(fl, 2, 200.0, 0.05, 0.1);
solve(fl, 1, 1);
obj = createSolution(obj,fl,f);


function obj = createSolution(obj,fl,f)
        % we should first clear the old values, because we indexing the new ones by
        % a "row set call"
        mm = molarMasses(obj.gas);
        obj.temp = solution(fl,'flow','T');
        obj.space = z(f);
        sn = speciesNames(obj.gas);
        
        obj = obj.clearFractions;
        massFrac = zeros(nSpecies(obj.gas) ,length(solution(fl,'flow',sn{1})));
        for kk = 1:nSpecies(obj.gas)   
            massFrac(kk,:) = solution(fl,'flow',sn{kk});    
        end
        D = sparse(1:length(mm),1:length(mm),1./mm); 
        M = 1./(1./mm'*massFrac);
        obj = obj.setMassFractions(massFrac);
        obj = obj.setMoleFractions((D*massFrac)*sparse(1:length(M),1:length(M),M));
        obj.u = solution(fl,'flow','u');
        obj.v = solution(fl,'flow','V').*z(f);
end


end

