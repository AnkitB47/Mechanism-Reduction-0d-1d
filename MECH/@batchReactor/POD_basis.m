function basis = POD_basis(obj,varargin)
% obj.POD_basis
% computes the POD basis wrt. the data in massFractions
switch nargin
    case 1        
        % no extra arguments
        nOfSnapshots = 30;         
    case 2
        % nOfSnapshots
        nOfSnapshots = round(varargin{1});
%     case
    otherwise
        % 
        ME = MException('batchReactor:POD_basis:TooManyArguments',...
            'Too many arguments, cf. help POD_basis:');
        help POD_basis;
        throw(ME);
end
 steps = round(length(obj.times)/nOfSnapshots);
Y = [obj.temp(:,1:steps:end);...
            obj.massFractions(:,1:steps:end)];
 
R = Y'*Y;
[EV,EW] = eig(R,'nobalance');
EW = diag(EW);
indx = EW>1e-3;
n = length(EW);

 
% basis = EV;

for k= 1:nOfSnapshots
     basis(:,k) = 1/sqrt(EW(k))*Y*EV(:,k);
end
 

end