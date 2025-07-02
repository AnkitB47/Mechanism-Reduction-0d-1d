function varargout = speciesProductionRates(obj,varargin)
% Computes the production rates of all or of given
% species from data at time point given by the vector
% times.
%
% rates = obj.speciesProductionRates
% rates = obj.speciesProductionRates(times)
% rates = obj.speciesProductionRates(times,speciesList)
%
% (c) 2011 U.P. for VIRTUHCON

timeIndex = 1:length(obj.times);
speciesIndex = 1:obj.nSpecies;

switch nargin   
    case 1
        % use standard values
    case 2
        timeIndex =  obj.getTimesIndex(varargin{1});       
    otherwise  % species list given  
        if ~isempty(varargin{1})
            timeIndex =  obj.getTimesIndex(varargin{1});       
        end
        speciesIndex = obj.getSpeciesIndexByName(varargin{2});
        % ignore the rest
end

rates = zeros(length(speciesIndex),length(timeIndex));

j = 1;
for k =  timeIndex
    moleFractions = vector2CanteraString(obj.speciesNames,...
        obj.moleFractions(:,k));
    set(obj.gas,'X',moleFractions,'T',obj.temp(k),'P',obj.pres(k));
    r = netProdRates(obj.gas); 
    rates(:,j) = r(speciesIndex);
    j = j+1;
end

% output
switch nargout
    case 0 % graphical output
        figure
        ratesIndex = 1:length(rates(:,1));
        if nargin<3 % plot only lelevant species rates
            ratesIndex = ratesIndex(( max(abs(rates'))>1e-8));          
        end
        bar3(rates(ratesIndex,:)) 
        if nargin<3 % the annotation must be adaptet case by case...
            annotation('textbox',[0,0.0,0,1],'String',...
                ['Species: ' obj.speciesNames(ratesIndex)],'FontSize',12');
        else
            annotation('textbox',[0,0.0,0,1],'String',...
                ['Species: ' ,varargin{2}],'FontSize',12');
        end
        xlabel('time t');
        ylabel('Species Number');
        title('Net production rates  over time','FontSize',14);
    case 1
        varargout{1} = rates;
    otherwise
        %
        error('batchReactor:speciesProductionRates:tooManyOutputs','The number of output must be zero or one.');
end
end