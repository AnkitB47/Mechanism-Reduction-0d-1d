function plotPress(obj,optional)
% plots the temperarture over time
switch nargin
    case 1
        plot(obj.times,obj.pres/oneatm);
    case 2
        try
            plot(obj.times,obj.pres/oneatm,optional);
        catch ME
            fprintf(['The argument of plotTemp may be not correct,',...
                ' check the result...\n The excemtions was:\n ',ME.message,'\n'])            
            try
                plot(obj.times,obj.pres/oneatm);
            catch ME
                throw ME
            end
        end
    otherwise
        error('batchReactor:plotPress:tooManyArguments',...
            'Too many arguments. Only one plot argument is allowed')
end
xlabel('t [s]','FontSize',16);
ylabel('p [Atm]','FontSize',16);
end
