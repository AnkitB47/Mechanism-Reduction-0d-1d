function plotTemp(obj,optional)
% plots the temperarture over time
switch nargin
    case 1
        plot(obj.times,obj.temp);
    case 2
        try
            plot(obj.times,obj.temp,optional);
        catch ME
            fprintf(['The argument of plotTemp may be not correct,',...
                ' check the result...\n The excemtions was:\n ',ME.message,'\n'])            
            plot(obj.times,obj.temp);
        end
    otherwise
end
xlabel('time [s]','FontSize',14);
ylabel('temperature [K]','FontSize',14);
end