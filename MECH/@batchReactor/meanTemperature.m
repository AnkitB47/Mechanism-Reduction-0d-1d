function varargout = meanTemperature(obj)
            % t_mean =  obj.meanTemperature 
            % computes the mean temperature of a flame1d object
            % by t_mean = \int_x_0^x_1 T dx / \int_x_0^x_1 1 dx.
            if isempty(obj.temp) || isempty(obj.times)
                error('batchReactor:meanTemperature:emptyTemp',...
                    'Before computing the location of the flame you must compute the temperature.')
            else
                h = obj.times(2:end)-obj.times(1:end-1);               
                mT = 0.5/(max(obj.times)-min(obj.times))*h*(obj.temp(1:end-1)+obj.temp(2:end))';
                switch nargout
                    case 0 
                        fprintf(['\t mean temperature = ',num2str(mT),' K \n']);
                    case 1                        
                        varargout{1} = mT;                    
                    otherwise
                        error('batchReactor:meanTemperature:tooManyOutputs',...
                            'Too many output arguments, try >>help meanTemperature.')
                end                       
            end 
        end