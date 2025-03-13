% Aidan Hunt
% 
% The BCBase class defines basic functions to be used by other classes that
% extend BCBase to implement blockage corrections.
classdef BCBase < matlab.mixin.Heterogeneous
    properties (Constant, Access = protected)
        % The correctorNames property establishes the field names that are
        % expected to be used by a class that extends BCBase to implement
        % a blockage correction.
            % CP   - performance coefficient
            % CT   - thrust coefficient
            % CL   - lateral force coefficient
            % CF   - resultant force coefficient
            % CQ   - torque coefficient
            % TSR  - Tip-speed ratio
            % Uinf   - Undisturbed freestream velocity (m/s)
            % h   - Undisturbed freestream depth (m)
            % beta - Blockage ratio (as a fraction)
        correctorNames = {'CP', 'CT', 'CL', 'CF', 'CQ', 'TSR', 'Uinf', 'h', 'beta'};
    end

    % Abstract methods: these must be implemented in order for a class to
    % extend the BCBase class
    methods (Abstract)
        % All subclasses of BCBase must have a predictUnconfined method for
        % converting confined performance data to unconfined performance
        % data
        unconf = predictUnconfined(bc, conf)
    end

    methods (Static, Access = public)

        %% Performing blockage correction
        function unconf = convertConfToUnconf(conf, scalingVel)
            % Converts confined turbine data to unconfined turbine data
            % using the unconfined freestream velocity, UinfPrime
            % Inputs:
            % conf        - A structure of turbine data at a confined state.
            % scalingVel  - The characteristic velocity (e.g., unconfined
            %               freestream velocity predicted from a linear
            %               momentum) that will be used to scale the data
            %               in conf.
            %
            % Outputs:
            % unconf      - Performance data in conf with a blockage correction
            %               applied 

            fieldList = fieldnames(conf);
            unconf = struct;

            % For each field of the confined input data
            for i = 1:length(fieldList)
                currField = fieldList{i};
                switch currField % Scale based on which field it is
                    case {'CP', 'CP_b'} % Efficiency
                        correctedData = BCBase.scalePowerMetric(conf.(currField), conf.Uinf, scalingVel);
                    case {'CT', 'CL', 'CF', 'CQ'} % Force/torque
                        correctedData = BCBase.scaleForcingMetric(conf.(currField), conf.Uinf, scalingVel);
                    case 'TSR' % Tip-speed ratio
                        correctedData = BCBase.scaleTSR(conf.(currField), conf.Uinf, scalingVel);
                    case 'Uinf'  % Dimensional velocity
                        correctedData = BCBase.scaleVel(scalingVel);
                    case 'theta' % Azimuthal position
                        correctedData = conf.(currField);
                    otherwise % If field not found, don't copy over.
                        correctedData = [];
                end

                % If scaled data produced, copy to unconfined structure.
                if (~isempty(correctedData))
                    unconf.(currField) = correctedData;
                end
            end
    
            % Add velocity ratio as final field
            unconf.velRatio = conf.Uinf ./ scalingVel;
        end

        function powerUnconf = scalePowerMetric(powerConf, Uinf, scalingVel)
            % Scales a power metric (e.g., CP) using the cubed ratio of
            % Uinf/scalingVel and returns the result.
            powerUnconf = powerConf .* ((Uinf ./ scalingVel).^3);
        end


        function loadUnconf = scaleForcingMetric(loadConf, Uinf, scalingVel)
            % Scales a forcing metric (e.g., CT) using the squared ratio of
            % Uinf/scalingVel and returns the result.
            loadUnconf = loadConf .* ((Uinf ./ scalingVel).^2);
        end

        function TSRUnconf = scaleTSR(TSRConf, Uinf, scalingVel)
            % Scales the tip speed ratio using the ratio of Uinf/scalingVel 
            % (Ross and Polagye (2020) EQ 14) and returns the result.
            TSRUnconf = TSRConf .* (Uinf ./ scalingVel);
        end

        function velUnconf = scaleVel(scalingVel)
            % Returns the unconfined velocity (this is a trivial method that
            % may be overwritten by subclasses as needed).
            velUnconf = scalingVel;
        end

        %% Input checking
        function confData = checkInputSizes(confData)
            % Checks that sizes of inputs to blockage corrections on
            % compatible (i.e., the fields of confData are vectors are all
            % the same length). If scalar values of h and beta are
            % provided, converts to vectors of appropriate length.
            for i = 1:length(confData)
                checkFields = BCBase.correctorNames(1:end-3);
                allSizes = [];
                for j = 1:length(checkFields)
                    if isfield(confData(i), checkFields{j})
                        allSizes  = [allSizes; size(confData(i).(checkFields{j}))];
                    end
                end
                sizeCheck = size(unique(allSizes, 'rows'), 1) == 1;
                if (~sizeCheck)
                    error('Sizes of input Uinf, CP, CT, and TSR must match.');
                else
                    targetSize = size(confData(i).Uinf);
                end
                
                % Check size of beta and h
                checkFields = {'beta', 'h'};
                for j = 1:length(checkFields)
                    if length(confData(i).(checkFields{j})) == 1
                        confData(i).(checkFields{j}) = repmat(confData(i).(checkFields{j}), targetSize(1), targetSize(2));
                    elseif ~all(size(confData(i).(checkFields{j})) == targetSize)
                        error('Size of input "%s" must either a) match size of Uinf, CP, CT, TSR, or b) be 1.', checkFields{j});
                    end
                end
            end
        end

        %% Diagnostics
        function solInfo = packageDiagnostics(convX, convY, convErr, exitFlag)
            % Packages metadata from an fzero or fminsearch iteration into
            % structure.
            % Inputs
                % convX    - The converged input value(s) of the function
                % convY    - The converged output value(s) of the function.
                % convErr  - The residual of the fzero optimization
                % exitFlag - Information from fzero on why iteration ended
            % Outputs
                % solInfo - A structure whose fields are the inputs

            solInfo = struct;
            solInfo.convX = convX;
            solInfo.convY = convY;
            solInfo.convErr = convErr;
            solInfo.exitFlag = exitFlag;
        end

        %% Checking physical validity

        function result = checkPhysicalValidity(conf)
            % Checks whether velocities produced by linear momentum are
            % physically valid as described below. True indicates physical
            % validity.
            %
            % Inputs
            % conf - confined performance data with values of uw, ub,
            %        and ut included in structure
            % Outputs
            % result - A structure with the the following fields
                % ubutuw - Whether the bypass velocity is faster than the 
                %          turbine velocity and the turbine velocity is
                %          faster than the wake velocity (true/false)
                % ubUinf   - Whether the bypass velocity is faster than the
                %          freestream velocity (true/false)
                % Uinfut   - Whether the freestream velocity is faster than
                %          the turbine velocity (true/false).
            for i = 1:size(conf, 1)
                for j = 1:size(conf, 2)
                    result = struct;
                    
                    % Check that bypass is faster than turbine velocity,
                    % and turbine velocity is faster than wake velocity
                    result.ubutuw = (conf(i,j).ub > conf(i,j).ut) & (conf(i,j).ut > conf(i,j).uw);

                    % Check that bypass is faster than freestream
                    result.ubUinf = conf(i,j).ub > conf(i,j).Uinf;

                    % Check that freestream is faster than velocity at
                    % turbine
                    result.Uinfut = conf(i,j).Uinf > conf(i,j).ut;

                    if any(~result.ubutuw) || any(~result.ubUinf) || any(~result.Uinfut)
                        warning('Non physical result detected.');
                    end
                end
            end
        end

        %% Quantifying accuracy of forecasts

        function [forecastErr] = quantifyForecastError(confPredict, confReference, options)
            % Quantifies error in each point of a blockage forecast
            % Inputs:
            %    confPredict   - Structure of predicted performance data at
            %                    a particular confined condition
            %    confReference - Structure of actual performance data at 
            %                    the same confined condition
            % Outputs:
            %   forecastErr    - A structure of the raw difference between
            %                    the elements of confPredict and
            %                    confReference. Error is evaluated for the
            %                    following fields, if present in both
            %                    structures: CP, CT, CL, CF, CQ, Uinf.
            arguments
                confPredict
                confReference
                options.eachPercentage = false;
                options.peakPercentage = false;
                options.absError = false;
            end

            forecastErr = struct;
            predictMetrics = {'CP', 'CT', 'CL', 'CF', 'CQ', 'Uinf'}; % Meaningful metrics for quantifying error

            for j = 1:size(confPredict, 2) % For each forecast blockage
                % Get the current TSR
                TSRInterp = confReference(j).TSR;
                
                % Get index of maximum performance
                [~, maxInd] = max(confReference(j).CP);

                for i = 1:size(confPredict, 1) % For each dataset that was used to make a forecast

                    % Save TSR
                    forecastErr(i,j).TSR = TSRInterp;

                    for k = 1:length(predictMetrics)
                        if isfield(confPredict, predictMetrics{k}) && isfield(confReference, predictMetrics{k})
                            % Interpolate onto basis TSR
                            currInterp = interp1(confPredict(i,j).TSR, confPredict(i,j).(predictMetrics{k}), TSRInterp, 'linear');
                            
                            % Compute error
                            forecastErr(i,j).(predictMetrics{k}) = currInterp - confReference(j).(predictMetrics{k});

                            if options.eachPercentage % Normalize each point by the reference point
                                forecastErr(i,j).(predictMetrics{k}) = forecastErr(i,j).(predictMetrics{k}) ./ confReference(j).(predictMetrics{k}) .* 100;

                            elseif options.peakPercentage % Normalize each point by max efficiency of reference

                                forecastErr(i,j).(predictMetrics{k}) = forecastErr(i,j).(predictMetrics{k}) ./ confReference(j).(predictMetrics{k})(maxInd) .* 100;
                            end

                            % Take absolute value if requested
                            if options.absError
                                forecastErr(i,j).(predictMetrics{k}) = abs(forecastErr(i,j).(predictMetrics{k}));
                            end
                        end
                    end

                end
            end
        end
    end
end