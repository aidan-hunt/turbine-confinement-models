% Aidan Hunt
% 
% The BCBase class defines basic functions to be used by other classes that
% extend BCBase to implement blockage corrections.
classdef BCBase < matlab.mixin.Heterogeneous
    properties (Constant, Access = public)
        % The correctorNames property establishes the field names that are
        % expected to be used by a class that extends BCBase to implement
        % a blockage correction.
            % CP   - performance coefficient
            % CT   - thrust coefficient
            % CL   - lateral force coefficient
            % CF   - resultant force coefficient
            % CQ   - torque coefficient
            % TSR  - Tip-speed ratio
            % V0   - Undisturbed freestream velocity (m/s)
            % d0   - Undisturbed freestream depth (m)
            % beta - Blockage ratio (as a fraction)
        correctorNames = {'CP', 'CT', 'CL', 'CF', 'CQ', 'TSR', 'V0', 'd0', 'beta'};
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
            % using the unconfined freestream velocity, V0Prime
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
                    case 'CP' % Efficiency
                        correctedData = BCBase.scalePowerMetric(conf.(currField), conf.V0, scalingVel);
                    case {'CT', 'CL', 'CF', 'CQ'} % Force/torque
                        correctedData = BCBase.scaleForcingMetric(conf.(currField), conf.V0, scalingVel);
                    case 'TSR' % Tip-speed ratio
                        correctedData = BCBase.scaleTSR(conf.(currField), conf.V0, scalingVel);
                    case 'V0'  % Dimensional velocity
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
            unconf.velRatio = conf.V0 ./ scalingVel;
        end

        function powerUnconf = scalePowerMetric(powerConf, V0, scalingVel)
            % Scales a power metric (e.g., CP) using the cubed ratio of
            % V0/V0Prime (Ross and Polagye (2020) EQ 12) and returns the
            % result.
            powerUnconf = powerConf .* ((V0 ./ scalingVel).^3);
        end


        function loadUnconf = scaleForcingMetric(loadConf, V0, scalingVel)
            % Scales a forcing metric (e.g., CT) using the squared ratio of
            % V0/V0Prime (Ross and Polagye (2020) EQ 13) and returns the
            % result.
            loadUnconf = loadConf .* ((V0 ./ scalingVel).^2);
        end

        function TSRUnconf = scaleTSR(TSRConf, V0, scalingVel)
            % Scales the tip speed ratio using the ratio of V0/V0Prime 
            % (Ross and Polagye (2020) EQ 14) and returns the result.
            TSRUnconf = TSRConf .* (V0 ./ scalingVel);
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
            % the same length). If scalar values of d0 and beta are
            % provided, converts to vectors of appropriate length.
            for i = 1:length(confData)
                checkFields = BCBase.correctorNames(1:end-2);
                allSizes = [];
                for j = 1:length(checkFields)
                    if isfield(confData(i), checkFields{j})
                        allSizes  = [allSizes; size(confData(i).(checkFields{j}))];
                    end
                end
                sizeCheck = size(unique(allSizes, 'rows'), 1) == 1;
                if (~sizeCheck)
                    error('Sizes of input V0, CP, CT, and TSR must match.');
                else
                    targetSize = size(confData(i).V0);
                end
                
                % Check size of beta and d0
                checkFields = BCBase.correctorNames(end-2:end);
                for j = 1:length(checkFields)
                    if length(confData(i).(checkFields{j})) == 1
                        confData(i).(checkFields{j}) = repmat(confData(i).(checkFields{j}), targetSize(1), targetSize(2));
                    elseif ~all(size(confData(i).(checkFields{j})) == targetSize)
                        error('Size of input "%s" must either a) match size of V0, CP, CT, TSR, or b) be 1.', checkFields{j});
                    end
                end
            end
        end

        %% Diagnostics
        function solInfo = packageDiagnostics(convX, convY, convErr, exitFlag)
            % Packages metadata from an fzero iteration into structure.
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

        function result = checkPhysicalValidity(conf)
            % Checks whether velocities produced by linear momentum are
            % physically valid as described below. True indicates physical
            % validity.
            %
            % Inputs
            % conf - confined performance data with values of u1, u2,
            %        and ut included in structure
            % Outputs
            % result - A structure with the the following fields
                % u2utu1 - Whether the bypass velocity is faster than the 
                %          turbine velocity and the turbine velocity is
                %          faster than the wake velocity (true/false)
                % u2V0   - Whether the bypass velocity is faster than the
                %          freestream velocity (true/false)
                % V0ut   - Whether the freestream velocity is faster than
                %          the turbine velocity (true/false).
            for i = 1:size(conf, 1)
                for j = 1:size(conf, 2)
                    result = struct;
                    
                    % Check that bypass is faster than turbine velocity,
                    % and turbine velocity is faster than wake velocity
                    result.u2utu1 = (conf(i,j).u2 > conf(i,j).ut) & (conf(i,j).ut > conf(i,j).u1);

                    % Check that bypass is faster than freestream
                    result.u2V0 = conf(i,j).u2 > conf(i,j).V0;

                    % Check that freestream is faster than velocity at
                    % turbine
                    result.V0ut = conf(i,j).V0 > conf(i,j).ut;
                end
            end
        end
    end
end