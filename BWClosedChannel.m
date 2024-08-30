% Aidan Hunt
%
% The BWClosedChannel class implements closed-channel linear momentum
% theory and an associated blockage correction following Barnsley and
% Wellicome (1990) as described by Ross and Polagye in "An experimental
% assessment of analytical blockage corrections" (2020;
% https://doi.org/10.1016/j.renene.2020.01.135).
%
% Linear momentum on an actuator disk in closed channel flow (parallel
% sided-tube) is employed to solve for the following flow properties:
%   uw      - Core wake velocity
%   ub      - Bypass velocity
%   ut      - Velocity at the turbine
%   UinfPrime - Unconfined freestream velocity
% These quantities can then be used for blockage correction or blockage
% forecasting.
%
% To use the BWClosedChannel class, construct a BWClosedChannel object
% using the following syntax:
%   bw = BWClosedChannel()
% and call methods using the dot notation (i.e., bw.solveLMAD(...))
% 
% BWClosedChannel Methods:
%   solveLMAD         - Solves for flow velocities in the vicinity of the rotor
%                       using closed-channel linear momentum on an actuator
%                       disk (LMAD) theory and confined performance data.
%                       No blockage correction is performed.
%   predictUnconfined - Uses closed-channel LMAD theory to predict unconfined 
%                       performance from confined performance data.
%                       "Standard" blockage corrections or "bluff body"
%                       blockage corrections may be used.
%   forecastConfined  - Uses closed-channel LMAD theory and a bluff body
%                       blockage correction basis to predict performance at
%                       one blockage ratio using data at a different
%                       blockage ratio.
% 
% The methods above expect that confined performance data is provided as an
% mxn structure array, conf, with the following fields:
%   beta (required)   - blockage ratio
%   Uinf   (required)   - undisturbed upstream freestream velocity (m/s)
%   CT   (required)   - thrust coefficient
%   CP   (optional)   - performance coefficient
%   CQ   (optional)   - torque coefficient
%   CL   (optional)   - lateral force coefficient
%   CF   (optional)   - resultant force coefficient
%   TSR  (optional)   - tip-speed ratio
% The fields of conf(i,j) must be vectors that are all the same size. You
% can use the methods above to apply corrections and forecasts to multiple
% datasets at once by specifying each dataset as an element of conf (e.g.,
% conf(i,j)).
%
% These methods utilize helper methods that implement the core equation set
% for the Barnsley and Wellicome correction as presented in Ross and
% Polagye (2020), which are also available for use. The BWClosedChannel
% class extends the BCBase class.
%
% See also: BCBase, HoulsbyOpenChannel
classdef BWClosedChannel < BCBase

    % The BWClosedChannel class inherits the properties of the BCBase
    % class, and defines additional properties
    properties (Constant, Access=protected)
        correctionModes = {'standard', 'bluff body'} % Available blockage correction modes
        expectedVel = {'UinfPrime', 'ub'}; % Default scaling velocity for each blockage correction mode
    end

    methods (Access = public)
        %% Main method for solving closed-channel linear momentum model
        function [conf] = solveLMAD(bw, conf, ubuwGuess)
            % Solves linear momentum on an actuator disk using a
            % closed-channel model. Velocities are estimated using
            % Equations 20-23 from Ross and Polagye (2020).
            %
            % Inputs
            %   conf      - A structure of confined performance data with fields
            %               as described in the BWClosedChannel class documentation.
            %   ubuwGuess - Initial guess for ub/uw for use in iteration
            %               (default: 1.4).
            % Outputs
            %   conf      - The input structure with the following fields
            %               added:
            %               uw       - Wake velocity estimated from closed-channel LMAD
            %               ub       - Bypass velocity estimated from closed-channel LMAD
            %               ut       - Velocity at the turbine estimated from closed-channel LMAD
            %               UinfPrime  - Unconfined freestream velocity estimated from closed-channel LMAD
            %               ubuwIter - Iteration diagnostics for ub/uw
            %               isPhys   - Results of physical validity checks on uw, ub, ut
            %
            % See also: BWClosedChannel, predictUnconfined, forecastConfined, checkPhysicalValidity

            % Define function input arguments, validation, and default
            % values.
            arguments
                bw
                conf
                ubuwGuess (1,:) = [1.4]
            end

            % Check input for correct sizing
            conf = bw.checkInputSizes(conf);

            for i = 1:size(conf, 1)
                for j = 1:size(conf, 2)

                    % Solve for ubuw via iteration
                    [ubuw, ubuwErr, ubuwExitFlag] = bw.convergeUbUw(ubuwGuess, conf(i,j));
    
                    % Compute channel velocities
                    conf(i,j).uw = bw.solveUw(ubuw, conf(i,j).CT, conf(i,j).Uinf);
                    conf(i,j).ub = bw.solveUb(ubuw, conf(i,j).CT, conf(i,j).Uinf);
                    conf(i,j).ut = bw.solveUt(ubuw, conf(i,j).beta, conf(i,j).CT, conf(i,j).Uinf);
                    conf(i,j).UinfPrime = bw.solveUinfPrime(conf(i,j).Uinf, conf(i,j).CT, conf(i,j).ut);
    
                    % Package diagnotics about ubuw iteration
                    [UinfUw_blockage, UinfUw_thrust] = bw.solveUinfUw_both(ubuw, conf(i,j).beta, conf(i,j).CT);
                    conf(i,j).ubuwIter = bw.packageDiagnostics(ubuw, [UinfUw_blockage UinfUw_thrust], ...
                                                                 ubuwErr, ubuwExitFlag);
    
                    % Check physical validity
                    conf(i,j).isPhys = bw.checkPhysicalValidity(conf(i,j));
                end
            end
        end

        %% Main method for performing blockage correction
        function [unconf, conf] = predictUnconfined(bw, conf, ubuwGuess, options)
            % Applies closed-channel linear momentum on an actuator disk to
            % predict unconfined performance from confined performance
            % data.
            %
            % Inputs (required)
            %   conf      - A structure of confined performance data with fields
            %               as described in the BWClosedChannel class documentation.
            %   ubuwGuess - Initial guess for ub/uw for use in iteration
            %               (default: 1.4)
            % Inputs (name-value pairs)
            %   correctionType - Type of blockage correction to apply: 
            %       "standard":   scales the confined data by the unconfined freestream velocity (UinfPrime, default)
            %       "bluff body": scales the confined data by the bypass velocity (ub)
            %   overrideScalingVel - Overrides the scaling velocity used in
            %                        the blockage correction to the one specified. Allowable
            %                        values are 'uw', 'ub', 'ut', or 'UinfPrime'
            % Outputs
            %   unconf    - A structure with the same size as the input
            %               conf structure and the following fields (if the
            %               field was not present in conf, it will not be
            %               present in unconf):
            %               Uinf  - unconfined freestream velocity (equal to specified scaling velocity)
            %               CT  - unconfined thrust coefficient
            %               CP  - unconfined performance coefficient
            %               CQ  - unconfined torque coefficient
            %               CL  - unconfined lateral force coefficient
            %               CF  - unconfined resultant force coefficient
            %               TSR - unconfined tip-speed ratio
            %               velRatio - Ratio between confined freestream
            %                          velocity and unconfined freestream velocity.
            %   conf      - The input structure with the following fields
            %               added:
            %               uw       - Wake velocity estimated from closed-channel LMAD
            %               ub       - Bypass velocity estimated from closed-channel LMAD
            %               ut       - Velocity at the turbine estimated from closed-channel LMAD
            %               UinfPrime  - Unconfined freestream velocity estimated from closed-channel LMAD
            %               ubuwIter - Iteration diagnostics for ub/uw
            %               isPhys   - Results of physical validity checks on uw, ub, ut
            %
            % See also: BWClosedChannel, solveLMAD, forecastConfined, checkPhysicalValidity

            % Define function input arguments, validation, and default
            % values.
            arguments
                bw
                conf
                ubuwGuess (1,:) = [1.4]
                options.correctionType {mustBeText, ismember(options.correctionType, {'standard', 'bluff body'})} = 'standard'
                options.scalingVelOverride {mustBeText, ismember(options.scalingVelOverride, {'', 'uw', 'ub', 'ut', 'UinfPrime'})} = ''
            end

            % Apply closed channel linear momentum model
            conf = bw.solveLMAD(conf, ubuwGuess);

            % Convert to unconfined using the appropriate (or requested)
            % scaling velocity
            scalingVelName = bw.checkScalingVel(options.correctionType, options.scalingVelOverride);

            for i = 1:size(conf,1)
                for j = 1:size(conf,2)
                    unconf(i,j) = bw.convertConfToUnconf(conf(i,j), conf(i,j).(scalingVelName));
                end
            end
        end

        function scalingVelName = checkScalingVel(bw, correctionType, scalingVelName)
            % Given the correction type and specified scaling velocity,
            % determines if the specified scaling velocity overrides the
            % default for that correction type. If so, a message is
            % printed. The name of the final scaling velocity to be used is
            % returned.

            % Compare model in use and expected scaling velocity
            matchInd = find(ismember(bw.correctionModes, correctionType));

            % If valid correction type supplied, check if the requested
            % scaling velocity matches
            if ~isempty(matchInd)
                expectedVelName = bw.expectedVel{matchInd};
            else
                error('Invalid correction type specified.');
            end

            % If requested scaling velocity is different than expected,
            % print a message
            if ~isempty(scalingVelName)
                check = strcmp(scalingVelName, expectedVelName);
                if ~check
                    fprintf('Overriding scaling velocity for %s %s blockage correction to "%s" (normally "%s").\n', ...
                            bw.getModelName, correctionType, scalingVelName, expectedVelName);
                end
            else
                scalingVelName = expectedVelName;
            end
        end

        %% Main method for performing a bluff-body analytical blockage forecasting
        function [conf_2, conf_1] = forecastConfined(bw, conf_1, beta_2, ubuwGuess)
            % Uses a closed-channel bluff-body blockage correction to
            % forecast performance at blockage 2 using performance data at
            % blockage 1. The forecast is performed using an analytical
            % relationship that is similar to that described by Hunt et al
            % (arxiv link TBD), but that instead uses the closed-channel
            % equations from Ross and Polagye.
            %
            % Inputs (required)
            %   conf_1    - A structure of confined performance data at blockage 1 with fields
            %               as described in the BWClosedChannel class documentation.
            %   beta_2    - The target blockage at which forecasted
            %               performance data is desired, as a fraction.
            %   ubuwGuess - Initial guess for ub/uw for use in iteration
            %               (default: 1.4)
            % Outputs
            %   conf_2    - A structure with the same size and fields as
            %               output conf_1, but with each field corresponding to
            %               the forecasted performance and velocities at beta_2. 
            %   conf_1    - The input structure with the following fields added:
            %               uw       - Wake velocity estimated from closed-channel LMAD
            %               ub       - Bypass velocity estimated from closed-channel LMAD
            %               ut       - Velocity at the turbine estimated from closed-channel LMAD
            %               UinfPrime  - Unconfined freestream velocity estimated from closed-channel LMAD
            %               ubuwIter - Iteration diagnostics for ub/uw
            %               isPhys   - Results of physical validity checks on uw, ub, ut
            %
            % See also: BWClosedChannel, solveLMAD, predictUnconfined, checkPhysicalValidity

            % Define function input arguments, validation, and default
            % values.
            arguments
                bw
                conf_1
                beta_2 (1,1)
                ubuwGuess (1,:) = [1.4]
            end

            % Solve closed-channel LMAD at confinement 1
            conf_1 = bw.solveLMAD(conf_1, ubuwGuess);

            for i = 1:size(conf_1, 1) % For each row of input data
                for j = 1:size(conf_1, 2) % For each column of input data

                    % Now, assume:
                    % Bypass velocity is the same at both blockages: ub_1 = ub_2
                    % Rotor thrust is the same at both blockages: CT_1*(Uinf_1)^2 = CT_2*(Uinf_2)^2
                    % Therefore, solve for uw using Equation 23
                    % This *should* be the same as conf_1(i,j).uw
                    uw_2 = bw.solveUw_direct(conf_1(i,j).ub, conf_1(i,j).CT, conf_1(i,j).Uinf);

                    % Solve for ut via known uw_2, ub_2, beta_2
                    ut_2 = bw.solveUt_direct(uw_2, conf_1(i,j).ub, beta_2);

                    % Solve for Uinf with known uw_2, ub_2, ut_2, beta_2
                    Uinf_2 = bw.solveUinf_direct(uw_2, conf_1(i,j).ub, ut_2, beta_2);

                    % Scale conf_1 by Uinf_2
                    currForecast = bw.convertConfToUnconf(conf_1(i,j), Uinf_2);

                    % Set h to avoid errors
                    currForecast.h = NaN .* ones(size(currForecast.Uinf));
                    currForecast.beta = beta_2 .* ones(size(currForecast.Uinf));
 
                    % Send this back through LMAD to get velocities and
                    % check
                    currForecast = bw.solveLMAD(currForecast, ubuwGuess);

                    % Save
                    conf_2(i,j) = currForecast;
                end
            end
        end
    end

    methods (Static, Access = public)
        %% Core Equations from Ross and Polagye
        
        % Ross and Polagye Equation 21
        function utuw = solveUtUw(ubuw, beta)
            % Ross and Polagye (2020) Equation 21
            % Solves for the ratio between the velocity at the turbine and the velocity
            % of the core flow using Barnsley and Wellicome's method.
            % Inputs:
            %   beta - blockage ratio
            %   ubuw - ratio between ub (bypass velocity) and uw (core wake velocity)
            % Outputs:
            %   utuw - ratio between ut (velocity at turbine) and uw (core flow 
            %          downstream of turbine)
            utuw = (-1 + sqrt(1 + beta .* (ubuw.^2 - 1))) ./ (beta .* (ubuw - 1));
        end
        
        % Ross and Polagye Equation 22
        function Uinfuw = solveUinfUw_blockage(ubuw, beta, utuw)
            % Ross and Polagye (2020) Equation 22
            % Solves for the ratio between the freestream velocity and the
            % core flow using Barnsley and Wellicome's method, via a
            % relationship with blockage.
            % Inputs:
            %   ubuw - ratio between ub (bypass velocity) and uw (core wake velocity)
            %   beta - blockage
            %   utuw - ratio between ut (velocity at turbine) and uw (core wake velocity)
            % Outputs:
            %   Uinfuw - ratio between Uinf (freestream velocity upstream of turbine) and 
            %          uw (core wake velocity)
            Uinfuw = ubuw - beta .* utuw .* (ubuw - 1);
        end
        
        % Ross and Polagye Equation 23
        function Uinfuw = solveUinfUw_thrust(ubuw, CT)
            % Ross and Polagye (2020) Equation 23
            % Solves for the ratio between the freestream velocity and the
            % core flow using Barnsley and Wellicome's method, via a
            % relationship with thrust
            % Inputs:
            %   ubuw - ratio between ub (bypass flow) and uw (core flow)
            %   CT   - thrust coefficient
            % Outputs:
            %   Uinfuw - ratio between Uinf (freestream velocity upstream of turbine) and 
            %          uw (core flow downstream of turbine)
            Uinfuw = sqrt((ubuw.^2 - 1) ./ CT);
        end

        %% Different forms of Ross and Polagye equations for solving for specific velocities

        function ut = solveUt(ubuw, beta, CT, Uinf)
            % Using Ross and Polagye (2020) Equations 21 and 23 and a known
            % ubuw, CT, and Uinf, solves for ut, the velocity at the turbine
            utuw = BWClosedChannel.solveUtUw(ubuw, beta);
            uw = BWClosedChannel.solveUw(ubuw, CT, Uinf);
            ut = utuw .* uw;
        end

        function uw = solveUw(ubuw, CT, Uinf)
            % Using Ross and Polagye Equation 23 and a known ubuw, CT, and
            % Uinf, solve for uw
            UinfUw = BWClosedChannel.solveUinfUw_thrust(ubuw, CT);
            uw = Uinf ./ UinfUw;
        end

        function ub = solveUb(ubuw, CT, Uinf)
            % Using Ross and Polagye Equation 23 and a known ubuw, CT, and Uinf,
            % solves for ub, the bypass velocity.
            UinfUw = BWClosedChannel.solveUinfUw_thrust(ubuw, CT);
            ub = ubuw ./ UinfUw .* Uinf;
        end

        % (Ross and Polagye (2020) EQ 20).
        function UinfPrime = solveUinfPrime(Uinf, CT, ut)
            % Solves for the unconfined freestream velocity using linear
            % momentum theory (Ross and Polagye (2020) EQ 20).
            % Inputs:
            % Uinf - Undisturbed freestream velocity at confined condition
            % CT - Thrust coefficient at confined condition
            % ut - Velocity at the turbine (common between confined and
            %      unconfined condition)
            % Outputs:
            % UinfPrime - The undisturbed freestream velocity at the unconfined
            %           condition
            UinfPrime = Uinf .* (CT./4+(ut./Uinf).^2)./(ut./Uinf);
        end

        %% Core equations rearranged for analytical forecasting

        function uw = solveUw_direct(ub, CT, Uinf)
            % Rearranged Ross and Polagye Eq 23: Given a known ub, CT and
            % Uinf, solve for uw.
            uw = sqrt(ub.^2 - CT .* Uinf.^2);
        end


        function ut = solveUt_direct(uw, ub, beta)
            % Rearranged Ross and Polagye Equation 21: Given a known uw,
            % ub, and beta, calculates ut.
            ubuw = ub ./ uw;
            ut = uw .* BWClosedChannel.solveUtUw(ubuw, beta);
        end

        function Uinf = solveUinf_direct(uw, ub, ut, beta)
            % Rearranged Ross and Polagye Equation 22: Given a known uw,
            % ub, ut, and beta, solve for Uinf.
            ubuw = ub ./ uw;
            utuw = ut ./ uw;
            Uinf = uw .* BWClosedChannel.solveUinfUw_blockage(ubuw, beta, utuw);
        end

        %% Iteration scheme for solving for ub/uw, ut, and ub via Ross and Polagye EQs 21-23
        
        function [ubuw, err, exitFlag] = convergeUbUw(ubuwGuess, conf)
            % Using Ross and Polagye Equations 21-23, iterates to find ub/uw that
            % satisfies both EQ 22 and EQ 23, and returns that ub/uw.
            % Inputs:
            %   ubuwGuess - Initial guess for ub/uw for use in iteration
            %               (default: 1.4)            
            %   conf      - Confined performance data as described in class
            %               documentation
            % Outputs:
            %   ubuw - Value of ubuw that satisfies both EQ 22 and EQ 23
            %           (scalar)
            %   err  - Error between EQ 22 and EQ 23 at convergence
            %   exitFlag - Exit condition at convergence

            % Preallocate
            nPoints = length(conf.CT);
            ubuw = zeros(size(conf.CT));
            err = zeros(size(conf.CT));
            exitFlag = zeros(size(conf.CT));

            % Loop through points
            for k = 1:nPoints

                % If good point, proceed
                if all(~isnan(ubuwGuess)) && (conf.CT(k) >= 0)
                    currFun = @(ubuw) BWClosedChannel.ubuwCompare(ubuw, conf.beta(k), conf.CT(k));
                    % [ubuw(k), err(k), exitFlag(k)] = fzero(currFun, ubuwGuess);
                    [ubuw(k), err(k), exitFlag(k)] = fminsearch(currFun, ubuwGuess);
                else
                    warning('Negative CT value or bad ubuwGuess: skipping application of LMAD for this point');
                    ubuw(k) = nan;
                    err(k) = nan;
                    exitFlag(k) = nan;
                end
            end
        end

        function [UinfUw_blockage, UinfUw_thrust] = solveUinfUw_both(ubuw, beta, CT)
            % Calculates the ratio between Uinf and uw via Ross and Polagye EQs
            % 21/22, as well as Ross and Polagye EQ 23, and returns both
            % values.
            utuw = BWClosedChannel.solveUtUw(ubuw, beta); % Use guess for ub/uw to solve for ut/uw
            UinfUw_blockage = BWClosedChannel.solveUinfUw_blockage(ubuw, beta, utuw); % Solve for Uinf/uw one way
            UinfUw_thrust = BWClosedChannel.solveUinfUw_thrust(ubuw, CT); % Solve for Uinf/uw another way
        end

        function err = ubuwCompare(ubuwGuess, beta, CT)
            % Calculates the ratio between Uinf and uw via Ross and Polagye EQs
            % 21/22, as well as Ross and Polagye EQ 23, and returns the error
            % between the two methods.
            [UinfUw_blockage, UinfUw_thrust] = BWClosedChannel.solveUinfUw_both(ubuwGuess, beta, CT);

            % Fzero error
            % err = UinfUw_blockage - UinfUw_thrust; % Compute error between those values

            % Fminsearch error
            err = abs(UinfUw_blockage - UinfUw_thrust); % Compute error between those values
        end

        %% Analytical forecasting - error analysis

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
                        end
                    end

                end
            end
        end

        %% Diagnostics

        function [err, ax] = plotUbUwConvergenceRegion(ubuwTest, beta, CT)
            % Plots the error between the two methods for calculating Uinf/uw
            % (used in ubuw iteration) for a specific beta and CT. Useful for
            % visually determining whether a given point can ever converge.
            if (isempty(ubuwTest))
                ubuwTest = 0:0.01:10;
            end
            % Generate test values for ubuw and evaluate error for each
            err = BWClosedChannel.ubuwCompare(ubuwTest, beta, CT);

            % Plot
            [fig] = figure();
            ax = axes(fig);
            grid(ax, 'on'); hold(ax, 'on');
            plot(ax, ubuwTest, zeros(size(ubuwTest)), '-k');
            plot(ax, ubuwTest, real(err), '-', 'marker', '.');
            plot(ax, ubuwTest, imag(err), '-', 'marker', '.');
            xlabel(ax, '$u_2$');
            ylabel(ax, '$u_1$ error');
        end

        %% Printing stuff out
        function modelName = getModelName()
            % Returns a label-friendly version of the blockage correction
            % model name.
            modelName = 'BW Closed-Channel';
        end
    end
end