% Aidan Hunt
%
% The BWClosedChannel class implements the Barnsley and Wellicome (1990)
% closed-channel blockage correction as described by Ross and Polagye in
% "An experimental assessment of analytical blockage corrections" (2020;
% https://doi.org/10.1016/j.renene.2020.01.135).
%
% Applies linear momentum to an actuator disk in closed channel flow
% (parallel sided-tube) to solve for the following flow properties:
%   u1      - Core wake velocity
%   u2      - Bypass velocity
%   ut      - Velocity at the turbine
%   V0Prime - Unconfined freestream velocity
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
%   V0   (required)   - undisturbed upstream freestream velocity (m/s)
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
        expectedVel = {'V0Prime', 'u2'}; % Default scaling velocity for each blockage correction mode
    end

    methods (Access = public)
        %% Main method for solving closed-channel linear momentum model
        function [conf] = solveLMAD(bw, conf, u2u1Guess)
            % Solves linear momentum on an actuator disk using a
            % closed-channel model. Velocities are estimated using
            % Equations 20-23 from Ross and Polagye (2020).
            %
            % Inputs
            %   conf      - A structure of confined performance data with fields
            %               as described in the BWClosedChannel class documentation.
            %   u2u1Guess - A two element vector containing the maximum and
            %               minimum guess values of u2/u1 for use in
            %               iteration (default: [1.04 20])
            % Outputs
            %   conf      - The input structure with the following fields
            %               added:
            %               u1       - Wake velocity estimated from closed-channel LMAD
            %               u2       - Bypass velocity estimated from closed-channel LMAD
            %               ut       - Velocity at the turbine estimated from closed-channel LMAD
            %               V0Prime  - Unconfined freestream velocity estimated from closed-channel LMAD
            %               u2u1Iter - Iteration diagnostics for u2/u1
            %               isPhys   - Results of physical validity checks on u1, u2, ut
            %
            % See also: BWClosedChannel, predictUnconfined, forecastConfined, checkPhysicalValidity

            % Define function input arguments, validation, and default
            % values.
            arguments
                bw
                conf
                u2u1Guess (1,2) = [1.4 20]
            end

            % Check input for correct sizing
            conf = bw.checkInputSizes(conf);

            for i = 1:size(conf, 1)
                for j = 1:size(conf, 2)

                    % Solve for u2u1 via iteration
                    [u2u1, u2u1Err, u2u1ExitFlag] = bw.convergeU2U1(u2u1Guess, conf(i,j));
    
                    % Compute channel velocities
                    conf(i,j).u1 = bw.solveU1(u2u1, conf(i,j).CT, conf(i,j).V0);
                    conf(i,j).u2 = bw.solveU2(u2u1, conf(i,j).CT, conf(i,j).V0);
                    conf(i,j).ut = bw.solveUt(u2u1, conf(i,j).beta, conf(i,j).CT, conf(i,j).V0);
                    conf(i,j).V0Prime = bw.solveV0Prime(conf(i,j).V0, conf(i,j).CT, conf(i,j).ut);
    
                    % Package diagnotics about u2u1 iteration
                    [V0U1_blockage, V0U1_thrust] = bw.solveV0U1_both(u2u1, conf(i,j).beta, conf(i,j).CT);
                    conf(i,j).u2u1Iter = bw.packageDiagnostics(u2u1, [V0U1_blockage V0U1_thrust], ...
                                                                 u2u1Err, u2u1ExitFlag);
    
                    % Check physical validity
                    conf(i,j).isPhys = bw.checkPhysicalValidity(conf(i,j));
                end
            end
        end

        %% Main method for performing blockage correction
        function [unconf, conf] = predictUnconfined(bw, conf, u2u1Guess, options)
            % Applies closed-channel linear momentum on an actuator disk to
            % predict unconfined performance from confined performance
            % data.
            %
            % Inputs (required)
            %   conf      - A structure of confined performance data with fields
            %               as described in the BWClosedChannel class documentation.
            %   u2u1Guess - A two element vector containing the maximum and
            %               minimum guess values of u2/u1 for use in
            %               iteration (default: [1.04 20])
            % Inputs (name-value pairs)
            %   correctionType - Type of blockage correction to apply: 
            %       "standard":   scales the confined data by the unconfined freestream velocity (V0Prime, default)
            %       "bluff body": scales the confined data by the bypass velocity (u2)
            %   overrideScalingVel - Overrides the scaling velocity used in
            %                        the blockage correction to the one specified. Allowable
            %                        values are 'u1', 'u2', 'ut', or 'V0Prime'
            % Outputs
            %   unconf    - A structure with the same size as the input
            %               conf structure and the following fields (if the
            %               field was not present in conf, it will not be
            %               present in unconf):
            %               V0  - unconfined freestream velocity (equal to specified scaling velocity)
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
            %               u1       - Wake velocity estimated from closed-channel LMAD
            %               u2       - Bypass velocity estimated from closed-channel LMAD
            %               ut       - Velocity at the turbine estimated from closed-channel LMAD
            %               V0Prime  - Unconfined freestream velocity estimated from closed-channel LMAD
            %               u2u1Iter - Iteration diagnostics for u2/u1
            %               isPhys   - Results of physical validity checks on u1, u2, ut
            %
            % See also: BWClosedChannel, solveLMAD, forecastConfined, checkPhysicalValidity

            % Define function input arguments, validation, and default
            % values.
            arguments
                bw
                conf
                u2u1Guess (1,2) = [1.4 20]
                options.correctionType {mustBeText, ismember(options.correctionType, {'standard', 'bluff body'})} = 'standard'
                options.scalingVelOverride {mustBeText, ismember(options.scalingVelOverride, {'', 'u1', 'u2', 'ut', 'V0Prime'})} = ''
            end

            % Apply closed channel linear momentum model
            conf = bw.solveLMAD(conf, u2u1Guess);

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
        function [conf_2, conf_1] = forecastConfined(bw, conf_1, beta_2, u2u1Guess)
            % Uses a closed-channel bluff-body blockage correction to
            % forecast performance at blockage 2 using performance data at
            % blockage 1. The forecast is performed using an analytical
            % relationship that will be described in an upcoming
            % publication.
            %
            % Inputs (required)
            %   conf_1    - A structure of confined performance data at blockage 1 with fields
            %               as described in the BWClosedChannel class documentation.
            %   beta_2    - The target blockage at which forecasted
            %               performance data is desired, as a fraction.
            %   u2u1Guess - A two element vector containing the maximum and
            %               minimum guess values of u2/u1 for use in
            %               iteration (default: [1.04 20])
            % Outputs
            %   conf_2    - A structure with the same size and fields as
            %               output conf_1, but with each field corresponding to
            %               the forecasted performance and velocities at beta_2. 
            %   conf_1    - The input structure with the following fields added:
            %               u1       - Wake velocity estimated from closed-channel LMAD
            %               u2       - Bypass velocity estimated from closed-channel LMAD
            %               ut       - Velocity at the turbine estimated from closed-channel LMAD
            %               V0Prime  - Unconfined freestream velocity estimated from closed-channel LMAD
            %               u2u1Iter - Iteration diagnostics for u2/u1
            %               isPhys   - Results of physical validity checks on u1, u2, ut
            %
            % See also: BWClosedChannel, solveLMAD, predictUnconfined, checkPhysicalValidity

            % Define function input arguments, validation, and default
            % values.
            arguments
                bw
                conf_1
                beta_2 (1,1)
                u2u1Guess (1,2) = [1.4 20]
            end

            % Solve closed-channel LMAD at confinement 1
            conf_1 = bw.solveLMAD(conf_1, u2u1Guess);

            for i = 1:size(conf_1, 1) % For each row of input data
                for j = 1:size(conf_1, 2) % For each column of input data

                    % Now, assume:
                    % Bypass velocity is the same at both blockages: u2_1 = u2_2
                    % Rotor thrust is the same at both blockages: CT_1*(V0_1)^2 = CT_2*(V0_2)^2
                    % Therefore, solve for u1 using Equation 23
                    % This *should* be the same as conf_1(i,j).u1
                    u1_2 = bw.solveU1_direct(conf_1(i,j).u2, conf_1(i,j).CT, conf_1(i,j).V0);

                    % Solve for ut via known u1_2, u2_2, beta_2
                    ut_2 = bw.solveUt_direct(u1_2, conf_1(i,j).u2, beta_2);

                    % Solve for V0 with known u1_2, u2_2, ut_2, beta_2
                    V0_2 = bw.solveV0_direct(u1_2, conf_1(i,j).u2, ut_2, beta_2);

                    % Scale conf_1 by V0_2
                    currForecast = bw.convertConfToUnconf(conf_1(i,j), V0_2);

                    % Set d0 to avoid errors
                    currForecast.d0 = NaN;
                    currForecast.beta = beta_2;
 
                    % Send this back through LMAD to get velocities and
                    % check
                    currForecast = bw.solveLMAD(currForecast, u2u1Guess);

                    % Save
                    conf_2(i,j) = currForecast;
                end
            end
        end
    end

    methods (Static, Access = public)
        %% Core Equations from Ross and Polagye
        
        % Ross and Polagye Equation 21
        function utu1 = solveUtU1(u2u1, beta)
            % Ross and Polagye (2020) Equation 21
            % Solves for the ratio between the velocity at the turbine and the velocity
            % of the core flow using Barnsley and Wellicome's method.
            % Inputs:
            %   beta - blockage ratio
            %   u2u1 - ratio between u2 (bypass velocity) and u1 (core wake velocity)
            % Outputs:
            %   utu1 - ratio between ut (velocity at turbine) and u1 (core flow 
            %          downstream of turbine)
            utu1 = (-1 + sqrt(1 + beta .* (u2u1.^2 - 1))) ./ (beta .* (u2u1 - 1));
        end
        
        % Ross and Polagye Equation 22
        function V0u1 = solveV0U1_blockage(u2u1, beta, utu1)
            % Ross and Polagye (2020) Equation 22
            % Solves for the ratio between the freestream velocity and the
            % core flow using Barnsley and Wellicome's method, via a
            % relationship with blockage.
            % Inputs:
            %   u2u1 - ratio between u2 (bypass velocity) and u1 (core wake velocity)
            %   beta - blockage
            %   utu1 - ratio between ut (velocity at turbine) and u1 (core wake velocity)
            % Outputs:
            %   V0u1 - ratio between V0 (freestream velocity upstream of turbine) and 
            %          u1 (core wake velocity)
            V0u1 = u2u1 - beta .* utu1 .* (u2u1 - 1);
        end
        
        % Ross and Polagye Equation 23
        function V0u1 = solveV0U1_thrust(u2u1, CT)
            % Ross and Polagye (2020) Equation 23
            % Solves for the ratio between the freestream velocity and the
            % core flow using Barnsley and Wellicome's method, via a
            % relationship with thrust
            % Inputs:
            %   u2u1 - ratio between u2 (bypass flow) and u1 (core flow)
            %   CT   - thrust coefficient
            % Outputs:
            %   V0u1 - ratio between V0 (freestream velocity upstream of turbine) and 
            %          u1 (core flow downstream of turbine)
            V0u1 = sqrt((u2u1.^2 - 1) ./ CT);
        end

        %% Different forms of Ross and Polagye equations for solving for specific velocities

        function ut = solveUt(u2u1, beta, CT, V0)
            % Using Ross and Polagye (2020) Equations 21 and 23 and a known
            % u2u1, CT, and V0, solves for ut, the velocity at the turbine
            utu1 = BWClosedChannel.solveUtU1(u2u1, beta);
            u1 = BWClosedChannel.solveU1(u2u1, CT, V0);
            ut = utu1 .* u1;
        end

        function u1 = solveU1(u2u1, CT, V0)
            % Using Ross and Polagye Equation 23 and a known u2u1, CT, and
            % V0, solve for u1
            V0U1 = BWClosedChannel.solveV0U1_thrust(u2u1, CT);
            u1 = V0 ./ V0U1;
        end

        function u2 = solveU2(u2u1, CT, V0)
            % Using Ross and Polagye Equation 23 and a known u2u1, CT, and V0,
            % solves for u2, the bypass velocity.
            V0U1 = BWClosedChannel.solveV0U1_thrust(u2u1, CT);
            u2 = u2u1 ./ V0U1 .* V0;
        end

        % (Ross and Polagye (2020) EQ 20).
        function V0Prime = solveV0Prime(V0, CT, ut)
            % Solves for the unconfined freestream velocity using linear
            % momentum theory (Ross and Polagye (2020) EQ 20).
            % Inputs:
            % V0 - Undisturbed freestream velocity at confined condition
            % CT - Thrust coefficient at confined condition
            % ut - Velocity at the turbine (common between confined and
            %      unconfined condition)
            % Outputs:
            % V0Prime - The undisturbed freestream velocity at the unconfined
            %           condition
            V0Prime = V0 .* (CT./4+(ut./V0).^2)./(ut./V0);
        end

        %% Core equations rearranged for analytical forecasting

        function u1 = solveU1_direct(u2, CT, V0)
            % Rearranged Ross and Polagye Eq 23: Given a known u2, CT and
            % V0, solve for u1.
            u1 = sqrt(u2.^2 - CT .* V0.^2);
        end


        function ut = solveUt_direct(u1, u2, beta)
            % Rearranged Ross and Polagye Equation 21: Given a known u1,
            % u2, and beta, calculates ut.
            u2u1 = u2 ./ u1;
            ut = u1 .* BWClosedChannel.solveUtU1(u2u1, beta);
        end

        function V0 = solveV0_direct(u1, u2, ut, beta)
            % Rearranged Ross and Polagye Equation 22: Given a known u1,
            % u2, ut, and beta, solve for V0.
            u2u1 = u2 ./ u1;
            utu1 = ut ./ u1;
            V0 = u1 .* BWClosedChannel.solveV0U1_blockage(u2u1, beta, utu1);
        end

        %% Iteration scheme for solving for u2/u1, ut, and u2 via Ross and Polagye EQs 21-23
        
        function [u2u1, err, exitFlag] = convergeU2U1(u2u1Guess, conf)
            % Using Ross and Polagye Equations 21-23, iterates to find u2/u1 that
            % satisfies both EQ 22 and EQ 23, and returns that u2/u1.
            % Inputs:
            %   u2u1Guess - Initial guess for u2/u1 (as a two element vector for fzero)
            %   conf      - Confined performance data as described in class
            %               documentation
            % Outputs:
            %   u2u1 - Value of u2u1 that satisfies both EQ 22 and EQ 23
            %           (scalar)
            %   err  - Error between EQ 22 and EQ 23 at convergence
            %   exitFlag - Exit condition at convergence

            % Preallocate
            nPoints = length(conf.CT);
            u2u1 = zeros(size(conf.CT));
            err = zeros(size(conf.CT));
            exitFlag = zeros(size(conf.CT));

            % Loop through points
            for k = 1:nPoints
                % Check if good value of u2u1 guess for this point
                u2u1Guess = BWClosedChannel.checkU2U1Guess(u2u1Guess, conf.beta(k), conf.CT(k));

                % If good point, proceed
                if all(~isnan(u2u1Guess)) && (conf.CT(k) >= 0)
                    currFun = @(u2u1) BWClosedChannel.u2u1Compare(u2u1, conf.beta(k), conf.CT(k));
                    [u2u1(k), err(k), exitFlag(k)] = fzero(currFun, u2u1Guess);
                else
                    warning('Negative CT value or bad u2u1Guess: skipping application of LMAD for this point');
                    u2u1(k) = nan;
                    err(k) = nan;
                    exitFlag(k) = nan;
                end
            end
        end

        function [V0U1_blockage, V0U1_thrust] = solveV0U1_both(u2u1, beta, CT)
            % Calculates the ratio between V0 and u1 via Ross and Polagye EQs
            % 21/22, as well as Ross and Polagye EQ 23, and returns both
            % values.
            utu1 = BWClosedChannel.solveUtU1(u2u1, beta); % Use guess for u2/u1 to solve for ut/u1
            V0U1_blockage = BWClosedChannel.solveV0U1_blockage(u2u1, beta, utu1); % Solve for V0/u1 one way
            V0U1_thrust = BWClosedChannel.solveV0U1_thrust(u2u1, CT); % Solve for V0/u1 another way
        end

        function err = u2u1Compare(u2u1Guess, beta, CT)
            % Calculates the ratio between V0 and u1 via Ross and Polagye EQs
            % 21/22, as well as Ross and Polagye EQ 23, and returns the error
            % between the two methods.
            [V0U1_blockage, V0U1_thrust] = BWClosedChannel.solveV0U1_both(u2u1Guess, beta, CT);
            err = V0U1_blockage - V0U1_thrust; % Compute error between those values
        end

        %% Analytical forecasting - error analysis

        function [forecastErr] = quantifyForecastError(confPredict, confReference)
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
            %                    structures: CP, CT, CL, CF, CQ, V0.

            forecastErr = struct;
            predictMetrics = {'CP', 'CT', 'CL', 'CF', 'CQ', 'V0'}; % Meaningful metrics for quantifying error

            for j = 1:size(confPredict, 2) % For each forecast blockage
                % Get the current TSR
                TSRInterp = confReference(j).TSR;

                for i = 1:size(confPredict, 1) % For each dataset that was used to make a forecast

                    % Save TSR
                    forecastErr(i,j).TSR = TSRInterp;

                    for k = 1:length(predictMetrics)
                        if isfield(confPredict, predictMetrics{k}) && isfield(confReference, predictMetrics{k})
                            % Interpolate onto basis TSR
                            currInterp = interp1(confPredict(i,j).TSR, confPredict(i,j).(predictMetrics{k}), TSRInterp, 'linear');
                            
                            % Compute error
                            forecastErr(i,j).(predictMetrics{k}) = currInterp - confReference(j).(predictMetrics{k});
                        end
                    end

                end
            end
        end

        %% Diagnostics

        function u2u1Guess = checkU2U1Guess(u2u1Guess, beta, CT)
            % Checks u1u2Guess to make sure that it returns a workable
            % value for iteration. Given an input guess for u2/u1, and
            % specified values of beta and CT, checks whether the guess
            % results in an interval with a sign change (necessary for
            % fzero) and attempts to adjust the guess if not.

            % Prime the loop
            badGuess = true;

            while badGuess
                % Evaluate current guess
                currErr = BWClosedChannel.u2u1Compare(u2u1Guess, beta, CT);

                if (isnan(currErr(1)) || ~isreal(currErr)) % If lower bound is too far left (and has become nan or complex), move back right
                    % u2u1Guess(1) = u2u1Guess(1) + 0.01;
                    warning('Initial u2u1Guess results in a complex value at an interval endpoint for beta = %g and CT = %g. This point may not converge. Skipping....', beta, CT);
                    u2u1Guess = nan;
                    badGuess = false;
                else % Check signs
                    errSigns = sign(currErr);
                    if (errSigns(1) == errSigns(2)) % If no sign change
                        warning('Initial u2u1Guess does not provide interval with sign change for beta = %g and CT = %g. Attempting to adjust....', beta, CT);
                        % Expand the interval and see what happens
                        u2u1Guess(1) = u2u1Guess(1) - 0.005;
                        u2u1Guess(2) = u2u1Guess(2) + 0.005;
%                         if (errSigns(1) < 0) % If left bound is negative, back off
%                             u2u1Guess(1) = u2u1Guess(1) - 0.005;
%                         end
%                         if (errSigns(2) > 0) % If right bound is positive, increase
%                             u2u1Guess(2) = u2u1Guess(2) + 1;
%                         end
                    else % If error is real and with sign change, exit loop.
                        badGuess = false;
                    end
                end
            end
        end

        %%  Visualization
        function [err, ax] = plotU2U1ConvergenceRegion(u2u1Test, beta, CT)
            % Plots the error between the two methods for calculating V0/u1
            % (used in u2u1 iteration) for a specific beta and CT. Useful for
            % visually determining whether a given point can ever converge.
            if (isempty(u2u1Test))
                u2u1Test = 0:0.01:10;
            end
            % Generate test values for u2u1 and evaluate error for each
            err = BWClosedChannel.u2u1Compare(u2u1Test, beta, CT);

            % Plot
            [fig] = figure();
            ax = axes(fig);
            grid(ax, 'on'); hold(ax, 'on');
            plot(ax, u2u1Test, zeros(size(u2u1Test)), '-k');
            plot(ax, u2u1Test, real(err), '-', 'marker', '.');
            plot(ax, u2u1Test, imag(err), '-', 'marker', '.');
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