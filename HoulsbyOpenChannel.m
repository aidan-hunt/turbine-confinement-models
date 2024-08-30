% Aidan Hunt
%
% The HoulsbyOpenChannel class implements the open-channel linear momentum
% theory of Houlsby et al. "Application of linear momentum actuator disc
% theory to open channel flow" (2008;
% https://ora.ox.ac.uk/objects/uuid:5576d575-7bac-44b6-ac79-f698edcda40e),
% see also "The Power Available to Tidal Turbines in an Open Channel Flow"
% by Houlsby and Vogel (2017; http://dx.doi.org/10.1680/jener.15.00035).
% Additionally, this code implements the Houlsby et al's linear momentum
% theory as a blockage correction following the implementation of Ross and
% Polagye in "An experimental assessment of analytical blockage
% corrections" (2020; https://doi.org/10.1016/j.renene.2020.01.135).
%
% Open-channel linear momentum on an actuator disk is used to solve for the
% following properties:
%   uw      - Core wake velocity
%   ub      - Bypass velocity
%   ut      - Velocity at the turbine
%   UinfPrime - Unconfined freestream velocity
%   Fr      - Depth-based Froude number
%   dhToh   - Free-surface drop across the turbine rotor normalized by
%             upstream depth
% These quantities can then be used for blockage correction or blockage
% forecasting.
%
% To use the HoulsbyOpenChannel class, construct a HoulsbyOpenChannel object
% using the following syntax:
%   hb = HoulsbyOpenChannel()
% and call methods using the dot notation (i.e., hb.solveLMAD(...))
% 
% HoulsbyOpenChannel Methods:
%   solveLMAD         - Solves for flow velocities in the vicinity of the rotor
%                       using open-channel linear momentum on an actuator
%                       disk (LMAD) theory and confined performance data.
%                       No blockage correction is performed.
%   predictUnconfined - Uses open-channel LMAD theory to predict unconfined 
%                       performance from confined performance data.
%                       "Standard" blockage corrections or "bluff body"
%                       blockage corrections may be used.
%   forecastConfined  - Uses open-channel LMAD theory and a bluff body
%                       blockage correction basis to predict performance at
%                       one blockage ratio using data at a different
%                       blockage ratio.
%   linearForecast    - Uses open-channel LMAD theory and an assumption of 
%                       a linear relationship between turbine performance
%                       (i.e., CP, CT) and blockage at constant TSR (e.g.,
%                       as observed by Kinsey and Dumas (2017)) to forecast
%                       performance across blockages.
% 
% The methods above expect that confined performance data is provided as an
% mxn structure array, conf, with the following fields:
%   beta (required)   - blockage ratio
%   Uinf   (required)   - undisturbed upstream freestream velocity (m/s)
%   h   (required)   - undisturbed upstream water depth (m)
%   CT   (required)   - thrust coefficient
%   CP   (optional)   - performance coefficient
%   CQ   (optional)   - torque coefficient
%   CL   (optional)   - lateral force coefficient
%   CF   (optional)   - resultant force coefficient
%   TSR  (optional)   - tip-speed ratio
% The fields of conf(i,j) must be vectors that are all the same size. The
% methods above may be used to apply corrections and forecasts to multiple
% datasets at once by specifying each dataset as an element of conf (e.g.,
% conf(i,j)).
%
% These methods utilize helper methods that implement the core equation set
% for the Houlsby-inspired correction as presented in Ross and
% Polagye (2020), which are also available for use. The HoulsbyOpenChannel
% class extends the BCBase and BWClosedChannel classes.
%
% See also: BCBase, BWClosedChannel

classdef HoulsbyOpenChannel < BWClosedChannel

    % The HoulsbyOpenChannel class inherits all properties of the BCBase
    % and BWClosedChannel classes.
    properties (Constant, Access=protected)
        g = 9.81; % Gravitational acceleration
    end

    methods (Access = public)
        %% Main method for solving open-channel linear momentum on an actuator disk
        function [conf] = solveLMAD(hb, conf, uGuess, options)
            % Solves linear momentum on an actuator disk using a
            % open-channel model. Velocities are estimated using
            % Equations 20, 43-45 from Ross and Polagye (2020).
            %
            % Inputs (required)
            %   conf      - A structure of confined performance data with fields
            %               as described in the HoulsbyOpenChannel class documentation.
            % Inputs (optional, positional)
            %   uGuess    - Guess for the value one of the linear momentum
            %               velocities (or ratio between velocities), as a
            %               scalar (default 1.4). See documentation for 
            %               'guessMode' input for valid guess types.
            % Inputs (optional, name-value pairs)
            %   guessMode - The linear momentum velocities (or ratio
            %               between velocities) represented by uGuess:
            %       'ubuw': ratio between ub and uw (default)
            %       'ubUinf': ratio between ub and Uinf
            %       'ub': the bypass velocity
            %   FrZeroLimit - Whether to evaluate linear momentum at the
            %                 closed-channel limit (Fr=0), ignoring the
            %                 provided Froude number. Default: false.
            % Outputs
            %   conf      - The input structure with the following fields
            %               added:
            %               Fr       - Depth-based Froude number (rectangular channel assumed)
            %               dhToh    - Normalized free-surface drop across the rotor
            %               uw       - Wake velocity estimated from closed-channel LMAD
            %               ub       - Bypass velocity estimated from closed-channel LMAD
            %               ut       - Velocity at the turbine estimated from closed-channel LMAD
            %               UinfPrime  - Unconfined freestream velocity estimated from closed-channel LMAD
            %               ubuwIter - Iteration diagnostics for ub/uw
            %               ubIter   - Iteration diagnostics for ub
            %               isPhys   - Results of physical validity checks on uw, ub, ut
            %
            % See also: HoulsbyOpenChannel, predictUnconfined, forecastConfined, checkPhysicalValidity

            % Define function input arguments, validation, and default
            % values.
            arguments
                hb
                conf
                uGuess (1,1) = [1.04]
                options.guessMode {mustBeText, mustBeMember(options.guessMode, {'ub', 'ubuw', 'ubUinf'})} = 'ubuw'
                options.FrZeroLimit (1,1) {mustBeNumericOrLogical} = false;
            end

            % Check input for correct sizing
            conf = hb.checkInputSizes(conf);

            % If considering limiting case of Fr -> 0, set depth to
            % high value.
            if options.FrZeroLimit
                conf = hb.setClosedChannelDepth(conf);
            end

            for i = 1:size(conf, 1) % For each row of the input structure
                for j = 1:size(conf, 2) % For each column of the input structure


                    % Calculate Froude numbers
                    conf(i,j).Fr = hb.calcFroude(conf(i,j).Uinf, conf(i,j).h);

                    % If ub guess is not provided, use closed channel LMAD to generate reasonable ub guess from a ubuw guess.
                    switch options.guessMode
                        case 'ubuw'
                            [ubuw, ubuwErr, ubuwExitFlag] = hb.convergeUbUw(uGuess, conf(i,j));
                            ubGuess = hb.solveUb(ubuw, conf(i,j).CT, conf(i,j).Uinf);
                        case 'ubUinf'
                            ubGuess = uGuess .* conf(i,j).Uinf;
                        case 'ub'
                            ubGuess = uGuess;
                    end
    
                    % Iterate on Ross and Polagye EQs 43 and 44 to find value of ub that solves for both.
                    [ub, ubErr, ubExitFlag] = hb.convergeUb(ubGuess, conf(i,j));
    
                    % Solve for velocities
                    conf(i,j).uw = hb.solveUw_Thrust(ub, conf(i,j).CT, conf(i,j).Uinf);
                    conf(i,j).ub = ub;
                    conf(i,j).ut = hb.solveUt(conf(i,j).uw, conf(i,j).ub, conf(i,j).h, conf(i,j).Uinf, conf(i,j).beta);
                    conf(i,j).UinfPrime = hb.solveUinfPrime(conf(i,j).Uinf, conf(i,j).CT, conf(i,j).ut);

                    % Solve for depths
                    conf(i,j).hBypass = hb.calcBypassDepth(conf(i,j).h, conf(i,j).Uinf, conf(i,j).ub);
                    conf(i,j).hdUp = hb.calcUpstreamDiskDepth(conf(i,j).h, conf(i,j).Uinf, conf(i,j).ut);
                    conf(i,j).hdDown = hb.calcDownstreamDiskDepth(conf(i,j).hBypass, conf(i,j).ut, conf(i,j).uw);
                    conf(i,j).dhDisk = hb.calcDiskDrop(conf(i,j).Uinf, conf(i,j).CT);
                    [conf(i,j).dhToh] = hb.calcTotalSurfaceDeformation(conf(i,j).CT, conf(i,j).beta, conf(i,j).Fr);
                    conf(i,j).hFinal = conf(i,j).h .* (1 - conf(i,j).dhToh);

                    % Package diagnostics about ubuw iteration
                    if strcmp(options.guessMode, 'ubuw')
                        [UinfUw_blockage, UinfUw_thrust] = hb.solveUinfUw_both(ubuw, conf(i,j).beta, conf(i,j).CT);
                        conf(i,j).ubuwIter = hb.packageDiagnostics(ubuw, [UinfUw_blockage, UinfUw_thrust], ...
                                                                     ubuwErr, ubuwExitFlag);
                    end

                    % Package diagnostics about ub iteration
                    [uw_Fr, uw_Thrust] = hb.solveUw_both(conf(i,j).ub, conf(i,j).beta, conf(i,j).CT, conf(i,j).Uinf, conf(i,j).Fr);
                    conf(i,j).ubIter = hb.packageDiagnostics(conf(i,j).ub, [uw_Fr, uw_Thrust], ...
                                                               ubErr, ubExitFlag);

                    % Check physical validity
                    conf(i,j).isPhys = hb.checkPhysicalValidity(conf(i,j));
                end
            end
        end

        %% Main method for performing blockage correction
        % "predictUnconfined" and "checkScalingVel" are directly inherited
        % from the BWClosedChannel class. predictUnconfined is repeated
        % here to update documenation.

        function [unconf, conf] = predictUnconfined(hb, conf, uGuess, options)
            % Applies open-channel linear momentum on an actuator disk to
            % predict unconfined performance from confined performance
            % data.
            %
            % Inputs (required)
            %   conf      - A structure of confined performance data with fields
            %               as described in the HoulsbyOpenChannel class documentation.
            % Inputs (optional, positional)
            %   uGuess    - Guess for the value one of the linear momentum
            %               velocities (or ratio between velocities), as a
            %               scalar (default 1.4). See documentation for 
            %               'guessMode' input for valid guess types.
            % Inputs (optional, name-value pairs)
            %   guessMode - The linear momentum velocities (or ratio
            %               between velocities) represented by uGuess:
            %       'ubuw': ratio between ub and uw (default)
            %       'ubUinf': ratio between ub and Uinf
            %       'ub': the bypass velocity
            %   FrZeroLimit - Whether to evaluate linear momentum at the
            %                 closed-channel limit (Fr=0), ignoring the
            %                 provided Froude number. Default: false.
            %   correctionType    - Type of blockage correction to apply: 
            %       "standard":     scales the confined data by the unconfined freestream velocity (UinfPrime, default)
            %       "bluff body":   scales the confined data by the bypass velocity (ub)
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
            %               Fr       - Depth-based Froude number (rectangular channel assumed)
            %               dhToh    - Normalized free-surface drop across the rotor
            %               uw       - Wake velocity estimated from closed-channel LMAD
            %               ub       - Bypass velocity estimated from closed-channel LMAD
            %               ut       - Velocity at the turbine estimated from closed-channel LMAD
            %               UinfPrime  - Unconfined freestream velocity estimated from closed-channel LMAD
            %               ubuwIter - Iteration diagnostics for ub/uw
            %               ubIter   - Iteration diagnostics for ub
            %               isPhys   - Results of physical validity checks on uw, ub, ut
            %
            % See also: HoulsbyOpenChannel, solveLMAD, forecastConfined, checkPhysicalValidity
    
            % Define function input arguments, validation, and default
            % values.
            arguments
                hb
                conf
                uGuess (1,1) = [1.4]
                options.guessMode {mustBeText, mustBeMember(options.guessMode, {'ub', 'ubuw', 'ubUinf'})} = 'ubuw'
                options.FrZeroLimit (1,1) {mustBeNumericOrLogical} = false;
                options.correctionType {mustBeText, ismember(options.correctionType, {'standard', 'bluff body'})} = 'standard'
                options.scalingVelOverride {mustBeText, ismember(options.scalingVelOverride, {'', 'uw', 'ub', 'ut', 'UinfPrime'})} = ''
            end

            % If considering the limiting closed-channel case
            if options.FrZeroLimit
                conf = hb.setClosedChannelDepth(conf);
            end

            % Apply closed channel linear momentum model
            conf = hb.solveLMAD(conf, uGuess, guessMode=options.guessMode, FrZeroLimit=options.FrZeroLimit);

            % Convert to unconfined using the appropriate (or requested)
            % scaling velocity
            scalingVelName = hb.checkScalingVel(options.correctionType, options.scalingVelOverride);

            for i = 1:size(conf,1)
                for j = 1:size(conf,2)
                    unconf(i,j) = hb.convertConfToUnconf(conf(i,j), conf(i,j).(scalingVelName));
                end
            end

            % Just call superclass method
            % [unconf, conf] = predictUnconfined@BWClosedChannel(hb, conf, uGuess, guessMode=options.guessMode, correctionType=options.correctionType, scalingVelOverride=options.scalingVelOverride);
        end

        %% Main method for performing bluff-body analytical blockage forecasting
        function [conf_2, conf_1] = forecastConfined(hb, conf_1, beta_2, uGuess, options)
            % Uses an open-channel bluff-body blockage correction to
            % forecast performance at blockage 2 using performance data at
            % blockage 1. The forecast is performed using an analytical
            % relationship that is described in Hunt et al (arxiv link TBD)
            %
            % Inputs (required)
            %   conf_1     - A structure of confined performance data at blockage 1 with fields
            %                as described in the BWClosedChannel class documentation.
            %   beta_2     - The target blockage at which forecasted
            %                performance data is desired, as a fraction.
            % Inputs (optional, positional)
            %   uGuess    - Guess for the value one of the linear momentum
            %               velocities (or ratio between velocities), as a
            %               scalar (default 1.4). See documentation for 
            %               'guessMode' input for valid guess types.
            % Inputs (optional, name-value pairs)
            %   guessMode - The linear momentum velocities (or ratio
            %               between velocities) represented by uGuess:
            %       'ubuw': ratio between ub and uw (default)
            %       'ubUinf': ratio between ub and Uinf
            %       'ub': the bypass velocity
            %   FrZeroLimit - Whether to evaluate linear momentum at the
            %                 closed-channel limit (Fr=0), ignoring the
            %                 provided Froude number. Default: false.
            %   constantFr - Whether to hold the Froude number constant
            %                between beta_1 and beta_2 (true, default) or allow to
            %                vary (false). If constantFr=true, h at beta_2
            %                is calculated using Fr at beta_1 and Uinf at
            %                beta_2. If constantFr=false, h at beta_2 is
            %                calculated using h at beta_1 and assuming
            %                constant channel width.
            % Outputs
            %   conf_2     - A structure with the same size and fields as
            %                output conf_1, but with each field corresponding to
            %                the forecasted performance and velocities at beta_2. 
            %   conf_1     - The input structure with the following fields
            %                added:
            %                Fr       - Depth-based Froude number (rectangular channel assumed)
            %                dhToh    - Normalized free-surface drop across the rotor
            %                uw       - Wake velocity estimated from closed-channel LMAD
            %                ub       - Bypass velocity estimated from closed-channel LMAD
            %                ut       - Velocity at the turbine estimated from closed-channel LMAD
            %                UinfPrime  - Unconfined freestream velocity estimated from closed-channel LMAD
            %                ubuwIter - Iteration diagnostics for ub/uw
            %                ubIter   - Iteration diagnostics for ub
            %                isPhys   - Results of physical validity checks on uw, ub, ut
            %
            % See also: HoulsbyOpenChannel, solveLMAD, predictUnconfined, checkPhysicalValidity

            % Define function input arguments, validation, and default
            % values.
            arguments
                hb
                conf_1
                beta_2 (1,1) {mustBeInRange(beta_2, 0, 1)}
                uGuess (1,1) = [1.4]
                options.guessMode {mustBeText, mustBeMember(options.guessMode, {'ub', 'ubuw', 'ubUinf'})} = 'ubuw'
                options.FrZeroLimit (1,1) {mustBeNumericOrLogical} = false;
                options.constantFr (1,1) {mustBeNumericOrLogical} = true;
            end

            if beta_2 == 0 % If 0 blockage requested, just predict unconfined
                [conf_2, conf_1] = hb.predictUnconfined(conf_1, uGuess, guessMode=options.guessMode, ...
                                                        correctionType='bluff body', FrZeroLimit=options.FrZeroLimit);
            else
                % Solve linear actuator disk
                conf_1 = hb.solveLMAD(conf_1, uGuess, guessMode=options.guessMode, FrZeroLimit=options.FrZeroLimit);
    
                % Now, assume that thrust on the rotor is the same at beta_1
                % and beta_2. If thrust is driven by the bypass velocity, then
                % the bypass velocity must be the same between the two cases as
                % well.
    
                for i = 1:size(conf_1, 1)
                    for j = 1:size(conf_1, 2)
    
                        % Generate UinfGuess based on beta1 and beta2
                        % Uinf_2Guess = conf_1(i,j).beta ./ beta_2 .* conf_1(i,j).Uinf;
                        Uinf_2Guess = conf_1(i,j).Uinf; % Try to constrain Uinf within reasonable values
    
                        [Uinf_2, Uinf_2Err, Uinf_2ExitFlag] = hb.convergeUinf2(Uinf_2Guess, beta_2, conf_1(i,j), options.constantFr);
    
                        % Re-scale data using converged velocity
                        currForecast = hb.convertConfToUnconf(conf_1(i,j), Uinf_2);
    
                        % Set h, beta to what was used.
                        if options.constantFr
                            % Calculate h_2 from Fr_1, Uinf_2
                            currForecast.h = hb.calcDepthFromFroude(currForecast.Uinf, conf_1(i,j).Fr);
                        else 
                            % Recalculate scaled depth that was used in iteration
                            currForecast.h = hb.calcDepthFromBlockage(conf_1(i,j).h, conf_1(i,j).beta, beta_2);
                        end
                        currForecast.Fr = hb.calcFroude(currForecast.Uinf, currForecast.h);
                        currForecast.beta = beta_2;
    
                        % Package diagnostics about Uinf_2 iteration
                        [uw_Fr, uw_Thrust] = hb.solveUw_both(conf_1(i,j).ub, beta_2, currForecast.CT, currForecast.Uinf, currForecast.Fr);
                        currForecast.Uinf_2Iter = hb.packageDiagnostics(currForecast.Uinf, [uw_Fr, uw_Thrust], ...
                                                                                      Uinf_2Err, Uinf_2ExitFlag);
    
                        % Send back through LMAD to get uw, ub
                        % NOTE: This should result in the same uw, ub as conf_1
                        currForecast = hb.solveLMAD(currForecast, uGuess, guessMode=options.guessMode);
    
                        % Save
                        conf_2(i,j) = currForecast;
                    end
                end
            end
        end

        function [conf_2, conf_1] = linearForecast(hb, conf_1, beta_2, uGuess, options) 
            % Performs a linear prediction of confined performance using
            % the method described by Kinsey and Dumas (2017;
            % 10.1016/j.renene.2016.11.021) and open-channel linear
            % momentum. Performance metrics are assumed to be linear with
            % blockage, such that performance at beta_2 is obtained by
            % linear interpolation/extrapolation between the provided data
            % at beta_1 and the zero-blockage condition obtained using an
            % open-channel blockage correction.
            % Inputs (required)
            %   conf_1     - A structure of confined performance data at blockage 1 with fields
            %                as described in the BWClosedChannel class documentation.
            %   beta_2     - The target blockage at which forecasted
            %                performance data is desired, as a fraction.
            % Inputs (optional, positional)
            %   uGuess    - Guess for the value one of the linear momentum
            %               velocities (or ratio between velocities), as a
            %               scalar (default 1.4). See documentation for 
            %               'guessMode' input for valid guess types.
            % Inputs (optional, name-value pairs)
            %   guessMode - The linear momentum velocities (or ratio
            %               between velocities) represented by uGuess:
            %       'ubuw': ratio between ub and uw (default)
            %       'ubUinf': ratio between ub and Uinf
            %       'ub': the bypass velocity
            %   FrZeroLimit - Whether to evaluate linear momentum at the
            %                 closed-channel limit (Fr=0), ignoring the
            %                 provided Froude number. Default: false.
            %   constantFr - Whether to hold the Froude number constant
            %                between beta_1 and beta_2 (true, default) or allow to
            %                vary (false). If constantFr=true, h at beta_2
            %                is calculated using Fr at beta_1 and Uinf at
            %                beta_2. If constantFr=false, h at beta_2 is
            %                calculated using h at beta_1 and assuming
            %                constant channel width.
            %   correctionType    - Type of blockage correction to apply: 
            %       "standard":     scales the confined data by the unconfined freestream velocity (UinfPrime, default)
            %       "bluff body":   scales the confined data by the bypass velocity (ub)
            % Outputs
            %   conf_2     - A structure with the same size and fields as
            %                output conf_1, but with each field corresponding to
            %                the forecasted performance at beta_2. 
            %   conf_1     - The input structure with the following fields
            %                added:
            %                Fr       - Depth-based Froude number (rectangular channel assumed)
            %                dhToh    - Normalized free-surface drop across the rotor
            %                uw       - Wake velocity estimated from closed-channel LMAD
            %                ub       - Bypass velocity estimated from closed-channel LMAD
            %                ut       - Velocity at the turbine estimated from closed-channel LMAD
            %                UinfPrime  - Unconfined freestream velocity estimated from closed-channel LMAD
            %                ubuwIter - Iteration diagnostics for ub/uw
            %                ubIter   - Iteration diagnostics for ub
            %                isPhys   - Results of physical validity checks on uw, ub, ut
            %
            % See also: HoulsbyOpenChannel, solveLMAD, predictUnconfined, checkPhysicalValidity
            arguments
                hb
                conf_1
                beta_2 (1,1) {mustBeInRange(beta_2, 0, 1)}
                uGuess (1,1) = [1.4];
                options.guessMode {mustBeText, mustBeMember(options.guessMode, {'ub', 'ubuw', 'ubUinf'})} = 'ubuw'
                options.correctionType {mustBeText, ismember(options.correctionType, {'standard', 'bluff body'})} = 'bluff body'
                options.FrZeroLimit (1,1) {mustBeNumericOrLogical} = false;
            end

            % Perform blockage correction
            unconf = hb.predictUnconfined(conf_1, uGuess, guessMode=options.guessMode, ...
                                                          correctionType=options.correctionType, ...
                                                          FrZeroLimit=options.FrZeroLimit);

            % Perform linear fits at each point
            conf_2 = struct;
            fitFields = {'TSR', 'CP', 'CT', 'CL', 'CQ', 'CF'};
            for i = 1:size(conf_1, 1) % For each row of conf_1
                for j = 1:size(conf_1, 2) % For each column of conf_1
                    % Calculate ratio between beta_2 and beta_1
                    betaRatio = (beta_2 ./ conf_1(i,j).beta);

                    % Set blockage in prediction structure
                    conf_2(i,j).beta = beta_2 .* ones(size(conf_1(i,j).beta));

                    for k = 1:length(fitFields)
                        % Calculate difference in performance between
                        % confined and unconfined performance
                        currData = conf_1(i,j).(fitFields{k});
                        currDataPrime = unconf(i,j).(fitFields{k});

                        % Perform linear interpolation/extrapolation
                        conf_2(i,j).(fitFields{k}) = (currData - currDataPrime) .* betaRatio + currDataPrime;
                    end

                end
            end

        end
    end

    methods (Static, Access = public)

        %% Core Equation set: Ross and Polagye EQs 43-45

        % Ross and Polagye Equation 43
        function [uw, num, den] = solveUw_Froude(ub, beta, CT, Uinf, Fr)
            % Ross and Polagye Equation 43: solves for uw, the core wake
            % velocity, using a blockage-Froude relationship.
            % Inputs:
            %   ub - the bypass velocity (m/s)
            %   beta - the channel blockage (decimal)
            %   CT - thrust coefficient at confined condition
            %   Uinf - undisturbed freestream velocity at the confined condition (m/s)
            %   Fr - depth-based Froude number at the confined condition
            % Ouputs:
            %   uw - The core wake velocity (m/s)
            %   num - Numerator of equation 43 (for convergence diagnostics)
            %   den - Denominator of equation 43 (for convergence diagnostics)
            num = (Fr.^2 .* ub.^4) - (4 + 2.*Fr.^2).*(Uinf.^2 .* ub.^2) + (8.*Uinf.^3 .* ub) - (4.*Uinf.^4) + (4.* beta .* CT .* Uinf.^4) + (Fr.^2 .* Uinf.^4);
            den = (-4 .* Fr.^2 .* ub.^3) + (4.*Fr.^2 + 8).*(Uinf.^2 .* ub) - (8.*Uinf.^3);
            uw = num ./ den;
        end

        % Ross and Polagye Equation 44:
        function uw = solveUw_Thrust(ub, CT, Uinf)
            % Ross and Polagye Equation 44: solves for uw, the core wake
            % velocity, using a thrust relation. Note that this is the same as
            % Ross and Polagye Equation 23, but rearranged.
            % Inputs:
            %   ub - the bypass velocity (m/s)
            %   CT - thrust coefficient at confined condition
            %   Uinf - undisturbed freestream velocity at the confined condition (m/s)
            % Ouputs:
            %   uw - The core wake velocity (m/s)
            uw = HoulsbyOpenChannel.solveUw_direct(ub, CT, Uinf);
        end
        
        % Ross and Polagye Equation 45:
        function ut = solveUt(uw, ub, h, Uinf, beta)
            % Ross and Polagye Equation 45: Solves for the velocity at the
            % turbine at the confined condition.
            % Inputs:
            %   uw - the core wake velocity (m/s)
            %   ub - the bypass velocity (m/s)
            %   h - undisturbed dynamic depth at the confined condition (m)
            %   Uinf - undisturbed freestream velocity at the confined condition (m/s)
            %   beta - the channel blockage (decimal)
            % Ouputs:
            %   ut - The velocity at the turbine at the confined condition (m/s)
            num = uw .* (ub - Uinf) .* (2.*HoulsbyOpenChannel.g.*h - ub.^2 - ub.*Uinf);
            den = 2 .* beta .* HoulsbyOpenChannel.g .* h .* (ub - uw);
            ut = num ./ den;
        end

        %% Froude number and dynamic depth calculations

        function Fr = calcFroude(Uinf, h)
            % Calculates the depth based Froude number from the input channel
            % depth (h) and freestream velocity (Uinf). A rectangular
            % channel is assumed.
            Fr = Uinf ./ sqrt(h .* HoulsbyOpenChannel.g);
        end
    
        function [h] = calcDepthFromFroude(UinfTarget, FrTarget)
            % Given a target Froude number and freestream velocity, calculates
            % the corresponding depth required to meet that condition.
            h = UinfTarget.^2 ./ FrTarget.^2 ./ HoulsbyOpenChannel.g;
        end

        function [h_2] = calcDepthFromBlockage(h_1, beta_1, beta_2)
            % Given two blockage states and the depth at one blockage state,
            % computes the depth of the other blockage state assuming constant
            % channel width and turbine area.
            h_2 = beta_1 ./ beta_2 .* h_1;
        end

        %% Linear mometum channel depth calculations

        function [dhToh] = calcTotalSurfaceDeformation(CT, beta, Fr)
            % Calculates the normalized free surface drop (dh/h) across
            % turbine rotor using Houlsby et al Equation 4e in Section 9.4,
            % given input CT, beta, and Fr vectors.
            dhToh = zeros(size(CT));
            for i = 1:length(CT)
                dhToh(i) = fzero(@(x) surfaceDefFun(x, CT(i), beta(i), Fr(i)), 0);
            end

            function err = surfaceDefFun(dhToh, CT, beta, Fr)
                % Calculates the normalized free surface drop (dh/h) across
                % turbine rotor using Houlsby et al Equation 4e in Section 9.4,
                % given a known CT, blockage ratio, and Froude number.
                err = 1/2*(dhToh.^3) - 3/2*(dhToh.^2) + (1 - Fr.^2 + (CT.*beta.*Fr.^2)/2)*(dhToh) - (CT .* beta .* Fr.^2)./2;
            end
        end

        function [hdUp] = calcUpstreamDiskDepth(h, Uinf, ut)
            % Calculates the depth in the core flow just upstream of the
            % actuator disk via Bernoulli, given the upstream undisturbed
            % upstream depth (h), undisturbed inflow velocity (Uinf), and
            % velocity at the turbine (ut).
            hdUp = h + 1/(2*HoulsbyOpenChannel.g) .* (Uinf.^2 - ut.^2);
        end

        function [hdDown] = calcDownstreamDiskDepth(h4, ut, uw)
            % Calculates the depth in the core flow just upstream of the
            % actuator disk via Bernoulli, given the depth in the bypass
            % (h4), velocity at the turbine (ut), and velocity in core wake (uw).
            hdDown = h4 + 1/(2*HoulsbyOpenChannel.g) .* (uw.^2 - ut.^2);
        end

        function [dhDisk] = calcDiskDrop(Uinf, CT)
            % Calculates the free surface drop across the disk given the
            % velcocity and thrust coefficient
            dhDisk = 1/(2*HoulsbyOpenChannel.g) .* CT .* Uinf.^2;
        end

        function [hBypass] = calcBypassDepth(h, Uinf, ub)
            % Calculates the depth in the bypass flow/core wake downstream
            % of the actuator disk via Bernoulli, given the undisturbed
            % upstream depth (h), velocity at the turbine (ut), and bypass
            % velocity (ub).
            hBypass = h + 1/(2*HoulsbyOpenChannel.g) .* (Uinf.^2 - ub.^2);
        end

        % function [conf] = calcLMADDepths(conf)
        %     for i = 1:size(conf,1)
        %         for j = 1:size(conf,2)
        %             % Calculate all depths
        %             conf(i,j).hBypass = calcBypassDepth(conf(i,j).h, conf(i,j).Uinf, conf(i,j).ub);
        %             conf(i,j).hdUp = calcUpstreamDiskDepth(conf(i,j).h, conf(i,j).Uinf, conf(i,j).ut);
        %             conf(i,j).hdDown = calcUpstreamDiskDepth(conf(i,j).hBypass, conf(i,j).ut, conf(i,j).uw);
        % 
        %             % Calculate drop across the disk predicted by thrust
        %             conf(i,j).dhDisk = calcDiskDrop(conf(i,j).Uinf, conf(i,j).CT);
        % 
        %             % Calculate total free surface deformation
        %             conf(i,j).dhToh = hb.calcTotalSurfaceDeformation(conf(i,j).CT, conf(i,j).beta, conf(i,j).Fr);
        %         end
        %     end
        % 
        % end

        function conf = setClosedChannelDepth(conf)
            % Overrides the depth in each element of conf to 1e6 to ensure
            % a very low Froude number, such that the limiting case of
            % Fr->0 can be approximated.
            depthVal = 1e6;
            fprintf('Setting depth to %g for open-channel limiting case of Fr -> 0 (closed-channel).\n', depthVal);
            for i = 1:numel(conf)
                conf(i).h = depthVal .* ones(size(conf(i).h));
            end
        end

        %% Iteration scheme for solving for ub using Ross and Polagye equations 43 and 44
        % CT, Uinf, Fr are known
        % ub, uw are unknown

        function [ub, err, exitFlag] = convergeUb(ubGuess, conf)
            % Using Ross and Polagye Equations 43-44, iterates to find ub
            % satisfies both EQ 43 and EQ 44, and returns that ub
            % Inputs:
            %   ubGuess - Initial guess for ub (as a two element vector for fzero)
            %   conf    - Confined performance data as described in the
            %             HoulsbyOpenChannel class documentation.
            % Outputs:
            %   ub      - Converged values of bypass velocity (m/s)
            %   err     - Error between EQ 43 and EQ 44 at convergence
            %   exitFlag - fzero exit condition

            % Preallocate
            nPoints = length(conf.TSR);
            ub = zeros(size(conf.TSR));
            err = zeros(size(conf.TSR));
            exitFlag = zeros(size(conf.TSR));

            % Iterate for each point
            for k = 1:nPoints
                % If good point, proceed
                if ~isnan(ubGuess(k)) && (conf.CT(k) >= 0)
                    currFun = @(ub) HoulsbyOpenChannel.ubCompare(ub, conf.beta(k), conf.CT(k), conf.Uinf(k), conf.Fr(k));

                    % [ub(k), err(k), exitFlag(k)] = fzero(currFun, ubGuess(k,:));
                    [ub(k), err(k), exitFlag(k)] = fminsearch(currFun, ubGuess(k,:));
                else
                    warning('Negative CT value or bad ubGuess: skipping application of LMAD for this point');
                    ub(k) = nan;
                    err(k) = nan;
                    exitFlag(k) = nan;
                end
            end
        end

        function [uw_Fr, uw_Thrust] = solveUw_both(ub, beta, CT, Uinf, Fr)
            % Calculates uw, the core wake velocity, using Ross and Polagye
            % EQs 43 and 44, and returns the values from each equation
            % Inputs:
            %   ubGuess - Bypass velocity (m/s)
            %   beta - channel blockage ratio
            %   CT - thrust coefficient
            %   Uinf - undisturbed freestream velocity (m/s)
            %   Fr - depth-based Froude number
            % Outputs:
            %   uw_Fr - uw calculated via Ross and Polagye equation 43
            %   uw_Thrust - uw calculated via Ross and Polagye equation 44
            uw_Fr = HoulsbyOpenChannel.solveUw_Froude(ub, beta, CT, Uinf, Fr);
            uw_Thrust = HoulsbyOpenChannel.solveUw_Thrust(ub, CT, Uinf);
        end

        function err = ubCompare(ubGuess, beta, CT, Uinf, Fr)
            % Calculates uw, the core wake velocity, using Ross and Polagye
            % EQs 43 and 44, and returns the error between the values yielded
            % by each method.
            % Inputs:
            %   ubGuess - Bypass velocity (m/s)
            %   beta - channel blockage ratio
            %   CT - thrust coefficient
            %   Uinf - undisturbed freestream velocity (m/s)
            %   Fr - depth-based Froude number
            % Outputs:
            %   err - Error between uw calculated via Ross and Polagye EQ 43 and
            %         uw calculated via Ross and Polagye EQ 44
            %
            % ## NOTE: Real part of solution is used to assess convergence
            %          to avoid issues with fminsearch and complex values.

            % Check if physical. If not, make error large
            if ubGuess / Uinf <= 1
                err = 1e6;
            else
                % Solve for uw both ways
                [uw_Fr, uw_Thrust] = HoulsbyOpenChannel.solveUw_both(ubGuess, beta, CT, Uinf, Fr);
    
                % Error for fzero:
                % err = real(uw_Fr - uw_Thrust); % Take only error of real parts to nudge away from complex solutions
    
                % Error for fminsearch
                err = abs(real(uw_Fr - uw_Thrust));
    
                % Compute error between values
                % err = uw_1 - uw_2; % Compute error between those values
                % err = abs(uw_1) - abs(uw_2); % Take error of magnitudes
            end
        end

        % Visualization
        function [uwErr, ax] = plotUbConvergenceRegion(ubTest, beta, CT, Uinf, Fr)
            % Plots the error between the two methods for calculating uw (used
            % in ub iteration) for a specific beta, CT, Uinf, and Fr Useful for
            % visually determining whether a given point can ever converge.

            % Generate test values for ub and evaluate uw error for each
            uwErr = HoulsbyOpenChannel.ubCompare(ubTest, beta, CT, Uinf, Fr);

            % Plot
            [fig] = figure();
            ax = axes(fig);
            grid(ax, 'on'); hold(ax, 'on');
            plot(ax, ubTest, zeros(size(ubTest)), '-k');
            plot(ax, ubTest, uwErr, '-', 'marker', '.');
            xlabel(ax, '$u_2$');
            ylabel(ax, '$u_1$ error');
        end

        %% Iteration scheme for solving for uw with known ub, unknown Uinf
        function [Uinf_2, err, exitFlag] = convergeUinf2(Uinf_2Guess, beta_2, conf_1, constantFr)
            % Using Ross and Polagye Equations 43-44 iterates to find Uinf 
            % that satisfies both EQ 43 and EQ 44, and returns that Uinf.
            % Used to forecast performance from blockage 1 to blockage 2.
            % Inputs:
            %   Uinf_2Guess  - A guess for the freestream velocity at beta_2
            %   beta_2     - Target blockage for forecasting (scalar, decimal)
            %   conf_1     - Performance data at starting blockage (beta_1)
            %                formatted as described in HoulsbyOpenChannel
            %                class documentation
            %   constantFr - Whether to hold the Froude number constant
            %                between beta_1 and beta_2 (true) or allow to
            %                vary (false). If constantFr=true, h at beta_2
            %                is calculated using Fr at beta_1 and Uinf at
            %                beta_2. If constantFr=false, h at beta_2 is
            %                calculated using h at beta_1 and assuming
            %                constant channel width.
            % Outputs:
            %   Uinf_2     - Converged value of freestream velocity at beta_2 (m/s)
            %   err      - Error between EQ 43 and EQ 44 at convergence
            %   exitFlag - fzero exit condition

            % Preallocate
            nPoints = length(conf_1.TSR);
            Uinf_2 = zeros(size(conf_1.TSR));
            err = zeros(size(conf_1.TSR));
            exitFlag = zeros(size(conf_1.TSR));

            % Compute h_2guess based on Froude number assumption
            if constantFr
                % Let h_2 vary to hold Fr_1 = Fr_2 as Uinf_2 varies
                h_2Guess = NaN .* ones(size(conf_1.h));
            else 
                % Fix h_2 based on target blockage, intentionally let Fr vary with Uinf_2
                h_2Guess = HoulsbyOpenChannel.calcDepthFromBlockage(conf_1.h, conf_1.beta, beta_2);
            end

            % Set options
            options = optimset('TolX', 1e-12);

            % Iterate for each point
            for k = 1:nPoints
                currFun = @(Uinf_2) HoulsbyOpenChannel.Uinf2Compare(Uinf_2, ...
                                                                h_2Guess(k), ...
                                                                beta_2, ...
                                                                conf_1.ub(k), ...
                                                                conf_1.CT(k), ...
                                                                conf_1.Uinf(k), ...
                                                                conf_1.Fr(k));

                % [Uinf_2(k), err(k), exitFlag] = fzero(currFun, Uinf_2Guess(k,:), options);
                [Uinf_2(k), err(k), exitFlag] = fminsearch(currFun, Uinf_2Guess(k,:), options);
                % [Uinf_2(k), err(k), exitFlag] = fminbnd(currFun, Uinf_2Guess(k,:), options);
            end
        end


        function err = Uinf2Compare(Uinf_2Guess, h_2Guess, beta_2, ub, CT_1, Uinf_1, Fr_1)
            % Calculates Uinf_2, the freestream velocity at blockage 2, using
            % Ross and Polagye EQs 43 and 44, and returns the error between
            % the values yielded by each method.
            % Inputs:
            %   Uinf_2Guess - Guess for freestream velocity at blockage 2 (m/s)
            %   h_2Guess - Guess for water depth at blockage 2 (m)
            %   beta_2    - channel blockage ratio at blockage 2
            %   CT_1      - thrust coefficient at blockage 1
            %   Uinf_1      - undisturbed freestream velocity at blockage 1 (m/s)
            %   Fr_1      - depth-based Froude number at blockage 1
            % Outputs:
            % err - Error between Uinf_2 calculated via Ross and Polagye EQ 43 and
            %       Uinf_2 calculated via Ross and Polagye EQ 44
            %
            % ## NOTE: Real part of solution is used to assess convergence
            %          to avoid issues with fminsearch and complex values.
            
            
            % Check that this guess is physical (bypass faster than
            % freestream). If not, give larger error to nudge away from
            % this spot.
            if (ub / Uinf_2Guess) <= 1
                err = 1e6; % Give large error to move away from this spot
            else

                % Compute CT
                CT_2 = HoulsbyOpenChannel.scaleForcingMetric(CT_1, Uinf_1, Uinf_2Guess);
    
                % Check Fr case
                if ~isnan(h_2Guess) 
                    % h_2 is given, then calculate Fr_2 from Uinf_2 and h_2
                    Fr_2 = HoulsbyOpenChannel.calcFroude(Uinf_2Guess, h_2Guess);
                else 
                    % Otherwise, assume Fr_2 = Fr_1
                    Fr_2 = Fr_1;
                end
                
                % Compute uw from equation for uw
                [uw_Fr, uw_Thrust] = HoulsbyOpenChannel.solveUw_both(ub, beta_2, CT_2, Uinf_2Guess, Fr_2);

                % Once again, check if physical (freestream faster than
                % wake). If not, assign large error.
                % if (uw_Fr / Uinf_2Guess >= 1)
                %     err = 1e6;
                % else
                % Error for fzero
                %err = real(uw_Fr - uw_Thrust); % Error approach for fzero. Take only error of real parts to nudge away from complex solutions
    
                % Error for fminsearch
                err = abs(uw_Fr - uw_Thrust); % Error approach for fminsearch
    
                % Alternative errors
                % err = uw_1 - uw_2; % Compute error between those values
                % err = abs(uw_1) - abs(uw_2); % Take error of magnitudes
                % end
            end
        end


        function [diag] = assessUinf2ConvergenceRegion(conf, ubUinfTest, beta_2)
            % Assesses convergence region for Uinf_2 in analytical blockage
            % forecasting by decomposing Ross and Polagye Equations 43 and
            % 44. Returns a structure of various intermediate values that
            % arise from the evaluation of these equations, which may be
            % used to assess conditions under which convergence succeeds or
            % fails.
            %
            % Inputs 
            %   conf     - Structure of confined performance data as 
            %              described in the class documentation
            %   ubUinfTest - The value of the ratio between ub and Uinf to
            %              evaluate the equations at
            %   beta_2   -
            % Outputs
            %   diag     - A structure of various intermediate values that
            %              arise from the evaluation of Equations 43 and 44,
            %              which may be used to assess conditions under
            %              which convergence succeeds or fails.

            diag = struct;
            for i = 1:size(conf,1)
                for j = 1:size(conf,2)
                    % Compute Uinf_2 test values from ubUinfTest
                    Uinf_2Test = conf(i,j).ub ./ ubUinfTest;

                    % Determine depth to hold Fr constant in freestream
                    h_2Test = HoulsbyOpenChannel.calcDepthFromFroude(Uinf_2Test, conf(i,j).Fr);

                    % Calculate corresponding Froude number in bypass
                    diag(i,j).FrBypass = conf(i,j).ub ./ sqrt(9.81 * h_2Test);

                    % Preallocate
                    nPoints = length(conf(i,j).CT);
                    diag(i,j).beta_2 = beta_2;
                    diag(i,j).Fr_2 = conf(i,j).Fr;
                    diag(i,j).Uinf_2 = Uinf_2Test;
                    diag(i,j).uw_Fr = zeros(nPoints, length(ubUinfTest));
                    diag(i,j).uw_Fr_num = zeros(nPoints, length(ubUinfTest));
                    diag(i,j).uw_Fr_den = zeros(nPoints, length(ubUinfTest));
                    diag(i,j).uw_Thrust = zeros(nPoints, length(ubUinfTest));
                    diag(i,j).CTBetaTerm = zeros(nPoints, length(ubUinfTest));
                    diag(i,j).CT = zeros(nPoints, length(ubUinfTest));

                    for k = 1:length(conf(i,j).CT) % For each point

                        % Get CT_2
                        CT_2 = HoulsbyOpenChannel.scaleForcingMetric(conf(i,j).CT(k), conf(i,j).Uinf(k), Uinf_2Test(k,:));

                        % Solve both equations for wake velocity
                        [uw_Fr, uw_Fr_num, uw_Fr_den] = HoulsbyOpenChannel.solveUw_Froude(conf(i,j).ub(k), beta_2, CT_2, Uinf_2Test(k,:), conf(i,j).Fr(k));
                        [uw_Thrust] = HoulsbyOpenChannel.solveUw_Thrust(conf(i,j).ub(k), CT_2, Uinf_2Test(k,:));

                        % Save results
                        diag(i,j).uw_Fr(k,:) = uw_Fr;
                        diag(i,j).uw_Fr_num(k,:) = uw_Fr_num;
                        diag(i,j).uw_Fr_den(k,:) = uw_Fr_den;
                        diag(i,j).uw_Thrust(k,:) = uw_Thrust;

                        % Also compute beta term
                        diag(i,j).CTBetaTerm(k,:) = 4 .* CT_2 .* beta_2 .* Uinf_2Test(k,:).^4;
                        diag(i,j).CT(k,:) = CT_2;
                    end
                end
            end
        end


        %% Printing stuff out
        function modelName = getModelName()
            % Returns a label-friendly version of the blockage correction
            % model name.
            modelName = 'Houlsby Open-Channel';
        end
    end
end