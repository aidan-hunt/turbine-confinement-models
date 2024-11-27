% Aidan Hunt
%
% Implements a two-scale closed-channel linear momentum model for turbine
% arrays as described by Nishino and Willden in "The efficiency of an array
% of tidal turbines partially blocking a wide channel" (2012;
% https://doi.org/10.1017/jfm.2012.349) and further expanded upon by
% Dehtyriov et al in "Fractal-like actuator disk theory for optimal energy
% extraction" (2021; https://doi.org/10.1017/jfm.2021.766). This
% implementation follows that of Dehtyriov et al. in "A two-scale blockage
% correction for an array of tidal turbines" (2023;
% https://doi.org/10.36688/ewtec-2023-366), who provided a procedure for
% estimating induction factors by using the measured blockages, thrust
% coefficients, and freestream velocities as inputs to the two-scale model
% (Figure 2 in Dehtyriov et al 2023). The resulting velocities may be used
% to perform a blockage correction.
%
% IMPORTANT NOTE: This code uses the convention of Dehtyriov et al 2021 for
% induction factors (alpha = ut/Uinf), which is different than that of
% Nishino and Willden 2012 and Dehtyriov et al 2023 (alpha = 1 - ut/Uinf).
% As we are following the implementation of Dehtyriov et al 2023, the
% appropriate equations were updated to match the former convention above.
%
% To use the NWTwoScale class, construct a NWTwoScale object
% using the following syntax:
%   nw = NWTwoScale()
% and call methods using the dot notation (i.e., nw.predictUnconfined(...)).
% 
% NWTwoScale Methods:
%   predictUnconfined - Uses a two-scale blockage correction to predict 
%                       unconfined performance from confined performance
%                       data.
%
% The methods above expect that confined performance data is provided as an
% mxn structure array, conf, with the following fields:
%   Uinf (required)   - undisturbed upstream freestream velocity (m/s)
%   CT   (required)   - thrust coefficient
%   CP   (optional)   - performance coefficient
%   CQ   (optional)   - torque coefficient
%   CL   (optional)   - lateral force coefficient
%   CF   (optional)   - resultant force coefficient
%   TSR  (optional)   - tip-speed ratio
%   h    (3D only)    - channel depth
% 
% The fields of conf(i,j) must be vectors that are all the same size. You
% can use the methods above to apply corrections and forecasts to multiple
% datasets at once by specifying each dataset as an element of conf (e.g.,
% conf(i,j)).
%
% See also: BCBase, BWClosedChannel, HoulsbyOpenChannel

classdef NWTwoScale < BCBase
    % The NWTwoScale class inherits all properties of the
    % BCBase superclass, and defines additional properties
    properties (Constant, Access=protected)
        correctionModes = {'standard', 'bluff body (array)', 'bluff body (device)'}
    end

    methods (Access = public)
        %% Main method for blockage correction
        function conf = solveLMAD(nw, conf, geom, gammaGuess, mode)
            % Applies the two-scale LMADT model of Nishino and Willden 2012
            % following the implementation of Dehtyriov et al 2023 (Figure
            % 2) to compute the turbine, wake, and bypass velocities at
            % each scale.
            % 
            % Inputs (required)
            %   conf       - A structure of confined performance data with fields
            %                as described in the NWTwoScale class
            %                documentation.
            %   geom       - A structure of array and channel geometry
            %                information with the following fields:
            %                devArea  - Turbine projected area [m^2]
            %                devWidth - Turbine width [m]
            %                n        - Number of turbines in the array
            %                w        -  Channel width [m]
            %                s        - Spacing between turbines [m]
            % 
            %   gammaGuess - An initial guess for the value of the local
            %                and array-scale wake velocity induction
            %                factor. A single guess is used to iterate on
            %                both values.
            %   mode       - Geometry evaluation mode: 2D or 3D. If 2D,
            %                device area is treated as device width, and
            %                channel depth is ignored. Default: 3D.
            % Outputs
            %   conf  - The input structure, but with the following
            %           fields added:
            %       UinfPrime - Unconfined freestream velocity [m/s] for
            %                   use in blockage correction.
            %       array     - A structure of flow quantities, thrust, and
            %                   blockage at the array scale.
            %       device    - A structure of flow quantities, thrust, and
            %                   blockage at the device scale. 
            %       The "array" and "device" substructures each have the
            %       following fields:
            %           beta  - Blockage ratio at this scale
            %           CT    - Thrust coefficient at this scale
            %           Uinf  - Freestream velocity at this scale
            %           ut    - Velocity through the turbine at this scale
            %           uw    - Wake velocity at this scale
            %           ub    - Bypass velocity at this scale
            %           alpha - Turbine induction factor at this scale
            %           gamma - Wake induction factor at this scale
            %           gammaIter - Information of the gamma iteration
            %                       procedure conducted at this scale

            arguments
                nw
                conf
                geom
                gammaGuess (1,1) = 0.8;
                mode {ismember(mode, {'2D', '3D'})} = '3D'
            end

            % Check input for correct sizing
            conf = nw.checkInputSizes(conf);

            % Check if scalar struct for geometry provided. If so, repmat
            if isequal(size(geom), [1,1])
                geom = repmat(geom, size(conf,1), size(conf,2));
            elseif all(size(conf) ~= size(geom))
                error('Size of conf and geom structures must match, or geom must be a 1x1 structure that applies to all entries of conf.');
            end

            % Loop through all cases
            for i = 1:size(conf, 1)
                for j = 1:size(conf, 2)
                    % Initialize structures for each scale
                    device = nw.initializeScaleStruct();
                    array = nw.initializeScaleStruct();

                    % Calculate blockage scales from the provided geometry data
                    switch mode
                        case '2D'
                            device.beta = nw.calcLocalBlockage(geom(i,j).devWidth, geom(i,j).devWidth, geom(i,j).s, ones(size(conf(i,j).CT)));
                            array.beta = nw.calcArrayBlockage(geom(i,j).devWidth, geom(i,j).n, geom(i,j).s, geom(i,j).w, ones(size(conf(i,j).CT)));
                        case '3D'
                            device.beta = nw.calcLocalBlockage(geom(i,j).devArea, geom(i,j).devWidth, geom(i,j).s, conf(i,j).h);
                            array.beta = nw.calcArrayBlockage(geom(i,j).devWidth, geom(i,j).n, geom(i,j).s, geom(i,j).w, conf(i,j).h);
                    end

                    % Check for BA = 1 (singularity)
                    if any(array.beta == 1)
                        warning('Array blockage = 1 (equally spaced turbines spanning entire channel) results in non-physical results. Please use "BWClosedChannel" or "HoulsbyOpenChannel" for these cases instead.')
                    end

                    % Iterate to solve for array-scale induction factors
                    % (Dehtyriov et al 2023 equations 28 and 9)
                    [array.gamma, gammaAErr, gammaAExit] = nw.convergeGammaA(gammaGuess, conf(i,j).CT, array.beta, device.beta);
                    array.alpha = nw.calcAlpha(array.gamma, array.beta);
            
                    % Solve for device-scale thrust
                    % (Dehtyriov et al 2023 equation 5) 
                    device.CT = nw.calcLocalThrustFromGlobal(conf(i,j).CT, array.alpha);
            
                    % Iterate to solve for device-scale induction factors
                    % (Dehtyriov et al 2023 equations 6 and 8)
                    [device.gamma, gammaLErr, gammaLExit] = nw.convergeGammaL(gammaGuess, device.CT, device.beta);
                    device.alpha = nw.calcAlpha(device.gamma, device.beta);
            
                    % Solve for array-scale thrust
                    % (Dehtryiov et al 2023 equation 4)
                    array.CT = nw.calcArrayThrust(array.alpha, device.beta, device.CT);

                    % Solve for unconfined freestream velocity for blockage
                    % correction
                    % (Dehtyriov et al 2023 equation 24)
                    velRatio = nw.calcVelocityRatio(array.alpha, array.CT);
                    conf(i,j).UinfPrime = conf(i,j).Uinf ./ velRatio;

                    % Compute turbine, wake, and bypass velocities at each scale
                    array.Uinf = conf(i,j).Uinf;
                    [array.ut, array.uw, array.ub] = nw.calcVelocitiesFromInduction(array.Uinf, array.CT, array.alpha, array.gamma);
                    device.Uinf = array.ut;
                    [device.ut, device.uw, device.ub] = nw.calcVelocitiesFromInduction(device.Uinf, device.CT, device.alpha, device.gamma);

                    % Package iteration diagnostics at array scale
                    array.gammaIter = nw.packageDiagnostics(array.gamma, nw.calcGlobalThrust(array.gamma, array.alpha, array.beta, device.beta), ... 
                                                            gammaAErr, gammaAExit);

                    % Package iteration diagnostics at device scale
                    device.gammaIter = nw.packageDiagnostics(device.gamma, nw.calcThrustFromInduction(device.gamma, device.alpha, device.beta), ...
                                                             gammaLErr, gammaLExit);

                    % Finally, add array and device structures to conf
                    conf(i,j).array = array;
                    conf(i,j).device = device;      
                end
            end
        end

        function [unconf, conf] = predictUnconfined(nw, conf, geom, gammaGuess, mode, options)
            % Applies Dehytriov et al (2023)'s two-scale blockage
            % correction to the given array performance data and array
            % layout.
            % 
            % Inputs (required)
            %   conf       - A structure of confined performance data with fields
            %                as described in the NWTwoScale class
            %                documentation.
            %   geom       - A structure of array and channel geometry
            %                information with the following fields:
            %                devArea  - Turbine projected area [m^2]
            %                devWidth - Turbine width [m]
            %                n        - Number of turbines in the array
            %                w        -  Channel width [m]
            %                s        - Spacing between turbines [m]
            % 
            %   gammaGuess - An initial guess for the value of the local
            %                and array-scale wake velocity induction
            %                factor. A single guess is used to iterate on
            %                both values.
            %   mode       - Geometry evaluation mode: 2D or 3D. If 2D,
            %                device area is treated as device width, and
            %                channel depth is ignored. Default: 3D.
            % Outputs
            %   unconf    - A structure with the same size as the input
            %               conf structure and the following fields (if the
            %               field was not present in conf, it will not be
            %               present in unconf):
            %               Uinf  - unconfined freestream velocity
            %               CT  - unconfined thrust coefficient
            %               CP  - unconfined performance coefficient
            %               CQ  - unconfined torque coefficient
            %               CL  - unconfined lateral force coefficient
            %               CF  - unconfined resultant force coefficient
            %               TSR - unconfined tip-speed ratio
            %               velRatio   - Ratio between confined freestream
            %                            velocity and unconfined freestream velocity.
            %   conf  - The input structure, but with the following
            %           fields added from solving LMADT:
            %       UinfPrime - Unconfined freestream velocity [m/s] for
            %                   use in blockage correction.
            %       array     - A structure of flow quantities, thrust, and
            %                   blockage at the array scale.
            %       device    - A structure of flow quantities, thrust, and
            %                   blockage at the device scale. 
            %       The "array" and "device" substructures each have the
            %       following fields:
            %           beta  - Blockage ratio at this scale
            %           CT    - Thrust coefficient at this scale
            %           Uinf  - Freestream velocity at this scale
            %           ut    - Velocity through the turbine at this scale
            %           uw    - Wake velocity at this scale
            %           ub    - Bypass velocity at this scale
            %           alpha - Turbine induction factor at this scale
            %           gamma - Wake induction factor at this scale
            %           gammaIter - Information of the gamma iteration
            %                       procedure conducted at this scale
            arguments
                nw
                conf
                geom
                gammaGuess (1,1) = 0.8
                mode {ismember(mode, {'2D', '3D'})} = '3D'
                options.correctionType {mustBeText, ismember(options.correctionType, {'standard', 'bluff body (array)', 'bluff body (device)'})} = 'standard';
            end

            % Solve LMAD
            conf = nw.solveLMAD(conf, geom, gammaGuess, mode);

            % Convert to unconfined using the appropriate scaling velocity
            for i = 1:size(conf,1)
                for j = 1:size(conf,2)
                    switch options.correctionType
                        case 'standard'
                            scalingVel = conf(i,j).UinfPrime;
                        case 'bluff body (array)'
                            scalingVel = conf(i,j).array.ub;
                        case 'bluff body (device)'
                            scalingVel = conf(i,j).device.ub;
                    end
                    unconf(i,j) = nw.convertConfToUnconf(conf(i,j), scalingVel);
                end
            end

        end
    end
   

    methods (Static, Access = public)

        %% Calculating blockage ratios at different scales

        % Nishino and Willden Equation 2.2, adapted to arbitrary turbine
        % geometry.
        function BL = calcLocalBlockage(devArea, devWidth, s, h)
            % Nishino and Willden 2.2, but adapated for an arbitrary
            % turbine geometry.
            % Calculates the local scale blockage, defined as the device area
            % divided by the local passage cross-sectional area
            % Inputs:
                % devArea  - Turbine projected area [m^2]
                % devWidth - Device width [m]
                % h        - Channel depth [m]
                % s        - Spacing between turbines [m]
            % Outputs:
                % BL       - Local scale blockage ratio

            BL = devArea ./ (h .* (devWidth + s));
        end

        % Nishino and Willden Equation 2.3, adapted to arbitrary turbine
        % geometry.
        function BA = calcArrayBlockage(devWidth, n, s, w, h)
            % Nishino and Willden Equation 2.3, but adapated for cross-flow
            % turbines. Calculates the array scale blockage, defined as the
            % array area (including spacing between turbines) divided by
            % the channel area.
            % Inputs:
                % devArea  - Turbine projected area [m^2]
                % devWidth - Turbine width [m]
                % n        - Number of turbines in the array
                % h        - Channel depth [m]
                % w        - Channel width [m]
                % s        - Spacing between turbines [m]
            % Outputs:
                % BA - Array scale blockage ratio

            BA = (h .*n .* (devWidth+s)) ./ (h.*w);
        end

        %% Multi-scale equations

        % Dehtryiov et al Equations 8 and 9, rearranged to solve for alphaA
        % or alphaL directly
        function alpha = calcAlpha(gamma, B)
            % Dehtryiov et al Equations 8 and 9, rearranged to solve for alphaA
            % or alphaL directly
            % Solves for the induction factor at a particular scale, given the
            % wake velocity induction factor and blockage ratio at that scale.
            % Inputs:
            %   gamma - Velocity induction factor at the scale of interest
            %   B  - Blockage ratio at the scale of interest
            % Outputs:
            %   alpha - Velocity induction factor at the same scale as gamma
            %           and B
            num = (1 + gamma);
            den = (1 + B) + sqrt((1 - B).^2 + B.*(1 - 1./gamma).^2);
            alpha = num ./ den;
        end

        % Dehtryiov et al Equations 6 and 7
        function CT = calcThrustFromInduction(gamma, alpha, B)
            % Dehtryiov et al Equations 6 and 7
            % Solves for the thrust coefficient at a scale of interest using
            % induction factors
            % Inputs:
            %   gamma - Wake velocity induction factor at the scale of interest
            %   alpha - Velocity induction factor at the scale of interest
            %   B     - Blockage ratio at the scale of interest
            % Outputs:
            %   CTL    - Thrust coefficient at the scale of interest
            num = (1 + gamma) - 2.*B.*alpha;
            den = (1 - B.*alpha./gamma).^2;
            CT = (1 - gamma) .* (num./den);
        end

        % Computing velocities from induction factors
        % NOTE: alpha is defined here as ut/Uinf (matching Dehtyriov et al
        %       2021) in contrast to 1 - ut/Uinf (as found in Nishino and
        %       Willden 2012 and Dehtryiov et al 2023)
        % NOTE: While an expression for the bypass velocity is not included
        %       in Nishino and Willden 2012 or Dehtyriov et al 2023, this
        %       relationship is derived in multiple closed-channel and
        %       open-channel LMADT models (as a consequence of the pressure
        %       drop across the actuator disk) and is noted by Dehytriov et
        %       al 2021 to apply at any given scale.
        function [ut, uw, ub] = calcVelocitiesFromInduction(Uinf, CT, alpha, gamma)
            % Uses the freestream velocity, thrust coefficient, and
            % induction factors alpha and gamma to solve for the velocity
            % at the turbine (ut), wake velocity (uw) and bypass velocity
            % (ub) at a particular scale.

            % Calculate turbine and wake velocities directly
            ut = Uinf .* alpha; 
            uw = Uinf .* gamma;

            % Calculate bypass induction factor
            bypassInduction = sqrt(CT + gamma.^2);

            % Compute bypass velocity from this induction factor
            ub = Uinf .* bypassInduction;
        end

        %% Global-scale equations

        % Dehtryiov et al Equation 28
        function CTG = calcGlobalThrust(gammaA, alphaA, BA, BL)
            % Dehtryiov et al Equation 28
            % Solves for the global thrust coefficient:
            % Inputs:
            %   gammaA - Array-scale wake velocity induction factor
            %   alphaA - Array-scale velocity induction factor
            %   BA  - Array-scale blockage
            %   BL  - Local-scale blockage
            % Outputs:
            %   CTG    - Global thrust coefficient
            num = (1 + gammaA) - 2.*BA.*alphaA;
            den = (1 - BA.*alphaA./gammaA).^2;
            CTG = ((1 - gammaA)./BL) .* (num./den);
        end

        %% Array-Scale Equations

        % Dehtryiov et al Equation 4
        function CTA = calcArrayThrust(alphaA, BL, CTL)
            % Dehtryiov et al Equation 4
            % Solves for the array thrust coefficient 
            % Inputs:
            %   alphaA - Array-scale velocity induction factor
            %   BL  - Local blockage ratio
            %   CTL    - Local thrust coefficient
            % Outputs:
            %   CTA    - Array-scale thrust coefficient
            CTA = (alphaA).^2 .* BL .* CTL;
        end

        %% Local-Scale Equations

        % Dehtryiov et al Equation 5 rearranged to solve for CTL
        function CTL = calcLocalThrustFromGlobal(CTG, alphaA)
            % Dehtryiov et al Equation 5 rearranged to solve for CTL
            % Solves for the local thrust coefficient using the global thrust
            % coefficient and array-scale induction
            % Inputs:
            %   CTG    - Global-scale thrust coefficient
            %   alphaA - Array-scale velocity induction factor
            % Outputs:
            %   CTL    - Local thrust coefficient
            CTL = CTG ./ (alphaA.^2);
        end
        
        %% Velocity correction factor

        % Dehtryiov et al Equation 24
        function velRatio = calcVelocityRatio(alphaA, CTA)
            % Dehtryiov et al Equation 24
            % Solves for the velocity correction factor Uc/Uc'
            % Inputs:
            %   alphaA - array-scale induction factor
            %   CTA    - array-scale thrust coefficient
            % Outputs:
            %   velRatio - Velocity correction factor: U/U'
            num = 4.*(alphaA);
            den = 4.*(alphaA).^2 + CTA;
            velRatio = num ./ den;
        end

        %% Iteration on GammaA
        function err = compareCTG(gammaAGuess, CTGActual, BA, BL)
            % Calculates the global thrust coefficient using the given guess
            % for gammaA and known BA and BL, and compares to the measured
            % global thrust coefficient. Returns the absolute error between the
            % two.

            % Use guess to compute alphaA
            alphaA = NWTwoScale.calcAlpha(gammaAGuess, BA);

            % Compute CTG via induction
            CTG = NWTwoScale.calcGlobalThrust(gammaAGuess, alphaA, BA, BL);

            % Compute error
            err = abs(CTGActual - CTG);
        end

        function [gammaA, gammaAErr, gammaAExit] = convergeGammaA(gammaAGuess, CTGActual, BA, BL)
            % Iterates on gammaA using Dehtryiov et al Equations 28 and 9,
            % comparing the calculated global thrust coefficient to the
            % measured global thrust coefficient. Returns the converged value
            % of gammaA, as well as the resulting error and exit flag of
            % iteration.

            % Preallocate
            nPoints = length(CTGActual);
            gammaA = zeros(nPoints, 1);
            gammaAErr = zeros(nPoints, 1);
            gammaAExit = zeros(nPoints, 1);

            % Set options
            opts = optimset('TolFun', 1e-10, 'TolX', 1e-10);

            % Iterate
            for k = 1:nPoints
                currFun = @(gammaA) NWTwoScale.compareCTG(gammaA, CTGActual(k), BA(k), BL(k));
                [gammaA(k), gammaAErr(k), gammaAExit(k)] = fminsearch(currFun, gammaAGuess, opts);
            end
        end

        %% Iteration on GammaL
        function err = compareCTL(gammaLGuess, CTLActual, BL)
            % Calculates the local thrust coefficient using the given guess
            % for gammaL and known BL, and compares to the local thrust
            % coefficient obtained from the Equation 5. Returns the error
            % between the two thrust calculations.

            % Use guess to compute alphaA
            alphaL = NWTwoScale.calcAlpha(gammaLGuess, BL);

            % Compute CTL via induction
            CTL = NWTwoScale.calcThrustFromInduction(gammaLGuess, alphaL, BL);

            % Compute error
            err = abs(CTLActual - CTL);
        end

        function [gammaL, gammaLErr, gammaLExit] = convergeGammaL(gammaLGuess, CTL, BL)
            % Iterates on gammaL using Dehtryiov et al Equations 6 and 8,
            % comparing the local thrust coefficient calculated from induction
            % parameters to that calculated from the global thrust. Returns
            % the converged value of gammaL, as well as the resulting error 
            % and exit flag of iteration.

            % Preallocate
            nPoints = length(CTL);
            gammaL = zeros(nPoints, 1);
            gammaLErr = zeros(nPoints, 1);
            gammaLExit = zeros(nPoints, 1);

            % Set options
            opts = optimset('TolFun', 1e-10, 'TolX', 1e-10);

            % Iterate
            for k = 1:nPoints
                currFun = @(gammaL) NWTwoScale.compareCTL(gammaL, CTL(k), BL(k));
                [gammaL(k), gammaLErr(k), gammaLExit(k)] = fminsearch(currFun, gammaLGuess);
            end
        end

        %% Initialization/organization
        function [s] = initializeScaleStruct()
            % Preallocates structures for storing information about a given
            % scale.
            s = struct;
            s.beta = [];
            s.CT = [];
            s.Uinf = [];
            s.ut = [];
            s.uw = [];
            s.ub = [];
            s.alpha = [];
            s.gamma = [];
        end

        %% Print stuff out
        function modelName = getModelName()
            % Returns a label-friendly version of the blockage correction
            % model name.
            modelName = 'NW Two-Scale';
        end
    end
end