% Aidan Hunt
%
% Implements a two-scale analytical blockage correction for turbine arrays
% as described by Dehtyriov et al. in "A two-scale blockage correction for
% an array of tidal turbines" (2023; https://doi.org/10.36688/ewtec-2023-366).
%
% To use the DehtyriovCorrector class, construct a DehtyriovCorrector object
% using the following syntax:
%   dc = DehtyriovCorrector()
% and call methods using the dot notation (i.e., dc.predictUnconfined(...)).
% 
% DehtyriovCorrector Methods:
%   predictUnconfined - Uses a two-scale blockage correction to predict 
%                       unconfined performance from confined performance
%                       data.
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
% The predictUnconfined method utilizes helper methods that implement the
% equation set described in Dehtyriov et al (2023). The DehtyriovCorrector
% class extends the BCBase class.
%
% See also: BCBase, BWClosedChannel, HoulsbyOpenChannel

classdef DehtyriovCorrector < BCBase
    % The DehtryiovCorrector class inherits all properties of the
    % BCBase superclass.
    properties (Access = public)
        % Calculated array-scale blockage is multiplied by this value. 
        % This can be used to artificially alter the array-scale blockage
        % used in the correction and assess sensitivity of the results.
        debugMult = 1; 
    end

    methods (Access = public)
        %% Main method for blockage correction
        function unconf = predictUnconfined(dc, conf, geom, gammaGuess)
            % Applies the two-scale blockage correction of Dehtyriov et al.
            % to confined performance data to predict performance in an
            % unconfined flow.
            %
            % Inputs (required)
            %   conf       - A structure of confined performance data with fields
            %                as described in the DehtyriovCorrector class
            %                documentation.
            %   geom       - A structure of array and channel geometry
            %                information with the following fields:
            %                 n - Number of turbines in the array
            %                 D - Diameter of each turbine (assumed equal
            %                     across all turbines)
            %                 H - Blade span of each turbine (assumed equal
            %                     across all turbines)
            %                 w - Channel width
            %                 s - Spacing between adjacent turbines (assumed
            %                     equal across all turbines)
            %   gammaGuess - An initial guess for the value of the local
            %                and array-scale wake velocity induction
            %                factor. A single guess is used to iterate on
            %                both values.
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
            %               BA         - Calculated array-scale blockage
            %               gammaAIter - Iteration information and
            %                            converged value of array-scale
            %                            wake velocity induction factor
            %               alphaA     - Array-scale velocity induction factor
            %                            at each point
            %               BL         - Calculated local-scale blockage
            %               gammaLIter - Iteration information and
            %                            converged value of local-scale
            %                            wake velocity induction factor
            %               alphaA     - Local-scale velocity induction factor
            %                            at each point
            %
            % See also: DehtyriovCorrector

            % Check input for correct sizing
            conf = dc.checkInputSizes(conf);

            % Loop through all cases
            for i = 1:size(conf, 1)
                for j = 1:size(conf, 2)
                    % Calculate blockage scales from the provided geometry data
                    BL = dc.calcLocalBlockage(geom(i,j).D, geom(i,j).H, conf(i,j).d0, geom(i,j).s);
                    BA = dc.calcArrayBlockage(geom(i,j).D, geom(i,j).n, conf(i,j).d0, geom(i,j).w, geom(i,j).s);

                    % Debug: set BA to a value that is not 1
                    BA = BA .* dc.debugMult;

                    % Preallocate
                    nPoints = length(conf(i,j).TSR);

                    gammaA = zeros(nPoints, 1);
                    gammaAErr = zeros(nPoints, 1);
                    gammaAExit = zeros(nPoints, 1);
                    alphaA = zeros(nPoints, 1);

                    gammaL = zeros(nPoints, 1);
                    gammaLErr = zeros(nPoints, 1);
                    gammaLExit = zeros(nPoints, 1);
                    alphaL = zeros(nPoints, 1);
                    CTL = zeros(nPoints, 1);

                    % Loop through all set-points to get gammas, alphas
                    for k = 1:nPoints
                        % Iterate to solve for gammaA, alphaA
                        [gammaA(k), gammaAErr(k), gammaAExit(k)] = dc.convergeGammaA(gammaGuess, conf(i,j).CT(k), BA(k), BL(k));
                        alphaA(k) = dc.calcAlpha(gammaA(k), BA(k));
            
                        % Solve for CTL
                        CTL(k) = dc.calcLocalThrustFromGlobal(conf(i,j).CT(k), alphaA(k));
            
                        % Iterate to solve for gammaL, alphaL
                        [gammaL(k), gammaLErr(k), gammaLExit(k)] = dc.convergeGammaL(gammaGuess, CTL(k), BL(k));
                        alphaL(k) = dc.calcAlpha(gammaL(k), BL(k));
                    end
            
                    % Solve for CTA
                    CTA = dc.calcArrayThrust(alphaA(k), BL, CTL);
            
                    % Solve for velocity correction factor
                    velRatio = dc.calcVelocityRatio(alphaA, CTA);
                    UinfPrime = conf(i,j).Uinf ./ velRatio;
        
                    % Perform correction
                    unconf_temp = dc.convertConfToUnconf(conf(i,j), UinfPrime);
        
                    % Package results
                    unconf_temp.BA = BA;
                    unconf_temp.gammaA_iter = dc.packageDiagnostics(gammaA, dc.calcGlobalThrust(gammaA, alphaA, BA, BL), ... 
                                                                    gammaAErr, gammaAExit);
                    unconf_temp.alphaA = alphaA;


                    unconf_temp.BL = BL;
                    unconf_temp.gammaL_iter = dc.packageDiagnostics(gammaL, dc.calcThrustFromInduction(gammaL, alphaL, BL), ...
                                                                    gammaLErr, gammaLExit);
                    unconf_temp.alphaL = alphaL;

                    unconf(i,j) = unconf_temp;
                end
            end
        end
    end

    methods (Static, Access = public)

        %% Calculating blockage ratios at different scales

        % Dehtryiov et al Equation 2, but adapated for cross-flow turbines.
        function BL = calcLocalBlockage(D, H, h, s)
            % Dehtryiov et al Equation 2, but adapated for cross-flow turbines.
            % Calculates the local scale blockage, defined as the device area
            % divided by the local passage cross-sectional area
            % Inputs:
                % D - Turbine diameter [m]
                % H - Turbine blade span [m]
                % h - Channel detph [m]
                % s - Spacing between turbines [m]
            % Outputs:
                % BL - Local scale blockage ratio
            A = D.*H;
            BL = A ./ (h.*(D+s));
        end

        % Dehtryiov et al Equation 1, but adapated for cross-flow turbines.
        function BA = calcArrayBlockage(D, n, h, w, s)
            % Dehtryiov et al Equation 1, but adapated for cross-flow turbines.
            % Calculates the array scale blockage, defined as the array area
            % (including spacing between turbines) divided by the channel area.
            % Inputs:
                % D - Turbine diameter [m]
                % n - Number of turbines in the array
                % h - Channel detph [m]
                % w - Channel width [m]
                % s - Spacing between turbines [m]
            % Outputs:
                % BL - Local scale blockage ratio
            BA = (h.*n.*(D+s)) ./ (h.*w);
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
            %   BA  - Blockage ratio at the scale of interest
            % Outputs:
            %   alpha - Velocity induction factor at the same scale as gamma
            %           and B
            num = (1 + gamma);
            den = (1 + B) + sqrt((1 - B).^2 + B.*(1 - 1./gamma).^2);
            alpha = 1 - (num./den);
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
            num = (1 + gamma) - 2.*B.*(1 - alpha);
            den = (1 - B.*(1 - alpha)./gamma).^2;
            CT = (1 - gamma) .* (num./den);
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
            num = (1 + gammaA) - 2.*BA.*(1 - alphaA);
            den = (1 - BA.*(1-alphaA)./gammaA).^2;
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
            CTA = (1-alphaA).^2 .* BL .* CTL;
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
            CTL = CTG ./ (1 - alphaA).^2;
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
            num = 4.*(1 - alphaA);
            den = 4.*(1 - alphaA).^2 + CTA;
            velRatio = num ./ den;
        end

        %% Iteration on GammaA

        function err = compareCTG(gammaGuess,CTGActual, BA, BL)
            % Calculates the global thrust coefficient using the given guess
            % for gammaA and known BA and BL, and compares to the measured
            % global thrust coefficient. Returns the absolute error between the
            % two.

            % Use guess to compute alphaA
            alphaA = DehtyriovCorrector.calcAlpha(gammaGuess, BA);

            % Compute CTG via induction
            CTG = DehtyriovCorrector.calcGlobalThrust(gammaGuess, alphaA, BA, BL);

            % Compute error
            err = CTGActual - CTG;
        end

        function [gammaA, err, exitFlag] = convergeGammaA(gammaGuess, CTGActual, BA, BL)
            % Iterates on gammaA using Dehtryiov et al Equations 28 and 9,
            % comparing the calculated global thrust coefficient to the
            % measured global thrust coefficient. Returns the converged value
            % of gammaA, as well as the resulting error and exit flag of
            % iteration.
            [gammaA, err, exitFlag] = fzero(@(gammaA) DehtyriovCorrector.compareCTG(gammaA, CTGActual, BA, BL), gammaGuess);
        end

        %% Iteration on GammaL
        function err = compareCTL(gammaGuess, CTLActual, BL)
            % Calculates the local thrust coefficient using the given guess
            % for gammaL and known BL, and compares to the local thrust
            % coefficient obtained from the Equation 5. Returns the error
            % between the two thrust calculations.

            % Use guess to compute alphaA
            alphaL = DehtyriovCorrector.calcAlpha(gammaGuess, BL);

            % Compute CTG via induction
            CTL = DehtyriovCorrector.calcThrustFromInduction(gammaGuess, alphaL, BL);

            % Compute error
            err = CTLActual - CTL;
        end

        function [gammaL, err, exitFlag] = convergeGammaL(gammaGuess, CLActual, BL)
            % Iterates on gammaL using Dehtryiov et al Equations 6 and 8,
            % comparing the local thrust coefficient calculated from induction
            % parameters to that calculated from the global thrust. Returns
            % the converged value of gammaL, as well as the resulting error 
            % and exit flag of iteration.
            [gammaL, err, exitFlag] = fzero(@(gammaL) DehtyriovCorrector.compareCTL(gammaL, CLActual, BL), gammaGuess);
        end
    end
end