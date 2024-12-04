% Aidan Hunt
%

classdef SteirosPotentialFlow < BCBase


    methods (Access=public)

        function [conf] = solvePotentialFlow(st, conf, utUinfGuess)
            arguments
                st
                conf
                utUinfGuess = 0.5
            end

            % Check input sizes
            conf = st.checkInputSizes(conf);

            for i = 1:numel(conf) % For each element of the confined flow
                % Solve for the velocity through the turbine as a function
                % of thrust and blockage
                [utUinf, err, exitFlag] = st.convergeUtUinf(utUinfGuess, conf(i));

                % Calculate all other flow quantities
                ubUinf = st.solveBypassVel(utUinf, conf(i).beta);
                uwUinf = st.solveWakeVel(utUinf, conf(i).beta);
                betaWake = st.solveWakeWidth(utUinf, conf(i).beta);

                % Package for return
                conf(i).ut = utUinf .* conf(i).Uinf;
                conf(i).uw = uwUinf .* conf(i).Uinf;
                conf(i).ub = ubUinf .* conf(i).Uinf;
                conf(i).betaWake = betaWake;

                conf(i).utUinfIter = st.packageDiagnostics(utUinf, [st.solveThrust(utUinf, conf(i).beta), conf(i).CT], err, exitFlag);
            end

        end

        function [conf2, conf1] = forecastConfined(st, conf1, beta2, uGuess, options)
            arguments
                st
                conf1
                beta2
                uGuess (1,1) double = 0.5
                options.mode {ismember(options.mode, {'standard', 'bluff body'})} = 'standard';
            end

            % Solve potential flow for provided data
            conf1 = st.solvePotentialFlow(conf1, 0.5);

            for i = 1:size(conf1,1)
                for j = 1:size(conf1,2)
        
                    currForecast = struct;
                    switch options.mode
                        case 'standard'
                            utUinf2 = st.betaMatchUtUinf(uGuess, conf1(i,j), beta2);
                            Uinf2 = conf1(i,j).ut ./ utUinf2;
                            currForecast = st.convertConfToUnconf(conf1(i,j), Uinf2);
    
                        case 'bluff body'
                            if beta2 == 0
                                ubUinf2 = ones(size(conf1(i,j).CT));
                            else
                                ubUinf2 = st.betaMatchUbUinf(uGuess, conf1(i,j), beta2);
                            end
                            Uinf2 = conf1(i,j).ub ./ ubUinf2;
                            currForecast = st.convertConfToUnconf(conf1(i,j), Uinf2);
                    end

                    % Solve potential flow at conf2
                    currForecast.beta = beta2;
                    currForecast.h = NaN .* ones(size(currForecast.Uinf));
                    currForecast = st.solvePotentialFlow(currForecast, 0.5);
                
                    % Save
                    conf2(i,j) = currForecast;
                end
            end
        end

        function [unconf, conf] = predictUnconfined(st, conf, uGuess, options)
            arguments
                st
                conf
                uGuess (1,1) double = 0.5
                options.mode {ismember(options.mode, {'standard', 'bluff body'})} = 'standard';
            end
            [unconf, conf] = forecastConfined(st, conf, 0, uGuess, mode=options.mode);
        end

    end

    methods (Static, Access=public)

        %% Steiros Equation Set

        function [CT] = solveThrust(utUinf, beta)
        % Solves for the thrust coefficient as a function of ut, Uinf, and
        % beta using Steiros et al. Equation 2.16.
            
            % Compute thrust
            num = 4 .* (utUinf.*beta - 1) .* (1 - utUinf);
            den = (1 - beta) .* (2 - utUinf - utUinf.*beta);

            CT = num./den .* ( (1-utUinf)./3 - (1-2.*utUinf.*beta + beta)./(1-beta) );
        end

        function [ubUinf] = solveBypassVel(utUinf, beta)
        % Solves for the normalized bypass velocity as a function of utUinf
        % and beta using Steiros et al. Equation 2.17


            % Compute bypass velocity
            num = 1 - 2.*utUinf.*beta + beta;
            den = 1 - beta;

            ubUinf = num./den;
        end

        function [utUinf] = solveTurbVelFromBypass(ubUinf, beta)
        % Solves for the normalized velocity through the turbine as a function
        % of ubUinf and beta using Steiros et al. Equation 2.17
        % Note, this equation cannot be used if beta == 0, as the bypass
        % velocity is equal to the freestream velocity

            num = ubUinf.*(1-beta) - 1 - beta;
            den = -2.*beta;
            utUinf = num ./ den;
        end

        function [betaWake] = solveWakeWidth(utUinf, beta)
        % Calculates the wake width normalized by the channel width as a
        % function of Uinf, beta, and ut using Steiros et al. Equation 2.18

            % Compute wake width
            num = (2 - utUinf - utUinf.*beta) .* beta;
            den = 1 - 2.*utUinf.*beta + beta;

            betaWake = num./den;
        end

        function [uwUinf] = solveWakeVel(utUinf, beta)
        % Calculates the wake velocity as a function of utUinf and beta,
        % using Steiros equation 2.7

            %uwUinf = utUinf .* beta ./ betaWake;
            num = utUinf .* (1 - 2.*utUinf.*beta + beta);
            den = 2 - utUinf - utUinf.*beta;

            uwUinf = num ./ den;
        end


        %% Manipulating equations for applying to experimental data

        function [utUinf, err, exitFlag] = convergeUtUinf(utUinfGuess, conf)
        % Solves Steiros et al. Equation 2.16 for the velocity through the
        % actuator disk, ut, as a function of the freestream velocity,
        % blockage ratio, and thrust coefficient.

            % Preallocate
            nPoints = length(conf.CT);
            utUinf = zeros(size(conf.CT));
            err = zeros(size(conf.CT));
            exitFlag = zeros(size(conf.CT));

            % Iterate for each point
            for i = 1:nPoints
                % Define function to optimize
                currFun = @(utUinf) SteirosPotentialFlow.utUinfCompare(utUinf, conf.beta(i), conf.CT(i));

                % Search for minimum error
                [utUinf(i), err(i), exitFlag(i)] = fminsearch(currFun, utUinfGuess);
            end


        end

        function err = utUinfCompare(utUinfGuess, beta, CT)

            % Solve for CT using Steiros equation 2.16
            CTGuess = SteirosPotentialFlow.solveThrust(utUinfGuess, beta);

            % Compute error
            err = abs(real(CTGuess - CT));
        end

        %% Manipulating equations for standard blockage correction/forecast

        function [utUinf2, err, exitFlag] = betaMatchUtUinf(utUinf2Guess, conf1, beta2)

            % Preallocate
            nPoints = length(conf1.CT);
            utUinf2 = zeros(size(conf1.CT));
            err = zeros(size(conf1.CT));
            exitFlag = zeros(size(conf1.CT));

            utUinf1 = conf1.ut ./ conf1.Uinf;

            % Iterate for each point
            for i = 1:nPoints
                % Define function to optimize
                currFun = @(utUinf2) SteirosPotentialFlow.utUinf2Compare(utUinf2, beta2, utUinf1(i), conf1.CT(i));

                % Search for minimum error
                [utUinf2(i), err(i), exitFlag(i)] = fminsearch(currFun, utUinf2Guess);
            end

        end

        function err = utUinf2Compare(utUinf2Guess, beta2, utUinf1, CT1)

            % Compute CT2 Steiros et al. Equation 2.22 using CT1, utUinf1, utUinf2
            CT2_fromCT1 = CT1 .* (utUinf2Guess ./ utUinf1).^2;

            % Compute CT2 from Steiros et al. Equation 2.16
            CT2_fromBeta2 = SteirosPotentialFlow.solveThrust(utUinf2Guess, beta2);

            err = abs(real(CT2_fromCT1 - CT2_fromBeta2));
        end

        %% Manipulating equations for bluff-body correction/forecast

        function [ubUinf2, err, exitFlag] = betaMatchUbUinf(ubUinf2Guess, conf1, beta2)

            % Preallocate
            nPoints = length(conf1.CT);
            ubUinf2 = zeros(size(conf1.CT));
            err = zeros(size(conf1.CT));
            exitFlag = zeros(size(conf1.CT));

            ubUinf1 = conf1.ub ./ conf1.Uinf;

            % Iterate for each point
            for i = 1:nPoints
                % Define function to optimize
                currFun = @(ubUinf2) SteirosPotentialFlow.ubUinf2Compare(ubUinf2, beta2, ubUinf1(i), conf1.CT(i));

                % Search for minimum error
                [ubUinf2(i), err(i), exitFlag(i)] = fminsearch(currFun, ubUinf2Guess);
            end

        end

        function err = ubUinf2Compare(ubUinf2Guess, beta2, ubUinf1, CT1)

            % Compute CT2 Steiros et al. Equation 2.22 using CT1, ubUinf1, ubUinf2
            CT2_fromCT1 = CT1 .* (ubUinf2Guess./ubUinf1).^2;

            % Compute utUinf2 from Steiros et al. Equation 2.17
            utUinf2 = SteirosPotentialFlow.solveTurbVelFromBypass(ubUinf2Guess, beta2);

            % Compute CT2 from Steiros et al. Equation 2.16
            CT2_fromBeta2 = SteirosPotentialFlow.solveThrust(utUinf2, beta2);

            err = abs(real(CT2_fromCT1 - CT2_fromBeta2));
        end
    end
end