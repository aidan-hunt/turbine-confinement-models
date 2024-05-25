% Aidan Hunt
%
% Translates data from UW cross-flow turbine data acquisition into
% the format used by the blockage correction classes. The BCTranslator
% class extends the BCBase class.

classdef BCTranslator < BCBase
    properties (Constant, Access = public)
        % Naming conventions of important quantities from data acquisition
        % system (daqNames) and their corresponding names from Ross and
        % Polagye 2020 (BCBase.correctorNames).
        % eff/CP              : coefficient of performance
        % cThrust/CT          : coefficient of thrust
        % cLat/CL             : coefficient of lateral force
        % cForce/CF           : coefficient of resultant horizontal force
        % cTorque/CQ          : coefficient of torque
        % TSR/TSR             : tip-speed ratio
        % Uinf/V0             : "Undisturbed" freestream velocity (m/s)
        % FST(:,d0FSTInd)/d0  : "Undisturbed" dynamic depth (m)
        % B(:,betaFSTInd)/beta: Channel blockage (decimal)
        daqNames = {'eff', 'cThrust', 'cLat', 'cForce', 'cTorque', 'TSR', 'Uinf', 'FST', 'B', 'temp'};
        flipList = {'FY', 'MX', 'MZ', 'cTorque', 'cLat', 'turbPos', 'turbVel', 'TSR'};
    end

    methods (Static, Access = public)
        %% Naming convention translation
        function confData = translateTimeAveConfData(tds, fds, d0FSTInd, betaFSTInd, useAveDepth)
            % Given time-average cross-flow turbine experimental data as
            % processed via v4PostProcess, parses quantities relevant to
            % blockage corrections and returns as a structure whose fields
            % use the variable naming conventions of Ross and Polagye 2020
            % (see documentation for the 'daqNames' class property.
            % Quantities relevant to blockage corrections are averaged
            % across all turbines (i.e., array average quantities are
            % computed for use in correction).
            %
            % Inputs
            % tds         - Time-averaged turbine data from v4PostProcess
            % fds         - Time-averaged flow data from v4PostProcess
            % d0FSTInd    -  Column of fds.FST to use as representative of the
            %               undisturbed free surface, d0, for blockage
            %               correction. If using old data format
            %               (fds.FST_single), use [] as input.
            %               Default: 1
            % betaFSTInd  - Column of fds.FST to use as representative of the
            %               channel blockage, beta, for blockage correction.
            %               Default: 1
            % useAveDepth - Whether to use the average depth over the entire
            %               test as representative of d0 and beta (true), 
            %               or to use the measured d0 for each point (false). 
            %               Default: false
            %
            % Ouputs:
            % confData    - A structure of array-average turbine data prepared
            %               for blockage correction. See the "correctorNames"
            %               property of the BCTranslator class for
            %               possible fields.
            %
            % See also: daqNames, correctorNames

            % Check inputs
            arguments
                tds
                fds
                d0FSTInd (1,1) = 1
                betaFSTInd (1,1) = 1
                useAveDepth = false
            end

            % Get naming conventions
            oldNames = BCTranslator.daqNames;
            newNames = BCTranslator.correctorNames;

            % Loop through tds and fds, computing array averages when
            % appropriate.
            confData = struct;

            % Get mean velocities
            meanVels = mean([tds.turbVel], 1);

            for i = 1:length(oldNames)

                % The below only works for single turbines...?
                % if isfield(tds, oldNames{i}) && isrow(tds.(oldNames{i})) % Catch old data case
                %     tds.(oldNames{i}) = tds.(oldNames{i})';
                % end
                % if isfield(fds, oldNames{i}) && isrow(fds.(oldNames{i})) % Catch old data case
                %     fds.(oldNames{i}) = fds.(oldNames{i})';
                % end
                switch oldNames{i}
                    case {'eff', 'cThrust', 'cForce', 'cLat', 'cTorque', 'TSR'}
                        % Get data
                        currData = [tds.(oldNames{i})];

                        % Flip data
                        currData = BCTranslator.flipMetric(currData, oldNames{i}, meanVels);

                        % Average
                        confData.(newNames{i}) = mean(currData, 2);
                    case 'Uinf'
                        confData.(newNames{i}) = fds.(oldNames{i});
                    case 'FST'
                        if isfield(fds, 'FST')
                            confData.(newNames{i}) = fds.(oldNames{i})(:,d0FSTInd);
                        elseif isfield(fds, 'FST_single') % Check for old format
                            confData.(newNames{i}) = fds.FST_single;
                        else
                            error('Free surface data not found as either fds.FST or fds.FST_single. Please check input data.');
                        end
                    case 'B'
                        if isfield(fds, 'B')
                            confData.(newNames{i}) = fds.(oldNames{i})(:,betaFSTInd);
                        else % Check for old data format
                            error(['Field "B" (blockage) not found in fds structure. '...
                                   'Input data may be from earlier processing code ' ...
                                   'where this was not calculated automatically. Please ' ...
                                   'calculate blockage from FST readings and add to fds as fds.B.'])
                        end
                    % case 'temp'
                    %     % Import density and kinematic viscosity
                    %     confData.temp = fds.temp;
                    %     confData.rho = zeros(size(fds.temp));
                    %     confData.nu = zeros(size(fds.temp));
                    %     for k = 1:length(confData.temp)
                    %         [confData.rho(k), confData.nu(k)] = getWaterProps(confData.temp(k));
                    %     end
                end

                % Force to be column vector
                if isrow(confData.(newNames{i}))
                    confData.(newNames{i}) = confData.(newNames{i})';
                end
            end
        
            % If user wants to use average depth for corrections, average
            % the d0 and beta fields.
            if (useAveDepth)
                confData.d0 = mean(confData.d0);
                confData.beta = mean(confData.beta);
            end
        end

        function confData = quickTranslate(fileList)
            % Given a cell array of cross-flow turbine data files, loads
            % each and converts to a structure of data in the blockage
            % correction format. The resulting structure is the same size
            % as fileList. 
            %
            % It is assumed that that the first FST in fds.FST corresponds
            % to the undisturbed upstream water depth and blockge.
            %
            % See also: daqNames, correctorNames

            % Set quick defaults
            d0FSTInd = 1;
            betaFSTInd = 1;
            useAveDepth = false;

            % Loop through files and translate each
            for i = 1:size(fileList,1)
                for j = 1:size(fileList, 2)
                    load(fileList{i,j}, 'tds', 'fds');
                    confData(i,j) = BCTranslator.translateTimeAveConfData(tds, fds, d0FSTInd, betaFSTInd, useAveDepth);
                end
            end
        end

        function confData = translatePhaseAveConfData(tdp, fds, d0FSTInd, betaFSTInd, phasePercInd)
            % Given phase-averaged cross-flow turbine experimental data as
            % processed via v4PostProcess, parses quantities relevant to
            % blockage corrections and returns as a structure whose fields
            % use the variable naming conventions of Ross and Polagye 2020
            % (see documentation for BCBase.daqNames). Quantities relevant
            % to blockage corrections are averaged across all turbines
            % (i.e., array average quantities are computed for use in
            % correction).
            % ## NOTE: RIGHT NOW THIS ONLY WORKS FOR ARRAYS WITH DTHETA = 0
            % Inputs
            % tdp          - Phase-averaged turbine data from v4PostProcess
            % fds          - Time-averaged flow data from v4PostProcess
            % d0FSTInd     - Column of fds.FST to use as representative of the
            %                undisturbed free surface, d0, for blockage correction
            %                Default: 1
            % betaFSTInd   - Column of fds.FST to use as representative of the
            %                channel blockage, beta, for blockage correction
            %                Default: 1
            % phasePercInd - Column of each field of tdp to use, which
            %                corresponds to a particular percentile (default: 1)
            %
            % Ouputs:
            % confData    - A structure of array-average turbine data prepared
            %               for blockage correction. See the "correctorNames"
            %               property of the BCTranslator class for
            %               possible fields.
            %
            % See also: daqNames, correctorNames

            % Check inputs
            arguments
                tdp
                fds
                d0FSTInd (1,1) = 1
                betaFSTInd (1,1) = 1
                phasePercInd (1,1) = 1
            end

            % Define inputs
            oldNames = BCTranslator.daqNames;
            newNames = BCTranslator.correctorNames;
            if (~exist('phasePercInd', 'var') || isempty(phasePercInd))
                phasePercInd = 1;
            end

            confData = struct;
            for i = 1:size(tdp,1) % For each set-point in TDP
                nThetas = length(tdp(i,1).theta);
                confData(i,1).theta = mean(abs([tdp(i,:).theta]), 2);
                for j = 1:length(oldNames)
                    currField = oldNames{j};
                    switch currField
                        case {'eff', 'cThrust', 'cForce'}
                            currPhaseData = cat(3,tdp(i,:).(oldNames{j}));
                            currPhaseData = mean(currPhaseData, 3);
                            confData(i,1).(newNames{j}) = currPhaseData(:,phasePercInd);
                        case {'cLat', 'cTorque'}
                            currPhaseData = cat(3,tdp(i,:).(oldNames{j}));
                            currPhaseData = mean(currPhaseData, 3);
                            confData(i,1).(newNames{j}) = currPhaseData(:,phasePercInd);
                        case {'TSR'}
                            currPhaseData = cat(3,tdp(i,:).(oldNames{j}));
                            currPhaseData = mean(abs(currPhaseData), 3);
                            confData(i,1).(newNames{j}) = currPhaseData(:,phasePercInd);
                        case 'Uinf'
                            confData(i,1).(newNames{j}) = repmat(fds.(oldNames{j})(i), nThetas, 1);
                        case 'FST'
                            confData(i,1).(newNames{j}) = repmat(fds.(oldNames{j})(i,d0FSTInd), nThetas, 1);
                        case 'B'
                            confData(i,1).(newNames{j}) = repmat(fds.(oldNames{j})(i,betaFSTInd), nThetas, 1);
                    end
                end
            end
        end

        function unconfDataOut = translateUnconfData(unconfDataIn)
            % Translates unconfined data using Ross and Polagye naming
            % conventions back to the naming convention produced by
            % v4PostProcess. In other words, converts a structure of
            % unconfined data from BCTranslator.correctorNames to
            % BCTranslator.daqNames.
            %
            % See also: daqNames, correctorNames

            % Load turbine data and extract important quantities
            newNames = BCTranslator.daqNames;
            oldNames = BCTranslator.correctorNames;

            unconfDataOut = struct;
            for i = 1:length(unconfDataIn)
                for j = 1:length(oldNames)
                    unconfDataOut(i).(newNames{j}) = unconfDataIn(i).(oldNames{j});
                end
            end
        end

        function geom = translateGeometryData(turbines, flow)
            % Translates geometry data from cross-flow turbine HWconfig
            % file into form usable in blockage corrections. Turbines are
            % assumed to have identical geometry and equivalently spaced.
            %
            % Inputs:
            %   turbines - A structure of cross-flow turbine geometry information as
            %              produced by the HWConfigGUI.
            %   flow     - A structure of channel properties and flow
            %              information as produced by the HWConfigGUI
            %
            % Outputs:
            %   geom     - A structure of turbine geometry information that is
            %              compatible with blockage corrector classes, with
            %              naming conventions following that of Dehtyriov et al
            %              (2023). This structure has the following fields:
            %                 n - Number of turbines in an array
            %                 D - Diameter of each turbine (assumed equal
            %                     across all turbines)
            %                 H - Blade span of each turbine (assumed equal
            %                     across all turbines)
            %                 w - Channel width
            %                 s - Spacing between adjacent turbines (assumed
            %                     equal across all turbines)

            geom = struct;

            % Get geometric parameters of the turbines
            activeInd = find([turbines.isActive]);
            geom.n = length(activeInd);
            geom.D = zeros(geom.n, 1);
            geom.H = zeros(geom.n, 1);
            for i = 1:length(activeInd)
                geom.D(i) = turbines(activeInd(i)).geometry.D;
                geom.H(i) = turbines(activeInd(i)).geometry.H;
            end
            geom.D = unique(geom.D);
            geom.H = unique(geom.H);
            if length(geom.D) ~= 1 || length(geom.H) ~= 1
                error('Turbines must have identical geometry for blockage corrections to apply.');
            end

            % Get geometric parameters of the channel
            geom.w = flow.width;
            geom.s = (geom.w - (geom.D .* geom.n)) ./ geom.n;
        end

        function dataOut = flipMetric(dataIn, metricName, meanVel)
            % Flips the sign of dataIn if it is a flippable metric (i.e., a
            % metric whose sign depends on the direction of rotation, see
            % BCTranslator.flipList).
            %
            % Inputs
            %   dataIn     - Input data
            %   metricName - The name of the performance metric that dataIn
            %                represents
            %   meanVel    - The mean angular velocity corresponding to
            %                dataIn
            % Outputs
            %   dataOut    - dataIn with its sign flipped, if appropriate

            flipSW = ones(1, size(dataIn, 2));
            if (ismember(lower(metricName), lower(BCTranslator.flipList)))
                flipSW(meanVel < 0) = -1;
            end
            dataOut = dataIn .* flipSW;
        end
    end
end