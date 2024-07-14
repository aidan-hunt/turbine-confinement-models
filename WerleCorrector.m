% Aidan Hunt
%
% NOTE: The method of Werle is included in this repository for reference
% and for completeness, but is not recommended for use as a blockage
% correction. It is a recommended that a physics-informed, linear momentum
% based blockage correction (e.g., those implemented in the BWClosedChannel
% and HoulsbyOpenChannel classes) be applied instead.
%
% Implements a closed-channel blockage correction based on the theory of
% Werle introduced in "Wind Turbine Wall-Blockage Performance Corrections"
% (2010; https://doi.org/10.2514/1.44602) and as implemented by Ross and
% Polagye in "An experimental assessment of analytical blockage
% corrections" (2020; https://doi.org/10.1016/j.renene.2020.01.135).
%
% To use the WerleCorrector class, construct a WerleCorrector object
% using the following syntax:
%   wl = WerleCorrector()
% and call methods using the dot notation (i.e., wlpredictUnconfined(...)).
% 
% WerleCorrector Methods:
%   predictUnconfined - Uses the closed-channel Werle correction to predict 
%                       unconfined performance from confined performance
%                       data.
%
% The methods above expect that confined performance data is provided as an
% mxn structure array, conf, with the following fields:
%   beta (required)   - blockage ratio
%   CT   (optional)   - thrust coefficient
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
% The WerleCorrector class extends the BCBase class
%
% See also: BCBase, BWClosedChannel, HoulsbyOpenChannel

classdef WerleCorrector < BCBase
    % The WerleCorrector class inherits all properties of the
    % BlockageCorrector superclass.

    methods (Static, Access = public)
        %% Main Method for Blockage Correction

        function [unconf] = predictUnconfined(conf)
            % Given a structure of confined data, converts to unconfined via
            % the Werle Correction.
            % Inputs:
            %   conf      - Performance data at a confined condition
            %               formatted as described in the WerleCorrector
            %               class definition.
            % Outputs:
            %   unconf    - Performance data at the unconfined condition. 
            %               unconf is the same size as conf and has the
            %               following fields, if present in conf:
            %                   CT  - unconfined thrust coefficient
            %                   CP  - unconfined performance coefficient
            %                   CQ  - unconfined torque coefficient
            %                   CL  - unconfined lateral force coefficient
            %                   CF  - unconfined resultant force coefficient
            %                   TSR - unconfined tip-speed ratio
            %
            % See also: WerleCorrector

            % Check input for correct sizing
            conf = WerleCorrector.checkInputSizes(conf);
            
            % Apply correction to each data set
            for i = 1:size(conf,1)
                for j = 1:size(conf,2)
                    unconf(i,j) = WerleCorrector.convertConfToUnconf(conf(i,j));
                end
            end
        end
    end

    methods (Static, Access = public)

        % Overwrite superclass method with Werle-style correction
        function unconfData = convertConfToUnconf(confData)
            % Given a 1x1 structure of confined data, converts to
            % unconfined via the Werle Correction.
            % Inputs:
            %   conf      - Performance data at a confined condition
            %               formatted as described in the WerleCorrector
            %               class definition.
            % Outputs:
            %   unconf    - Performance data at the unconfined condition. unconf
            %               is the same size as conf and has the following
            %               fields, if present in conf:
            %                   CT  - unconfined thrust coefficient
            %                   CP  - unconfined performance coefficient
            %                   CQ  - unconfined torque coefficient
            %                   CL  - unconfined lateral force coefficient
            %                   CF  - unconfined resultant force coefficient
            %                   TSR - unconfined tip-speed ratio
            %
            % See also: WerleCorrector, predictUnconfined
            fieldList = fieldnames(confData);
            unconfData = struct;
            for i = 1:length(fieldList)
                currField = fieldList{i};
                switch currField
                    case 'CP'
                        unconfData.(currField) = confData.(currField) .* ((1 - confData.beta).^2);
                    case {'CT', 'CL', 'CF', 'CQ'}
                        unconfData.(currField) = confData.(currField) .* ((1 - confData.beta).^2) ./ (1 + confData.beta);
                    case 'TSR'
                        unconfData.(currField) = confData.(currField) .* (1 - confData.beta);
                        unconfData.([currField '_correctSign']) = confData.(currField) .* (1 + confData.beta);
                end
            end
        end
    end
end