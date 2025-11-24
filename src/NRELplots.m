classdef NRELplots
% NRELplots Class for generating property comparison plots for NREL paper
%
%   Usage:
%       np = NRELplots()
%       np = NRELplots('save',0,'high_fidelity',0,'zoom',0,'titles_on',1);
%
%   Inputs:
%       save            - (0/1) whether to save images (default: 0)
%       high_fidelity   - (0/1) whether to use high-fidelity data (default: 0)
%       zoom            - (0/1) whether to use zoomed in data (default: 0)
%       titles_on       - (0/1) whether to display titles on plots (default: 1)
%
%  Methods:
%       print_params()                   - prints parameter files for Pele runs
%       collect_data()                   - runs Pele simulations to collect data
%       plot_density()                   - plots density comparisons
%       plot_density_error()             - plots density error comparisons
%       plot_viscosity()                 - plots viscosity comparisons
%       plot_viscosity_error()           - plots viscosity error comparisons
%       plot_dmudrho()                   - plots dmu/drho comparisons
%       plot_conductivity()              - plots conductivity comparisons
%       plot_conductivity_error()        - plots conductivity error comparisons
%       plot_specificHeat()              - plots specific heat comparisons
%       plot_specificHeat_error()        - plots specific heat error comparisons
%       plot_dlamdrho()                  - plots dlam/drho comparisons
%       plot_density_floor()             - plots density floor comparisons
%       plot_kinematic_viscosity_ratio() - plots kinematic viscosity ratio
%
%  Author: Tyler Fara                  Date: November 23, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        opts                % options       
        data                % data storage
        F                   % computed fields
    end

    properties (Dependent)
        key                 % key for accessing fields in F
        F_name              % names for F
        F_filename          % filenames for F
        F_unit              % units for F
        F_contourLines      % contour lines for F
        C_min               % color limits
        C_max               % color limits
        C_scheme            % color scheme
        C_label             % colorbar label
        C_manual_levels     % manual contour levels
    end

    methods
        function self = NRELplots(args)

            arguments
                args.save = 0;
                args.high_fidelity = 0;
                args.zoom = 0;
                args.titles_on = 1;
            end

            % SET OPTIONS 
            self.opts.save_images = args.save;
            self.opts.high_fidelity = args.high_fidelity;
            self.opts.titles_on = args.titles_on;
            self.opts.zoom = args.zoom;

            % if zoomed in dat requestd
            if args.zoom

                % if high fidelity requested, do high fidelity zoom
                if self.opts.high_fidelity
                    self.opts.high_fidelity_zoom = 1;
                    self.opts.low_fidelity_zoom = 0;

                % else do low fidelity zoom
                else
                    self.opts.high_fidelity_zoom = 0;
                    self.opts.low_fidelity_zoom = 1;
                end

            % if no zoom requested
            else
                self.opts.high_fidelity_zoom = 0;
                self.opts.low_fidelity_zoom = 0;
            end


            % STORE KEY INDICES FOR DATA ACCESS
            % Define species data
            self.data.species = {'ch4','co2','n2'};
            self.data.T_crit  = [190.564, 304.1282, 126.192];
            self.data.P_crit  = [45.391, 72.808, 33.514];

            % define data labels
            self.data.species_str = {'CH$_4$','CO$_2$','N$_2$'};
            self.data.property = {'density','viscosity','conductivity','specificHeat','dmudrho','dlamdrho'};
            self.data.property_sym = {'\rho','\mu','\lambda','C_p','\partial \mu / \partial \rho','\partial \lambda / \partial \rho'};
            self.data.source_str = {'PelePhysics','CoolProp','Error'};

            % LOAD DATA
            self = self.read_data();
            self = self.compute_F();

        end

        function print_params(self)

            % Define parameters
            T_step_size = 0;
            T_step_num  = 100;
            T_start     = self.data.T_crit;
            T_end       = 3 * self.data.T_crit;
            P_step_size = 0;
            P_step_num  = 100;
            P_start     = self.data.P_crit;
            P_end       = 3 * self.data.P_crit;
            P_unit      = 'atm';
            R_unit      = 'kg/m3';
            mu_unit     = 'uPa_s';
            lam_unit    = 'W/m_K';
            isotherm_num         = 0;

            % Open file for writing
            for i = 1:length(self.data.species)
                f_str = sprintf('run_pelephysics/parameters/%s_2d_parms.txt', self.data.species{i});
                fid = fopen(f_str, 'w');
                if fid == -1
                    error('Failed to open file for writing.');
                end

                % Write each parameter line using the variables
                fprintf(fid, 'chemical_species          %s\n', upper(self.data.species{i}));
                fprintf(fid, 'temperature_step_size     %d\n', T_step_size);
                fprintf(fid, 'temperature_step_number   %d\n', T_step_num);
                fprintf(fid, 'temperature_start         %f\n', T_start(i));
                fprintf(fid, 'temperature_end           %f\n', T_end(i));
                fprintf(fid, 'pressure_step_size        %d\n', P_step_size);
                fprintf(fid, 'pressure_step_number      %d\n', P_step_num);
                fprintf(fid, 'pressure_start            %f\n', P_start(i));
                fprintf(fid, 'pressure_end              %f\n', P_end(i));
                fprintf(fid, 'pressure_unit             %s\n', P_unit);
                fprintf(fid, 'density_unit              %s\n', R_unit);
                fprintf(fid, 'viscosity_unit            %s\n', mu_unit);
                fprintf(fid, 'lam_unit                  %s\n', lam_unit);
                fprintf(fid, 'isotherm_num              %d\n', isotherm_num);

                % Close the file
                fclose(fid);
            end
        end

        function collect_data(self)
            for i = 1:length(self.data.species)
                fid = sprintf('run_pelephysics/parameters/%s_2d_parms.txt', self.data.species{i});
                outPath = sprintf('results/');
                modifier = sprintf('2d_');
                cmd = sprintf('./run_pelephysics/Pele2d.gnu.ex %s %s %s', fid, outPath, modifier);
                system(cmd);
            end
        end

        function self = read_data(self)
            
            dir = 'data/low_fidelity';

            if self.opts.high_fidelity
                dir = 'data/high_fidelity';
            end

            if self.opts.high_fidelity_zoom
                dir = 'data/high_fidelity_zoom';
            end

            if self.opts.low_fidelity_zoom
                dir = 'data/low_fidelity_zoom';
            end

            % read independent variables (T,P)
            T = self.read_temp(self.data,dir);
            P = self.read_pres(self.data,dir);

            % get pele data
            for spec = 1:length(self.data.species)
                for prop = 1:4
                    D_data_str = sprintf('%s/%s_2d_pelephys_%s.txt', dir, upper(self.data.species{spec}), self.data.property{prop});
                    D_tab = readtable(D_data_str,'ReadVariableNames',false);
                    D_pele{spec}{prop} = table2array(D_tab);
                end
            end

            % get pele derivative data
            for spec = 1:length(self.data.species)
                for prop = 5:6
                    D_data_str = sprintf('%s/%s_2d_mixed_%s.txt', dir, upper(self.data.species{spec}), self.data.property{prop});
                    D_tab = readtable(D_data_str,'ReadVariableNames',false);
                    D_pele{spec}{prop} = table2array(D_tab);
                end
            end

            % get coolprop data
            for spec = 1:length(self.data.species)
                for prop = 1:4
                    D_data_str = sprintf('%s/%s_2d_coolprop_%s.txt', dir, upper(self.data.species{spec}), self.data.property{prop});
                    D_tab = readtable(D_data_str,'ReadVariableNames',false);
                    D_cool{spec}{prop} = table2array(D_tab);
                end
            end

            % get cool derivative data
            %    NOTE: no derivative data from CoolProp, so set to NaN
            for spec = 1:length(self.data.species)
                for prop = 5:6
                    D_cool{spec}{prop} = NaN;
                end
            end

            % get mixed data
            for spec = 1:length(self.data.species)
                for prop = 1:4
                    D_data_str = sprintf('%s/%s_2d_mixed_%s.txt', dir, upper(self.data.species{spec}), self.data.property{prop});
                    D_tab = readtable(D_data_str,'ReadVariableNames',false);
                    D_test{spec}{prop} = table2array(D_tab);
                end
            end

            % get mixed derivative data
            for spec = 1:length(self.data.species)
                for prop = 5:6
                    D_data_str = sprintf('%s/%s_2d_mixed_%s.txt', dir, upper(self.data.species{spec}), self.data.property{prop});
                    D_tab = readtable(D_data_str,'ReadVariableNames',false);
                    D_test{spec}{prop} = table2array(D_tab);
                end
            end

            % get critical line data
            key = self.key;
            Cp_max = {};
            Cp_idx = {};
            T_crit_line = {};
            P_crit_line = {};
            for i = 1:length(self.data.species)
                for j = 1:2

                    % get max Cp values along isobars
                    %[Cp_max,Cp_idx] = max(D{i}{Cp}{j},[],2);
                    [Cp_max,Cp_idx] = max(D_cool{i}{key.Cp},[],2);

                    % lookup corresponding temperatures
                    T_crit_line_temp = T{i}(Cp_idx);
                    
                    % get corresponding pressure data
                    P_crit_line_temp = P{i};

                    % get indices of pressures greater than critical pressure
                    P_crit_idx = find(P_crit_line_temp > self.data.P_crit(i));

                    % keep only those data points corresponding to larger pressures
                    P_crit_line_temp = P_crit_line_temp(P_crit_idx);
                    T_crit_line_temp = T_crit_line_temp(P_crit_idx);

                    % get unique values
                    [T_crit_line_temp_unique, unique_idx] = unique(T_crit_line_temp);
                    P_crit_line_temp_unique = P_crit_line_temp(unique_idx);

                    % replace first data point with critical point
                    T_crit_line_temp_unique(1) = self.data.T_crit(i);
                    P_crit_line_temp_unique(1) = self.data.P_crit(i);

                    % append end point
                    T_crit_line_temp_unique = [T_crit_line_temp_unique; T_crit_line_temp(end)];
                    P_crit_line_temp_unique = [P_crit_line_temp_unique; P_crit_line_temp(end)];

                    % keep only those temperature and pressure data corresponding to these idx
                    T_crit_line{i}{j} = T_crit_line_temp_unique;
                    P_crit_line{i}{j} = P_crit_line_temp_unique;

                end
            end

            % store some results
            self.data.T = T;
            self.data.P = P;
            self.data.D_pele = D_pele;
            self.data.D_cool = D_cool;
            self.data.D_test = D_test;
            self.data.T_crit_line = T_crit_line;
            self.data.P_crit_line = P_crit_line;

        end


        function C_min = get.C_min(self)
            [C_min,~] = self.compute_color_limits();
        end

        function C_max = get.C_max(self)
            [~,C_max] = self.compute_color_limits();
        end

        function [C_min,C_max] = compute_color_limits(self)

            % store variables
            key = self.key();
            F = self.F;
            nSpecies = 3;
            nFields  = numel(F{1});

            % INITIAL PASS: compute min/max over all computed fields
            C_max = cell(nSpecies, 1);
            C_min = cell(nSpecies, 1);

            for i = 1:nSpecies
                C_max{i} = cell(1, nFields);
                C_min{i} = cell(1, nFields);
                for j = 1:nFields
                    Fij = F{i}{j};
                    C_max{i}{j} = max(Fij(:));
                    C_min{i}{j} = min(Fij(:));
                end
            end

            % STANDARDIZE COLOR LIMITS FOR DENSITY (gen vs ref)
            densPair = [key.dens_gen, key.dens_ref];

            for i = 1:nSpecies
                valsMax = [C_max{i}{densPair}];
                valsMin = [C_min{i}{densPair}];

                temp_max = max(valsMax);
                temp_min = min(valsMin);

                [C_max{i}, C_min{i}] = self.setGroupLimits(C_max{i}, C_min{i}, ...
                                                    densPair, temp_min, temp_max);
            end

            % STANDARDIZE COLOR LIMITS FOR VISCOSITY (gen, ref, hyb)
            viscGroup = [key.visc_gen, key.visc_ref, key.visc_hyb];

            for i = 1:nSpecies
                valsMax = [C_max{i}{viscGroup}];
                valsMin = [C_min{i}{viscGroup}];

                temp_max = max(valsMax);
                temp_min = min(valsMin);

                [C_max{i}, C_min{i}] = self.setGroupLimits(C_max{i}, C_min{i}, ...
                                                    viscGroup, temp_min, temp_max);
            end

            % DENSITY ABSOLUTE ERROR (symmetric about 0)
            for i = 1:nSpecies
                temp_max = C_max{i}{key.dens_ERR};
                temp_min = C_min{i}{key.dens_ERR};
                [temp_min, temp_max] = self.makeSymmetric(temp_min, temp_max);

                [C_max{i}, C_min{i}] = self.setGroupLimits(C_max{i}, C_min{i}, ...
                                                    key.dens_ERR, temp_min, temp_max);
            end

            % DENSITY RELATIVE ERROR (symmetric about 0, then unified across species)
            % first: per-species symmetric limits
            for i = 1:nSpecies
                temp_max = C_max{i}{key.dens_err};
                temp_min = C_min{i}{key.dens_err};
                [temp_min, temp_max] = self.makeSymmetric(temp_min, temp_max);

                [C_max{i}, C_min{i}] = self.setGroupLimits(C_max{i}, C_min{i}, ...
                                                    key.dens_err, temp_min, temp_max);
            end

            % then: override so every species has same limits (use last computed)
            temp_min = C_min{nSpecies}{key.dens_err};
            temp_max = C_max{nSpecies}{key.dens_err};
            for i = 1:nSpecies
                [C_max{i}, C_min{i}] = self.setGroupLimits(C_max{i}, C_min{i}, ...
                                                    key.dens_err, temp_min, temp_max);
            end

            % VISCOSITY ABSOLUTE ERROR (group of 5, symmetric)
            viscErrAbs = [ ...
                key.visc_ERR_tot, ...
                key.visc_ERR_corr, ...
                key.visc_ERR_eos, ...
                key.visc_ERR_assembled, ...
                key.visc_ERR_secondorder];

            for i = 1:nSpecies
                valsMax = [C_max{i}{viscErrAbs}];
                valsMin = [C_min{i}{viscErrAbs}];
                temp_max = max(valsMax);
                temp_min = min(valsMin);
                [temp_min, temp_max] = self.makeSymmetric(temp_min, temp_max);

                [C_max{i}, C_min{i}] = self.setGroupLimits(C_max{i}, C_min{i}, ...
                                                    viscErrAbs, temp_min, temp_max);
            end

            % VISCOSITY RELATIVE ERROR (you manually clamp to [-0.2, 0.2])
            viscErrRel = [ ...
                key.visc_err_tot, ...
                key.visc_err_corr, ...
                key.visc_err_eos, ...
                key.visc_err_assembled, ...
                key.visc_err_secondorder];

            temp_Min = -0.2;
            temp_Max =  0.2;

            for i = 1:nSpecies
                [C_max{i}, C_min{i}] = self.setGroupLimits(C_max{i}, C_min{i}, ...
                                                    viscErrRel, temp_Min, temp_Max);
            end

            % KINEMATIC VISCOSITY RATIO
            % (note: preserves your original formula, just formatted)
            for i = 1:nSpecies
                maxVal = abs(1 - max(self.max([key.visc_kinematic_ratio])));
                minVal = abs(1 - min(self.min([key.visc_kinematic_ratio])));
                val    = max(maxVal, minVal);

                C_max{i}{key.visc_kinematic_ratio} = 1 + val;
                C_min{i}{key.visc_kinematic_ratio} = 1 - val;
            end

            % CONDUCTIVITY (gen, ref, hyb) – with your manual overrides
            condGroup = [key.cond_gen, key.cond_ref, key.cond_hyb];

            for i = 1:nSpecies
                % original auto-scaling (kept but simplified, though overridden below)
                temp_max = C_max{i}{key.cond_gen};
                temp_min = C_min{i}{key.cond_gen};
                temp_min = temp_min - 0.05 * temp_min;

                % manual species-specific limits (unchanged)
                if i == 1
                    temp_max = 0.099;
                    temp_min = 0.025;
                elseif i == 2
                    temp_max = 0.11;
                    temp_min = 0.03;
                elseif i == 3
                    temp_max = 0.07;
                    temp_min = 0.01;
                end

                [C_max{i}, C_min{i}] = self.setGroupLimits(C_max{i}, C_min{i}, ...
                                                    condGroup, temp_min, temp_max);
            end

            % CONDUCTIVITY ABSOLUTE ERROR (manual species-specific limits)
            condErrAbs = [ ...
                key.cond_ERR_tot, ...
                key.cond_ERR_corr, ...
                key.cond_ERR_eos, ...
                key.cond_ERR_assembled, ...
                key.cond_ERR_secondorder];

            for i = 1:nSpecies
                % NOTE: you had some commented logic using sorted errors; I leave
                % that out and keep the hand-tuned species limits.

                if i == 1
                    temp_max = 20;   temp_min = -20;
                elseif i == 2
                    temp_max = 50;   temp_min = -50;
                else % i == 3
                    temp_max = 8;    temp_min = -8;
                end

                [C_max{i}, C_min{i}] = self.setGroupLimits(C_max{i}, C_min{i}, ...
                                                    condErrAbs, temp_min, temp_max);
            end

            % CONDUCTIVITY RELATIVE ERROR (global, then CO2 override)
            condErrRel = [ ...
                key.cond_err_tot, ...
                key.cond_err_corr, ...
                key.cond_err_eos, ...
                key.cond_err_assembled, ...
                key.cond_err_secondorder];

            % global manual clamp
            temp_Min = -0.25;
            temp_Max =  0.25;

            for i = 1:nSpecies
                [C_max{i}, C_min{i}] = self.setGroupLimits(C_max{i}, C_min{i}, ...
                                                    condErrRel, temp_Min, temp_Max);
            end

            % override for CO2 (species 2)
            temp_Min = -0.5;
            temp_Max =  0.5;
            i = 2;
            [C_max{i}, C_min{i}] = self.setGroupLimits(C_max{i}, C_min{i}, ...
                                                condErrRel, temp_Min, temp_Max);

            % DENSITY FLOOR LIMITS
            for i = 1:nSpecies
                C_max{i}{key.visc_density_FLOOR} = -150;
                C_min{i}{key.visc_density_FLOOR} = -1300;
                C_max{i}{key.visc_density_floor} = 0;
                C_min{i}{key.visc_density_floor} = -50;

                C_max{i}{key.cond_density_FLOOR} = -0.1;
                C_min{i}{key.cond_density_FLOOR} = -1000;
                C_max{i}{key.cond_density_floor} = -1;
                C_min{i}{key.cond_density_floor} = -50;
            end

            % species-specific overrides for cond_density_FLOOR
            C_min{1}{key.cond_density_FLOOR} = -700;
            C_max{1}{key.cond_density_FLOOR} = -100;
            C_min{2}{key.cond_density_FLOOR} = -900;
            C_max{2}{key.cond_density_FLOOR} = -300;
            C_min{3}{key.cond_density_FLOOR} = -800;
            C_max{3}{key.cond_density_FLOOR} = -300;

        end

        function [minVal, maxVal] = makeSymmetric(self, minVal, maxVal)
        % Force symmetric limits around 0 using the larger magnitude.
            if abs(minVal) > abs(maxVal)
                maxVal = abs(minVal);
            else
                minVal = -abs(maxVal);
            end
        end

        function [CmaxRow, CminRow] = setGroupLimits(self, CmaxRow, CminRow, idxGroup, minVal, maxVal)
        % Set the same min/max for a group of field indices.
            for idx = idxGroup
                CmaxRow{idx} = maxVal;
                CminRow{idx} = minVal;
            end
        end

 

        % PLOTTERS 
        function plot_density(self,spec2plot)
        % plots generalized and reference densities
            key = self.key;
            if nargin < 2, spec2plot = [key.ch4, key.co2, key.n2]; end
            props2plot = [key.dens_gen, key.dens_ref];
            self.plot2d(spec2plot,props2plot);
        end

        function plot_density_error(self,spec2plot)
        % plots absolute and relative density errors
            key = self.key;
            if nargin < 2, spec2plot = [key.ch4, key.co2, key.n2]; end
            props2plot = [key.dens_ERR, key.dens_err];
            self.plot2d(spec2plot,props2plot);
        end

        function plot_viscosity(self,spec2plot)
        % plots generalized, reference, and hybrid viscosities
            key = self.key;
            if nargin < 2, spec2plot = [key.ch4, key.co2, key.n2]; end
            props2plot = [key.visc_gen, key.visc_ref, key.visc_hyb];
            self.plot2d(spec2plot,props2plot);
        end

        function plot_viscosity_error(self,spec2plot)
        % plots both absolute and relative viscosity errors
            if nargin < 2, spec2plot = [self.key.ch4, self.key.co2, self.key.n2]; end
            self.plot_viscosity_error_absolute(spec2plot);
            self.plot_viscosity_error_relative(spec2plot);
        end

        function plot_viscosity_error_absolute(self,spec2plot)
        % plots absolute viscosity errors
            key = self.key;
            if nargin < 2, spec2plot = [key.ch4, key.co2, key.n2]; end
            props2plot = [key.visc_ERR_tot, ...
                            key.visc_ERR_corr, ...
                            key.visc_ERR_eos, ...
                            key.visc_ERR_assembled, ...
                            key.visc_ERR_secondorder];
            

            self.plot2d(spec2plot,props2plot);
        end

        function plot_viscosity_error_relative(self,spec2plot)
        % plots relative viscosity errors
            key = self.key;
            if nargin < 2, spec2plot = [key.ch4, key.co2, key.n2]; end
            props2plot = [key.visc_err_tot, ...
                            key.visc_err_corr, ...
                            key.visc_err_eos, ...
                            key.visc_err_assembled, ...
                            key.visc_err_secondorder];

            self.plot2d(spec2plot,props2plot);
        end

        function plot_dmudrho(self,spec2plot)
            key = self.key;
            if nargin < 2, spec2plot = [key.ch4, key.co2, key.n2]; end
            props2plot = [key.visc_dmudrho];
            self.plot2d(spec2plot,props2plot);
        end

        function plot_density_floor_viscosity(self,spec2plot)
            key = self.key;
            if nargin < 2, spec2plot = [key.ch4, key.co2, key.n2]; end
            props2plot = [key.visc_density_FLOOR, key.visc_density_floor];

            self.plot2d(spec2plot,props2plot);
        end

        function plot_kinematic_viscosity_ratio(self,spec2plot)
            key = self.key;
            if nargin < 2, spec2plot = [key.ch4, key.co2, key.n2]; end
            props2plot = [key.visc_kinematic_ratio];

            self.plot2d(spec2plot,props2plot);
        end


        % CONDUCTIVITY PLOTTERS
        function plot_conductivity(self,spec2plot)
        % plots generalized, reference, and hybrid conductivities
            key = self.key;
            if nargin < 2, spec2plot = [key.ch4, key.co2, key.n2]; end
            props2plot = [key.cond_gen, key.cond_ref, key.cond_hyb];
            self.plot2d(spec2plot,props2plot);
        end

        function plot_conductivity_error(self,spec2plot,prop2plot)
        % plots both absolute and relative conductivity errors
            key = self.key;
            if nargin < 2, spec2plot = [self.key.ch4, self.key.co2, self.key.n2]; end
            if nargin < 3, prop2plot = [self.key.cond_ERR_tot, ...
                            key.cond_ERR_corr, ...
                            key.cond_ERR_eos, ...
                            key.cond_ERR_assembled, ...
                            key.cond_ERR_secondorder]; end
            self.plot_conductivity_error_absolute(spec2plot,prop2plot);
            self.plot_conductivity_error_relative(spec2plot,prop2plot);
        end

        function plot_conductivity_error_absolute(self,spec2plot,prop2plot)
        % plots absolute conductivity errors
            key = self.key;
            if nargin < 2, spec2plot = [key.ch4, key.co2, key.n2]; end
            if nargin < 3, prop2plot = [key.cond_ERR_tot, ...
                            key.cond_ERR_corr, ...
                            key.cond_ERR_eos, ...
                            key.cond_ERR_assembled, ...
                            key.cond_ERR_secondorder]; end
            self.plot2d(spec2plot,prop2plot);
        end

        function plot_conductivity_error_relative(self,spec2plot,prop2plot)
        % plots relative conductivity errors
            key = self.key;
            if nargin < 2, spec2plot = [key.ch4, key.co2, key.n2]; end
            if nargin < 3, prop2plot = [key.cond_err_tot, ...
                            key.cond_err_corr, ...
                            key.cond_err_eos, ...
                            key.cond_err_assembled, ...
                            key.cond_err_secondorder]; end
            self.plot2d(spec2plot,prop2plot);
        end

        function plot_dlamdrho(self,spec2plot)
            key = self.key;
            if nargin < 2, spec2plot = [key.ch4, key.co2, key.n2]; end
            props2plot = [key.cond_dlamdrho];
            self.plot2d(spec2plot,props2plot);
        end

        function plot_density_floor_conductivity(self,spec2plot)
            key = self.key;
            if nargin < 2, spec2plot = [key.ch4, key.co2, key.n2]; end
            props2plot = [key.cond_density_FLOOR, key.cond_density_floor];
            self.plot2d(spec2plot,props2plot);
        end

        function plot_density_floor(self,spec2plot)
            key = self.key;
            if nargin < 2, spec2plot = [key.ch4, key.co2, key.n2]; end
            self.plot_density_floor_viscosity(spec2plot);
            self.plot_density_floor_conductivity(spec2plot);
        end


        % PLOTTERS (Auxiliary Functions)
        function plot2d(self,spec2plot,props2plot)

            data = self.data;

            for spec = spec2plot, for prop = props2plot,

                % store data
                T_val = data.T{spec};
                P_val = data.P{spec};
                D_val = self.F{spec}{prop};

                % plot data
                h = pcolor(T_val,P_val,D_val);
                c = colorbar;

                % get handles to figure; h is surface
                g = h.Parent; % axes
                f = g.Parent; % figure
                f.Position(1) = 100;

                % run formatting functions
                self.format_axis(g,spec,prop);
                self.format_plot(h,spec,prop);
                self.format_colorbar(c,spec,prop);

                % plot Widom line
                key = self.key;
                hold on;
                    % using less accurate pele data:
                    %plot(data.T_crit_line{spec}{key.pele}, data.P_crit_line{spec}{key.pele}, 'w--', 'LineWidth', 3);

                    % using more accurate coolprop data:
                    plot(data.T_crit_line{spec}{key.cool}, data.P_crit_line{spec}{key.cool}, 'w--', 'LineWidth', 3);
                    %f.Children(2).Children(1).Color = [0.6 0.6 0.6 0.9]; % make line semi-transparent
                    %f.Children(2).Children(1).Color = [1 1 1 0.9]; % make line semi-transparent
                    f.Children(2).Children(1).Color = [0 0 0 0.9]; % make line semi-transparent
                hold off;

                % overlay contour plots
                %{
                %   NOTE: umute the following for AUTOMATIC contour line labels
                hold on;
                [C,hh] = contour(T_val,P_val,D_val,self.F_contourLines{prop}{spec},...
                            'k-', ...
                            'LineWidth',1.5, ...
                            'ShowText',true, ...
                            'EdgeAlpha',0.2);
                clabel(C, hh, '', ...
                        'FontSize',24, ...
                        'Color','k', ...
                        'LabelSpacing',1000);
                hold off;
                %}

                %{
                %   NOTE: umute the following for MANUAL contour line labels
                hold on;
                [C,hh] = contour(T_val,P_val,D_val,self.F_contourLines{prop}{spec},...
                            'k-', ...
                            'LineWidth',1.5, ...
                            'ShowText',false, ...
                            'EdgeAlpha',0.2);
                clabel(C, hh, 'manual', ...
                        'FontSize',24, ...
                        'Color','k', ...
                        'LabelSpacing',1000);
                hold off;
                %}

                hold on;
                [C,hh] = contour(T_val,P_val,D_val,self.F_contourLines{prop}{spec},...
                            'k-', ...
                            'LineWidth',1.5, ...
                            'ShowText',~self.C_manual_levels{prop}, ...
                            'EdgeAlpha',0.2);
                if self.C_manual_levels{prop} == true
                    clabel(C, hh, 'manual', ...
                            'FontSize',24, ...
                            'Color','k', ...
                            'LabelSpacing',1000);
                else
                    clabel(C, hh, ...
                            'FontSize',24, ...
                            'Color','k', ...
                            'LabelSpacing',1000);
                end
                hold off;

                % save and close figure
                pause();
                if self.opts.save_images
                    dir = 'images/';
                    saveas(f, sprintf('%s/%s_2d_%s.png', dir, upper(data.species{spec}), self.F_filename{prop}), 'png');
                end

            end, end
        end

        function plot2d_zoom(self,spec2plot,props2plot)

            data = self.data;

            for spec = spec2plot, for prop = props2plot,

                % store data
                T_val = data.T{spec};
                P_val = data.P{spec};
                D_val = self.F{spec}{prop};

                % plot data
                h = pcolor(T_val,P_val,D_val);
                c = colorbar;

                % get handles to figure; h is surface
                g = h.Parent; % axes
                f = g.Parent; % figure
                f.Position(1) = 100;

                % run formatting functions
                self.format_axis(g,spec,prop);
                self.format_plot(h,spec,prop);
                self.format_colorbar(c,spec,prop);

                % turn off axis square
                % axis auto 

                % change ticks and tick labels
                g.XTickLabel = {'1', '1.025'};
                g.YTickLabel = {'1', '1.05', '1.1'};
                g.XTick = [g.XTick(1), g.XTick(1)*1.0249];
                g.YTick = [g.YTick(1), g.YTick(1)*1.05, g.YTick(1)*1.0999];

                % change figure size
                g.Position(3) = g.Position(4) / 2;
                c.Position(1) = 0.54;

                % plot Widom line
                key = self.key;
                hold on;
                    % using more accurate coolprop data:
                    plot(data.T_crit_line{spec}{key.cool}, data.P_crit_line{spec}{key.cool}, 'w--', 'LineWidth', 3);
                    f.Children(2).Children(1).Color = [0 0 0 0.9]; % make line semi-transparent
                hold off;

                % overlay contour plots
                hold on;
                %[C,hh] = contour(T_val,P_val,D_val,self.F_contourLines{prop}{spec},...
                
                % contours for cond_ref
                if prop == self.key.cond_ref
                    [C,hh] = contour(T_val,P_val,D_val,[0.05,0.06,0.08,0.1],...
                                'k-', ...
                                'LineWidth',1.5, ...
                                'ShowText',~self.C_manual_levels{prop}, ...
                                'EdgeAlpha',0.2);
                    if self.C_manual_levels{prop} == true
                        clabel(C, hh, 'manual', ...
                                'FontSize',24, ...
                                'Color','k', ...
                                'LabelSpacing',1000);
                    else
                        clabel(C, hh, ...
                                'FontSize',24, ...
                                'Color','k', ...
                                'LabelSpacing',1000);
                    end
                    hold off;
                end

                
                % conttours for cond_gen
                if prop == self.key.cond_gen
                    %[C,hh] = contour(T_val,P_val,D_val,[0.05,0.06,0.08,0.1],...
                    [C,hh] = contour(T_val,P_val,D_val,[0.035,0.04,0.05,0.06,0.065],...
                                'k-', ...
                                'LineWidth',1.5, ...
                                'ShowText',~self.C_manual_levels{prop}, ...
                                'EdgeAlpha',0.2);
                    if self.C_manual_levels{prop} == true
                        clabel(C, hh, 'manual', ...
                                'FontSize',24, ...
                                'Color','k', ...
                                'LabelSpacing',1000);
                    else
                        clabel(C, hh, ...
                                'FontSize',24, ...
                                'Color','k', ...
                                'LabelSpacing',1000);
                    end
                    hold off;
                end

               %[C,hh] = contour(T_val,P_val,D_val,[0.035,0.04,0.05,0.06,0.065],...

                % contours for cond_err
                if prop == key.cond_err_tot
                    [C,hh] = contour(T_val,P_val,D_val,[-.1,-.2,-.3,-.4,-.5,-.6],...
                                'k-', ...
                                'LineWidth',1.5, ...
                                'ShowText',~self.C_manual_levels{prop}, ...
                                'EdgeAlpha',0.2);
                    if self.C_manual_levels{prop} == true
                        clabel(C, hh, 'manual', ...
                                'FontSize',24, ...
                                'Color','k', ...
                                'LabelSpacing',1000);
                    else
                        clabel(C, hh, ...
                                'FontSize',24, ...
                                'Color','k', ...
                                'LabelSpacing',1000);
                    end
                    hold off;
                end

                % save and close figure
                pause();
                if self.opts.save_images
                    dir = 'images/';
                    saveas(f, sprintf('%s/%s_2d_%s.png', dir, upper(data.species{spec}), self.F_filename{prop}), 'png');
                end

            end, end
        end

        function format_axis(self,ax,spec,prop)

            % make axis square
            axis square

            % format color scheme
            ax.Colormap = colormap_creator(self.C_scheme{prop});
            %ax.CLim = [self.data.C_min{spec}{prop}, self.data.C_max{spec}{prop}];
            ax.CLim = [self.C_min{spec}{prop}, self.C_max{spec}{prop}];

            % format font
            ax.FontName = 'Roboto';
            ax.FontSize = 32;

            % format title
            if self.opts.titles_on
            ax.TitleFontSizeMultiplier = 1.3;
            ax.Title.String = sprintf('%s, $%s$', self.data.species_str{spec}, self.F_name{prop});
            ax.Title.Interpreter = 'latex';
            end

            % format x and y axes
            ax.LineWidth = 5;
            ax.LabelFontSizeMultiplier = 1.2;
            ax.XTick = round(self.data.T_crit(spec) * [1:1:3],3);
            ax.XTickLabel = [1:1:3];
            ax.XLabel.String = '$T^*$';
            ax.XLabel.Interpreter = 'latex';
            ax.YTick = round(self.data.P_crit(spec) * [1:1:3],3);
            ax.YTickLabel = [1:1:3];
            ax.YLabel.String = '$P^*$';
            ax.YLabel.Interpreter = 'latex';

            % KLUDGE: CO2's critical point was mistakenly recorded as 73.8080
            %   atm instead of the correct 72.8080 atm. The critical point has
            %   been corrected in this NRELplots object. But the data for CO2
            %   were generated with this incorrect value, so we have to adjust
            %   the tick marks here to match the data. If ever the data are
            %   regenerated, they will be generated with the correct critical
            %   pressure, and this kludge should be removed. Also, note that the
            %   data for the "zoomed in" plots were generated with the correct
            %   critical pressure, so this kludge is not applied in those plots.
            if spec == self.key.co2 && self.opts.zoom == 0
                ax.YTick = round((self.data.P_crit(spec) + 1) * [1:1:3],3);
            end

        end

        function format_plot(self,h,spec,prop)
            
            % adjust plot 
            h.EdgeAlpha = 0;
            %h.FaceColor = 'interp';

        end

        function format_colorbar(self,c,spec,prop)

            % adjust colorbar
            c.LineWidth = 1;
            c.Position(1) = c.Position(1) + 0.06;

            c.Label.String = self.C_label{prop};
            c.Label.Interpreter = 'latex';
            c.Label.FontSize = 32;
            c.Label.FontWeight = 'normal';
            c.Label.Rotation = 90;
            c.Label.Position(1) = 5.7;
            %c.Label.Position(1) = 3.2;
            %c.Label.Position(2) = -0.184;
            %c.Label.Position(2) = self.data.C_min{spec}{prop};

        end


        % DATA FUNCTIONS
        function val = max(self,prop,spec)

            % process inputs (default to all species)
            if nargin < 3, spec = 1:3; end
            
            % instantiate array
            val = zeros(length(prop),length(spec));

            % store values
            for i = 1:length(prop), for j = 1:length(spec),
                val(i,j) = max(max(self.F{spec(j)}{prop(i)}));
            end, end
        end

        function val = min(self,prop,spec)

            % process inputs (default to all species)
            if nargin < 3, spec = 1:3; end
            
            % instantiate array
            val = zeros(length(prop),length(spec));

            % store values
            for i = 1:length(prop), for j = 1:length(spec),
                val(i,j) = min(min(self.F{spec(j)}{prop(i)}));
            end, end
        end

        function val = max_excludeCritical(self,prop,spec,r)

            % process inputs (default to all species)
            if nargin < 3, spec = 1:3; end
            if nargin < 4, r = 1; end
            
            % instantiate array
            val = zeros(length(prop),length(spec));

            % store values
            for i = 1:length(prop), for j = 1:length(spec),

                % exclude critical region
                temp = self.F{spec(j)}{prop(i)};
                temp = self.excludeCritical(temp,r);
                val(i,j) = max(max(temp));
            end, end
        end

        function val = min_excludeCritical(self,prop,spec,r)

            % process inputs (default to all species)
            if nargin < 3, spec = 1:3; end
            if nargin < 4, r = 1; end
            
            % instantiate array
            val = zeros(length(prop),length(spec));

            % store values
            for i = 1:length(prop), for j = 1:length(spec),
                % exclude critical region
                temp = self.F{spec(j)}{prop(i)};
                temp = self.excludeCritical(temp,r);
                val(i,j) = min(min(temp));
            end, end
        end

        function val = excludeCritical(self,A,r)

            % get size of input array
            [n1,n2] = size(A);
            assert(n1 == n2, 'A must be n×n');
            n = n1;

            % Convert normalized radius r (in T*,P* space) to index-space radius
            r_idx = ((n - 1) * r) / 2;

            % get indices to exclude
            [i, j] = ndgrid(1:n, 1:n);
            d2_idx = (i - 1).^2 + (j - 1).^2;
            mask = d2_idx < r_idx^2;   % use < if you want to keep the boundary

            % exclude values
            val = A;
            val(mask) = NaN;

        end

        function val = TP_at_min(self,prop,spec)

            % process inputs (default to all species)
            if nargin < 3, spec = 1:3; end
            
            % instantiate array
            val = zeros(4,length(spec));

            % store values
            for i = 1:length(prop), for j = 1:length(spec),

                % store main values
                [min_val, min_idx] = min(min(self.F{spec(j)}{prop(i)}));
                [row, col] = find(self.F{spec(j)}{prop(i)} == min_val);
                T_val = self.data.T{spec(j)}(col(1));
                P_val = self.data.P{spec(j)}(row(1));
                val(1:2,j) = [T_val; P_val];

                % normalize by critical values
                T_star = T_val / self.data.T_crit(spec(j));
                P_star = P_val / self.data.P_crit(spec(j));
                val(3:4,j) = [T_star; P_star];

            end, end

            % convert to table
            val = array2table(val, ...
                    'VariableNames', self.data.species(spec), ...
                    'RowNames', {'T_min', 'P_min', '(T_min)^*', '(P_min)^*'});
        end

        function val = TP_at_max(self,prop,spec)

            % process inputs (default to all species)
            if nargin < 3, spec = 1:3; end
            
            % instantiate array
            val = zeros(4,length(spec));

            % store values
            for i = 1:length(prop), for j = 1:length(spec),

                % store main values
                [min_val, min_idx] = max(max(self.F{spec(j)}{prop(i)}));
                [row, col] = find(self.F{spec(j)}{prop(i)} == min_val);
                T_val = self.data.T{spec(j)}(col(1));
                P_val = self.data.P{spec(j)}(row(1));
                val(1:2,j) = [T_val; P_val];

                % normalize by critical values
                T_star = T_val / self.data.T_crit(spec(j));
                P_star = P_val / self.data.P_crit(spec(j));
                val(3:4,j) = [T_star; P_star];

            end, end

            % convert to table
            val = array2table(val, ...
                    'VariableNames', self.data.species(spec), ...
                    'RowNames', {'T_max', 'P_max', '(T_max)^*', '(P_max)^*'});
        end




        % GETTERS
        function self = compute_F(self)

            D_pele = self.data.D_pele;
            D_cool = self.data.D_cool;
            D_test = self.data.D_test;
            key = self.key;

            % store data
            for spec = 1:3, 

                % compute density properties
                F{spec}{key.dens_gen} = D_pele{spec}{key.rho};
                F{spec}{key.dens_ref} = D_cool{spec}{key.rho};
                F{spec}{key.dens_ERR} = F{spec}{key.dens_gen} - F{spec}{key.dens_ref};
                F{spec}{key.dens_err} = F{spec}{key.dens_ERR} ./ F{spec}{key.dens_ref};

                % compute viscosity properties
                F{spec}{key.visc_gen} = D_pele{spec}{key.mu};
                F{spec}{key.visc_ref} = D_cool{spec}{key.mu};
                F{spec}{key.visc_hyb} = D_test{spec}{key.mu};
                F{spec}{key.visc_ERR_tot} = F{spec}{key.visc_gen} - F{spec}{key.visc_ref};
                F{spec}{key.visc_err_tot} = F{spec}{key.visc_ERR_tot} ./ F{spec}{key.visc_ref};
                F{spec}{key.visc_ERR_corr} = F{spec}{key.visc_hyb} - F{spec}{key.visc_ref};
                F{spec}{key.visc_err_corr} = F{spec}{key.visc_ERR_corr} ./ F{spec}{key.visc_ref};
                F{spec}{key.visc_dmudrho} = D_test{spec}{key.dmu};
                F{spec}{key.visc_mismatch} = F{spec}{key.visc_hyb} ./ F{spec}{key.visc_ref};
                F{spec}{key.visc_ERR_eos} = F{spec}{key.dens_ERR} .* F{spec}{key.visc_dmudrho};
                F{spec}{key.visc_err_eos} = F{spec}{key.dens_err} .* F{spec}{key.visc_mismatch} .* F{spec}{key.dens_ref} .* F{spec}{key.visc_dmudrho} ./ F{spec}{key.visc_hyb};
                F{spec}{key.visc_ERR_assembled} = F{spec}{key.visc_ERR_corr} + F{spec}{key.visc_ERR_eos};
                F{spec}{key.visc_err_assembled} = F{spec}{key.visc_err_corr} + F{spec}{key.visc_err_eos};
                F{spec}{key.visc_ERR_secondorder} = (F{spec}{key.visc_gen} - F{spec}{key.visc_hyb}) - F{spec}{key.visc_ERR_eos};
                F{spec}{key.visc_err_secondorder} = (F{spec}{key.visc_gen} - F{spec}{key.visc_hyb}) ./ F{spec}{key.visc_ref} - F{spec}{key.visc_err_eos};
                F{spec}{key.visc_density_FLOOR} = -F{spec}{key.visc_ref} ./ D_pele{spec}{key.dmu};
                F{spec}{key.visc_density_floor} = F{spec}{key.visc_density_FLOOR} ./ F{spec}{key.dens_ref};
                F{spec}{key.visc_kinematic_ratio} = (F{spec}{key.visc_ref} ./ F{spec}{key.visc_gen}) ./ (F{spec}{key.dens_ref} ./ F{spec}{key.dens_gen});

                % compute conductivity properties
                F{spec}{key.cond_gen} = D_pele{spec}{key.lam};
                F{spec}{key.cond_ref} = D_cool{spec}{key.lam};
                F{spec}{key.cond_hyb} = D_test{spec}{key.lam};
                F{spec}{key.cond_ERR_tot} = F{spec}{key.cond_gen} - F{spec}{key.cond_ref};
                F{spec}{key.cond_err_tot} = F{spec}{key.cond_ERR_tot} ./ F{spec}{key.cond_ref};
                F{spec}{key.cond_ERR_corr} = F{spec}{key.cond_hyb} - F{spec}{key.cond_ref};
                F{spec}{key.cond_err_corr} = F{spec}{key.cond_ERR_corr} ./ F{spec}{key.cond_ref};
                F{spec}{key.cond_dlamdrho} = D_test{spec}{key.dlam};
                F{spec}{key.cond_mismatch} = F{spec}{key.cond_hyb} ./ F{spec}{key.cond_ref};
                F{spec}{key.cond_ERR_eos} = F{spec}{key.dens_ERR} .* F{spec}{key.cond_dlamdrho};
                F{spec}{key.cond_err_eos} = F{spec}{key.dens_err} .* F{spec}{key.cond_mismatch} .* F{spec}{key.dens_ref} .* F{spec}{key.cond_dlamdrho} ./ F{spec}{key.cond_hyb};
                F{spec}{key.cond_ERR_assembled} = F{spec}{key.cond_ERR_corr} + F{spec}{key.cond_ERR_eos};
                F{spec}{key.cond_err_assembled} = F{spec}{key.cond_err_corr} + F{spec}{key.cond_err_eos};
                F{spec}{key.cond_ERR_secondorder} = (F{spec}{key.cond_gen} - F{spec}{key.cond_hyb}) - F{spec}{key.cond_ERR_eos};
                F{spec}{key.cond_err_secondorder} = (F{spec}{key.cond_gen} - F{spec}{key.cond_hyb}) ./ F{spec}{key.cond_ref} - F{spec}{key.cond_err_eos};
                F{spec}{key.cond_density_FLOOR} = -F{spec}{key.cond_ref} ./ D_pele{spec}{key.dlam};
                F{spec}{key.cond_density_floor} = F{spec}{key.cond_density_FLOOR} ./ F{spec}{key.dens_ref};

                % convert absolute conductivity errors to mW/m_K
                F{spec}{key.cond_ERR_tot} = 1000 * F{spec}{key.cond_ERR_tot};
                F{spec}{key.cond_ERR_corr} = 1000 * F{spec}{key.cond_ERR_corr};
                F{spec}{key.cond_ERR_eos} = 1000 * F{spec}{key.cond_ERR_eos};
                F{spec}{key.cond_ERR_assembled} = 1000 * F{spec}{key.cond_ERR_assembled};
                F{spec}{key.cond_ERR_secondorder} = 1000 * F{spec}{key.cond_ERR_secondorder};   

            end
            
            % store result
            self.F = F;

        end

        function key = get.key(self)

            % species
            i = 1;
            key.ch4 = i;    i = i + 1;
            key.co2 = i;    i = i + 1;
            key.n2 = i;     i = i + 1;

            % legacy properties
            i = 1;
            key.rho = i;    i = i + 1;
            key.mu = i;     i = i + 1;
            key.lam = i;    i = i + 1;
            key.Cp = i;     i = i + 1;
            key.dmu = i;    i = i + 1;
            key.dlam = i;   i = i + 1;

            % sources
            i = 1;
            key.pele = i;                   i = i + 1;
            key.cool = i;                   i = i + 1;
            key.err = i;                    i = i + 1;
            key.re = i;                     i = i + 1;

            % density property key
            i = 1;
            key.dens_gen = i;               i = i + 1;
            key.dens_ref = i;               i = i + 1;
            key.dens_ERR = i;               i = i + 1;
            key.dens_err = i;               i = i + 1;

            % viscosity property key
            key.visc_gen = i;               i = i + 1;
            key.visc_ref = i;               i = i + 1;
            key.visc_hyb = i;               i = i + 1;
            key.visc_ERR_tot = i;           i = i + 1;
            key.visc_err_tot = i;           i = i + 1;
            key.visc_ERR_corr = i;          i = i + 1;
            key.visc_err_corr = i;          i = i + 1;
            key.visc_dmudrho = i;           i = i + 1;
            key.visc_mismatch = i;          i = i + 1;
            key.visc_ERR_eos = i;           i = i + 1;
            key.visc_err_eos = i;           i = i + 1;
            key.visc_ERR_assembled = i;     i = i + 1;
            key.visc_err_assembled = i;     i = i + 1;
            key.visc_ERR_secondorder = i;   i = i + 1;
            key.visc_err_secondorder = i;   i = i + 1;
            key.visc_density_FLOOR = i;     i = i + 1;
            key.visc_density_floor = i;     i = i + 1;
            key.visc_kinematic_ratio = i;   i = i + 1;

            % conductivity property key
            key.cond_gen = i;               i = i + 1;
            key.cond_ref = i;               i = i + 1;
            key.cond_hyb = i;               i = i + 1;
            key.cond_ERR_tot = i;           i = i + 1;
            key.cond_err_tot = i;           i = i + 1;
            key.cond_ERR_corr = i;          i = i + 1;
            key.cond_err_corr = i;          i = i + 1;
            key.cond_dlamdrho = i;          i = i + 1;
            key.cond_mismatch = i;          i = i + 1;
            key.cond_ERR_eos = i;           i = i + 1;    
            key.cond_err_eos = i;           i = i + 1;
            key.cond_ERR_assembled = i;     i = i + 1;
            key.cond_err_assembled = i;     i = i + 1;
            key.cond_ERR_secondorder = i;   i = i + 1;
            key.cond_err_secondorder = i;   i = i + 1;
            key.cond_density_FLOOR = i;     i = i + 1;
            key.cond_density_floor = i;     i = i + 1;

        end

        function F_name = get.F_name(self)

            % setup
            F_name = {};
            key = self.key;

            % density
            F_name{key.dens_gen} = '\rho_\textrm{gen}';
            F_name{key.dens_ref} = '\rho_\textrm{ref}';
            F_name{key.dens_ERR} = '\Delta \rho_\textrm{tot}';
            F_name{key.dens_err} = '\delta \rho_\textrm{tot}';

            % viscosity
            F_name{key.visc_gen} = '\eta_\textrm{gen}';
            F_name{key.visc_ref} = '\eta_\textrm{ref}';
            F_name{key.visc_hyb} = '\eta_\textrm{test}';
            F_name{key.visc_ERR_tot} = '\Delta \eta_\textrm{tot}';
            F_name{key.visc_err_tot} = '\delta \eta_\textrm{tot}';
            F_name{key.visc_ERR_corr} = '\Delta \eta_\textrm{corr}';
            F_name{key.visc_err_corr} = '\delta \eta_\textrm{corr}';
            %F_name{key.visc_dmudrho} = '\partial \hat{\eta}_\textrm{gen} / \partial \rho';
            F_name{key.visc_dmudrho} = 'S_1';
            F_name{key.visc_mismatch} = '\eta_\textrm{hyb} / \eta_\textrm{ref}';
            F_name{key.visc_ERR_eos} = '\Delta \eta_\textrm{eos}^{(1)}';
            F_name{key.visc_err_eos} = '\delta \eta_\textrm{eos}^{(1)}';
            F_name{key.visc_ERR_assembled} = '\Delta \eta_\textrm{corr} + \Delta \eta_\textrm{eos}^{(1)}';
            F_name{key.visc_err_assembled} = '\delta \eta_\textrm{corr} + \delta \eta_\textrm{eos}^{(1)}';
            F_name{key.visc_ERR_secondorder} = '\Delta \eta_\textrm{eos}^{(2)}';
            F_name{key.visc_err_secondorder} = '\delta \eta_\textrm{eos}^{(2)}';
            F_name{key.visc_density_FLOOR} = '\Delta \rho_*^{(\eta)}';
            F_name{key.visc_density_floor} = '\delta \rho_*^{(\eta)}';
            F_name{key.visc_kinematic_ratio} = '\nu_\textrm{ref} / \nu_\textrm{gen}';

            % conductivity
            F_name{key.cond_gen} = '\lambda_\textrm{gen}';
            F_name{key.cond_ref} = '\lambda_\textrm{ref}';
            F_name{key.cond_hyb} = '\lambda_\textrm{test}';
            F_name{key.cond_ERR_tot} = '\Delta \lambda_\textrm{tot}';
            F_name{key.cond_err_tot} = '\delta \lambda_\textrm{tot}';
            F_name{key.cond_ERR_corr} = '\Delta \lambda_\textrm{corr}';
            F_name{key.cond_err_corr} = '\delta \lambda_\textrm{corr}';
            %F_name{key.cond_dlamdrho} = '\partial \lambda / \partial \rho';
            F_name{key.cond_dlamdrho} = 'S_2';
            F_name{key.cond_mismatch} = '\lambda_\textrm{hyb} / \lambda_\textrm{ref}';
            F_name{key.cond_ERR_eos} = '\Delta \lambda_\textrm{eos}^{(1)}';
            F_name{key.cond_err_eos} = '\delta \lambda_\textrm{eos}^{(1)}';
            F_name{key.cond_ERR_assembled} = '\Delta \lambda_\textrm{corr} + \Delta \lambda_\textrm{eos}^{(1)}';
            F_name{key.cond_err_assembled} = '\delta \lambda_\textrm{corr} + \delta \lambda_\textrm{eos}^{(1)}';
            F_name{key.cond_ERR_secondorder} = '\Delta \lambda_\textrm{eos}^{(2)}';
            F_name{key.cond_err_secondorder} = '\delta \lambda_\textrm{eos}^{(2)}';
            F_name{key.cond_density_FLOOR} = '\Delta \rho_*^{(\lambda)}';
            F_name{key.cond_density_floor} = '\delta \rho_*^{(\lambda)}';

        end

        function F_filename = get.F_filename(self)

            % setup
            F_filename = {};
            key = self.key;

            % density
            F_filename{key.dens_gen} = 'dens_gen';
            F_filename{key.dens_ref} = 'dens_ref';
            F_filename{key.dens_ERR} = 'dens_ERR';
            F_filename{key.dens_err} = 'dens_err';

            % viscosity
            F_filename{key.visc_gen} = 'visc_gen';
            F_filename{key.visc_ref} = 'visc_ref';
            F_filename{key.visc_hyb} = 'visc_hyb';
            F_filename{key.visc_ERR_tot} = 'visc_ERR_tot';
            F_filename{key.visc_err_tot} = 'visc_err_tot';
            F_filename{key.visc_ERR_corr} = 'visc_ERR_corr';
            F_filename{key.visc_err_corr} = 'visc_err_corr';
            F_filename{key.visc_dmudrho} = 'visc_dmudrho';
            F_filename{key.visc_mismatch} = 'visc_mismatch';
            F_filename{key.visc_ERR_eos} = 'visc_ERR_eos';
            F_filename{key.visc_err_eos} = 'visc_err_eos';
            F_filename{key.visc_ERR_assembled} = 'visc_ERR_assembled';
            F_filename{key.visc_err_assembled} = 'visc_err_assembled';
            F_filename{key.visc_ERR_secondorder} = 'visc_ERR_secondorder';
            F_filename{key.visc_err_secondorder} = 'visc_err_secondorder';
            F_filename{key.visc_density_FLOOR} = 'visc_density_FLOOR';
            F_filename{key.visc_density_floor} = 'visc_density_floor';
            F_filename{key.visc_kinematic_ratio} = 'visc_kinematic_ratio';

            % conductivity
            F_filename{key.cond_gen} = 'cond_gen';
            F_filename{key.cond_ref} = 'cond_ref';
            F_filename{key.cond_hyb} = 'cond_hyb';
            F_filename{key.cond_ERR_tot} = 'cond_ERR_tot';
            F_filename{key.cond_err_tot} = 'cond_err_tot';
            F_filename{key.cond_ERR_corr} = 'cond_ERR_corr';
            F_filename{key.cond_err_corr} = 'cond_err_corr';
            F_filename{key.cond_dlamdrho} = 'cond_dlamdrho';
            F_filename{key.cond_mismatch} = 'cond_mismatch';
            F_filename{key.cond_ERR_eos} = 'cond_ERR_eos';  
            F_filename{key.cond_err_eos} = 'cond_err_eos';
            F_filename{key.cond_ERR_assembled} = 'cond_ERR_assembled';
            F_filename{key.cond_err_assembled} = 'cond_err_assembled';
            F_filename{key.cond_ERR_secondorder} = 'cond_ERR_secondorder';
            F_filename{key.cond_err_secondorder} = 'cond_err_secondorder';
            F_filename{key.cond_density_FLOOR} = 'cond_density_FLOOR';
            F_filename{key.cond_density_floor} = 'cond_density_floor';

        end

        function C_scheme = get.C_scheme(self)

            % setup
            C_scheme = {};
            key = self.key;

            % color scheme names
            fast = sprintf('fast%s','1024');
            bluered = sprintf('bluered%s','1024');
            kindlemann = sprintf('kindlmann%s','1024');
            kindlemannE = sprintf('kindlmannextended%s','1024');
            inferno = sprintf('inferno%s','1024');
            blackbody = sprintf('blackbody%s','1024');
            plasma = sprintf('plasma%s','1024');

            kindlemann = blackbody;

            % density
            C_scheme{key.dens_gen} = fast;
            C_scheme{key.dens_ref} = fast;
            C_scheme{key.dens_ERR} = bluered;
            C_scheme{key.dens_err} = bluered;

            % viscosity
            C_scheme{key.visc_gen} = fast;
            C_scheme{key.visc_ref} = fast;
            C_scheme{key.visc_hyb} = fast;
            C_scheme{key.visc_ERR_tot} = bluered;
            C_scheme{key.visc_err_tot} = bluered;
            C_scheme{key.visc_ERR_corr} = bluered;
            C_scheme{key.visc_err_corr} = bluered;
            C_scheme{key.visc_dmudrho} = fast;
            C_scheme{key.visc_mismatch} = bluered;
            C_scheme{key.visc_ERR_eos} = bluered;
            C_scheme{key.visc_err_eos} = bluered;
            C_scheme{key.visc_ERR_assembled} = bluered;
            C_scheme{key.visc_err_assembled} = bluered;
            C_scheme{key.visc_ERR_secondorder} = bluered;
            C_scheme{key.visc_err_secondorder} = bluered;
            C_scheme{key.visc_density_FLOOR} = kindlemann;
            C_scheme{key.visc_density_floor} = kindlemann;
            C_scheme{key.visc_kinematic_ratio} = bluered;

            % conductivity
            C_scheme{key.cond_gen} = fast;
            C_scheme{key.cond_ref} = fast;
            C_scheme{key.cond_hyb} = fast;
            C_scheme{key.cond_ERR_tot} = bluered;
            C_scheme{key.cond_err_tot} = bluered;
            C_scheme{key.cond_ERR_corr} = bluered;
            C_scheme{key.cond_err_corr} = bluered;
            C_scheme{key.cond_dlamdrho} = fast;
            C_scheme{key.cond_mismatch} = bluered;
            C_scheme{key.cond_ERR_eos} = bluered;
            C_scheme{key.cond_err_eos} = bluered;
            C_scheme{key.cond_ERR_assembled} = bluered;
            C_scheme{key.cond_err_assembled} = bluered;
            C_scheme{key.cond_ERR_secondorder} = bluered;
            C_scheme{key.cond_err_secondorder} = bluered;
            C_scheme{key.cond_density_FLOOR} = kindlemann;
            C_scheme{key.cond_density_floor} = kindlemann;

        end

        function C_label = get.C_label(self)

            % setup
            C_label = {};
            key = self.key;

            % density
            C_label{key.dens_gen} = 'Density [kg/m$^3$]';
            C_label{key.dens_ref} = 'Density [kg/m$^3$]';
            C_label{key.dens_ERR} = 'Error [kg/m$^3$]';
            C_label{key.dens_err} = 'Relative Error';

            % viscosity
            C_label{key.visc_gen} = 'Viscosity [$\mu$Pa/s]';
            C_label{key.visc_ref} = 'Viscosity [$\mu$Pa/s]';
            C_label{key.visc_hyb} = 'Viscosity [$\mu$Pa/s]';
            C_label{key.visc_ERR_tot} = 'Error [$\mu$Pa/s]';
            C_label{key.visc_err_tot} = 'Relative Error';
            C_label{key.visc_ERR_corr} = 'Error [$\mu$Pa/s]';
            C_label{key.visc_err_corr} = 'Relative Error';
            C_label{key.visc_dmudrho} = '';
            C_label{key.visc_mismatch} = '';
            C_label{key.visc_ERR_eos} = 'Error [$\mu$Pa/s]';
            C_label{key.visc_err_eos} = 'Relative Error';
            C_label{key.visc_ERR_assembled} = 'Error [$\mu$Pa/s]';
            C_label{key.visc_err_assembled} = 'Relative Error';
            C_label{key.visc_ERR_secondorder} = 'Error [$\mu$Pa/s]';
            C_label{key.visc_err_secondorder} = 'Relative Error';
            C_label{key.visc_density_FLOOR} = 'Error [kg/m$^3$]';
            C_label{key.visc_density_floor} = 'Relative Error';
            C_label{key.visc_kinematic_ratio} = 'Ratio';

            % conductivity
            C_label{key.cond_gen} = 'Conductivity [W/m-K]';
            C_label{key.cond_ref} = 'Conductivity [W/m-K]';
            C_label{key.cond_hyb} = 'Conductivity [W/m-K]';
            %C_label{key.cond_ERR_tot} = 'Error [W/m-K]';
            C_label{key.cond_ERR_tot} = 'Error [mW/m-K]';
            C_label{key.cond_err_tot} = 'Relative Error';
            %C_label{key.cond_ERR_corr} = 'Error [W/m-K]';
            C_label{key.cond_ERR_corr} = 'Error [mW/m-K]';
            C_label{key.cond_err_corr} = 'Relative Error';
            C_label{key.cond_dlamdrho} = '';
            C_label{key.cond_mismatch} = '';
            %C_label{key.cond_ERR_eos} = 'Error [W/m-K]';
            C_label{key.cond_ERR_eos} = 'Error [mW/m-K]';
            C_label{key.cond_err_eos} = 'Relative Error';
            %C_label{key.cond_ERR_assembled} = 'Error [W/m-K]';
            C_label{key.cond_ERR_assembled} = 'Error [mW/m-K]';
            C_label{key.cond_err_assembled} = 'Relative Error';
            %C_label{key.cond_ERR_secondorder} = 'Error [W/m-K]';
            C_label{key.cond_ERR_secondorder} = 'Error [mW/m-K]';
            C_label{key.cond_err_secondorder} = 'Relative Error';
            C_label{key.cond_density_FLOOR} = 'Error [kg/m$^3$]';
            C_label{key.cond_density_floor} = 'Relative Error';

        end

        function unit = get.F_unit(self)

            % setup
            unit = {};
            key = self.key;

            % density
            unit{key.dens_gen} = 'kg/m$^3$';
            unit{key.dens_ref} = 'kg/m$^3$';
            unit{key.dens_ERR} = 'kg/m$^3$';
            unit{key.dens_err} = '';

            % viscosity
            unit{key.visc_gen} = '$\mu$Pa/s';
            unit{key.visc_ref} = '$\mu$Pa/s';
            unit{key.visc_hyb} = '$\mu$Pa/s';
            unit{key.visc_ERR_tot} = '$\mu$Pa/s';
            unit{key.visc_err_tot} = '';
            unit{key.visc_ERR_corr} = '$\mu$Pa/s';
            unit{key.visc_err_corr} = '';
            unit{key.visc_dmudrho} = '';
            unit{key.visc_mismatch} = '';
            unit{key.visc_ERR_eos} = '$\mu$Pa/s';
            unit{key.visc_err_eos} = '';
            unit{key.visc_ERR_assembled} = '$\mu$Pa/s';
            unit{key.visc_err_assembled} = '';
            unit{key.visc_ERR_secondorder} = '$\mu$Pa/s';
            unit{key.visc_err_secondorder} = '';
            unit{key.visc_density_FLOOR} = 'kg/m$^3$';
            unit{key.visc_density_floor} = '';
            unit{key.visc_kinematic_ratio} = '';

            % conductivity
            unit{key.cond_gen} = 'W/m-K';
            unit{key.cond_ref} = 'W/m-K';
            unit{key.cond_hyb} = 'W/m-K';
            unit{key.cond_ERR_tot} = 'mW/m-K';
            unit{key.cond_err_tot} = '';
            unit{key.cond_ERR_corr} = 'mW/m-K';
            unit{key.cond_err_corr} = '';
            unit{key.cond_dlamdrho} = '';
            unit{key.cond_mismatch} = '';
            unit{key.cond_ERR_eos} = 'mW/m-K';
            unit{key.cond_err_eos} = '';
            unit{key.cond_ERR_assembled} = 'mW/m-K';
            unit{key.cond_err_assembled} = '';
            unit{key.cond_ERR_secondorder} = 'mW/m-K';
            unit{key.cond_err_secondorder} = '';
            unit{key.cond_density_FLOOR} = 'kg/m$^3$';
            unit{key.cond_density_floor} = '';

            % store result
            F_unit = unit;

        end

        function val = get.F_contourLines(self)
        % NOTE: contour line values are manually specified to give the "best
        % looking" plots. Values for gen, ref, and hyb are matched. 
        
            % setup
            val = {};
            key = self.key;

            % DENSITY
            % density, values
            val{key.dens_gen}{key.ch4} = [25,50,100,150,250];
            val{key.dens_gen}{key.co2} = [50,100,150,250,400,700];
            val{key.dens_gen}{key.n2}  = [50,100,150,250,500];
            val{key.dens_ref}{key.ch4} = val{key.dens_gen}{key.ch4};
            val{key.dens_ref}{key.co2} = val{key.dens_gen}{key.co2};
            val{key.dens_ref}{key.n2}  = val{key.dens_gen}{key.n2};

            % density, absolute error
            val{key.dens_ERR}{key.ch4} = [-0.1,-0.25,-0.5,-1,-2.5,-22];
            val{key.dens_ERR}{key.co2} = [-0.1,-0.25,-0.5,-1,-2.5,-20,-100];
            val{key.dens_ERR}{key.n2}  = [-0.3,-0.5,-1,-3,-40];

            % density, relative error
            val{key.dens_err}{key.ch4} = [-0.005,-0.01,-0.02,-0.10];
            val{key.dens_err}{key.co2} = [-0.015,-0.02,-0.03,-0.05,-0.15];
            val{key.dens_err}{key.n2}  = [-0.015,-0.02,-0.03,-0.05,-0.15];


            % VISCOSITY
            % viscosity, values
            val{key.visc_gen}{key.ch4} = [10 12.5 15 17.5 30];
            val{key.visc_ref}{key.ch4} = val{key.visc_gen}{key.ch4};
            val{key.visc_hyb}{key.ch4} = val{key.visc_gen}{key.ch4};
            val{key.visc_gen}{key.co2} = [25 30 35 70];
            val{key.visc_ref}{key.co2} = val{key.visc_gen}{key.co2};
            val{key.visc_hyb}{key.co2} = val{key.visc_gen}{key.co2};
            val{key.visc_gen}{key.n2} = [15 17.5 20 30 40];
            val{key.visc_ref}{key.n2} = val{key.visc_gen}{key.n2};
            val{key.visc_hyb}{key.n2} = val{key.visc_gen}{key.n2};

            % viscosity, absolute error
            val{key.visc_ERR_tot}{key.ch4} = [-2 0 .4 .6];
            val{key.visc_ERR_tot}{key.co2} = [-10 1 2];
            val{key.visc_ERR_tot}{key.n2}  = [-5 .1 .2];

            % viscosity, relative error
            val{key.visc_err_tot}{key.ch4} = [-.1 0 .01];
            val{key.visc_err_tot}{key.co2} = [-.12 0.04 0.12];
            val{key.visc_err_tot}{key.n2}  = [-0.1 .01];

            % viscosity, correlation absolute error
            val{key.visc_ERR_corr}{key.ch4} = [-.1 0 .1 0.3];
            val{key.visc_ERR_corr}{key.co2} = [.5 1 2 3 6];
            val{key.visc_ERR_corr}{key.n2}  = [-1.5 .1];

            % viscosity, correlation relative error
            val{key.visc_err_corr}{key.ch4} = [0 .01 0.04 0.08];
            val{key.visc_err_corr}{key.co2} = [0.03 0.06 0.09 0.12];
            val{key.visc_err_corr}{key.n2}  = [-0.03 0.01 0.04];

            % viscosity, dmudrho
            val{key.visc_dmudrho}{key.ch4} = [0.04 0.05 0.07 0.1 0.15 0.2];
            val{key.visc_dmudrho}{key.co2} = [0.04 0.05 0.06 0.2];
            val{key.visc_dmudrho}{key.n2}  = [0.02 0.03 0.04 0.06 .12];

            % viscosity, mismatch
            val{key.visc_mismatch}{key.ch4} = [5];
            val{key.visc_mismatch}{key.co2} = [5];
            val{key.visc_mismatch}{key.n2}  = [5];

            % viscosity equation of state absolute error
            val{key.visc_ERR_eos}{key.ch4} = [-4, -.4 -.05, -0.01];
            val{key.visc_ERR_eos}{key.co2} = [-.05 -.1 -1 -15];
            val{key.visc_ERR_eos}{key.n2}  = [-4 -.1 -.01];

            % viscosity equation of state relative error
            val{key.visc_err_eos}{key.ch4} = [-.12 -.06 -.01 -.001];
            val{key.visc_err_eos}{key.co2} = [-.25 -.1 -.01 -.001];
            val{key.visc_err_eos}{key.n2}  = [-.12 -.06 -.01 -.001];

            % viscosity assembled absolute error
            val{key.visc_ERR_assembled}{key.ch4} = [-3 0 .4];
            val{key.visc_ERR_assembled}{key.co2} = [-10 1 2];
            val{key.visc_ERR_assembled}{key.n2}  = [-4 0.1 .2];

            % viscosity assembled relative error
            val{key.visc_err_assembled}{key.ch4} = [-.1 0 0.01];
            val{key.visc_err_assembled}{key.co2} = [-.2 0.02 0.05 .1];
            val{key.visc_err_assembled}{key.n2}  = [-.12 .005 0.01 0.02];

            % viscosity second order absolute error
            val{key.visc_ERR_secondorder}{key.ch4} = [0.01 .1 .3];
            val{key.visc_ERR_secondorder}{key.co2} = [0.01 .1 2];
            val{key.visc_ERR_secondorder}{key.n2}  = [.01 .1 .3];

            % viscosity second order relative error
            val{key.visc_err_secondorder}{key.ch4} = [0.001 .01];
            val{key.visc_err_secondorder}{key.co2} = [0.001 0.01 0.04];
            val{key.visc_err_secondorder}{key.n2}  = [0.001 0.01 0.04];

            % viscosity density floor absolute error
            val{key.visc_density_FLOOR}{key.ch4} = -[180 250 350 450 550];
            val{key.visc_density_FLOOR}{key.co2} = -[360 500 600 800 1000 1200];
            val{key.visc_density_FLOOR}{key.n2}  = -[340 400 500 600 800 1000 1200];

            % viscosity density floor relative error
            val{key.visc_density_floor}{key.ch4} = -[0.65 1 2 3 5 10 20 40];
            val{key.visc_density_floor}{key.co2} = -[0.5 1 2 3 5 10 20];
            val{key.visc_density_floor}{key.n2}  = -[0.65 1 2 3 5 10 20];
             
            % viscosity kinematic ratio
            val{key.visc_kinematic_ratio}{key.ch4} = [0.92 1 1.05];
            val{key.visc_kinematic_ratio}{key.co2} = [0.85 0.9 0.95 1];
            val{key.visc_kinematic_ratio}{key.n2}  = [0.96 0.99 1 1.1];


            % CONDUCTIVITY
            % conductivity, value
            val{key.cond_gen}{key.ch4} = [0.03 0.04 0.05 0.06 0.08];
            val{key.cond_ref}{key.ch4} = val{key.cond_gen}{key.ch4};
            val{key.cond_hyb}{key.ch4} = val{key.cond_gen}{key.ch4};
            val{key.cond_gen}{key.co2} = [0.03 0.04 0.05 0.06 0.09 0.13];
            val{key.cond_ref}{key.co2} = val{key.cond_gen}{key.co2};
            val{key.cond_hyb}{key.co2} = val{key.cond_gen}{key.co2};
            val{key.cond_gen}{key.n2} = [0.025 0.03 0.035 0.06];
            val{key.cond_ref}{key.n2} = val{key.cond_gen}{key.n2};
            val{key.cond_hyb}{key.n2} = val{key.cond_gen}{key.n2};

            % conductivity, absolute error
            val{key.cond_ERR_tot}{key.ch4} = [-15 -5 -1 0 1];
            val{key.cond_ERR_tot}{key.co2} = [-5 0 1 3 5];
            val{key.cond_ERR_tot}{key.n2} = [-6 -4 0 0.5];

            % conductivity, relative error
            val{key.cond_err_tot}{key.ch4} = [-.2 -.1 0 .02];
            val{key.cond_err_tot}{key.co2} = [-.1 0 .05 .1 .15];
            val{key.cond_err_tot}{key.n2} = [-.1 0 .05];

            % conductivity, correlation absolute error
            val{key.cond_ERR_corr}{key.ch4} = [-15 -5 -1 0 1];
            val{key.cond_ERR_corr}{key.co2} = [-5 0 1 5 10 40];
            val{key.cond_ERR_corr}{key.n2} = [-4 0 1 4 8];

            % conductivity, correlation relative error
            val{key.cond_err_corr}{key.ch4} = [-.15 -.05 0 .02];
            val{key.cond_err_corr}{key.co2} = [-.05 .05 .1 .2 0.4];
            val{key.cond_err_corr}{key.n2} = [-.1 0 .05];

            % conductivity, dlamdrho
            val{key.cond_dlamdrho}{key.ch4} = [.00015 0.0002 .00045];
            val{key.cond_dlamdrho}{key.co2} = [0.0001 0.0002 0.00035];
            val{key.cond_dlamdrho}{key.n2} = [.00005 .00007 0.00015 ];

            % conductivity, mismatch
            val{key.cond_mismatch}{key.ch4} = [5];
            val{key.cond_mismatch}{key.co2} = [5];
            val{key.cond_mismatch}{key.n2} = [5];

            % conductivity equation of state absolute error
            val{key.cond_ERR_eos}{key.ch4} = [0 -.05 -.1 -1 -10];
            val{key.cond_ERR_eos}{key.co2} = [0 -0.05 -0.1 -1 -22];
            val{key.cond_ERR_eos}{key.n2} = [0 -0.01 -0.1 -1 -5];

            % conductivity equation of state relative error
            val{key.cond_err_eos}{key.ch4} = [-.1 -.01 -0.005 0];
            val{key.cond_err_eos}{key.co2} = [-.1 -.01 -0.005 0];
            val{key.cond_err_eos}{key.n2} = [-.1 -.01 -0.005 0];

            % conductivity assembled absolute error
            val{key.cond_ERR_assembled}{key.ch4} = [0 1 2 -1 -3 -10];
            val{key.cond_ERR_assembled}{key.co2} = [2 5 10 0 -4];
            val{key.cond_ERR_assembled}{key.n2} = [-7 -4 0 1 1.5];

            % conductivity assembled relative error
            val{key.cond_err_assembled}{key.ch4} = [-.18 -.1 0 .02];
            val{key.cond_err_assembled}{key.co2} = [-.1 0 .05 .1];
            val{key.cond_err_assembled}{key.n2} = [-.1 0 .05];

            % conductivity second order absolute error
            val{key.cond_ERR_secondorder}{key.ch4} = [0 -0.1 -0.01 -0.005];
            val{key.cond_ERR_secondorder}{key.co2} = [-0.1 -0.01];
            val{key.cond_ERR_secondorder}{key.n2} = [0 -0.45 -0.1 -0.01 -0.005];

            % conductivity second order relative error
            val{key.cond_err_secondorder}{key.ch4} = [-.001 -.004 -.0005 0];
            val{key.cond_err_secondorder}{key.co2} = [-.001 -.004];
            val{key.cond_err_secondorder}{key.n2} = [-.001 -.004 -.0005 0];            

            % conductivity density floor absolute error
            val{key.cond_density_FLOOR}{key.ch4} = [-200 -250 -300 -400 -500];
            val{key.cond_density_FLOOR}{key.co2} = [-400 -500 -600 -700 -800];
            val{key.cond_density_FLOOR}{key.n2} = [-300 -400 -500 -600 -700];

            % conductivity density floor relative error
            val{key.cond_density_floor}{key.ch4} = [-0.5 -1 -2 -3 -5 -10 -20];
            val{key.cond_density_floor}{key.co2} = [-0.5 -1 -2 -3 -5 -10 -20];
            val{key.cond_density_floor}{key.n2} = [-0.5 -1 -2 -3 -5 -10 -20];

        end

        function val = get.C_manual_levels(self)

            % setup
            unit = {};
            key = self.key;

            % values
            %   note: 'false' means 'don't use manual levels', i.e. 'use automatic levels'
            auto = false;
            manual = true;

            % density
            unit{key.dens_gen} = auto;
            unit{key.dens_ref} = auto;
            unit{key.dens_ERR} = manual;
            unit{key.dens_err} = manual;

            % viscosity
            unit{key.visc_gen} = auto;
            unit{key.visc_ref} = auto;
            unit{key.visc_hyb} = auto;
            unit{key.visc_ERR_tot} = manual;
            unit{key.visc_err_tot} = manual;
            unit{key.visc_ERR_corr} = manual;
            unit{key.visc_err_corr} = manual;
            unit{key.visc_dmudrho} = auto;
            unit{key.visc_mismatch} = auto;
            unit{key.visc_ERR_eos} = manual;
            unit{key.visc_err_eos} = manual;
            unit{key.visc_ERR_assembled} = manual;
            unit{key.visc_err_assembled} = manual;
            unit{key.visc_ERR_secondorder} = manual;
            unit{key.visc_err_secondorder} = manual;
            unit{key.visc_density_FLOOR} = auto;
            unit{key.visc_density_floor} = auto;
            unit{key.visc_kinematic_ratio} = manual;

            % conductivity
            unit{key.cond_gen} = auto;
            unit{key.cond_ref} = auto;
            unit{key.cond_hyb} = auto;
            unit{key.cond_ERR_tot} = manual;
            unit{key.cond_err_tot} = manual;
            unit{key.cond_ERR_corr} = manual;
            unit{key.cond_err_corr} = manual;
            unit{key.cond_dlamdrho} = auto;
            unit{key.cond_mismatch} = manual;
            unit{key.cond_ERR_eos} = manual;
            unit{key.cond_err_eos} = manual;
            unit{key.cond_ERR_assembled} = manual;
            unit{key.cond_err_assembled} = manual;
            unit{key.cond_ERR_secondorder} = manual;
            unit{key.cond_err_secondorder} = manual;
            unit{key.cond_density_FLOOR} = auto;
            unit{key.cond_density_floor} = auto;
            %{
            unit{key.cond_ERR_tot} = auto;
            unit{key.cond_err_tot} = auto;
            unit{key.cond_ERR_corr} = auto;
            unit{key.cond_err_corr} = auto;
            unit{key.cond_dlamdrho} = auto;
            unit{key.cond_mismatch} = auto;
            unit{key.cond_ERR_eos} = auto;
            unit{key.cond_err_eos} = auto;
            unit{key.cond_ERR_assembled} = auto;
            unit{key.cond_err_assembled} = auto;
            unit{key.cond_ERR_secondorder} = auto;
            unit{key.cond_err_secondorder} = auto;
            unit{key.cond_density_FLOOR} = auto;
            unit{key.cond_density_floor} = auto;
            %}

            % store result
            val = unit;

        end

    end % end methods


    methods (Static)

        function P = read_pres(data,dir)

            P = {};

            for i = 1:length(data.species)
                P_data_str = sprintf('%s/%s_2d_pressure.txt', dir, upper(data.species{i}));
                P_tab = readtable(P_data_str,'ReadVariableNames',false);
                P{i} = table2array(P_tab)';
            end

        end

        function T = read_temp(data,dir)

            T = {};

            for i = 1:length(data.species)
                T_data_str = sprintf('%s/%s_2d_temperature.txt', dir, upper(data.species{i}));
                T_tab = readtable(T_data_str,'ReadVariableNames',false);
                T{i} = table2array(T_tab)';
            end

        end

    end
end
