classdef MrProtocol
% MRPROTOCOL Represents an MRI diffusion protocol with computed q-space parameters.
%
% This class encapsulates the diffusion MRI acquisition parameters, including 
% gradient strengths, durations, separations, and computed q-space values.
% It provides methods for constructing an instance from input parameters 
% or loading from a JSON file.
%
% PROPERTIES:
%   gradient_amplitude_per_shell     - (Nx1 vector) Gradient amplitudes per shell.
%   gradient_duration_per_shell      - (Nx1 vector) Gradient durations per shell.
%   gradient_separation_per_shell    - (Nx1 vector) Gradient separations per shell.
%   q_per_shell                      - (Nx1 vector) Computed q-values per shell.
%   bval_per_shell                   - (Nx1 vector) b-values per shell.
%   bvec_per_shell                   - (Nx1 cell) Normalized gradient directions per shell.
%   gradient_amplitude_per_bvec      - (Mx1 vector) Gradient amplitude per direction.
%   gradient_duration_per_bvec       - (Mx1 vector) Gradient duration per direction.
%   gradient_separation_per_bvec     - (Mx1 vector) Gradient separation per direction.
%   q_per_bvec                       - (Mx1 vector) Computed q-values per direction.
%   bval_per_bvec                    - (Mx1 vector) b-values per direction.
%   bvec                             - (Mx3 matrix) Gradient directions for all acquisitions.
%   num_shells                       - (scalar) Number of unique shells.
%   num_gradient_directions_per_shell - (Nx1 vector) Number of directions per shell.
%   bvec_idx_per_shell               - (Nx1 cell) Indices of bvecs for each shell.
%   t_e                              - (scalar) Echo time (TE).
%   snr                              - (scalar) Signal-to-noise ratio (if available).
%
%  USAGE:
%   % Create a protocol manually
%   protocol = MrProtocol(gradient_amplitudes, durations, separations, ...
%                         bvals, gradient_dirs, t_e);
%
%   % Load protocol from json file
%   protocol = MrProtocol.fromJsonFile(json_file)
%       - Static method to load protocol parameters from a JSON file.

    
    properties(SetAccess=public, GetAccess=public)
        gradient_amplitude_per_shell
        gradient_duration_per_shell
        gradient_separation_per_shell
        q_per_shell
        bval_per_shell
        bvec_per_shell
        gradient_amplitude_per_bvec
        gradient_duration_per_bvec
        gradient_separation_per_bvec
        q_per_bvec
        bval_per_bvec
        bvec
        num_shells
        num_gradient_directions_per_shell
        bvec_idx_per_shell
        t_e
        snr
    end
    
    methods
        
        function obj = MrProtocol(gradient_amplitude_per_shell, gradient_duration_per_shell, gradient_separation_per_shell, bval_per_shell, gradient_directions_per_shell, t_e)
            %MRPROTOCOL Construct an instance of this class
            %   Detailed explanation goes here
            obj.gradient_amplitude_per_shell = gradient_amplitude_per_shell;
            obj.gradient_duration_per_shell = gradient_duration_per_shell;
            obj.gradient_separation_per_shell = gradient_separation_per_shell;
            obj.bval_per_shell = bval_per_shell;      
            obj.num_shells = length(obj.bval_per_shell);
            obj.t_e = t_e;
            if ~iscell(gradient_directions_per_shell)
                cell_gradient_directions_per_shell = cell(1, size(gradient_directions_per_shell,1));
                for i = 1:size(gradient_directions_per_shell, 1)
                    cell_gradient_directions_per_shell{i} = squeeze(gradient_directions_per_shell(i, :, :));  % Extract the i-th slice along the first dimension
                end
                gradient_directions_per_shell = cell_gradient_directions_per_shell;
            end
            obj.num_gradient_directions_per_shell = cellfun(@(x) size(x, 1), gradient_directions_per_shell);
         

            bvec_per_bval = cell(1, length(bval_per_shell));
            bvec = [];
            bval_per_bvec = [];
            for gradient_dir_index = 1:length(obj.num_gradient_directions_per_shell)
                bvec_per_bval{gradient_dir_index} = gradient_directions_per_shell{gradient_dir_index}';
                bvec = [bvec bvec_per_bval{gradient_dir_index}];
                bval_per_bvec = [bval_per_bvec repmat(bval_per_shell(gradient_dir_index), 1, obj.num_gradient_directions_per_shell(gradient_dir_index))];
            end   
            bvec = bvec';
            bval_per_bvec = bval_per_bvec';

            obj.gradient_amplitude_per_bvec = zeros(length(bvec),1);
            obj.gradient_duration_per_bvec = zeros(length(bvec),1);
            obj.gradient_separation_per_bvec = zeros(length(bvec),1);
            obj.bvec_per_shell = cell(obj.num_shells, 1);
            obj.bvec_idx_per_shell = cell(obj.num_shells, 1);
            for shell_index = 1:obj.num_shells
                shell_mask_idx = find(bval_per_bvec == obj.bval_per_shell(shell_index));
                obj.gradient_amplitude_per_bvec(shell_mask_idx) = obj.gradient_amplitude_per_shell(shell_index);
                obj.gradient_duration_per_bvec(shell_mask_idx) = obj.gradient_duration_per_shell(shell_index);
                obj.gradient_separation_per_bvec(shell_mask_idx) = obj.gradient_separation_per_shell(shell_index);
                obj.bvec_per_shell{shell_index} = bvec(shell_mask_idx,:);
                obj.bvec_idx_per_shell{shell_index} = shell_mask_idx;
            end

            obj.bval_per_bvec = bval_per_bvec;
            obj.bvec = bvec;

            gyroMagnRatio =  267.513*10^(-6);
            obj.q_per_shell = obj.gradient_amplitude_per_shell*gyroMagnRatio;
            obj.q_per_bvec = obj.gradient_amplitude_per_bvec*gyroMagnRatio;
        end
    end

    methods (Static)
        function obj = fromJsonFile(json_file)

            json_data = jsondecode(fileread(json_file));
            obj = MrProtocol( ...
                json_data.diffusion_gradient_amplitude, ...
                json_data.diffusion_gradient_separation, ...
                json_data.diffusion_gradient_duration, ...
                json_data.b_value, ...
                json_data.diffusion_gradient_direction, ...
                json_data.echo_time);
            obj.snr = json_data.signal_to_noise_ratio;
        end
    end


end

