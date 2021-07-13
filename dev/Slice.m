classdef Slice < handle
    properties
        id              % Slice id
        pixel_positions % Pixel positions of the interference image        
        pixel_values    % Interference image pixel values
        source_image    % Handle to source image
        radius          % Radius of the slice
        phase_shift_vals % Measured phase shift values
        phase_shift_pos % Measured phase shift positions
        wavelength      % Wavelength of reference beam
        n0              % Oustide medium refractive index
        reference_ray_half_shift % halfshift (0 or lambda/2) for the zero value 
        n               % aproximated refractive index for each shell
        nShells         % Number of shells
        shellR          % Radius for each shell
    end
    
    methods
        function this = Slice(id,pixel_values,pixel_positions,radius,n0,wavelength,reference_ray_half_shift)
            % Input variables
            this.id = id;
            this.pixel_positions = pixel_positions;
            this.pixel_values = pixel_values;
            this.radius = radius;
            this.n0 = n0;
            this.wavelength = wavelength;
            this.reference_ray_half_shift = reference_ray_half_shift;
        end
        
        % simple unweighted function for extracting line information from
        % image pixels
        function get_peak_positions(this)
            v = 
            [peakVal,peakPos,nrPeaks] = findPeaks(coords,vals,threshold);
        end
        
    end
end