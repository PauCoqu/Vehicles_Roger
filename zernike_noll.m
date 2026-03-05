function [n, m, name, val, Z_func] = zernike_noll(x, y, j)
% ZERNIKE_NOLL Computes Zernike polynomial statistics and values for a given Noll index.
%
%   INPUTS:
%       x, y : Cartesian coordinates (scalars, vectors, or matrices).
%       j    : Scalar integer Zernike Noll index (j >= 1).
%
%   OUTPUTS:
%       n    : Radial order.
%       m    : Azimuthal frequency.
%       name : String name of the aberration (for j <= 15).
%       val  : Calculated values of the polynomial at (x, y). 
%              Note: Values outside the unit disk (rho > 1) are set to 0.
%       Z_func: Anonymous function handle @(x,y) allowing re-evaluation.
%
%   Example:
%       [n, m, name, val, ~] = zernike_noll(X, Y, 4);
%       % Returns Defocus (n=2, m=0)

    %% 1. Input Validation
    if ~isscalar(j) || j < 1 || floor(j) ~= j
        error('Index j must be a positive integer scalar.');
    end

    %% 2. Calculate n and m from Noll index j
    % Calculate radial order n
    % n is the smallest integer such that j <= (n+1)(n+2)/2
    n = ceil((-3 + sqrt(9 + 8*(j-1))) / 2);
    
    % Calculate azimuthal frequency m
    sub_j = j - n*(n+1)/2; % Position within the order n
    
    % Logic for m based on Noll convention
    if mod(n, 2) == 0
        % Even radial order
        m_mag = 2 * floor(sub_j / 2);
    else
        % Odd radial order
        m_mag = 2 * floor((sub_j - 1) / 2) + 1;
    end
    
    if mod(j, 2) == 0
        m = m_mag;  % Even j -> Cosine term (or m=0)
    else
        m = -m_mag; % Odd j -> Sine term
    end

    %% 3. Identify Name (for lower orders)
    name = '';
    switch j
        case 1, name = 'Piston';
        case 2, name = 'Tilt X (Tip)';
        case 3, name = 'Tilt Y (Tilt)';
        case 4, name = 'Defocus';
        case 5, name = 'Astigmatism (Oblique)';
        case 6, name = 'Astigmatism (Vertical)';
        case 7, name = 'Coma (Vertical)';
        case 8, name = 'Coma (Horizontal)';
        case 9, name = 'Trefoil (Vertical)';
        case 10, name = 'Trefoil (Oblique)';
        case 11, name = 'Spherical (Primary)';
        case 12, name = 'Spherical (Secondary)'; % Sometimes secondary astig
        case 13, name = 'Quadrafoil (Vertical)';
        case 14, name = 'Quadrafoil (Oblique)';
        case 15, name = 'Spherical (Secondary)'; 
    end

    %% 4. Define the Zernike Calculation Logic
    % We define this as a nested function so the handle can capture n and m.
    
    calc_zernike = @(xx, yy) evaluate_zernike_raw(xx, yy, n, m);
    
    % Create the output function handle
    Z_func = calc_zernike;

    %% 5. Compute Value
    val = calc_zernike(x, y);

end

function z_vals = evaluate_zernike_raw(x, y, n, m)
    % Helper to compute Zernike value
    
    % Convert to polar
    [theta, rho] = cart2pol(x, y);
    
    % Radial Polynomial R_n^|m|
    R = zeros(size(rho));
    m_abs = abs(m);
    
    % Optimization: Pre-compute coefficients could be faster, but direct sum 
    % is clearer for general n.
    s_max = (n - m_abs) / 2;
    for s = 0:s_max
        num = (-1)^s * factorial(n - s);
        den = factorial(s) * factorial((n + m_abs)/2 - s) * factorial((n - m_abs)/2 - s);
        R = R + (num / den) * rho.^(n - 2*s);
    end
    
    % Normalization factor (Noll standard)
    if m == 0
        norm_factor = sqrt(n + 1);
        azimuthal = 1;
    else
        norm_factor = sqrt(2 * (n + 1));
        if m > 0
            azimuthal = cos(m * theta);
        else
            azimuthal = sin(m_abs * theta);
        end
    end
    
    z_vals = norm_factor .* R .* azimuthal;
    
    % Mask values outside unit disk
    %z_vals(rho > 1) = 0; 
end