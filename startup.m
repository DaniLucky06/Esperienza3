addpath(genpath(pwd))

% Some useful variables
    g = 9.806;                    % Gravitational acceleration
   
    c     = 299792458;            % Speed of light in vacuum
    G     = 6.6742867e-11;        % Newtonian constant of gravitation
    h     = 6.6260689633e-34;     % Planck constant
    hred  = h/(2*pi);             % Reduced Planck constant
    
    me    = 9.1093821545e-31;     % Electron mass
    mp    = 1.67262163783e-27;    % Proton mass
    mn    = 1.674927351e-27;      % Neutron mass
    amu   = 1.660538921e-27;      % Unified atomic mass unit
    
    u0    = 4*pi*1e-7;            % Magnetic constant
    e0    = 1/(u0*c^2);           % Electric constant
    Z0    = u0*c;                 % Impedance of vacuum
    ke    = 1/(4*pi*e0);          % Coulomb's constant
    e     = 1.60217648740e-19;    % Elementary charge
    uB    = (e*hred)/(2*me);      % Bohr magneton
    G0    = (2*e^2)/h;            % Conductance quantum
    G0inv = h/(2*e^2);            % Inverse conductance quantum
    KJ    = (2*e)/h;              % Josephson constant
    O0    = h/(2*e);              % Magnetic flux quantum
    uN    = (e*hred)/(2*mp);      % Nuclear magneton
    RK    = h / e^2;              % von Klitzing constant

    alpha = (u0*e^2*c)/(2*h);     % fine-structure constant
    Rinf  = (alpha^2*me*c)/(2*h); % Rydberg constant
    a0    = alpha/(4*pi*Rinf);    % Bohr radius
    re    = e^2/(4*pi*e0*me*c^2); % Classical electron radius
    Eh    = 2*Rinf*h*c;           % Hartree energy
    
    kB    = 1.3806504e-23;        % Boltzmann constant
    Na    = 6.02214129e23;        % Avogadro's number
    R     = 8.3144621;            % Gas constant