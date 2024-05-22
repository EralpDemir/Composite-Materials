% This script defines the fiber-matrix properties

%_______________________ Fiber and matrix properties _____________

% Fiber properties (S-Glass)
% Ref [Ever J. Barbero, page 29, table 2.1]
%vol_f = 0.45;
%E_f = 85e3; % MPa
%nu_f = 0.22;
%G_f = E_f/2/(1+nu_f); % MPa

% Matrix properties (Epoxy 9310/9360 @23degC)
% Ref [Ever J. Barbero, page 29, table 2.4]
%vol_m = 1 - vol_f;
%E_m = 3.12e3;
%nu_m = 0.38;
%G_m = E_m/2/(1+nu_m);

% Fiber-Matrix homogenization
%E11 = vol_f*E_f + vol_m*E_m;
%E22 = 1/(vol_f/E_f + vol_m/E_m);
%nu12 = vol_f*nu_f + vol_m*nu_m;
%G12 = 1/(vol_f/G_f + vol_m/G_m);


% mass density
rho = 1950;


% Longitudinal modulus
% E11 = 181e3;
E11 = 25e9;
% Lateral modulus
% E22 = 10.273e3;
E22 = 1e9;
% Shear modulus
% G12 = 7.1705e3;
G12 = 0.5e9;
% G13=2.6E3;
% G23=2.6E3;
% Poissons ratio
% nu12 = 0.28;
nu12 = 0.25;



% Calculation of stiffness parameters
Q11 = E11^2/(E11-(nu12^2*E22));
Q12 = nu12*E11*E22/(E11-(nu12^2*E22));
Q22 = E11*E22/(E11-(nu12^2*E22));
Q66 = G12;

%         % Shear stiffness parameters
%        Q44 = G23;
%        Q55 = G13;

%____________________ Calculation of lamination parameters _______
% Lamina invariants in the principal material coordinate system
U1 = (3*Q11+3*Q22+2*Q12+4*Q66)/8;
U2 = (Q11-Q22)/2;
U3 = (Q11+Q22-2*Q12-4*Q66)/8;
U4 = (Q11+Q22+6*Q12-4*Q66)/8;
% U5 = Q44+Q55;
% U6 = Q44-Q55;

U = [    U1   U2     0    U3    0;
    U4    0     0   -U3    0;
    U1  -U2     0    U3    0;
    0    0  U2/2     0   U3;
    0    0  U2/2     0  -U3;
    (U1-U4)/2    0     0   -U3    0 ];

% U_ = [   U5   U6     0     0    0;
%     U5  -U6     0     0    0;
%     0    0   -U6     0    0 ];




% In Abdalla et al. "b" refers to shorter length, but it is not!
% Original expression in the reference
% nondim =  (Ly^2 / t) * sqrt (rho/E22) ;
% Use the following expression instead!!! Note Lx is used instead of Ly.
 nondim =  (Lx^2 / t) * sqrt (rho/E22) ;

 
% % The other normalization
% nu21=nu12*E22/E11;
% D0 = E22*t^3/12/(1-nu12*nu21);
% nondim = Lx^2*sqrt(rho/D0);







