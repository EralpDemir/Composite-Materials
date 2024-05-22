% This script defines the fiber-matrix properties

%_______________________ Fiber and matrix properties _____________

% Fiber properties (S-Glass)
% Ref [Ever J. Barbero, page 29, table 2.1]
f_f = 0.45;
E_f = 85e3; % MPa
nu_f = 0.22;
G_f = E_f/2/(1+nu_f); % MPa

% Matrix properties (Epoxy 9310/9360 @23degC)
% Ref [Ever J. Barbero, page 29, table 2.4]
f_m = 1 - f_f;
E_m = 3.12e3;
nu_m = 0.38;
G_m = E_m/2/(1+nu_m);

% Fiber-Matrix homogenization
E11 = f_f*E_f + f_m*E_m;
E22 = 1/(f_f/E_f + f_m/E_m);
nu12 = f_f*nu_f + f_m*nu_m;
G12 = 1/(f_f/G_f + f_m/G_m);

%         % Out of plane Shear Moduli - counted as ZERO !!!
%         G13 = 5/6*G12;
%         G23 = 5/6*G12;

% Calculation of stiffness parameters
Q11 = E11^2/(E11-(nu12^2*E22));
Q12 = nu12*E11*E22/(E11-(nu12^2*E22));
Q22 = E11*E22/(E11-(nu12^2*E22));
Q66 = G12;

%         % Shear stiffness parameters
%         Q44 = G23;
%         Q55 = G13;

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
