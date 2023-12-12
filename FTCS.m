L = 0.4;                       % length of domain in x direction [m]

I0 = 0.0000829*24*60*60;        % Bq m^-2 day^-1
R1 = (0.693/53.2);              % Decay constant day^-1
tmax = 35;                      % end time [day]

nx = 41;                        % number of nodes in x direction
nt = 50401;                     % number of time steps

dx = L/(nx-1);
dt = tmax/(nt-1);

alpha= 3.5447*10^-9*24*3600;   % Mass advection [m^2/day]
r = alpha*dt/dx^2; rr = 1 - 2*r-R1*dt;

T = 25+273;                            % Absolute temperature [K]
mu_oil = 51*24*3600;              % Oil Viscosity [kg/(m.day)]
mu_water = 0.001*24*3600;         % Water Viscosity [kg/(m.day)]

dEdx = 150;                       % Voltage gradient [V/m]
V = dEdx*L;                       % Voltage [V]
zeta = -0.0027;                   % Zeta potential [V]
u_e = 0;                          % Electromigration mobility [m^2/V.s] 
e_0 = 8.854*10^-12;               % permittivity of free space [F/m]
e_r = 7.5;                        % relative permittivity of clay [F/m]
e_oil = 2.3;                      % Crude oil permittivity [F/m]
e_clay = e_0*e_r;                 % Clay permittivity [F/m]
e_w = 5;                          % Water permittivity [F/m]
n=0.64;                           % Porosity
tau = 0.44;                       % Tortuosity
z = 1;                            % Valency
rho = 1620;                       % Clay Density [kg/m^3]

% Calculated constants
k_eo1 = (e_oil*zeta)*n/(mu_oil);                 % Electroosmotic mobility Alshawabkeh [m^2/V.s]
k_eo2 = (e_oil*zeta)*n/(mu_oil*(tau^2));         % Electroosmotic mobility Vane [m^2/V.s]
k_eo3 = (e_oil*zeta*dEdx)/(mu_oil);              % Electroosmotic mobility Shapiro [S/m]

v0 = 2.5*10^-6; % Convection velocity m day^-1

v = u_e + k_eo1;
beta = v*dEdx*10;
rr1 = beta/(2*dx); rrr = (R1*dt);


J0 = I0/sqrt(R1*alpha); % total inventory of Be-7 in soil
J0 = J0/100;

% --- Create arrays to save data for export
x = linspace(0,L,nx);
t = linspace(0,tmax,nt);
U = zeros(nx,nt);
Q = zeros(nx,nt);
E = zeros(nx,nt);
Y = zeros(nx,nt);
G = zeros(nx,nt);

W = zeros(nx,nt);
X = zeros(nx,nt);
J = zeros(nx,nt);
O = zeros(nx,nt);
B = zeros(nx,nt);
% --- Set IC and BC
U(:,1)= 2000;
U(:,2)= 2000;
Q(:,1)= 4000;
Q(:,2)= 4000;
E(:,1)= 6000;
E(:,2)= 6000;
Y(:,1)= 8000;
Y(:,2)= 8000;
G(:,1)= 10000;
G(:,2)= 10000;

W(:,1)= 2000;
X(:,1)= 4000;
J(:,1)= 6000;
O(:,1)= 8000;
B(:,1)= 10000;

% --- Loop over time steps

for m= 2:nt-1
    U(1,m) =J0; %--- Upper boundary
    U(end,m) = 0; %--- Lower boundary
    Q(1,m) =J0; %--- Upper boundary
    Q(end,m) = 0; %--- Lower boundary
    E(1,m) =J0; %--- Upper boundary
    E(end,m) = 0; %--- Lower boundary
    Y(1,m) =J0; %--- Upper boundary
    Y(end,m) = 0; %--- Lower boundary
    G(1,m) =J0; %--- Upper boundary
    G(end,m) = 0; %--- Lower boundary

    for i= 2:nx-1
        
        
        U(i,m+1) = U(i,m) + r*(U(i+1,m) -2*U(i,m) + U(i-1,m)) + rr1*(U(i+1,m) - U(i-1,m)) + rrr; 
        Q(i,m+1) = Q(i,m) + r*(Q(i+1,m) -2*Q(i,m) + Q(i-1,m)) + rr1*(Q(i+1,m) - Q(i-1,m)) + rrr;
        E(i,m+1) = E(i,m) + r*(E(i+1,m) -2*E(i,m) + E(i-1,m)) + rr1*(E(i+1,m) - E(i-1,m)) + rrr;
        Y(i,m+1) = Y(i,m) + r*(Y(i+1,m) -2*Y(i,m) + Y(i-1,m)) + rr1*(Y(i+1,m) - Y(i-1,m)) + rrr;
        G(i,m+1) = G(i,m) + r*(G(i+1,m) -2*G(i,m) + G(i-1,m)) + rr1*(G(i+1,m) - G(i-1,m)) + rrr;

        %U(i,m) = r*U(i-1,m-1)+ rr*U(i,m-1)+ r*U(i+1,m-1);


        
        
    end
end

for m= 2:nt
    W(1,m) =J0; %--- Upper boundary
    W(end,m) = 0; %--- Lower boundary
    X(1,m) =J0; %--- Upper boundary
    X(end,m) = 0; %--- Lower boundary
    J(1,m) =J0; %--- Upper boundary
    J(end,m) = 0; %--- Lower boundary
    O(1,m) =J0; %--- Upper boundary
    O(end,m) = 0; %--- Lower boundary
    B(1,m) =J0; %--- Upper boundary
    B(end,m) = 0; %--- Lower boundary

    for i= 2:nx-1
        
         
        W(i,m) = r*W(i-1,m-1)+ rr*W(i,m-1)+ r*W(i+1,m-1);
        X(i,m) = r*X(i-1,m-1)+ rr*X(i,m-1)+ r*X(i+1,m-1);
        J(i,m) = r*J(i-1,m-1)+ rr*J(i,m-1)+ r*J(i+1,m-1);
        O(i,m) = r*O(i-1,m-1)+ rr*O(i,m-1)+ r*O(i+1,m-1);
        B(i,m) = r*B(i-1,m-1)+ rr*B(i,m-1)+ r*B(i+1,m-1);

        
    end
end

figure;
plot(t,W(10,:),'--','DisplayName', 'Original FDM 2000');
hold on;
plot(t,X(10,:),'--','DisplayName', 'Original FDM 4000');

plot(t,J(10,:),'--','DisplayName', 'Original FDM 6000');

plot(t,O(10,:),'--','DisplayName', 'Original FDM 8000');

plot(t,B(10,:),'--','DisplayName', 'Original FDM 10000');



plot(t,U(10,:),'-','DisplayName', 'Optimized FDM 2000');

plot(t,Q(10,:),'-','DisplayName', 'Optimized FDM 4000');

plot(t,E(10,:),'-','DisplayName', 'Optimized FDM 6000');

plot(t,Y(10,:),'-','DisplayName', 'Optimized FDM 8000');

plot(t,G(10,:),'-','DisplayName', 'Optimized FDM 10000');

xlabel('Time');
ylabel('Conc(mg/kg)');

legend();
