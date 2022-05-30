clc
clear


myred           = [216 30 49]/255;
myblue          = [27 99 157]/255;
myblack         = [0 0 0]/255;
mygreen         = [0 128 0]/255;
mycyan          = [2 169 226]/255;
myyellow        = [251 194 13]/255;
mygray          = [89 89 89]/255;

set(groot,'defaultAxesColorOrder',[myblack;myblue;myred;mygreen;myyellow;mycyan;mygray]);

alw             = 1;                        % AxesLineWidth
fsz             = 26;                       % Fontsize
lw              = 2;                        % LineWidth
msz             = 40;                       % MarkerSize

set(0,'defaultLineLineWidth',lw);           % set the default line width to lw
set(0,'defaultLineMarkerSize',msz);         % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);           % set the default line width to lw
set(0,'defaultLineMarkerSize',msz);         % set the default line marker size to msz


fn = 'Cambria Math';    

% Initialise the state space

A = [-0.04891 0.09804 0 ; -1.426*10^(-4) -2.85*10^(-4) 0 ; 0 126.7 0 ];
B = [0.03625; 3.8*10^(-5); 0];

% A = [-0.313 56.7 0; -0.0139 -0.426 0; 0 56.7 0];
% B = [0.232; 0.0203; 0];

C = eye(3);

%Reachability condition
W = [B A*B A^2*B];
Wrank = rank(W);

%Calculate the characteristic equation
syms s
char = det(s * eye(Wrank) - A);
coeffs = eval(coeffs(char));
coeffs = fliplr(coeffs(1:end));
coeffs(end+1) = 0;


%Set the system parameters
Mp = 0.05; %Overshoot percentage
ts = 4; %Settling time
tr = 3; %Rise time

sigma = 4/ts;
wn = 1.8/tr;
zeta = 0.707; %From 5% overshoot

%Calculate the position of the non dominant poles
pole3 = 10*sigma;
pole4 = pole3;

%Desired characteristic polynomial
% p = (s+pole3) * (s+pole4) * (s^2+2*zeta*wn*s+wn^2);
% p_coeffs = eval(coeffs(p));

p_coeffs = roots([1 2*zeta*wn wn^2]);
p_coeffs(end+1) = pole3;

%Find the desired characteristic polynomial
syms s
cp = vpa(expand((s-p_coeffs(1))*(s-p_coeffs(2))*(s-p_coeffs(3))));
cp_coeffs = sym2poly(cp);

%Calculate the gain array
Kt = cp_coeffs-coeffs;

%Arrange the gain array
Kt = Kt(2:4);

%Calculate the new reachability matrix
Bt = [1; 0; 0];
At = [-coeffs(2) -coeffs(3) -coeffs(4); 1 0 0 ; 0 1 0];
Wt = [Bt At*Bt At^2*Bt];

%Find the transformation matrix
Transf = Wt * inv(W);
Kgain = Kt * Transf;

%Feedforward gain design
kr = -1 / (C*inv(A-B*Kgain)*B);
kr = kr(3);

%Observers
%Determine the rank and observability
C = [0 0 1];
Wobs = [C; C*A; C*A^2];
Wobsrank = rank(Wobs);

%Find the characteristic polynomial
syms s
L = sym('m',[3 1]);
CPobs = det(s*eye(3) - A+ L*C);

%Use simultaneous equations for coefficient matching
syms x y z
eqn1 = z+0.0492 == 4;
eqn2 = 126.7*y + 0.0492*z + 2.792*10^(-5) == 12;
eqn3 = -0.01807*x + 6.197*y + 2.792*10^(-5) == 8;

[A1, B1] = equationsToMatrix([eqn1, eqn2, eqn3], [x,y,z]);
obsgain = linsolve(A1,B1);

Lgain = [-410.; 0.0932; 3.9508];


%Simulate the closed loop system with observers
T = 10;
r = 5;
x0 = [0.0363;3.8e-5;0];
xh0 = zeros(3,1);
[t,x] = ode45(@(t,xt)closedloop(t,xt,A,B,C,Kgain,kr,Lgain,r),[0,T],[x0; xh0]);
x = x'; y = C*x(1:3,:);

figure(1)
grid minor
hold on
plot([0,T],[r r],'k:')
plot(t,y)
xlabel('Time (s)')
ylabel('Pitch Angle (Degrees)')

set(gcf, 'Color', [1 1 1]);      
set(gca, 'Color', [1 1 1]);  
set(gca,'GridLineStyle','-')                            % set gridline and font properties
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(gca,'FontSize',13)

%Simulate the closed loop system without observers
[t,x] = ode45(@(t,x)dynam(t,x,A,B,Kgain,kr,r),[0 T],x0);
x = x'; y = C*x;
