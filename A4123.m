%% q 1 to 3

clearvars; clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle','docked')  % 'docked' 'normal'
set(0,'DefaultLineLineWidth',1)


% dV = 10;
% dI = 10;
% V = I * R;
% I = C* dV;
% V2 = L*dI;


R1 = 1;
R2 = 2;
R3 = 10 ; %%changed due to slope of current vs voltage sweep
R4 = 0.1;
Ro = 1000;
C1 = 0.25;
L = 0.2;
alpha = 100;

G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
Go = 1/Ro;


%equations
% V1 = Vin 
% I1 = (V1 - V2)*G1; +  C*d(V1-V2)/dt;
% (V1 - V2)*G1; +  C*d(V1-V2)/dt - V2*G2 - IL = 0 ;
% IL - V3*G3 = 0
% V2-V3 = L*d(IL)/dt;
% V4 = alpha*IL;
% I4 = (V4-Vo)*G4;
% 0 =  (V4-Vo)*G4 + (Vo*Go);

%unknowns:  V = [V1, I1, V2, IL, V3, V4, I4, Vo]

G = [1   0   0     0      0  0  0  0;
     G1 -1  -G1    0      0  0  0  0;
    -G1  0  G1+G2  0      0  0  0  0;
      0  0    0    1     G3 0  0  0;
      0  0   -1    0      1  0  0  0;
      0  0    0   -alpha  0  -1 0  0;
      0  0    0     0     0  G4 -1 -G4;
      0  0    0     0     0  -G4  0 Go+G4]


F = [1;0;0;0;0;0;0;0];

C = [ 0 0  0 0 0 0 0 0;
      C1 0 -C1 0 0 0 0 0;
     -C1 0  C1 0 0 0 0 0;
      0 0  0 0 0 0 0 0;
      0 0  0 L 0 0 0 0
      0 0  0 0 0 0 0 0;
      0 0  0 0 0 0 0 0;
      0 0  0 0 0 0 0 0;]

Vin = [];
V3 = [];
Vo = [];
 %% time domaain stuff 
for i = 1:21
    Vin(i) = i - 11;
    F = [Vin(i);0;0;0;0;0;0;0];

    %A*x = B
    %G*V = F

    V = G\F;

    V3(i) = V(5);
    Vo(i) = V(8);
    
end
figure(1)
plot (Vin, V3);
hold on;
plot (Vin, Vo);
title('DC Sweep Plot')
xlabel('Vin')
ylabel('Voltage')
legend('V3', 'Vo')



% %conductance = G = 1/R
% V1 = Vin 
Vin = 10; %%sweep 
% I1 = (V1 - V2)*G1; +  C*d(V1-V2)/dt;
% 
% %node 2
% 
% % I(n1->n2) + I(n0->n2) + I(n3->n2) = 0
% (V1 - V2)*G1; +  C*d(V1-V2)/dt - V2*G2 -IL = 0 ;
% 
% %node 3
% IL - V3*G3 = 0
% V2-V3 = L*d(IL)/dt;
% % IL = I3;
% %node 4
% 
% V4 = alpha*IL;
% I4 = (V4-Vo)*G4;
% 
% %node 5
% 
% 0 =  (V4-Vo)*G4 + (Vo*Go);

%% freqency response 
F = [1;0;0;0;0;0;0;0];

omegaV = [];
for omega = 1:100

    omegaV(omega) = omega;

H = G + 1j*omega*C;

V = H\F;
Vo(omega) = V(8);

end

% plot (omegaV, abs(Vo));
% title('Frequency Response')
% xlabel('Omega')
% ylabel('Voltage')


%dB plot
figure (2);
plot (omegaV, 20*log10(abs(Vo)));
title('Frequency Response')
xlabel('Omega')
ylabel(' Gain (dB)')
legend('omega', 'abs(Vo)')



% figure (3)
% G = tf((Vo/Vin),1);
% bode(G);


std = .05;
omega = pi;

C1_vary = std*randn(1000,1) + C1;
figure(100)
hist(C1_vary,50)


for i = 1:length(C1_vary)
    C(2,1) = C1_vary(i);
    C(2,3) = -C1_vary(i);
    C(3,1) = -C1_vary(i);
    C(3,3) = C1_vary(i);

    %same as line 154 to 157
%     C = [ 0 0  0 0 0 0 0 0;
%         C1_vary(i) 0 -C1_vary(i) 0 0 0 0 0;
%         -C1_vary(i) 0  C1_vary(i) 0 0 0 0 0;
%         0 0  0 0 0 0 0 0;
%         0 0  0 L 0 0 0 0
%         0 0  0 0 0 0 0 0;
%         0 0  0 0 0 0 0 0;
%         0 0  0 0 0 0 0 0;]

    H = G + 1j*omega*C;

    V = H\F;
    
%     Vin(i) = V(1);
    Vo(i) = abs(V(8));

end

figure(5)
subplot(2,1,1);
hist(C1_vary,15)

subplot(2,1,2);
hist(Vo,50)
