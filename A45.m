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
Cn = 0.00001; 
L = 0.2;
alpha = 100;
In = 0.001;

G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
Go = 1/Ro;


%equations
% V1 = Vin 
% I1 = (V1 - V2)*G1; +  C*d(V1-V2)/dt;
% (V1 - V2)*G1; +  C*d(V1-V2)/dt - V2*G2 - IL = 0 ;
% IL - V3*G3 + Cn*d(V3)/dt + In = 0
% V2-V3 = L*d(IL)/dt;
% V4 = alpha*IL;
% I4 = (V4-Vo)*G4;
% 0 =  (V4-Vo)*G4 + (Vo*Go);

%unknowns:  V = [V1, I1, V2, IL, V3, V4, I4,Vo]

G = [1   0   0     0      0  0  0  0;
     G1 -1  -G1    0      0  0  0  0;
    -G1  0  G1+G2  0      0  0  0  0;
      0  0    0    1     G3 0  0  0;
      0  0   -1    0      1  0  0  0;
      0  0    0   -alpha  0  -1 0  0;
      0  0    0     0     0  G4 -1 -G4;
      0  0    0     0     0  -G4  0 Go+G4];  


F = [1;0;0;0;0;0;0;0];

%%add CN
C = [ 0 0  0 0 0 0 0 0;
      C1 0 -C1 0 0 0 0 0;
     -C1 0  C1 0 0 0 0 0;
      0 0  0 0 Cn 0 0 0;
      0 0  0 -L 0 0 0 0
      0 0  0 0 0 0 0 0;
      0 0  0 0 0 0 0 0;
      0 0  0 0 0 0 0 0;]


%%setting up
dt = 0.001;
numit= 1000;
f = 1/(0.03);
gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;

x = linspace(0,1,1000);
mu = 0.2; %%middle
sig = 0.03;  %standard dev (in mili secs)
amp = 1; %amplitude
vo = 0;
plot(x,gaus(x,mu,sig,amp,vo))


Vin = [];
V3 = [];
Vo = [];
Vp = zeros(8,1);


for i = 1:numit
    t = (i-1)*dt;
    ts(i) = (i-1)*dt;

    Vin(i) = amp*exp(-(((t-mu)^2)/(2*sig^2)))+vo;
    F = [Vin(i);0;0;0;In*randn();0;0;0];  %voltage source and random curr source

  
    Vi(i) = Vin(i) ;

    F(1) = Vin(i);

    H = C/dt+ G;

    V = H\(F+(C/dt*Vp)) ;

    Vo(i) = V(8);
    Vp = V;
    
end

% %conductance = G = 1/R
% V1 = Vin 
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


figure(1)
leg = {};
plot (ts, Vo); leg{end+1}='Vo'; hold on
plot (ts, Vi); leg{end+1}='Vin'; hold on
hold off
title('Gaussian Pulse Plot With Noise')
xlabel('Time (s)')
ylabel('Voltage (V)' )
legend(leg)



%%fft stuff

VoFFT = abs(fftshift(fft(Vo)));
ViFFT = abs(fftshift(fft(Vi)));
freqV = (0:1000-1) -(1/(2*dt));


figure (2)
hold on
plot (ts,  20*log10(VoFFT));
plot (ts,  20*log10(ViFFT));
title('Gaussian Pulse Frequency Response')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
legend('Output', 'Input')

