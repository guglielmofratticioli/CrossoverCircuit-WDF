% Fratticioli Pirrello
clear all
close all
clc

%% Import Input Audio Signal
[Vin,~] = audioread('ExpSweep.wav');

%% LTSpice Files for Ground-Truth
[OutLowSpice,~]=audioread('outlowsweep.wav');
[OutMidSpice,~]=audioread('outmidsweep.wav');
[OutHighSpice,FsLTSpice]=audioread('outhighsweep.wav');
TsLTSpice=1/FsLTSpice;

%% Sampling frequency (to be varied: FsLTSpice/downSampFact, with downSampFact={4,3,2})
downSampFact=4;
Fs =FsLTSpice/downSampFact; 

%% Downsample Input Signal
Vin=Vin([1:downSampFact:end]);

%% Sampling Period
Ts=1/Fs;
%% Number of Samples
Nsamp=length(Vin);
%% Simulated time
tstop=Nsamp*Ts;
%% Parameters of Dynamic Element
L1=0.35*10^(-3);
L2=0.35*10^(-3);
L3=3.5*10^(-3);
L4=3.5*10^(-3);
C1= 2.8*10^(-6);
C2= 2.8*10^(-6);
C3= 28*10^(-6);
C4= 4.7*10^(-6);
C5=28*10^(-6);
C6=47*10^(-6);
%% Resistive Parameters
R1=10;
RspkLow=8;
R2=10;
RspkMid=8;
RspkHigh=8;

%% Variables Initialization


a_h = zeros(6,1); 
b_h = zeros(6,1);
Z_h = zeros(6,1);

%  HIGH BAND CIRCUIT                                                      
%                                                        
%              IN                                        
%           a1    b1                                     
%          +--------+       +--------+                   
%       a2 |      * | a3=b4 |        |  a5                 
%   C1     |   S1   |       |   P2   |     L2              
%      *b2 |       *| b3=a4 |        | *b5                 
%          +--------+       +--------+                   
%                            a6    *b6                    
%                               Rh                       
%                                 -->OUT                         
%                                                        
                                                        
            
a_m = zeros(18,1); 
b_m = zeros(18,1);
Z_m = zeros(18,1);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
%  MID BAND CIRCUIT
%
%              IN                                 C3                                                    R1             
%           b1     a1                          a8    *b8                                             a17   *b17        
%          +---------+       +---------+      +---------+       +---------+       +---------+       +---------+        
%      a2  |         | a3=b4 |         |a6=b7 |         |a9=b10 |         |a12=b13|         |a15=b16|         | a18    
%  L2      |   S1    |       |   P2    |      |   S3    |       |   P4    |       |   P5    |       |   S6    |      C4 
%     *b2  |        *| b3=a4 |        *|b6=a7 |        *|b9=a10 |        *|b12=a13|        *|b15=a16|         | *b18    
%          +---------+       +---------+      +---------+       +---------+       +---------+       +---------+        
%                             a5    *b5                          a11   *b11        a14  *b14                           
%                                C2                                 L3                Rm                               
%                                                                                       -->OUT                          
                                                                                                                      
                                                                                                                                                                                                                                           
a_l = zeros(12,1); 
b_l = zeros(12,1);
Z_l = zeros(12,1);
                                                                                                                                                                                                                        
%  HIGH BAND CIRCUIT
%                                                                                                                                                                                                                                                                                                                                                                                                             
%              IN                                 RL                                                   
%           a1     b1                          a8    *b8                                               
%          +---------+       +---------+      +---------+       +---------+                            
%      a2  |       * | a3=b4 |         |a6=b7 |         |a9=b10 |         |  a12                       
%  L4      |   S1    |       |   P2    |      |   P3    |       |   S4    |       R2                     
%     *b2  |        *| b3=a4 |        *|b6=a7 |        *|b9=a10 |         | *b12                       
%          +---------+       +---------+      +---------+*      +---------+                            
%                            a5    *b5                          a11   *b11                            
%                                C5                                 C6                                 
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            


%% WDF setting of free parameters (adaptation conditions)

    % High  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .     
    
    Z_h(6) = RspkHigh;
    Z_h(5) = 2*L1/Ts; 
    Z_h(4) = (Z_h(6)^-1 + Z_h(5)^-1)^-1; 
    Z_h(3) = Z_h(4); 
    Z_h(2) = Ts/(2*C1);
    Z_h(1) = Z_h(2) + Z_h(3);
    
    % Mid   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   . 

    Z_m(18) = Ts/(2*C4); 
    Z_m(17) = R1;
    Z_m(16) = Z_m(17) + Z_m(18); 
    Z_m(15) = Z_m(16);
    Z_m(14) = RspkMid; 
    Z_m(13) = (Z_m(14)^-1 + Z_m(15)^-1)^-1;
    Z_m(12) = Z_m(13);
    Z_m(11) = 2*L3/Ts;
    Z_m(10) = (Z_m(11)^-1 + Z_m(12)^-1)^-1;
    Z_m(9) = Z_m(10);
    Z_m(8) = Ts/(2*C3);
    Z_m(7) = Z_m(8) + Z_m(9); 
    Z_m(6) = Z_m(7); 
    Z_m(5) = Ts/(2*C2); 
    Z_m(4) = (Z_m(5)^-1 + Z_m(6)^-1)^-1;
    Z_m(3) = Z_m(4); 
    Z_m(2) = 2*L2/Ts; 
    Z_m(1) = Z_m(2)+Z_m(3);

    
    % Low   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .       
    
    Z_l(12) = R2; 
    Z_l(11) = Ts/(2*C6);
    Z_l(10) = Z_l(11)+Z_l(12);
    Z_l(9) = Z_l(10);
    Z_l(8) = RspkLow; 
    Z_l(7) = (Z_l(8)^-1+Z_l(9)^-1)^-1;
    Z_l(6) = Z_l(7);
    Z_l(5) = Ts/(2*C5); 
    Z_l(4) = (Z_l(5)^-1+Z_l(6)^-1)^-1; 
    Z_l(3) = Z_l(4); 
    Z_l(2) = 2*L4/Ts; 
    Z_l(1) = Z_l(2) + Z_l(3); 

%% Computation of Scattering Matrices

    % High  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   . 
    
    gamma_2 = Z_h(5)/(Z_h(5)+Z_h(6)); 
    gamma_1 = Z_h(2)/(Z_h(2)+Z_h(3));
    
    Z_Sh1 = [Z_h(1); Z_h(2); Z_h(3)];  
    G_Ph2 = [Z_h(4)^-1; Z_h(5)^-1; Z_h(6)^-1]; 
    
    I = eye(3);  
    one = ones(3, 1);  
    
    % Calculate the coefficent value
    K_Sh1 = 2 / sum(Z_Sh1);
    K_Ph2 = 2 / sum(G_Ph2);
    
    % Calculate the Scattering matrices
    S_h1 = I - K_Sh1 * Z_Sh1 * one';
    P_h2 = K_Ph2 * G_Ph2 * one' - I;
    P_h2 = P_h2.';

    
    % Mid   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  
    
    Z_Sm1 = [Z_m(1); Z_m(2) ;Z_m(3) ];
    G_Pm2 = [Z_m(4)^-1; Z_m(5)^-1; Z_m(6)^-1]; 
    Z_Sm3 = [Z_m(7); Z_m(8) ;Z_m(9) ]; 
    G_Pm4 = [Z_m(10)^-1; Z_m(11)^-1; Z_m(12)^-1];
    G_Pm5 = [Z_m(13)^-1; Z_m(14)^-1; Z_m(15)^-1];
    Z_Sm6 = [Z_m(16); Z_m(17) ;Z_m(18) ]; 

    I = eye(3);  
    one = ones(3, 1);

    % Calculate the coefficent value
    K_Sm1 = 2 / sum(Z_Sm1);
    K_Pm2 = 2 / sum(G_Pm2);
    K_Sm3 = 2 / sum(Z_Sm3);
    K_Pm4 = 2 / sum(G_Pm4);
    K_Pm5 = 2 / sum(G_Pm5);
    K_Sm6 = 2 / sum(Z_Sm6);

    % Calculate the Scattering matrices
    S_m1 = I - K_Sm1 * Z_Sm1 * one';
    P_m2 = K_Pm2 * G_Pm2 * one' - I;
    P_m2 = P_m2.';
    S_m3 = I - K_Sm3 * Z_Sm3 * one';
    P_m4 = K_Pm4 * G_Pm4 * one' - I;
    P_m4 = P_m4.';
    P_m5 = K_Pm5 * G_Pm5 * one' - I;
    P_m5 = P_m5.';
    S_m6 = I - K_Sm6 * Z_Sm6 * one';
    
    % Low   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .

    Z_Sl1 = [Z_l(1); Z_l(2); Z_l(3); ]; 
    G_Pl2 = [Z_l(4)^-1; Z_l(5)^-1; Z_l(6)^-1; ]; 
    G_Pl3 = [Z_l(7)^-1; Z_l(8)^-1; Z_l(9)^-1; ];
    Z_Sl4 = [Z_l(10); Z_l(11); Z_l(12); ]; 

    % Calculate the coefficent value 
    K_Sl1 = 2 / sum(Z_Sl1);
    K_Pl2 = 2 / sum(G_Pl2);
    K_Pl3 = 2 / sum(G_Pl3);
    K_Sl4 = 2 / sum(Z_Sl4);

    % Calculate the Scattering matrices  
    S_l1 = I - K_Sl1 * Z_Sl1 * one';
    P_l2 = (K_Pl2 * G_Pl2 * one'-I).';
    P_l3 = (K_Pl3 * G_Pl3 * one'-I).';
    S_l4 = I - K_Sl4 * Z_Sl4 * one';

%% Initialize Output Signals
    % Low
    VoutLow=zeros(size(Vin));
    % Mid
    VoutMid=zeros(size(Vin));
    % High
    VoutHigh=zeros(size(Vin));

ii=0;

while (ii<Nsamp)
    ii=ii+1;

    %% Manage Dynamic Elements

        % High  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   . 
        
        a_h(2) = b_h(2);        
        a_h(5) = -b_h(5);       
    
        % Mid   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  
        
        a_m(2) = -b_m(2); 
        a_m(5) = b_m(5);
        a_m(8) = b_m(8);
        a_m(11) = -b_m(11); 
        a_m(18) = b_m(18);

        % Low   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .        
        
        a_l(2) = -b_l(2); 
        a_l(5) = b_l(5); 
        a_l(11) = b_l(11); 

    %% Forward Scan

        % High  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   . 
        
        % P_h2 Junction 
        b_h(4:6) = P_h2 * a_h(4:6); 
        a_h(3) = b_h(4); 
        % S_h1 Junction
        b_h(1:3) = S_h1* a_h(1:3);
    
        % Mid   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  
        
        % S_m6 Junction
        b_m(16:18) = S_m6 * a_m(16:18);
        a_m(15) = b_m(16);
        % P_m5 Junction
        b_m(13:15) = P_m5 * a_m(13:15);
        a_m(12) = b_m(13);
        % P_m4 Junction
        b_m(10:12) = P_m4 * a_m(10:12);      
        a_m(9) = b_m(10);
        % S_m3 Junction
        b_m(7:9) = S_m3 * a_m(7:9);     
        a_m(6) = b_m(7);
        % P_m2 Junction
        b_m(4:6) = P_m2 * a_m(4:6);
        a_m(3) = b_m(4);
        % S_m1 Junction
        b_m(1:3) = S_m1 * a_m(1:3);

        % Low   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .

        % S_l4 Junction
        b_l(10:12) = S_l4 * a_l(10:12); 
        a_l(9) = b_l(10);
        % P_l3 Junction
        b_l(7:9) = P_l3 * a_l(7:9); 
        a_l(6) = b_l(7);
        % P_l2 Junction
        b_l(4:6) = P_l2 * a_l(4:6); 
        a_l(3) = b_l(4);
        % S_l1 Junction
        b_l(1:3) = S_l1 * a_l(1:3); 


    %% Local Root Scattering
    
        % High  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   . 

        a_h(1) = 2 * Vin(ii) - b_h(1);
    
        % Mid   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   . 

        a_m(1) = 2 * Vin(ii) - b_m(1);
        
        % Low   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .

        a_l(1) = 2 * Vin(ii) - b_l(1);

    %% Backward Scan

        % High  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   . 

        % S_h1 Junction
        b_h(1:3) = S_h1* a_h(1:3);
        a_h(4) = b_h(3);
        % P_h2 Junction
        b_h(4:6) = P_h2 * a_h(4:6); 

        % Mid   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  

        % S_m1 Junction
        b_m(1:3) = S_m1 * a_m(1:3); 
        a_m(4) = b_m(3);
        % P_m2 Junction
        b_m(4:6) = P_m2 * a_m(4:6);
        a_m(7) = b_m(6);
        % S_m3 Junction
        b_m(7:9) = S_m3 * a_m(7:9);
        a_m(10) = b_m(9);
        % P_m4 Junction
        b_m(10:12) = P_m4 * a_m(10:12);
        a_m(13) = b_m(12);
        % P_m5 Junction
        b_m(13:15) = P_m5 * a_m(13:15);
        a_m(16) = b_m(15);
        % S_m6 Junction
        b_m(16:18) = S_m6 * a_m(16:18);

        % Low   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .

        % S_l1 Junction
        b_l(1:3) = S_l1 * a_l(1:3);
        a_l(4) = b_l(3);
        % P_l2 Junction
        b_l(4:6) = P_l2 * a_l(4:6); 
        a_l(7) = b_l(6);
        % P_l3 Junction
        b_l(7:9) = P_l3 * a_l(7:9);
        a_l(10) = b_l(9);
        % S_l4 Junction
        b_l(10:12) = S_l4 * a_l(10:12);

    %% Read Output

        % High  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   . 

        VoutHigh(ii) = -(a_h(6) + b_h(6))/2;
    
        % Mid   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  

        VoutMid(ii) = (a_m(14) + b_m(14))/2;
        
        % Low   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
        
        VoutLow(ii) = -(a_l(8) + b_l(8))/2;
    
end


%% Output Plots
figure
set(gcf, 'Color', 'w');
subplot(311)
plot(TsLTSpice*[1:length(OutLowSpice)],OutLowSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutLow,'b--','Linewidth',1); grid on; xlim([0,tstop]); 
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
legend('LTspice','WDF','Fontsize',16,'interpreter','latex');
title('Output Signals','Fontsize',18,'interpreter','latex');
subplot(312)
plot(TsLTSpice*[1:length(OutMidSpice)],OutMidSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutMid,'b--','Linewidth',1); grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
subplot(313)
plot(TsLTSpice*[1:length(OutHighSpice)],OutHighSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutHigh,'b--','Linewidth',1); grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');

%% Error Plots
figure
set(gcf, 'Color', 'w');
subplot(311)
plot(Ts*[1:Nsamp],OutLowSpice([1:downSampFact:end])-VoutLow,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
title(['Error Signals. $F_{\mathrm{s}}=$ ',num2str(Fs),' Hz'],'Fontsize',18,'interpreter','latex');
subplot(312)
plot(Ts*[1:Nsamp],OutMidSpice([1:downSampFact:end])-VoutMid,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
subplot(313)
plot(Ts*[1:Nsamp],OutHighSpice([1:downSampFact:end])-VoutHigh,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');

