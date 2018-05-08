function fitnessFcn = minimizeFun(param)

%% Multiple assessment points
expect = [0 36 19.11 20.4 16.68 12.95; 0 15.2 8.46 6.9 4.072 3.984; 0 4.424 2.088 1.316 1.864 0.86; %RP
          0 48.88 82.9 107.14 18.38 9.77; 0 10.95 50.44 100.88 26.19 7.04; 0 0 12.12 26.6 13.3 4.3; %MERKPP
          0 6.77 83.43 150.4 910.2 69.9; 0 20.1 26.87 157.07 473.94 53.33; 0 6.77 13.54 60.51 177.58 25.56; %ERKPP
          0 1.59 3.488 2.94 2.53 2.13; 0 1.03 1.61 1.41 1.06 0.76; 0 0.3 0.51 0.66 0.43 0.33];    %AktPIPP
time1 = [0;30;60;120;300;600];
time2 = [0;60;120;300;600;1800];

%% Time span
tspan1 = [min(time1),max(time1)];
tspan2 = [min(time2),max(time2)];

%% Assigning weights
w1 = 10;
w2 = 10;
w3 = 5;
w4 = 7;
w5 = 7;
w6 = 7;
w7 = 5;
w8 = 10;
w9 = 10;
w10 = 10;
w11 = 10;
w12 = 10;


%% Define parameters
prm(1) = 0.0012;%k1
prm(2) = 0.01;%k2
prm(3) = 1.0;%k3
prm(4) = 10^param(1);%V4 %62.5 in paper %450 in ref[5]
prm(5) = 0.09;%k5 %0.1 in paper %not match with ref[29] %0.09 in [5]
prm(6) = 10^param(2);%k6 %20 in paper %not match with ref[29]
prm(7) = 10^param(3);%k7 %Estimate
prm(8) = 10^param(4);%k8 %Estimate
prm(9) = 10^param(5);%k9 %Estimate
prm(10) = 10^param(6);%V10 %Estimate
prm(11) = 10^param(7);%k11 %Estimate
prm(12) = 10^param(8);%V12 %Estimate
prm(13) = 10^param(9);%k13 %Estimate
prm(14) = 10^param(10);%k14 %Estimate
prm(15) = 10^param(11);%k15 %Estimate
prm(16) = 0.058;%k16
prm(17) = 2.9;%k17
prm(18) = 0.058;%k18
prm(19) = 16;%k19 %9.5 in paper %16 from ref[4]
prm(20) = 0.3;%k20
prm(21) = 5.7;%k21 %16 in paper %5.7 in ref[4]
prm(22) = 0.27;%k22
prm(23) = 10^param(12);%K23 %not in ref[29] %0.1 in paper
prm(24) = 10^param(13);%k24  %not in ref[29] %0.0985 in paper
prm(25) = 10^param(14);%k25 %not in ref[29] %45.8 in paper
prm(26) = 10^param(15);%k26 %2620 in paper %Estimate
prm(27) = 10^param(16);%k27 %16.9 in paper %Estimate
prm(28) = 10^param(17);%k28 %17000 in paper %Estimate
prm(29) = 10^param(18);%k29 %507 in paper %Estimate
prm(30) = 20000;%k30;
prm(31) = 10^param(19);%k31 %0.107 in paper %Estimate
prm(32) = 20000;%v32
prm(33) = 10^param(20);%k33 %0.211 in paper %Estimate
prm(34) = 0.001;%k34
prm(35) = 7.6E-4;%k-1 
prm(36) = 0.1;%k-2
prm(37) = 0.01;%k-3
prm(38) = 50;%K4
prm(39) = 10^param(21);%k-5 %1 in paper &not match ref[5]
prm(40) = 10^param(22);%k-6 %5 in paper &not match ref[5]
prm(41) = 10^param(23);%k-7 %546 in paper %Estimate
prm(42) = 10^param(24);%k-8 %15700 in paper %Estimate
prm(43) = 10^param(25);%k-9 %0 in paper %Estimate
prm(44) = 340;%K10
prm(45) = 10^param(26);%K11 %0.181 in paper %Estimate
prm(46) = 10^param(27);%K12 %0.0571 in paper %Estimate
prm(47) = 10^param(28);%K13 %11.7 in paper %Estimate
prm(48) = 10^param(29);%K14 %8.07 in paper %Estimate
prm(49) = 317;%K15
prm(50) = 2200;%K16
prm(51) = 263;%K17 %317 in paper
prm(52) = 60;%K18
prm(53) = 146000;%K19
prm(54) = 160;%K20
prm(55) = 5.2E4;%K21 %1.46E5 in paper
prm(56) = 60;%K22 
prm(57) = 10^param(30);%k-23 %2 in paper %not found in ref[29]
prm(58) = 10^param(31);%k-24 %0.0985 in paper %not found in ref[29]
prm(59) = 10^param(32);%k-25 %0.047 in paper %not found in ref[29]
prm(60) = 10^param(33);%K26 %3680 in paper %Estimate
prm(61) = 10^param(34);%K27 %39.1 in paper %Estimate
prm(62) = 10^param(35);%k28 %9.02 in paper %Estimate
prm(63) = 10^param(36);%k29 %234 in paper %Estimate
prm(64) = 80000;%K30
prm(65) = 10^param(37);%K31 %4.35 in paper %Estimate
prm(66) = 80000;%K32
prm(67) = 10^param(38);%K33 %12 in paper %Estimate
prm(68) = 0;%k-34
prm(69) = 7.4;%Enzyme E
prm(70) = 11.4;%PP2A
prm(71) = 2.4;%MKP3
prm(72) = 1;%dummy for when rate equations do not fit particular form

%% Define initial conditions
%% Akt type
Akt_0 = 10;%initial Akt concentration (from Hatakeyama=10nm) x1 
AktPIP3_0 = 0;
AktPIP_0 = 0;
AktPIPP_0 = 0;
%% ERK type
ERK_0 = 1000;%initial ERK concentration (from Hatakeyama=1000nm)
ERKP_0 = 0;
ERKPP_0 = 0;
%% GS
GS_0 = 10;%initial GS concentration (from Hatakeyama=10nm) 
%%HRG
HRG_01 = 10;%initial HRG concentration (from Fig4 can be varied 0.1 1 10nM)
HRG_02 = 1;
HRG_03 = 0.1;
%% internalization
internalization_0 = 0;
%% MEK
MEK_0 = 120;%initial MEK concentration (from Hatakeyama=10nm)
MEKP_0 = 0;
MEKPP_0 = 0;
%% PI
PI_0 = 800;%initial PI concentration (from Hatakeyama=800nm)
PIP3_0 = 0;
%% PI3K
PI3K_0 = 10;%initial PI3k concentration (from Hatakeyama=10nm)
PI3KA_0 = 0;
%% R
R_0 = 80;%initial PI3k concentration (from Hatakeyama=80nm)
RP_0 = 0;
%% R plus other substrates
RHRG_0 = 0;
RHRG2_0 = 0;
RPI3K_0 = 0;
RPI3KA_0 = 0;
RShGS_0 = 0;
RShP_0 = 0;
RShc_0 = 0;
%% Raf
Raf_0 = 100;%initialRaf concentration (from Hatakeyama=100nm)
RafA_0 = 0;
%%Ras GdP
RasGDP_0 = 120;%initial RasGDPconcentration (from Hatakeyama=120nm)
RasGTP_0 = 0;
%% ShGS
ShGS_0 = 0;
%%SHP
ShP_0 = 0;
%%Shc
Shc_0 = 1000; %initial shc concentration (from Hatakeyama=1000nm)

%% Calculate fitness function
X_01 = [Akt_0,AktPIP3_0,AktPIP_0,AktPIPP_0,ERK_0,ERKP_0,ERKPP_0, GS_0,HRG_01,internalization_0,MEK_0,MEKP_0,MEKPP_0,PI_0,PI3K_0,PI3KA_0,PIP3_0,R_0,...
RP_0,RHRG_0,RHRG2_0,RPI3K_0,RPI3KA_0,RShGS_0,RShP_0,RShc_0,Raf_0,RafA_0,RasGDP_0,RasGTP_0,ShGS_0,ShP_0,Shc_0];

X_02 = [Akt_0,AktPIP3_0,AktPIP_0,AktPIPP_0,ERK_0,ERKP_0,ERKPP_0, GS_0,HRG_02,internalization_0,MEK_0,MEKP_0,MEKPP_0,PI_0,PI3K_0,PI3KA_0,PIP3_0,R_0,...
RP_0,RHRG_0,RHRG2_0,RPI3K_0,RPI3KA_0,RShGS_0,RShP_0,RShc_0,Raf_0,RafA_0,RasGDP_0,RasGTP_0,ShGS_0,ShP_0,Shc_0];

X_03 = [Akt_0,AktPIP3_0,AktPIP_0,AktPIPP_0,ERK_0,ERKP_0,ERKPP_0, GS_0,HRG_03,internalization_0,MEK_0,MEKP_0,MEKPP_0,PI_0,PI3K_0,PI3KA_0,PIP3_0,R_0,...
RP_0,RHRG_0,RHRG2_0,RPI3K_0,RPI3KA_0,RShGS_0,RShP_0,RShc_0,Raf_0,RafA_0,RasGDP_0,RasGTP_0,ShGS_0,ShP_0,Shc_0];

options=odeset('RelTol',1e-6,'AbsTol',1e-8);

%%Create solutions
sol1 = ode15s(@Hatakeyama_ODE,tspan1,X_01,options,prm);
y1 = deval(sol1, time1');
fitnessFcn(1) = sum((y1(22,:) - expect(1,:)).^2) + sum((y1(23,:)-expect(1,:)).^2) + sum((y1(26,:)-expect(1,:)).^2) + sum((y1(24,:)- expect(1,:)).^2) + sum((y1(25,:)- expect(1,:)).^2) + sum((y1(19,:) - expect(1,:)).^2); %RPI3K+RPI3KA+RShc+RShGS+RShP+RP
fitnessFcn(1) = w1*fitnessFcn(1);

sol2 = ode15s(@Hatakeyama_ODE,tspan1,X_02,options,prm);
y2 = deval(sol2, time1');
fitnessFcn(2) = sum((y2(22,:) - expect(2,:)).^2) + sum((y2(23,:)-expect(2,:)).^2) + sum((y2(26,:)-expect(2,:)).^2) + sum((y2(24,:)- expect(2,:)).^2) + sum((y2(25,:)- expect(2,:)).^2) + sum((y2(19,:) - expect(2,:)).^2);
fitnessFcn(2) = w2*fitnessFcn(2);

sol3 = ode15s(@Hatakeyama_ODE,tspan1,X_03,options,prm);
y3 = deval(sol3, time1');
fitnessFcn(3) = sum((y3(22,:) - expect(3,:)).^2) + sum((y3(23,:)-expect(3,:)).^2) + sum((y3(26,:)-expect(3,:)).^2) + sum((y3(24,:)- expect(3,:)).^2) + sum((y3(25,:)- expect(3,:)).^2) + sum((y3(19,:) - expect(3,:)).^2);
fitnessFcn(3) = w3*fitnessFcn(3);

sol4 = ode15s(@Hatakeyama_ODE,tspan2,X_01,options,prm);
y4 = deval(sol4, time2');
fitnessFcn(4) = sum((y4(13,:) - expect(4,:)).^2);
fitnessFcn(4) = w4*fitnessFcn(4);

sol5 = ode15s(@Hatakeyama_ODE,tspan2,X_02,options,prm);
y5 = deval(sol5, time2');
fitnessFcn(5) = sum((y5(13,:) - expect(5,:)).^2);
fitnessFcn(5) = w5*fitnessFcn(5);

sol6 = ode15s(@Hatakeyama_ODE,tspan2,X_03,options,prm);
y6 = deval(sol6, time2');
fitnessFcn(6) = sum((y6(13,:) - expect(6,:)).^2);
fitnessFcn(6) = w6*fitnessFcn(6);

sol7 = ode15s(@Hatakeyama_ODE,tspan2,X_01,options,prm);
y7 = deval(sol7, time2');
fitnessFcn(7) = sum((y7(7,:) - expect(7,:)).^2);
fitnessFcn(7) = w7*fitnessFcn(7);

sol8 = ode15s(@Hatakeyama_ODE,tspan2,X_02,options,prm);
y8 = deval(sol8, time2');
fitnessFcn(8) = sum((y8(7,:) - expect(8,:)).^2);
fitnessFcn(8) = w8*fitnessFcn(8);

sol9 = ode15s(@Hatakeyama_ODE,tspan2,X_03,options,prm);
y9 = deval(sol9, time2');
fitnessFcn(9) = sum((y9(7,:) - expect(9,:)).^2);
fitnessFcn(9) = w9*fitnessFcn(9);

sol10 = ode15s(@Hatakeyama_ODE,tspan2,X_01,options,prm);
y10 = deval(sol10, time2');
fitnessFcn(10) = sum((y10(4,:) - expect(10,:)).^2);
fitnessFcn(10) = w10*fitnessFcn(10);

sol11 = ode15s(@Hatakeyama_ODE,tspan2,X_02,options,prm);
y11 = deval(sol11, time2');
fitnessFcn(11) = sum((y11(4,:) - expect(11,:)).^2);
fitnessFcn(11) = w11*fitnessFcn(11);

sol12 = ode15s(@Hatakeyama_ODE,tspan2,X_03,options,prm);
y12 = deval(sol12, time2');
fitnessFcn(12) = sum((y12(4,:) - expect(12,:)).^2);
fitnessFcn(12) = w12*fitnessFcn(12);
end