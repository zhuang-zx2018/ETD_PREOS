clc
clear all

%%%% 常数
PI = 3.14159265359;
kb = 1.38E-23;

%%%% 状态方程系数（甲烷CH4）
ms=1; sigma=3.7412; epslon_k=153.36;
m=12.65; n=6;
%%%%%%%%%指数

T  = 101;

Eta1=0.001:0.001:0.1;
Eta2=0.4425:0.00001:0.4426;
Eta3=0.46017;%%%%%0.45064
Eta=[Eta1,Eta2,Eta3];

[Beta,eps,C,d]=dH(T,m,n,sigma,epslon_k);

rhos=Eta/(PI*d^3/6);
dE_rho=PI*d^3/6;
drho_E=6/(PI*d^3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hard sphere contribution  equation (11)      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ahs= ms*(4*Eta-3*Eta.^2)./(1-Eta).^2;             % (eq-11)  
 dahs= ms*(4-2*Eta)./(1-Eta).^3;
d2ahs= ms*(10-4*Eta)./(1-Eta).^4;
d3ahs= ms*(36-12*Eta)./(1-Eta).^5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dispersion contribution                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the effective packing fraction          (Eq-41)
 x=  [m,n,2*m,m+n,2*n];
c1=  0.81096 +1.7888./x -37.578./x.^2 +92.284./x.^3;
c2=  1.0205  -19.341./x +151.26./x.^2 -463.50./x.^3;
c3= -1.9057  +22.845./x -228.14./x.^2 +973.92./x.^3;
c4=  1.0885  -6.1962./x +106.98./x.^2 -677.64./x.^3;  

% effective packing fraction                     (Eq-40)
  Eta_e(1,:) =c1(1)*Eta +c2(1)*Eta.^2 +c3(1)*Eta.^3 +c4(1)*Eta.^4;
  Eta_e(2,:) =c1(2)*Eta +c2(2)*Eta.^2 +c3(2)*Eta.^3 +c4(2)*Eta.^4;
  Eta_e(3,:) =c1(3)*Eta +c2(3)*Eta.^2 +c3(3)*Eta.^3 +c4(3)*Eta.^4;
  Eta_e(4,:) =c1(4)*Eta +c2(4)*Eta.^2 +c3(4)*Eta.^3 +c4(4)*Eta.^4;
  Eta_e(5,:) =c1(5)*Eta +c2(5)*Eta.^2 +c3(5)*Eta.^3 +c4(5)*Eta.^4;

 dEta_e(1,:) =c1(1) +2*c2(1)*Eta +3*c3(1)*Eta.^2 +4*c4(1)*Eta.^3;
 dEta_e(2,:) =c1(2) +2*c2(2)*Eta +3*c3(2)*Eta.^2 +4*c4(2)*Eta.^3;
 dEta_e(3,:) =c1(3) +2*c2(3)*Eta +3*c3(3)*Eta.^2 +4*c4(3)*Eta.^3;
 dEta_e(4,:) =c1(4) +2*c2(4)*Eta +3*c3(4)*Eta.^2 +4*c4(4)*Eta.^3;
 dEta_e(5,:) =c1(5) +2*c2(5)*Eta +3*c3(5)*Eta.^2 +4*c4(5)*Eta.^3;

 d2Eta_e(1,:) =  2*c2(1) +6*c3(1)*Eta  +12*c4(1)*Eta.^2;
 d2Eta_e(2,:) =  2*c2(2) +6*c3(2)*Eta  +12*c4(2)*Eta.^2;
 d2Eta_e(3,:) =  2*c2(3) +6*c3(3)*Eta  +12*c4(3)*Eta.^2;
 d2Eta_e(4,:) =  2*c2(4) +6*c3(4)*Eta  +12*c4(4)*Eta.^2;
 d2Eta_e(5,:) =  2*c2(5) +6*c3(5)*Eta  +12*c4(5)*Eta.^2;

d3Eta_e(1,:) = 6*c3(1) +24*c4(1)*Eta;
d3Eta_e(2,:) = 6*c3(2) +24*c4(2)*Eta;
d3Eta_e(3,:) = 6*c3(3) +24*c4(3)*Eta;
d3Eta_e(4,:) = 6*c3(4) +24*c4(4)*Eta;
d3Eta_e(5,:) = 6*c3(5) +24*c4(5)*Eta;
% The expression corresponding to the mean attractive energy (Eq-39)
  t1=  1-Eta_e/2; 
 dt1= -1/2.*dEta_e;      
d2t1= -1/2.*d2Eta_e;
d3t1= -1/2.*d3Eta_e;

  t0= 1./(1-Eta_e);
  t2= t0.^3;
 dt2= 3*t0.^4.*dEta_e;
d2t2= 12*t0.^5.*dEta_e.^2 + 3*t0.^4.*d2Eta_e;
d3t2= 60*t0.^6.*dEta_e.^3 + 36*t0.^5.*dEta_e.*d2Eta_e+ 3*t0.^4.*d3Eta_e;

  t12= t1.*t2;
 dt12= dt1.*t2 + t1.*dt2;
d2t12= d2t1.*t2+ 2*dt1.*dt2+ t1.*d2t2;
d3t12= d3t1.*t2+ 3*d2t1.*dt2+ 3*dt1.*d2t2+ t1.*d3t2;

 c0 = -12*eps./(x-3);
 a1Sx(1,:) = c0(1).*Eta.*t12(1,:);         % Eq-39
 a1Sx(2,:) = c0(2).*Eta.*t12(2,:);  
 a1Sx(3,:) = c0(3).*Eta.*t12(3,:);  
 a1Sx(4,:) = c0(4).*Eta.*t12(4,:);  
 a1Sx(5,:) = c0(5).*Eta.*t12(5,:);  
 
da1Sx(1,:) = c0(1).*(t12(1,:)+Eta.*dt12(1,:));
da1Sx(2,:) = c0(2).*(t12(2,:)+Eta.*dt12(2,:));
da1Sx(3,:) = c0(3).*(t12(3,:)+Eta.*dt12(3,:));
da1Sx(4,:) = c0(4).*(t12(4,:)+Eta.*dt12(4,:));
da1Sx(5,:) = c0(5).*(t12(5,:)+Eta.*dt12(5,:));

d2a1Sx(1,:)= c0(1).*(2*dt12(1,:)+Eta.*d2t12(1,:));
d2a1Sx(2,:)= c0(2).*(2*dt12(2,:)+Eta.*d2t12(2,:));
d2a1Sx(3,:)= c0(3).*(2*dt12(3,:)+Eta.*d2t12(3,:));
d2a1Sx(4,:)= c0(4).*(2*dt12(4,:)+Eta.*d2t12(4,:));
d2a1Sx(5,:)= c0(5).*(2*dt12(5,:)+Eta.*d2t12(5,:));

d3a1Sx(1,:)= c0(1).*(3*d2t12(1,:)+Eta.*d3t12(1,:));
d3a1Sx(2,:)= c0(2).*(3*d2t12(2,:)+Eta.*d3t12(2,:));
d3a1Sx(3,:)= c0(3).*(3*d2t12(3,:)+Eta.*d3t12(3,:));
d3a1Sx(4,:)= c0(4).*(3*d2t12(4,:)+Eta.*d3t12(4,:));
d3a1Sx(5,:)= c0(5).*(3*d2t12(5,:)+Eta.*d3t12(5,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The first perturbation term
% The computation of the first integral I1A Eq-24
  x0=sigma/d;                                 %(Eq-17)
  I1A= C*(x0^n*a1Sx(2,:)-x0^m*a1Sx(1,:));         %(Eq-24)
 dI1A= C*(x0^n*da1Sx(2,:)-x0^m*da1Sx(1,:));
d2I1A= C*(x0^n*d2a1Sx(2,:)-x0^m*d2a1Sx(1,:));
d3I1A= C*(x0^n*d3a1Sx(2,:)-x0^m*d3a1Sx(1,:));

% The second integral I1B can then be written as Eq 28-29,33
Ix=-(x0.^(3-x)-1)./(x-3);
Jx=-(x0.^(4-x).*(x-3)-x0.^(3-x).*(x-4)-1)./(x-3)./(x-4);

  t1(1,:)=12.0*eps*Ix(1)*(Eta.*(1-Eta/2)./(1-Eta).^3);
  t1(2,:)=12.0*eps*Ix(2)*(Eta.*(1-Eta/2)./(1-Eta).^3);
  t1(3,:)=12.0*eps*Ix(3)*(Eta.*(1-Eta/2)./(1-Eta).^3);
  t1(4,:)=12.0*eps*Ix(4)*(Eta.*(1-Eta/2)./(1-Eta).^3);
  t1(5,:)=12.0*eps*Ix(5)*(Eta.*(1-Eta/2)./(1-Eta).^3);

 dt1(1,:)=12.0*eps*Ix(1)*(2+2*Eta-Eta.^2)./(1-Eta).^4/2;
 dt1(2,:)=12.0*eps*Ix(2)*(2+2*Eta-Eta.^2)./(1-Eta).^4/2;
 dt1(3,:)=12.0*eps*Ix(3)*(2+2*Eta-Eta.^2)./(1-Eta).^4/2;
 dt1(4,:)=12.0*eps*Ix(4)*(2+2*Eta-Eta.^2)./(1-Eta).^4/2;
 dt1(5,:)=12.0*eps*Ix(5)*(2+2*Eta-Eta.^2)./(1-Eta).^4/2;
 
d2t1(1,:)=12.0*eps*Ix(1)*(5+2*Eta-Eta.^2)./(1-Eta).^5;
d2t1(2,:)=12.0*eps*Ix(2)*(5+2*Eta-Eta.^2)./(1-Eta).^5;
d2t1(3,:)=12.0*eps*Ix(3)*(5+2*Eta-Eta.^2)./(1-Eta).^5;
d2t1(4,:)=12.0*eps*Ix(4)*(5+2*Eta-Eta.^2)./(1-Eta).^5;
d2t1(5,:)=12.0*eps*Ix(5)*(5+2*Eta-Eta.^2)./(1-Eta).^5;

d3t1(1,:)=12.0*eps*Ix(1)*3*(9+2*Eta-Eta.^2)./(1-Eta).^6;
d3t1(2,:)=12.0*eps*Ix(2)*3*(9+2*Eta-Eta.^2)./(1-Eta).^6;
d3t1(3,:)=12.0*eps*Ix(3)*3*(9+2*Eta-Eta.^2)./(1-Eta).^6;
d3t1(4,:)=12.0*eps*Ix(4)*3*(9+2*Eta-Eta.^2)./(1-Eta).^6;
d3t1(5,:)=12.0*eps*Ix(5)*3*(9+2*Eta-Eta.^2)./(1-Eta).^6;

  t2(1,:)=12.0*eps*Jx(1).*9.0/2*(Eta.^2.*(1+Eta)./(1-Eta).^3);
  t2(2,:)=12.0*eps*Jx(2).*9.0/2*(Eta.^2.*(1+Eta)./(1-Eta).^3);
  t2(3,:)=12.0*eps*Jx(3).*9.0/2*(Eta.^2.*(1+Eta)./(1-Eta).^3);  
  t2(4,:)=12.0*eps*Jx(4).*9.0/2*(Eta.^2.*(1+Eta)./(1-Eta).^3);
  t2(5,:)=12.0*eps*Jx(5).*9.0/2*(Eta.^2.*(1+Eta)./(1-Eta).^3);
  
  
 dt2(1,:)=12.0*eps*Jx(1)*9.0/2*(2*Eta+4*Eta.^2)./(1-Eta).^4;
 dt2(2,:)=12.0*eps*Jx(2)*9.0/2*(2*Eta+4*Eta.^2)./(1-Eta).^4; 
 dt2(3,:)=12.0*eps*Jx(3)*9.0/2*(2*Eta+4*Eta.^2)./(1-Eta).^4;
 dt2(4,:)=12.0*eps*Jx(4)*9.0/2*(2*Eta+4*Eta.^2)./(1-Eta).^4; 
 dt2(5,:)=12.0*eps*Jx(5)*9.0/2*(2*Eta+4*Eta.^2)./(1-Eta).^4;
 
d2t2(1,:)=12.0*eps*Jx(1)*9.0/2*2*(1+7*Eta+4*Eta.^2)./(1-Eta).^5;
d2t2(2,:)=12.0*eps*Jx(2)*9.0/2*2*(1+7*Eta+4*Eta.^2)./(1-Eta).^5;
d2t2(3,:)=12.0*eps*Jx(3)*9.0/2*2*(1+7*Eta+4*Eta.^2)./(1-Eta).^5;
d2t2(4,:)=12.0*eps*Jx(4)*9.0/2*2*(1+7*Eta+4*Eta.^2)./(1-Eta).^5;
d2t2(5,:)=12.0*eps*Jx(5)*9.0/2*2*(1+7*Eta+4*Eta.^2)./(1-Eta).^5;

d3t2(1,:)=12.0*eps*Jx(1)*9.0/2*24*(1+3*Eta+Eta.^2)./(1-Eta).^6;
d3t2(2,:)=12.0*eps*Jx(2)*9.0/2*24*(1+3*Eta+Eta.^2)./(1-Eta).^6;
d3t2(3,:)=12.0*eps*Jx(3)*9.0/2*24*(1+3*Eta+Eta.^2)./(1-Eta).^6;
d3t2(4,:)=12.0*eps*Jx(4)*9.0/2*24*(1+3*Eta+Eta.^2)./(1-Eta).^6;
d3t2(5,:)=12.0*eps*Jx(5)*9.0/2*24*(1+3*Eta+Eta.^2)./(1-Eta).^6;

 Bx=t1-t2;
dBx=dt1-dt2;
d2Bx=d2t1-d2t2;
d3Bx=d3t1-d3t2;

  I1B= C*(x0^n*Bx(2,:)-x0^m*Bx(1,:));            % Eq 32
 dI1B= C*(x0^n*dBx(2,:)-x0^m*dBx(1,:));  
d2I1B= C*(x0^n*d2Bx(2,:)-x0^m*d2Bx(1,:));
d3I1B= C*(x0^n*d3Bx(2,:)-x0^m*d3Bx(1,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the first-order term a1 for the Mie fluid
  a1=I1A+I1B;                                % Eq 21
 da1=dI1A+dI1B;
d2a1=d2I1A+d2I1B;
d3a1=d3I1A+d3I1B;         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % The second perturbation term  % Eq 16 
t1=(1-Eta).^4;  t2=(1+4*Eta+4*Eta.^2-4*Eta.^3+Eta.^4);

Khs=t1./t2;
dKhs=(-4*Eta.^5+32*Eta.^4-64*Eta.^3+40*Eta.^2+4*Eta-8)./t2.^2;
d2Khs=4*(3*Eta.^8-36*Eta.^7+140*Eta.^6-264*Eta.^5+276*Eta.^4-76*Eta.^3-108*Eta.^2+48*Eta+17)./t2.^3;
d3Khs=48*(-Eta.^11+16*Eta.^10-92*Eta.^9+287*Eta.^8-582*Eta.^7+724*Eta.^6-352*Eta.^5-238*Eta.^4+347*Eta.^3-12*Eta.^2-84*Eta-13)./t2.^4;


Qc=[7.5365557 , -359.44  , 1550.9 , -1.19932 ,  -1911.28 , 9236.9 ,  10
    -37.60463 ,  1825.6  ,-5070.1 , 9.063632 , 21390.175 ,-129430 ,  10
    71.745953 , -3168.0  , 6534.6 , -17.9482  , -51320.7 , 357230 ,0.57
    -46.83552 ,  1884.2  ,-3288.7 , 11.34027 ,  37064.54 ,-315530 ,-6.7
    -2.467982 ,-0.82376  ,-2.7171 , 20.52142 ,  1103.742 , 1390.2 ,  -8
     -0.50272 , -3.1935  , 2.0883 , -56.6377  , -3264.61 ,-4518.2 ,   0
    8.0956883,   3.7090  ,      0 , 40.53683 ,  2556.181 , 4241.6 ,   0];
 
a=C*(1/(n-3)-1/(m-3));   % Eq 18

% Eq 20
f1a=(Qc(1,1)+Qc(2,1)*a+Qc(3,1)*a^2+Qc(4,1)*a^3)/(1+Qc(5,1)*a+Qc(6,1)*a^2+Qc(7,1)*a^3);
f2a=(Qc(1,2)+Qc(2,2)*a+Qc(3,2)*a^2+Qc(4,2)*a^3)/(1+Qc(5,2)*a+Qc(6,2)*a^2+Qc(7,2)*a^3);
f3a=(Qc(1,3)+Qc(2,3)*a+Qc(3,3)*a^2+Qc(4,3)*a^3)/(1+Qc(5,3)*a+Qc(6,3)*a^2+Qc(7,3)*a^3);
f4a=(Qc(1,4)+Qc(2,4)*a+Qc(3,4)*a^2+Qc(4,4)*a^3)/(1+Qc(5,4)*a+Qc(6,4)*a^2+Qc(7,4)*a^3);
f5a=(Qc(1,5)+Qc(2,5)*a+Qc(3,5)*a^2+Qc(4,5)*a^3)/(1+Qc(5,5)*a+Qc(6,5)*a^2+Qc(7,5)*a^3);
f6a=(Qc(1,6)+Qc(2,6)*a+Qc(3,6)*a^2+Qc(4,6)*a^3)/(1+Qc(5,6)*a+Qc(6,6)*a^2+Qc(7,6)*a^3);


X=f1a*Eta*x0^3+f2a*(Eta*x0^3).^5+f3a*(Eta*x0^3).^8;   % Eq17
dX=f1a*x0^3+f2a*5*(Eta).^4*x0^15+f3a*8*(Eta).^7*x0^24;
d2X=f2a*20*(Eta).^3*x0^15+f3a*56*(Eta).^6*x0^24;
d3X=f2a*60*(Eta).^2*x0^15+f3a*56*6*(Eta).^5*x0^24;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1=eps/2*C^2*Khs.*(1+X);                          %Eq36
dt1=eps/2*C^2*(dKhs.*(1+X)+ Khs.*dX);
d2t1=eps/2*C^2*(d2Khs.*(1+X)+2*dKhs.*dX+Khs.*d2X);
d3t1=eps/2*C^2*(d3Khs.*(1+X)+3*d2Khs.*dX+3*dKhs.*d2X+Khs.*d3X);

t2=x0^x(5)*(a1Sx(5,:)+Bx(5,:))-2*x0^x(4)*(a1Sx(4,:)+Bx(4,:))+x0^x(3)*(a1Sx(3,:)+Bx(3,:));
dt2=x0^x(5)*(da1Sx(5,:)+dBx(5,:))-2*x0^x(4)*(da1Sx(4,:)+dBx(4,:))+x0^x(3)*(da1Sx(3,:)+dBx(3,:));
d2t2=x0^x(5)*(d2a1Sx(5,:)+d2Bx(5,:))-2*x0^x(4)*(d2a1Sx(4,:)+d2Bx(4,:))+x0^x(3)*(d2a1Sx(3,:)+d2Bx(3,:));
d3t2=x0^x(5)*(d3a1Sx(5,:)+d3Bx(5,:))-2*x0^x(4)*(d3a1Sx(4,:)+d3Bx(4,:))+x0^x(3)*(d3a1Sx(3,:)+d3Bx(3,:));

a2=t1.*t2;                         %Eq36
da2=dt1.*t2+t1.*dt2; 
d2a2=d2t1.*t2+2*dt1.*dt2+t1.*d2t2; 
d3a2=d3t1.*t2+3*d2t1.*dt2+3*dt1.*d2t2+t1.*d3t2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The three perturbation term        %Eq19
t1=f5a*Eta*x0^3+f6a*Eta.^2*x0^6;
dt1=f5a*x0^3+f6a*2*Eta*x0^6;
d2t1=f6a*2*x0^6;
d3t1=0;

t2=exp(t1);
dt2=t2.*dt1;
d2t2=dt2.*dt1+t2.*d2t1;
d3t2=d2t2.*dt1+2*dt2.*d2t1+t2.*d3t1;

a3=-eps^3*f4a*x0^3*Eta.*t2;                 %Eq19
da3=-eps^3*f4a*x0^3*(t2+Eta.*dt2);
d2a3=-eps^3*f4a*x0^3*(2*dt2+Eta.*d2t2);
d3a3=-eps^3*f4a*x0^3*(3*d2t2+Eta.*d3t2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SAFT-VR-Mie’s dispersion contribution   % Eq10
adsp=ms*(Beta*a1+Beta^2*a2+Beta^3*a3);         % Eq10-ahs
dadsp=ms*(Beta*da1+Beta^2*da2+Beta^3*da3);         
d2adsp=ms*(Beta*d2a1+Beta^2*d2a2+Beta^3*d2a3);   
d3adsp=ms*(Beta*d3a1+Beta^2*d3a2+Beta^3*d3a3);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% The SAFT-VR-Mie’s chain contribution (Eq-72)
% Eq 47
t1=42*Eta-39*Eta.^2+9*Eta.^3-2*Eta.^4;
dt1=42-39*2*Eta+27*Eta.^2-8*Eta.^3;
d2t1=-39*2+54*Eta-24*Eta.^2;
d3t1=54-48*Eta;

t0=1./(1-Eta);
t2=t0.^3/6;
dt2=t0.^4/2;
d2t2=t0.^5*2;
d3t2=t0.^6*10;

k0= -log(1-Eta)+t1.*t2;                            
dk0=t0+(dt1.*t2+t1.*dt2);
d2k0=t0.^2+(d2t1.*t2+2*dt1.*dt2+t1.*d2t2);
d3k0=2*t0.^3+(d3t1.*t2+3*d2t1.*dt2+3*dt1.*d2t2+t1.*d3t2);

% Eq 48
t1=Eta.^4+6*Eta.^2-12*Eta;
dt1=4*Eta.^3+12*Eta-12;
d2t1=12*Eta.^2+12;
d3t1=24*Eta;

k1=3*t1.*t2;                                          
dk1=3*(dt1.*t2+t1.*dt2); 
d2k1=3*(d2t1.*t2+2*dt1.*dt2+t1.*d2t2); 
d3k1=3*(d3t1.*t2+3*d2t1.*dt2+3*dt1.*d2t2+t1.*d3t2); 

% Eq 49
k2= -3/8*Eta.^2./(1-Eta).^2;                        
dk2=-3/4*Eta./(1-Eta).^3;
d2k2=-3/4*(1+2*Eta)./(1-Eta).^4;
d3k2=-9/2*(1+Eta)./(1-Eta).^5;

% Eq 50 
t1=-Eta.^4+3*Eta.^2+3*Eta;
dt1=-4*Eta.^3+6*Eta+3;
d2t1=-12*Eta.^2+6;
d3t1=-24*Eta;
       
 k3=t1.*t2;                                        
dk3=dt1.*t2+t1.*dt2; 
d2k3=d2t1.*t2+2*dt1.*dt2+t1.*d2t2; 
d3k3=d3t1.*t2+3*d2t1.*dt2+3*dt1.*d2t2+t1.*d3t2; 

t1=k0+k1*x0+k2*x0^2+k3*x0^3;
dt1=dk0+dk1*x0+dk2*x0^2+dk3*x0^3;
d2t1=d2k0+d2k1*x0+d2k2*x0^2+d2k3*x0^3;
d3t1=d3k0+d3k1*x0+d3k2*x0^2+d3k3*x0^3;

% Eq 46
gHs=exp(t1);                                      
dgHs=gHs.*dt1; 
d2gHs=dgHs.*dt1+gHs.*d2t1; 
d3gHs=d2gHs.*dt1+2*dgHs.*d2t1+gHs.*d3t1; 


% (the 1st term in Eq 64)
  a1dr=3*da1*dE_rho;  
 da1dr=3*d2a1*dE_rho; 
d2a1dr=3*d3a1*dE_rho; 

% rhos=Eta/(PI*d^3/6);
% dE_rho=PI*d^3/6;
% drho_E=6/(PI*d^3);
t0=1./rhos;
dt0=-t0.^2*drho_E;
d2t0=2*t0.^3*drho_E^2;

% (the rest term in Eq 64)
  t1=C*m*x0^m*(a1Sx(1,:)+Bx(1,:)).*t0;
 dt1=C*m*x0^m*((da1Sx(1,:)+dBx(1,:)).*t0 + (a1Sx(1,:)+Bx(1,:)).*dt0);
d2t1=C*m*x0^m*((d2a1Sx(1,:)+d2Bx(1,:)).*t0 +2*(da1Sx(1,:)+dBx(1,:)).*dt0 +(a1Sx(1,:)+Bx(1,:)).*d2t0);

  t2=C*n*x0^n*(a1Sx(2,:)+Bx(2,:)).*t0;
 dt2=C*n*x0^n*((da1Sx(2,:)+dBx(2,:)).*t0 + (a1Sx(2,:)+Bx(2,:)).*dt0);
d2t2=C*n*x0^n*((d2a1Sx(2,:)+d2Bx(2,:)).*t0 +2*(da1Sx(2,:)+dBx(2,:)).*dt0 +(a1Sx(2,:)+Bx(2,:)).*d2t0);

  gd1=1/(2*PI*eps*d^3)*(a1dr+t1-t2);
 dgd1=1/(2*PI*eps*d^3)*(da1dr+dt1-dt2);
d2gd1=1/(2*PI*eps*d^3)*(d2a1dr+d2t1-d2t2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Eq(66)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (the 1st term in Eq 66)
t1=1./(1+X);
dt1=-t1.^2.*dX;
d2t1=2*t1.^3.*dX.^2-t1.^2.*d2X;

  a2dr=3*(da2.*t1+a2.*dt1)*dE_rho;  
 da2dr=3*(d2a2.*t1+2*da2.*dt1+a2.*d2t1)*dE_rho;  
d2a2dr=3*(d3a2.*t1+3*d2a2.*dt1+3*da2.*d2t1+a2.*d3t1)*dE_rho;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   t1=Khs.*t0;
  dt1=dKhs.*t0+Khs.*dt0;
 d2t1=d2Khs.*t0+2*dKhs.*dt0+Khs.*d2t0;
 
  t11=eps*C^2*m*x0^x(3)*((a1Sx(3,:)+Bx(3,:)).*t1);
 dt11=eps*C^2*m*x0^x(3)*((da1Sx(3,:)+dBx(3,:)).*t1+(a1Sx(3,:)+Bx(3,:)).*dt1);
d2t11=eps*C^2*m*x0^x(3)*((d2a1Sx(3,:)+d2Bx(3,:)).*t1+2*(da1Sx(3,:)+dBx(3,:)).*dt1+(a1Sx(3,:)+Bx(3,:)).*d2t1);
 
  t12=eps*C^2*(m+n)*x0^x(4)*((a1Sx(4,:)+Bx(4,:)).*t1);
 dt12=eps*C^2*(m+n)*x0^x(4)*((da1Sx(4,:)+dBx(4,:)).*t1+(a1Sx(4,:)+Bx(4,:)).*dt1);
d2t12=eps*C^2*(m+n)*x0^x(4)*((d2a1Sx(4,:)+d2Bx(4,:)).*t1+2*(da1Sx(4,:)+dBx(4,:)).*dt1+(a1Sx(4,:)+Bx(4,:)).*d2t1);

  t13=eps*C^2*n*x0^x(5)*((a1Sx(5,:)+Bx(5,:)).*t1);
 dt13=eps*C^2*n*x0^x(5)*((da1Sx(5,:)+dBx(5,:)).*t1+(a1Sx(5,:)+Bx(5,:)).*dt1);
d2t13=eps*C^2*n*x0^x(5)*((d2a1Sx(5,:)+d2Bx(5,:)).*t1+2*(da1Sx(5,:)+dBx(5,:)).*dt1+(a1Sx(5,:)+Bx(5,:)).*d2t1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Eq 66)
  g2MCA= 1/(2*PI*eps^2*d^3)*(a2dr-t11+t12-t13);                 
 dg2MCA= 1/(2*PI*eps^2*d^3)*(da2dr-dt11+dt12-dt13);
d2g2MCA= 1/(2*PI*eps^2*d^3)*(d2a2dr-d2t11+d2t12-d2t13);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Eq 63)
t1=Qc(1,7)*(-tanh(Qc(2,7)*(Qc(3,7)-a))+1)*x0^3*(exp(Beta*eps)-1)*Eta;
dt1=t1./Eta;
d2t1=0;

t2=exp(Qc(4,7)*Eta*x0^3+Qc(5,7)*Eta.^2*x0^6);
dt2=t2.*(Qc(4,7)*x0^3+Qc(5,7)*2*Eta*x0^6);
d2t2=dt2.*(Qc(4,7)*x0^3+Qc(5,7)*2*Eta*x0^6)+t2*(Qc(5,7)*2*x0^6);

rc=t1.*t2;               
drc=dt1.*t2+t1.*dt2;
d2rc=d2t1.*t2+2*dt1.*dt2+t1.*d2t2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Eq 63)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Eq 62)
gd2=(1+rc).*g2MCA;           
dgd2=(1+rc).*dg2MCA+drc.*g2MCA;
d2gd2=(1+rc).*d2g2MCA+2*drc.*dg2MCA+d2rc.*g2MCA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Eq 62)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Eq 72)
  t0=1./gHs;
 dt0=-t0.^2.*dgHs;
d2t0=2*t0.^3.*dgHs.^2-t0.^2.*d2gHs;

  t1=(Beta*eps*gd1+(Beta*eps)^2*gd2).*t0;
 dt1=(Beta*eps*dgd1+(Beta*eps)^2*dgd2).*t0+(Beta*eps*gd1+(Beta*eps)^2*gd2).*dt0;
d2t1=(Beta*eps*d2gd1+(Beta*eps)^2*d2gd2).*t0+2*(Beta*eps*dgd1+(Beta*eps)^2*dgd2).*dt0+(Beta*eps*gd1+(Beta*eps)^2*gd2).*d2t0;

  gMie=gHs.*exp(t1);                                 
 dgMie=dgHs.*exp(t1)+gMie.*dt1;  
d2gMie=d2gHs.*exp(t1)+dgHs.*exp(t1).*dt1+dgMie.*dt1+gMie.*d2t1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Eq 72)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Eq 71)
  ach=(1-ms)*log(gMie);                             
 dach=(1-ms)*1./gMie.*dgMie;
d2ach=(1-ms)*(1./gMie.*d2gMie-1./gMie.^2.*dgMie.^2);

% residual Helmholtz’s energy
  Are=ahs+adsp+ach;
 dAre=dahs+dadsp+dach;
d2Are=d2ahs+d2adsp+d2ach;

Z=1+Eta.*dAre;
dZ=dAre+Eta.*d2Are;

P=Z*kb*T.*rhos*1e30;
dP=kb*T*1e30*(dZ.*rhos+Z*drho_E);

%%% Fugacity coefficients 
lnphi=Are+Z-1-log(abs(Z));
phi=exp(lnphi);
P0=[Eta;P;phi]';
% dphi=diff(phi,Eta);
% dPdp=dP/dphi;


% % 初始化
% Vtem=PI/6*ms*Vm*NA*dBH^3*1e-27;    % Eta reduced density
% dphi=1;                            % 汽液两相活度差
% Tol=1e-2;
% cof=1.0;
% m=0;
% time=1;
% while (dphi>Tol)
%     m=m+1;
%     for i=1:2
%         max=200;
%         j=0;  
%         V1=Vtem(i);
%         while(j<=max)
%             % 调节汽液相对比密度，使得汽液相压力相等
%             P1=subs(P,Eta,V1);
%             dP1=subs(dP,Eta,V1);
%             y1=(P1-P0)/P0;
%             V1=V1-(P1-P0)/dP1;
%             if (V1<0.0)
%                 V1=0.0;
%             end
%             if (abs(y1)<1e-10)
%                 j=max+1;
%             end
%             j=j+1;
%         end
%          phii(i)=subs(phi,Eta,V1); % 汽液相逸度
%          Vtem(i)=V1;               % 汽液相对比密度
%     end
%     % 调节压力使得汽液相逸度相等
%      if m>1
%          Dphi=abs(phii(1)-phii(2));
%          time=Dphi/abs(Dphi-dphi)
%          if time<1.0
%              time=1.0;
%          end
%          if Dphi-dphi>0
%              cof=-0.8*cof          % 压力调节调节系数
%          end
%      end
%          dphi=abs(phii(1)-phii(2))    % 平衡判据：汽液两相逸度之差
%          % 调节压力使得汽液两相逸度相等
%          if dphi>1e-2
%              Tol=1e-3;
%              P0=P0+500*cof*time;
%          elseif dphi>1e-3
%              Tol=1e-3;
%              P0=P0+200*cof*time;   
%          elseif dphi>1e-4
%              Tol=1e-4;
%              P0=P0+100*cof*time;
%          % elseif dphi>1e-4
%          %    Tol=1e-4;
%          %    P0=P0+50*cof;
%          end
% 		 P0
% end
% Vtem=Vtem/(PI/6*ms*dBH^3*NA*1e-27);
% % dPdp=subs(dPdp,Eta,V1);
