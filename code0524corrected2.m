
%function [x2] = LCVM_vaporliquid_ternary_corrected(x3)
T = 273.15+25;
R = 8.31451;
% 1 = CO2; 2 = Ethyl lactate; 3 = Sclareol
Tc(1) = 304.2; w(1) = 0.225; Pc(1) = 7.3815e6;
Tc(2) = 588; w(2) = 0.6; Pc(2) = 3.86e6;
Tc(3) = 845.5; w(3) = 1.023; Pc(3) = 1.58e6;

PB = 600*6894.757; z = 10;
% Sclareol physical properties
Ptp = 1e4; v3s = .300e-3; DeltaH_fus = 26225; Ttp = 376.46; DeltaCp = 0.0;
for i = 1:3
    K(i) = 0.37464 + 1.5422*w(i) - 0.26992*w(i)^2;
    alfa(i) = ( 1 + K(i)*( 1 - sqrt(T/Tc(i)) ) )^2;
    ai(i) = 0.45724*R^2*Tc(i)^2/Pc(i)*alfa(i);
    bi(i) = 0.07780*R*Tc(i)/Pc(i);
end
x(2) = 0.0000001; x2 = 0.01; x3 = 0.11
g1new = 50;
Ethyllactate_sat=40;
g2new = Ethyllactate_sat/PB;
Sclareol_sat = 4;
g3new = Sclareol_sat/PB;

while abs(x(2)-x2) > 10^(-6)
    x(3) = x3;
    x(2) = (1+(g1new-g3new)*x(3)-g1new)/(g2new-g1new);
    x(1) = 1-x(3)-x(2);
    
    y1 = g1new*x(1); y2 = g2new*x(2); y3 = g3new*x(3);
    m = y1 + y2 + y3;
    y(1) = y1/m; y(2) = y2/m; y(3) = y3/m;
    %%   % VAPOR CALCULATION
    bV = bi*y';
    % New mixing rule incorporating UNIFAC
    ulambda=0.36;
    Av=-0.623;
    Am=-0.52;
    %Inserting UNIFAC Excess Gibbs Energy function
    %C02 (1) : C02=1;
    %Ethyl-s-Lactate(2) : nCH2COO = 1; nCH3 = 2; nOH = 1; nCH = 1;
    %Sclareol(3) : nCH3 = 5; nCH2 = 7; nCH = 2; nC = 4; nCH2=CH2 = 1; nOH = 2;
    %r = Volume parameters
    %q = Area parameters
    %z = coordination number
    % PART E
    %CH3 = 1; CH2 = 2; CH = 3; C = 4; C=C = 5; OH = 6; CH2COO = 7; CO2 =8;
    r(1) = 0.9011; r(2) = 0.6744; r(3) = 0.4469; r(4) = 0.2195; r(5) = 1.3454; r(6) = 1.000; r(7) = 1.6764; r(8) =...
        1.2960;
    q(1) = 0.848; q(2) = 0.540; q(3) = 0.228; q(4) = 0.000; q(5) = 1.176; q(6) = 1.200; q(7) = 1.420; q(8) =...
        1.261;
    
    Rv(1) = 1*r(8);
    Rv(2) = 2*r(1)+1*r(3)+1*r(6)+1*r(7);
    Rv(3) = 5*r(1)+7*r(2)+2*r(3)+4*r(4)+1*r(5)+2*r(6);
    Q(1) = 1*q(8);
    Q(2) = 2*q(1)+1*q(3)+1*q(6)+1*q(7);
    Q(3) = 5*q(1)+7*q(2)+2*q(3)+4*q(4)+1*q(5)+2*q(6);
    Rsigma_V = Rv*y';
    Qsigma_V = Q*y';
    for i = 1:3
        l(i) = z/2*(Rv(i)-Q(i))-(Rv(i)-1);
        theta_V(i) = Rv(i)*y(i)/Rsigma_V;
        pi_V(i) = Q(i)*y(i)/Qsigma_V;
    end
    lsigma_V = l*y';
    for i=1:3
        ln_gamma_combinatorial_V(i)=log(theta_V(i)/y(i))+(z/2)*Q(i)*log(pi_V(i)/theta_V(i))+l(i)-theta_V(i)/y(i)*lsigma_V;
    end
    
    %a is in degree K
    a = [0, 0, 0, 0,-200.0, 986.5, 232.1, 116.7
        0, 0, 0, 0, -200.0, 986.5, 232.1, 116.7
        0, 0, 0, 0, -200.0, 986.5, 232.1, 116.7
        0, 0, 0, 0, -200.0, 986.5, 232.1, 116.7
        2520, 2520, 2520, 2520, 0, 693.9, 71.23, 48.57
        156.4, 156.4, 156.4, 156.4, 8694, 0, 101.1,471.83
        114.8, 114.8, 114.8, 114.8, 269.3, 245.4, 0, 102.75
        110.6,110.6,110.6,110.6,55.74,87.1,-126.9,0];
    
    for i=1:8
        for j=1:8
            w(i,j) = exp(-a(i,j)/T);
        end
    end
    % PART F
    %For pure CO2
    %p indicates pure
    %CO2 = 1
    Y_p1(1)=0; Y_p1(2)=0; Y_p1(3)=0; Y_p1(4)=0; Y_p1(5)=0; Y_p1(6)=0; Y_p1(7)=0;Y_p1(8) = 1;
    Area_thetaSigma_p1_V=q*Y_p1';
    
    for i=1:8
        Area_theta_p1_V(i) = q(i)*Y_p1(i)/Area_thetaSigma_p1_V;
    end
    b1_V(8) = log( Area_theta_p1_V(8)*w(8,8) );
    c1_V(8) = 1;
    ln_t_p1_V(8) = q(8) * (1 - b1_V(8) - c1_V(8));
    % PART G
    %For pure Ethyl lactate
    %p indicates pure
    Y_p2(1)=2/5; Y_p2(2)=0; Y_p2(3)=1/5; Y_p2(4)=0; Y_p2(5)=0; Y_p2(6)=1/5; Y_p2(7)=1/5;
    Y_p2(8) = 0;
    Area_thetaSigma_p2_V=q*Y_p2';
    for i=1:8
        Area_theta_p2_V(i) = q(i)*Y_p2(i)/Area_thetaSigma_p2_V;
    end
    for i=1:8
        b2_V(i)=log(Area_theta_p2_V(1)*w(1,i)+Area_theta_p2_V(3)*w(3,i)+ Area_theta_p2_V(6)*w(6,i) +...
            Area_theta_p2_V(7)*w(7,i) );
        c2_V(i)=...
            Area_theta_p2_V(1)*w(i,1)/(Area_theta_p2_V(1)*w(1,1)+Area_theta_p2_V(3)*w(3,1)+...
            Area_theta_p2_V(6)*w(6,1) + Area_theta_p2_V(7)*w(7,1)) +...
            Area_theta_p2_V(3)*w(i,3)/(Area_theta_p2_V(1)*w(1,3)+Area_theta_p2_V(3)*w(3,3)+...
            Area_theta_p2_V(6)*w(6,3) + Area_theta_p2_V(7)*w(7,3)) +...
            Area_theta_p2_V(6)*w(i,6)/(Area_theta_p2_V(1)*w(1,6)+Area_theta_p2_V(3)*w(3,6)+...
            Area_theta_p2_V(6)*w(6,6) + Area_theta_p2_V(7)*w(7,6)) +...
            Area_theta_p2_V(7)*w(i,7)/(Area_theta_p2_V(1)*w(1,7)+Area_theta_p2_V(3)*w(3,7)+...
            Area_theta_p2_V(6)*w(6,7) + Area_theta_p2_V(7)*w(7,7));
        ln_t_p2_V(i) = q(i) * (1 - b2_V(i) - c2_V(i));
    end
    % PART H
    %For pure Sclareol
    %p indicates pure
    Y_p3(1)=5/21; Y_p3(2)=7/21; Y_p3(3)=2/21; Y_p3(4)=4/21; Y_p3(5)=1/21;
    Y_p3(6)=2/21; Y_p3(7)=0; Y_p3(8)=0;
    Area_thetaSigma_p3_V = q*Y_p3';
    
    for i=1:8
        Area_theta_p3_V(i) = q(i)*Y_p3(i)/Area_thetaSigma_p3_V;
    end
    for i=1:8
        b3_V(i)=log(Area_theta_p3_V(1)*w(1,i)+Area_theta_p3_V(2)*w(2,i)+Area_theta_p3_V(3)*w(3,i)+...
            Area_theta_p3_V(4)*w(4,i) + Area_theta_p3_V(5)*w(5,i) + Area_theta_p3_V(6)*w(6,i) );
        c3_V(i)=...
            Area_theta_p3_V(1)*w(i,1)/(Area_theta_p3_V(1)*w(1,1)+Area_theta_p3_V(2)*w(2,1)+...
            Area_theta_p3_V(3)*w(3,1)+Area_theta_p3_V(4)*w(4,1)+Area_theta_p3_V(5)*w(5,1)+...
            Area_theta_p3_V(6)*w(6,1))+...
            Area_theta_p3_V(2)*w(i,2)/(Area_theta_p3_V(1)*w(1,2)+Area_theta_p3_V(2)*w(2,2)+...
            Area_theta_p3_V(3)*w(3,2)+Area_theta_p3_V(4)*w(4,2)+Area_theta_p3_V(5)*w(5,2)+...
            Area_theta_p3_V(6)*w(6,2))+...
            Area_theta_p3_V(3)*w(i,3)/(Area_theta_p3_V(1)*w(1,3)+Area_theta_p3_V(2)*w(2,3)+...
            Area_theta_p3_V(3)*w(3,3)+Area_theta_p3_V(4)*w(4,3)+Area_theta_p3_V(5)*w(5,3)+...
            Area_theta_p3_V(6)*w(6,3))+...
            Area_theta_p3_V(4)*w(i,4)/(Area_theta_p3_V(1)*w(1,4)+Area_theta_p3_V(2)*w(2,4)+...
            Area_theta_p3_V(3)*w(3,4)+Area_theta_p3_V(4)*w(4,4)+Area_theta_p3_V(5)*w(5,4)+...
            Area_theta_p3_V(6)*w(6,4))+...
            Area_theta_p3_V(5)*w(i,5)/(Area_theta_p3_V(1)*w(1,5)+Area_theta_p3_V(2)*w(2,5)+...
            Area_theta_p3_V(3)*w(3,5)+Area_theta_p3_V(4)*w(4,5)+Area_theta_p3_V(5)*w(5,5)+...
            Area_theta_p3_V(6)*w(6,5))+...
            Area_theta_p3_V(6)*w(i,6)/(Area_theta_p3_V(1)*w(1,6)+Area_theta_p3_V(2)*w(2,6)+...
            Area_theta_p3_V(3)*w(3,6)+Area_theta_p3_V(4)*w(4,6)+Area_theta_p3_V(5)*w(5,6)+...
            Area_theta_p3_V(6)*w(6,6));
        ln_t_p3_V(i) = q(i) * (1 - b3_V(i) - c3_V(i));
    end
    % PART I
    %For Ternary
    Ysigma = y(1)*1+y(2)*5+y(3)*21;
    Y(1) = (y(1)*0 + y(2)*2 + y(3)*5)/Ysigma;
    Y(2) = (y(1)*0 + y(2)*0 + y(3)*7)/Ysigma;
    Y(3) = (y(1)*0 + y(2)*1 + y(3)*2)/Ysigma;
    Y(4) = (y(1)*0 + y(2)*0 + y(3)*4)/Ysigma;
    Y(5) = (y(1)*0 + y(2)*0 + y(3)*1)/Ysigma;
    Y(6) = (y(1)*0 + y(2)*1 + y(3)*2)/Ysigma;
    Y(7) = (y(1)*0 + y(2)*1 + y(3)*0)/Ysigma;
    Y(8) = (y(1)*1 + y(2)*0 + y(3)*0)/Ysigma;
    
    Area_thetaSigma_V=q*Y';
    
    for i=1:8
        Area_theta_V(i) = q(i)*Y(i)/Area_thetaSigma_V;
    end
    for i=1:8
        b_V(i)=log(Area_theta_V(1)*w(1,i)+Area_theta_V(2)*w(2,i)+Area_theta_V(3)*w(3,i)+Area_theta_V...
            (4)*w(4,i)+Area_theta_V(5)*w(5,i)+Area_theta_V(6)*w(6,i)+Area_theta_V(7)*w(7,i)+...
            Area_theta_V(8)*w(8,i));
        c_V(i)=...
            Area_theta_V(1)*w(i,1)/(Area_theta_V(1)*w(1,1)+Area_theta_V(2)*w(2,1)+Area_theta_V(3)*w(3,1)...
            +Area_theta_V(4)*w(4,1)+Area_theta_V(5)*w(5,1)+Area_theta_V(6)*w(6,1)+Area_theta_V(7)*w(7,...
            1)+Area_theta_V(8)*w(8,1))+...
            Area_theta_V(2)*w(i,2)/(Area_theta_V(1)*w(1,2)+Area_theta_V(2)*w(2,2)+Area_theta_V(3)*w(3,2)...
            +Area_theta_V(4)*w(4,2)+Area_theta_V(5)*w(5,2)+Area_theta_V(6)*w(6,2)+Area_theta_V(7)*w(7,...
            2)+Area_theta_V(8)*w(8,2))+...
            Area_theta_V(3)*w(i,3)/(Area_theta_V(1)*w(1,3)+Area_theta_V(2)*w(2,3)+Area_theta_V(3)*w(3,3)...
            +Area_theta_V(4)*w(4,3)+Area_theta_V(5)*w(5,3)+Area_theta_V(6)*w(6,3)+Area_theta_V(7)*w(7,...
            3)+Area_theta_V(8)*w(8,3))+...
            Area_theta_V(4)*w(i,4)/(Area_theta_V(1)*w(1,4)+Area_theta_V(2)*w(2,4)+Area_theta_V(3)*w(3,4)...
            +Area_theta_V(4)*w(4,4)+Area_theta_V(5)*w(5,4)+Area_theta_V(6)*w(6,4)+Area_theta_V(7)*w(7,...
            4)+Area_theta_V(8)*w(8,4))+...
            Area_theta_V(5)*w(i,5)/(Area_theta_V(1)*w(1,5)+Area_theta_V(2)*w(2,5)+Area_theta_V(3)*w(3,5)...
            +Area_theta_V(4)*w(4,5)+Area_theta_V(5)*w(5,5)+Area_theta_V(6)*w(6,5)+Area_theta_V(7)*w(7,...
            5)+Area_theta_V(8)*w(8,5))+...
            Area_theta_V(6)*w(i,6)/(Area_theta_V(1)*w(1,6)+Area_theta_V(2)*w(2,6)+Area_theta_V(3)*w(3,6)...
            +Area_theta_V(4)*w(4,6)+Area_theta_V(5)*w(5,6)+Area_theta_V(6)*w(6,6)+Area_theta_V(7)*w(7,...
            6)+Area_theta_V(8)*w(8,6))+...
            Area_theta_V(7)*w(i,7)/(Area_theta_V(1)*w(1,7)+Area_theta_V(2)*w(2,7)+Area_theta_V(3)*w(3,7)...
            +Area_theta_V(4)*w(4,7)+Area_theta_V(5)*w(5,7)+Area_theta_V(6)*w(6,7)+Area_theta_V(7)*w(7,...
            7)+Area_theta_V(8)*w(8,7))+...
            Area_theta_V(8)*w(i,8)/(Area_theta_V(1)*w(1,8)+Area_theta_V(2)*w(2,8)+Area_theta_V(3)*w(3,8)...
            +Area_theta_V(4)*w(4,8)+Area_theta_V(5)*w(5,8)+Area_theta_V(6)*w(6,8)+Area_theta_V(7)*w(7,...
            8)+Area_theta_V(8)*w(8,8));
        ln_t_V(i) = q(i) * (1 - b_V(i) - c_V(i));
    end
    ln_gamma_residual_V(1)=1*( ln_t_V(8) - ln_t_p1_V(8) );
    ln_gamma_residual_V(2)=2*( ln_t_V(1) - ln_t_p2_V(1) )+1*( ln_t_V(3) - ln_t_p2_V(3) )+1*( ln_t_V(6)...
        - ln_t_p2_V(6) ) + 1*( ln_t_V(7) - ln_t_p2_V(7) );
    ln_gamma_residual_V(3)=5*( ln_t_V(1) - ln_t_p3_V(1) )+7*( ln_t_V(2) - ln_t_p3_V(2) )+2*( ln_t_V(3)...
        - ln_t_p3_V(3) ) + 4*( ln_t_V(4) - ln_t_p3_V(4) ) + 1*( ln_t_V(5) - ln_t_p3_V(5) ) + 2*( ln_t_V(6) -...
        ln_t_p3_V(6) );
    
    ln_gamma_total_V(1) = ln_gamma_combinatorial_V(1) + ln_gamma_residual_V(1);
    ln_gamma_total_V(2) = ln_gamma_combinatorial_V(2) + ln_gamma_residual_V(2);
    ln_gamma_total_V(3) = ln_gamma_combinatorial_V(3) + ln_gamma_residual_V(3);
    % Introduction of GE function
    
    GE_V=(y(1)*ln_gamma_total_V(1)+y(2)*ln_gamma_total_V(2)+y(3)*ln_gamma_total_V(3))*R*T;
    
    %Inserting the new mixing rule incorporating UNIFAC
    
    for i=1:3
        bratioV(i) = log(bV/bi(i));
        alphai(i) = ai(i)/bi(i);
    end
    
    aV1=(ulambda/Av + (1-ulambda)/Am)*GE_V/R/T;
    aV2=((1-ulambda)/Am)*(y*bratioV');
    aV3=y*alphai';
    jV= aV3*bV;
    wV= (aV1+aV2)*R*T*bV;
    aV= jV+wV;
    
    AV = aV*PB/(R^2*T^2);
    BV = bV*PB/(R*T);
    
    % V(1)*X^N + ... + V(N)*X + V(N+1)
    % Peng Robinson coefficients of the cubic eq.
    V(1) = 1;
    V(2) = BV - 1;
    V(3) = AV - 3*BV^2 - 2*BV;
    V(4) = BV^3 + BV^2 - AV*BV;
    % Solve for the roots.
    ZtempV = roots(V);
    % I will consider only the largest value. The smallest value for liquid phase.
    j = 1;
    for i = 1:3
        if imag(ZtempV(i))==0
            zV(j) = ZtempV(i);
            j = j + 1;
        end
    end
    ZV = max(zV);
    for i = 1:3
        alpha_prime_V(i)=(ulambda/Av+(1-ulambda)/Am)*ln_gamma_total_V(i)+(1-...
            ulambda)/Am*(log(bV/bi(i))+bi(i)/bV-1) + ai(i)/bi(i)/R/T;
        phiV(i)= exp( bi(i)/bV*(ZV-1) - log(ZV-BV) - alpha_prime_V(i)/2.828*log((ZV+2.414*BV)/(ZV-0.414*BV)) );
    end
    %%
    % LIQUID CALCULATION
    % the same calculation procedure as the vapor phase can be applied to obtain the fugacity coefficient in the liquid phase
    bL = bi*x';
    % New mixing rule incorporating UNIFAC
    ulambda=0.36;
    Av=-0.623;
    Am=-0.52;
    %Inserting UNIFAC Excess Gibbs Energy function
    %C02 (1) : C02=1;
    %Ethyl-s-Lactate(2) : nCH2COO = 1; nCH3 = 2; nOH = 1; nCH = 1;
    %Sclareol(3) : nCH3 = 5; nCH2 = 7; nCH = 2; nC = 4; nCH2=CH2 = 1; nOH = 2;
    %r = Lolume parameters
    %q = Area parameters
    %z = coordination number
    % PART E
    %CH3 = 1; CH2 = 2; CH = 3; C = 4; C=C = 5; OH = 6; CH2COO = 7; CO2 =8;
    r(1) = 0.9011; r(2) = 0.6744; r(3) = 0.4469; r(4) = 0.2195; r(5) = 1.3454; r(6) = 1.000; r(7) = 1.6764; r(8) =...
        1.2960;
    q(1) = 0.848; q(2) = 0.540; q(3) = 0.228; q(4) = 0.000; q(5) = 1.176; q(6) = 1.200; q(7) = 1.420; q(8) =...
        1.261;
    
    Rv(1) = 1*r(8);
    Rv(2) = 2*r(1)+1*r(3)+1*r(6)+1*r(7);
    Rv(3) = 5*r(1)+7*r(2)+2*r(3)+4*r(4)+1*r(5)+2*r(6);
    Q(1) = 1*q(8);
    Q(2) = 2*q(1)+1*q(3)+1*q(6)+1*q(7);
    Q(3) = 5*q(1)+7*q(2)+2*q(3)+4*q(4)+1*q(5)+2*q(6);
    Rsigma_L = Rv*x';
    Qsigma_L = Q*x';
    for i = 1:3
        l(i) = z/2*(Rv(i)-Q(i))-(Rv(i)-1);
        theta_L(i) = Rv(i)*x(i)/Rsigma_L;
        pi_L(i) = Q(i)*x(i)/Qsigma_L;
    end
    lsigma_L = l*x';
    for i=1:3
        ln_gamma_combinatorial_L(i)=log(theta_L(i)/x(i))+z/2*Q(i)*log(pi_L(i)/theta_L(i))+l(i)-theta_L(i)/x(i)*lsigma_L;
    end
    
    %a is in degree K
    a = [0, 0, 0, 0,-200.0, 986.5, 232.1, 116.7
        0, 0, 0, 0, -200.0, 986.5, 232.1, 116.7
        0, 0, 0, 0, -200.0, 986.5, 232.1, 116.7
        0, 0, 0, 0, -200.0, 986.5, 232.1, 116.7
        2520, 2520, 2520, 2520, 0, 693.9, 71.23, 48.57
        156.4, 156.4, 156.4, 156.4, 8694, 0, 101.1,471.83
        114.8, 114.8, 114.8, 114.8, 269.3, 245.4, 0, 102.75
        110.6,110.6,110.6,110.6,55.74,87.1,-126.9,0];
    
    for i=1:8
        for j=1:8
            w(i,j) = exp(-a(i,j)/T);
        end
    end
    % PART F
    %For pure CO2
    %p indicates pure
    %CO2 = 1
    X_p1(1)=0; X_p1(2)=0; X_p1(3)=0; X_p1(4)=0; X_p1(5)=0; X_p1(6)=0; X_p1(7)=0;X_p1(8) = 1;
    Area_thetaSigma_p1_L=q*X_p1';
    
    for i=1:8
        Area_theta_p1_L(i) = q(i)*X_p1(i)/Area_thetaSigma_p1_L;
    end
    b1_L(8) = log( Area_theta_p1_L(8)*w(8,8) );
    c1_L(8) = 1;
    ln_t_p1_L(8) = q(8) * (1 - b1_L(8) - c1_L(8));
    % PART G
    %For pure Ethyl lactate
    %p indicates pure
    X_p2(1)=2/5; X_p2(2)=0; X_p2(3)=1/5; X_p2(4)=0; X_p2(5)=0; X_p2(6)=1/5; X_p2(7)=1/5;
    X_p2(8) = 0;
    Area_thetaSigma_p2_L=q*X_p2';
    for i=1:8
        Area_theta_p2_L(i) = q(i)*X_p2(i)/Area_thetaSigma_p2_L;
    end
    for i=1:8
        b2_L(i)=log(Area_theta_p2_L(1)*w(1,i)+Area_theta_p2_L(3)*w(3,i)+ Area_theta_p2_L(6)*w(6,i) +...
            Area_theta_p2_L(7)*w(7,i) );
        c2_L(i)=...
            Area_theta_p2_L(1)*w(i,1)/(Area_theta_p2_L(1)*w(1,1)+Area_theta_p2_L(3)*w(3,1)+...
            Area_theta_p2_L(6)*w(6,1) + Area_theta_p2_L(7)*w(7,1)) +...
            Area_theta_p2_L(3)*w(i,3)/(Area_theta_p2_L(1)*w(1,3)+Area_theta_p2_L(3)*w(3,3)+...
            Area_theta_p2_L(6)*w(6,3) + Area_theta_p2_L(7)*w(7,3)) +...
            Area_theta_p2_L(6)*w(i,6)/(Area_theta_p2_L(1)*w(1,6)+Area_theta_p2_L(3)*w(3,6)+...
            Area_theta_p2_L(6)*w(6,6) + Area_theta_p2_L(7)*w(7,6)) +...
            Area_theta_p2_L(7)*w(i,7)/(Area_theta_p2_L(1)*w(1,7)+Area_theta_p2_L(3)*w(3,7)+...
            Area_theta_p2_L(6)*w(6,7) + Area_theta_p2_L(7)*w(7,7));
        ln_t_p2_L(i) = q(i) * (1 - b2_L(i) - c2_L(i));
    end
    % PART H
    %For pure Sclareol
    %p indicates pure
    X_p3(1)=5/21; X_p3(2)=7/21; X_p3(3)=2/21; X_p3(4)=4/21; X_p3(5)=1/21;
    X_p3(6)=2/21; X_p3(7)=0; X_p3(8)=0;
    Area_thetaSigma_p3_L = q*X_p3';
    
    for i=1:8
        Area_theta_p3_L(i) = q(i)*X_p3(i)/Area_thetaSigma_p3_L;
    end
    for i=1:8
        b3_L(i)=log(Area_theta_p3_L(1)*w(1,i)+Area_theta_p3_L(2)*w(2,i)+Area_theta_p3_L(3)*w(3,i)+...
            Area_theta_p3_L(4)*w(4,i) + Area_theta_p3_L(5)*w(5,i) + Area_theta_p3_L(6)*w(6,i) );
        c3_L(i)=...
            Area_theta_p3_L(1)*w(i,1)/(Area_theta_p3_L(1)*w(1,1)+Area_theta_p3_L(2)*w(2,1)+...
            Area_theta_p3_L(3)*w(3,1)+Area_theta_p3_L(4)*w(4,1)+Area_theta_p3_L(5)*w(5,1)+...
            Area_theta_p3_L(6)*w(6,1))+...
            Area_theta_p3_L(2)*w(i,2)/(Area_theta_p3_L(1)*w(1,2)+Area_theta_p3_L(2)*w(2,2)+...
            Area_theta_p3_L(3)*w(3,2)+Area_theta_p3_L(4)*w(4,2)+Area_theta_p3_L(5)*w(5,2)+...
            Area_theta_p3_L(6)*w(6,2))+...
            Area_theta_p3_L(3)*w(i,3)/(Area_theta_p3_L(1)*w(1,3)+Area_theta_p3_L(2)*w(2,3)+...
            Area_theta_p3_L(3)*w(3,3)+Area_theta_p3_L(4)*w(4,3)+Area_theta_p3_L(5)*w(5,3)+...
            Area_theta_p3_L(6)*w(6,3))+...
            Area_theta_p3_L(4)*w(i,4)/(Area_theta_p3_L(1)*w(1,4)+Area_theta_p3_L(2)*w(2,4)+...
            Area_theta_p3_L(3)*w(3,4)+Area_theta_p3_L(4)*w(4,4)+Area_theta_p3_L(5)*w(5,4)+...
            Area_theta_p3_L(6)*w(6,4))+...
            Area_theta_p3_L(5)*w(i,5)/(Area_theta_p3_L(1)*w(1,5)+Area_theta_p3_L(2)*w(2,5)+...
            Area_theta_p3_L(3)*w(3,5)+Area_theta_p3_L(4)*w(4,5)+Area_theta_p3_L(5)*w(5,5)+...
            Area_theta_p3_L(6)*w(6,5))+...
            Area_theta_p3_L(6)*w(i,6)/(Area_theta_p3_L(1)*w(1,6)+Area_theta_p3_L(2)*w(2,6)+...
            Area_theta_p3_L(3)*w(3,6)+Area_theta_p3_L(4)*w(4,6)+Area_theta_p3_L(5)*w(5,6)+...
            Area_theta_p3_L(6)*w(6,6));
        ln_t_p3_L(i) = q(i) * (1 - b3_L(i) - c3_L(i));
    end
    % PART I
    %For Ternary
    Xsigma = x(1)*1+x(2)*5+x(3)*21;
    X(1) = (x(1)*0 + x(2)*2 + x(3)*5)/Xsigma;
    X(2) = (x(1)*0 + x(2)*0 + x(3)*7)/Xsigma;
    X(3) = (x(1)*0 + x(2)*1 + x(3)*2)/Xsigma;
    X(4) = (x(1)*0 + x(2)*0 + x(3)*4)/Xsigma;
    X(5) = (x(1)*0 + x(2)*0 + x(3)*1)/Xsigma;
    X(6) = (x(1)*0 + x(2)*1 + x(3)*2)/Xsigma;
    X(7) = (x(1)*0 + x(2)*1 + x(3)*0)/Xsigma;
    X(8) = (x(1)*1 + x(2)*0 + x(3)*0)/Xsigma;
    
    Area_thetaSigma_L=q*X';
    
    for i=1:8
        Area_theta_L(i) = q(i)*X(i)/Area_thetaSigma_L;
    end
    for i=1:8
        b_L(i)=log(Area_theta_L(1)*w(1,i)+Area_theta_L(2)*w(2,i)+Area_theta_L(3)*w(3,i)+Area_theta_L...
            (4)*w(4,i)+Area_theta_L(5)*w(5,i)+Area_theta_L(6)*w(6,i)+Area_theta_L(7)*w(7,i)+...
            Area_theta_L(8)*w(8,i));
        c_L(i)=...
            Area_theta_L(1)*w(i,1)/(Area_theta_L(1)*w(1,1)+Area_theta_L(2)*w(2,1)+Area_theta_L(3)*w(3,1)...
            +Area_theta_L(4)*w(4,1)+Area_theta_L(5)*w(5,1)+Area_theta_L(6)*w(6,1)+Area_theta_L(7)*w(7,...
            1)+Area_theta_L(8)*w(8,1))+...
            Area_theta_L(2)*w(i,2)/(Area_theta_L(1)*w(1,2)+Area_theta_L(2)*w(2,2)+Area_theta_L(3)*w(3,2)...
            +Area_theta_L(4)*w(4,2)+Area_theta_L(5)*w(5,2)+Area_theta_L(6)*w(6,2)+Area_theta_L(7)*w(7,...
            2)+Area_theta_L(8)*w(8,2))+...
            Area_theta_L(3)*w(i,3)/(Area_theta_L(1)*w(1,3)+Area_theta_L(2)*w(2,3)+Area_theta_L(3)*w(3,3)...
            +Area_theta_L(4)*w(4,3)+Area_theta_L(5)*w(5,3)+Area_theta_L(6)*w(6,3)+Area_theta_L(7)*w(7,...
            3)+Area_theta_L(8)*w(8,3))+...
            Area_theta_L(4)*w(i,4)/(Area_theta_L(1)*w(1,4)+Area_theta_L(2)*w(2,4)+Area_theta_L(3)*w(3,4)...
            +Area_theta_L(4)*w(4,4)+Area_theta_L(5)*w(5,4)+Area_theta_L(6)*w(6,4)+Area_theta_L(7)*w(7,...
            4)+Area_theta_L(8)*w(8,4))+...
            Area_theta_L(5)*w(i,5)/(Area_theta_L(1)*w(1,5)+Area_theta_L(2)*w(2,5)+Area_theta_L(3)*w(3,5)...
            +Area_theta_L(4)*w(4,5)+Area_theta_L(5)*w(5,5)+Area_theta_L(6)*w(6,5)+Area_theta_L(7)*w(7,...
            5)+Area_theta_L(8)*w(8,5))+...
            Area_theta_L(6)*w(i,6)/(Area_theta_L(1)*w(1,6)+Area_theta_L(2)*w(2,6)+Area_theta_L(3)*w(3,6)...
            +Area_theta_L(4)*w(4,6)+Area_theta_L(5)*w(5,6)+Area_theta_L(6)*w(6,6)+Area_theta_L(7)*w(7,...
            6)+Area_theta_L(8)*w(8,6))+...
            Area_theta_L(7)*w(i,7)/(Area_theta_L(1)*w(1,7)+Area_theta_L(2)*w(2,7)+Area_theta_L(3)*w(3,7)...
            +Area_theta_L(4)*w(4,7)+Area_theta_L(5)*w(5,7)+Area_theta_L(6)*w(6,7)+Area_theta_L(7)*w(7,...
            7)+Area_theta_L(8)*w(8,7))+...
            Area_theta_L(8)*w(i,8)/(Area_theta_L(1)*w(1,8)+Area_theta_L(2)*w(2,8)+Area_theta_L(3)*w(3,8)...
            +Area_theta_L(4)*w(4,8)+Area_theta_L(5)*w(5,8)+Area_theta_L(6)*w(6,8)+Area_theta_L(7)*w(7,...
            8)+Area_theta_L(8)*w(8,8));
        ln_t_L(i) = q(i) * (1 - b_L(i) - c_L(i));
    end
    ln_gamma_residual_L(1)=1*( ln_t_L(8) - ln_t_p1_L(8) );
    ln_gamma_residual_L(2)=2*( ln_t_L(1) - ln_t_p2_L(1) )+1*( ln_t_L(3) - ln_t_p2_L(3) )+1*( ln_t_L(6)...
        - ln_t_p2_L(6) ) + 1*( ln_t_L(7) - ln_t_p2_L(7) );
    ln_gamma_residual_L(3)=5*( ln_t_L(1) - ln_t_p3_L(1) )+7*( ln_t_L(2) - ln_t_p3_L(2) )+2*( ln_t_L(3)...
        - ln_t_p3_L(3) ) + 4*( ln_t_L(4) - ln_t_p3_L(4) ) + 1*( ln_t_L(5) - ln_t_p3_L(5) ) + 2*( ln_t_L(6) -...
        ln_t_p3_L(6) );
    
    ln_gamma_total_L(1) = ln_gamma_combinatorial_L(1) + ln_gamma_residual_L(1);
    ln_gamma_total_L(2) = ln_gamma_combinatorial_L(2) + ln_gamma_residual_L(2);
    ln_gamma_total_L(3) = ln_gamma_combinatorial_L(3) + ln_gamma_residual_L(3);
    % Introduction of GE function
    
    GE_L=(x(1)*ln_gamma_total_L(1)+x(2)*ln_gamma_total_L(2)+x(3)*ln_gamma_total_L(3))*R*T;
    
    %Inserting the new mixing rule incorporating UNIFAC
    
    for i=1:3
        bratioL(i) = log(bL/bi(i));
        alphai(i) = ai(i)/bi(i);
    end
    
    aL1=(ulambda/Av + (1-ulambda)/Am)*GE_L/R/T;
    aL2=((1-ulambda)/Am)*(x*bratioL');
    aL3=x*alphai';
    jL= aL3*bL;
    wL= (aL1+aL2)*R*T*bL;
    aL= jL+wL;
    
    AL = aL*PB/(R^2*T^2);
    BL = bL*PB/R/T;
    
    % L(1)*X^N + ... + L(N)*X + L(N+1)
    % Peng Robinson coefficients of the cubic eq.
    L(1) = 1;
    L(2) = BL - 1;
    L(3) = AL - 3*BL^2 - 2*BL;
    L(4) = BL^3 + BL^2 - AL*BL;
    % Solve for the roots.
    ZtempL = roots(L);
    % I will consider only the largest value. The smallest value for liquid phase.
    j = 1;
    for i = 1:3
        if imag(ZtempL(i))==0
            zL(j) = ZtempL(i);
            j = j + 1;
        end
    end
    ZL = min(zL);
    for i = 1:3
        alpha_prime_L(i)=(ulambda/Av+(1-ulambda)/Am)*ln_gamma_total_L(i)+(1-...
            ulambda)/Am*(log(bL/bi(i))+bi(i)/bL-1) + ai(i)/bi(i)/R/T;
        phiL(i)= exp( bi(i)/bL*(ZL-1) - log(ZL-BL) - alpha_prime_L(i)/2.828*log((ZL+2.414*BL)/(ZL-0.414*BL)) );
    end
    g1new = phiL(1)/phiV(1); g2new = phiL(2)/phiV(2); g3new = phiL(3)/phiV(3)
    x2 = (1+(g1new-g3new)*x(3)-g1new)/(g2new-g1new)
end
