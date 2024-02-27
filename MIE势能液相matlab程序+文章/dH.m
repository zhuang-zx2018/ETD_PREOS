
function [beta,epslon,C,dHS]=dH(T,m,n,sigma,epslon_k)
kb=1.38E-23; 
beta=1/(kb*T);
epslon=epslon_k*kb;
a=0; b=sigma;

C=m/(m-n)*(m/n)^(n/(m-n));
         
    GauTen=[    
1	0.2955242247147529	-0.1488743389816312
2	0.2955242247147529	0.1488743389816312
3	0.2692667193099963	-0.4333953941292472
4	0.2692667193099963	0.4333953941292472
5	0.2190863625159820	-0.6794095682990244
6	0.2190863625159820	0.6794095682990244
7	0.1494513491505806	-0.8650633666889845
8	0.1494513491505806	0.8650633666889845
9	0.0666713443086881	-0.9739065285171717
10	0.0666713443086881	0.9739065285171717];      
             
    t=GauTen(:,3);
    r=(b+a)/2+(b-a)/2*t;
    u_mie=C*epslon*((sigma./r).^m-(sigma./r).^n);
    fx=1-exp(-beta*u_mie);
    dHS=(b-a)/2*dot(GauTen(:,2),fx);
end


