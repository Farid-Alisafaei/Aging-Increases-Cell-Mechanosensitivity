function [f] = f5( x , E_v1 , E_v2 , E_p )
%
% f4(x) = 0                  represents a system of 4 non-linear equations
% x = [x(1);x(2);x(3);x(4)]  unknown variables
% x(1) = e1                  strain in contractile element
% x(2) = e2                  strain in tensile element
% x(3) = ep                  strain in ECM
% x(4) = sigma               stress
% x(5) = encl                strain in nucleus element
% x(6) = efa                 strain in focal adhesion element
%
global E_MT
global E_a
global Beta
global Alpha
global Rho0
global e_v1
global e_v2
global e_MT
global e_a
global L_v1
global L_v2
global m_v1
global m_v2
global L_MT
global m_MT
global L_a
global m_a
global RhoBar0

global E_NC
global e_NC
global L_NC
global m_NC
global E_NL
global e_NL
global L_NL
global m_NL

global E_FA
global e_FA
global L_FA
global m_FA

%
f1 = x(1) + x(2) + 2*x(3) + x(5) + 2*x(6);
%
f2 = x(4) - E_p * x(3) ;
%
EBar = E_MT + E_v1 ;
if ( abs(x(1))>abs(e_v1) && abs(x(1))<=abs(e_MT) )
    EBar = E_MT                                                   + E_v1 * ( 1 + L_v1 * ((abs(x(1))-abs(e_v1))^(m_v1-0)) ) ;
elseif ( abs(x(1))>abs(e_MT) )
    EBar = E_MT * ( 1 + L_MT * ((abs(x(1))-abs(e_MT))^(m_MT-0)) ) + E_v1 * ( 1 + L_v1 * ((abs(x(1))-abs(e_v1))^(m_v1-0)) ) ;
end
KBar = (EBar*Beta-1)  / (Beta-Alpha) ;
f3 = x(4) - KBar * x(1) - RhoBar0 ;
%
E_av = E_a + E_v2 ;
if ( (x(2)>e_a) && (x(2)<=e_v2) )
    E_av = E_a * ( 1 + L_a * ((x(2)-e_a)^m_a) ) + E_v2 ;
elseif ( x(2)>e_v2 )
    E_av = E_a * ( 1 + L_a * ((x(2)-e_a)^m_a) ) + E_v2 * ( 1 + L_v2 * ((x(2)-e_v2)^m_v2) );
end
f4 = x(4) - E_av * x(2) ;
%

E_NCL = E_NC + E_NL ;
if ( (x(5)>e_NC) && (x(5)<=e_NL) )
    E_NCL = E_NC * ( 1 + L_NC * ((x(5)-e_NC)^m_NC) ) + E_NL ;
elseif ( x(5)>e_NL )
    E_NCL = E_NC * ( 1 + L_NC * ((x(5)-e_NC)^m_NC) ) + E_NL * ( 1 + L_NL * ((x(5)-e_NL)^m_NL) );
end
f5 = x(4) - E_NCL * x(5) ;

%

E_FAD = E_FA ;
if (x(6)>e_FA)
    E_FAD = E_FA * ( 1 + L_FA * ((x(6)-e_FA)^m_FA) ) ;
end
f6 = x(4) - E_FAD * x(6) ;

%
f = [f1;f2;f3;f4;f5;f6] ;







