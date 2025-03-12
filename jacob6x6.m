function [J] = jacob6x6( x , E_v1 , E_v2 , E_p )
%
% Evaluates the Jacobian of a 6x6 system of non-linear equations
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
J(1,1) = 1 ;
%
J(1,2) = 1 ;
%
J(1,3) = 2 ;
%
J(1,4) = 0 ;
%
J(1,5) = 1 ;
%
J(1,6) = 2 ;
%
J(2,1) = 0 ;
%
J(2,2) = 0 ;
%
J(2,3) = -E_p ;
%
J(2,4) = 1 ;
%
J(2,5) = 0 ;
%
J(2,6) = 0 ;
%
EBar = E_MT + E_v1 ;
if ( abs(x(1))>abs(e_v1) && abs(x(1))<=abs(e_MT) )
    EBar = E_MT                                                   + E_v1 * ( 1 + L_v1 * ((abs(x(1))-abs(e_v1))^(m_v1-0)) ) ;
elseif ( abs(x(1))>abs(e_MT) )
    EBar = E_MT * ( 1 + L_MT * ((abs(x(1))-abs(e_MT))^(m_MT-0)) ) + E_v1 * ( 1 + L_v1 * ((abs(x(1))-abs(e_v1))^(m_v1-0)) ) ;
end
KBar = (EBar*Beta-1)  / (Beta-Alpha) ; 
J(3,1) = - KBar ;
if ( abs(x(1))>abs(e_v1) && abs(x(1))<=abs(e_MT) )
    J(3,1) = - KBar - x(1)*Beta*E_v1*L_v1*m_v1*((abs(x(1))-abs(e_v1))^(m_v1-1))/(Beta-Alpha) ;
elseif ( abs(x(1))>abs(e_MT) )
    J(3,1) = - KBar - x(1)*Beta*E_v1*L_v1*m_v1*((abs(x(1))-abs(e_v1))^(m_v1-1))/(Beta-Alpha) - x(1)*Beta*E_MT*L_MT*m_MT*((abs(x(1))-abs(e_MT))^(m_MT-1))/(Beta-Alpha) ;
end
%
J(3,2) = 0 ;
%
J(3,3) = 0 ;
%
J(3,4) = 1 ;
%
J(3,5) = 0 ;
%
J(3,6) = 0 ;
%
J(4,1) = 0 ;
%
E_av   = E_a + E_v2 ;
J(4,2) = - E_av ;
if ( (x(2)>e_a) && (x(2)<=e_v2) )
    E_av   = E_a * ( 1 + L_a * ((x(2)-e_a)^m_a) ) + E_v2 ;
    J(4,2) = - E_av - x(2)*E_a*L_a*m_a*((x(2)-e_a)^(m_a-1)) ;
elseif ( x(2)>e_v2 )
    E_av   = E_a * ( 1 + L_a * ((x(2)-e_a)^m_a) ) + E_v2 * ( 1 + L_v2 * ((x(2)-e_v2)^m_v2) );
    J(4,2) = - E_av - x(2)*E_a*L_a*m_a*((x(2)-e_a)^(m_a-1)) - x(2)*E_v2*L_v2*m_v2*((x(2)-e_v2)^(m_v2-1)) ;
end
%
J(4,3) = 0 ;
%
J(4,4) = 1 ;
%
J(4,5) = 0 ;
%
J(4,6) = 0 ;
%
J(5,1) = 0 ;
%
J(5,2) = 0 ;
%
J(5,3) = 0 ;
%
J(5,4) = 1 ;
%
E_NCL   = E_NC + E_NL ;
J(5,5) = - E_NCL ;
if ( (x(5)>e_NC) && (x(5)<=e_NL) )
    E_NCL   = E_NC * ( 1 + L_NC * ((x(5)-e_NC)^m_NC) ) + E_NL ;
    J(5,5) = - E_NCL - x(5)*E_NC*L_NC*m_NC*((x(5)-e_NC)^(m_NC-1)) ;
    % J(5,5) = - E_NCL ;
elseif ( x(5)>e_NL )
    E_NCL   = E_NC * ( 1 + L_NC * ((x(5)-e_NC)^m_NC) ) + E_NL * ( 1 + L_NL * ((x(5)-e_NL)^m_NL) );
    J(5,5) = - E_NCL - x(5)*E_NC*L_NC*m_NC*((x(5)-e_NC)^(m_NC-1)) - x(5)*E_NL*L_NL*m_NL*((x(5)-e_NL)^(m_NL-1)) ;
    % J(5,5) = - E_NCL ;
    % J(5,5)
end

%
J(6,1) = 0 ;
%
J(6,2) = 0 ;
%
J(6,3) = 0 ;
%
J(6,4) = 1 ;
%
J(6,5) = 0 ;
%
E_FAD   = E_FA ;
J(6,6) = - E_FAD ;
if (x(6)>e_FA)
    E_FAD   = E_FA * ( 1 + L_FA * ((x(6)-e_FA)^m_FA) ) ;
    J(6,6) = - E_FAD - x(6)*E_FA*L_FA*m_FA*((x(6)-e_FA)^(m_FA-1)) ;
    % J(5,5) = - E_NCL ;
end








