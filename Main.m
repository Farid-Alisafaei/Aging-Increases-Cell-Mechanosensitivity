% Farid Alisafaei (New Jersey Institute of Technology)
% Manuscript: Aging increases mechanosensitivity of muscle progenitor cells to nanofibers
% 
clear all
clc
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
E_MT    = 9.000  ;         % stiffness of the microtubule element
E_a     = 0.500  ;         % initial stiffness of the actin element
Beta    = 1.200  ;         % chemical stiffness parameter
Alpha   = 0.800  ;         % chemo-mechanical feedback parameter
Rho0    = 10.00  ;         % motor density in the quiescent state

E_NC    = 5.000  ;         % initial stiffness of the chromatin element
E_NL    = 80.000 ;         % initial stiffness of the lamin element

E_FA    = 20.000 ;         % initial stiffness of the focal adhesion element

%
e_a     = 0.050  ;         % critical strain for strain-stiffening of the actin element
L_a     = 200.0  ;         % stiffening parameter for the actin element
m_a     = 15.00  ;         % stiffening parameter for the actin element
%
e_v2    = 0.000  ;         % critical strain for strain-stiffening of the intermediate filament element in tension
L_v2    = 0.000  ;         % stiffening parameter for the intermediate filament element in tension
m_v2    = 0.000  ;         % stiffening parameter for the intermediate filament element in tension
%
e_v1    = 0.000  ;         % critical strain for strain-stiffening of the intermediate filament element in compression
L_v1    = 0.000  ;         % stiffening parameter for the intermediate filament element in compression 
m_v1    = 0.000  ;         % stiffening parameter for the intermediate filament element in compression  
%
e_MT    = 0.000  ;         % critical strain for the microtubule element
L_MT    = 0.000  ;         % stiffening parameter for microtubule element
m_MT    = 0.000  ;         % stiffening parameter for microtubule element

e_NC    = 0.000  ;         % critical strain for the chromatin element
L_NC    = 0.000  ;         % stiffening parameter the chromatin element
m_NC    = 0.000  ;         % stiffening parameter the chromatin element

e_NL    = 0.000  ;         % critical strain for the lamin element
L_NL    = 0.000  ;         % stiffening parameter for the lamin element
m_NL    = 0.000  ;         % stiffening parameter for the lamin element

e_FA    = 0.020  ;         % critical strain for the focal adhesion element
L_FA    = 200.000;         % stiffening parameter for the focal adhesion element
m_FA    = 1.000  ;         % stiffening parameter for the focal adhesion element
%

% YAP translocation parameters
Kin0 = 2 ;
Kout = 1.5 ;
alphan = 0.8 ;
nn = 3 ;

%
RhoBar0 = Beta * Rho0 / (Beta-Alpha) ;
%
%
% check parameter
if ( E_MT*Alpha <= 1 )
    disp('Error')
    disp('E_MT*Alpha <= 1')
    return
end
if ( Beta <= Alpha )
    disp('Error')
    disp('Beta <= Alpha')
    return
end
%
% ECM stiffness
E_p = linspace(3,15,200) ;  
%
% stiffness of the intermediate filament element in compression
E_v1 = linspace(0,10,200) ;
%
% factor to indicate the difference between intermediate filament in compression and tension
f = 1.0 ;
%
% stiffness of the intermediate filament element in tension
E_v2 = E_v1 / f ;
%
% initial guess
x0 = [0;0;0;0;0;0] ;
%
% strain in the contractile element of the cytoskeleton
Epsilon_1    = zeros(length(E_v1),length(E_p)) ;
%
% strain in the tensile element of the cytoskeleton
Epsilon_2    = zeros(length(E_v1),length(E_p)) ;
%
% strain in the ECM element
Epsilon_p    = zeros(length(E_v1),length(E_p)) ;
%
% difference in the ECM strain
dEpsilon_p   = zeros(length(E_v1),length(E_p)) ;
%
% difference in total cell strain
dCellStrain  = zeros(length(E_v1),length(E_p)) ;
%
% stress
Sigma        = zeros(length(E_v1),length(E_p)) ;
%
% stress difference
dSigma       = zeros(length(E_v1),length(E_p)) ;
%
% contractility
Rho          = zeros(length(E_v1),length(E_p)) ;
%
% contractility difference
dRho         = zeros(length(E_v1),length(E_p)) ;
%
% stiffness of tensile element
E_av         = zeros(length(E_v1),length(E_p)) ;
%
% stiffness of contractile element
E_mv         = zeros(length(E_v1),length(E_p)) ;
%
% total stiffness 
E_tot        = zeros(length(E_v1),length(E_p)) ;
%
% total stiffness difference
dE_tot       = zeros(length(E_v1),length(E_p)) ;

% strain of nucleus element
Epsilon_3    = zeros(length(E_v1),length(E_p)) ;

% difference in strain of nucleus element
dEpsilon_3    = zeros(length(E_v1),length(E_p)) ;

% stiffness of nucleus element
E_NCL    = zeros(length(E_v1),length(E_p)) ;

% strain of focal adhesion element
Epsilon_4    = zeros(length(E_v1),length(E_p)) ;

% difference in strain of focal adhesion element
dEpsilon_4    = zeros(length(E_v1),length(E_p)) ;

% stiffness of nucleus element
E_FAD    = zeros(length(E_v1),length(E_p)) ;

%
% solving equations
for i = 1:length(E_v1)
    for j = 1:length(E_p)
        [x,iter] = newtonm( x0 , E_v1(i) , E_v2(i) , E_p(j) ) ;
        Epsilon_1(i,j)  = x(1) ;
        Epsilon_2(i,j)  = x(2) ;
        Epsilon_p(i,j)  = x(3) ;
        Sigma(i,j)      = x(4) ;
        Epsilon_3(i,j)  = x(5) ;
        Epsilon_4(i,j)  = x(6) ;
    end
end
%
% contractility and stiffness of contractile component
for i = 1:length(E_v1)
    for j = 1:length(E_p)
        EBar = E_MT + E_v1(i) ;
        if ( abs(Epsilon_1(i,j))>abs(e_v1) && abs(Epsilon_1(i,j))<=abs(e_MT) )
            EBar = E_MT                                                             + E_v1(i) * ( 1 + L_v1 * ((abs(Epsilon_1(i,j))-abs(e_v1))^(m_v1-0)) ) ;
        elseif ( abs(Epsilon_1(i,j))>abs(e_MT) )
            EBar = E_MT * ( 1 + L_MT * ((abs(Epsilon_1(i,j))-abs(e_MT))^(m_MT-0)) ) + E_v1(i) * ( 1 + L_v1 * ((abs(Epsilon_1(i,j))-abs(e_v1))^(m_v1-0)) ) ;
        end
        KpBar    = (EBar*Alpha-1)  / (Beta-Alpha)   ;
        KBar     = (EBar*Beta-1)   / (Beta-Alpha)   ;
        Rho(i,j) = KpBar * Epsilon_1(i,j) + RhoBar0 ;
        E_mv(i,j)= KBar                             ;
    end
end
%
% stiffness of tensile component
for i = 1:length(E_v1)
    for j = 1:length(E_p)
        E_av(i,j)       = E_a + E_v2(i) ;
        if ( (Epsilon_2(i,j)>e_a) && (Epsilon_2(i,j)<=e_v2) )
            E_av(i,j)   = E_a * ( 1 + L_a * ((Epsilon_2(i,j)-e_a)^m_a) ) + E_v2(i) ;
        elseif ( Epsilon_2(i,j)>e_v2 )
            E_av(i,j)   = E_a * ( 1 + L_a * ((Epsilon_2(i,j)-e_a)^m_a) ) + E_v2(i) * ( 1 + L_v2 * ((Epsilon_2(i,j)-e_v2)^m_v2) ) ;
        end
    end
end      

% stiffness of nucleus component
for i = 1:length(E_v1)
    for j = 1:length(E_p)
        E_NCL(i,j)       = E_NC + E_NL ;
        if ( (Epsilon_3(i,j)>e_NC) && (Epsilon_3(i,j)<=e_NL) )
            E_NCL(i,j)   = E_NC * ( 1 + L_NC * ((Epsilon_3(i,j)-e_NC)^m_NC) ) + E_NL ;
        elseif ( Epsilon_3(i,j)>e_NL )
            E_NCL(i,j)   = E_NC * ( 1 + L_NC * ((Epsilon_3(i,j)-e_NC)^m_NC) ) + E_NL * ( 1 + L_NL * ((Epsilon_3(i,j)-e_NL)^m_NL) ) ;
        end
    end
end      

% stiffness of focal adhesion component
for i = 1:length(E_v1)
    for j = 1:length(E_p)
        E_FAD(i,j)       = E_FA ; % both FA elements have same stiffness and are connected in series
        if (Epsilon_4(i,j)>e_FA)
            E_FAD(i,j)   = E_FA * ( 1 + L_FA * ((Epsilon_4(i,j)-e_FA)^m_FA) ) ;
        end
    end
end      

%
% total stiffness 
for i = 1:length(E_v1)
    for j = 1:length(E_p)
        E_tot(i,j) = E_av(i,j) * E_mv(i,j) * E_NCL(i,j) * E_FAD(i,j)/ ( ( E_av(i,j) * E_mv(i,j) * E_FAD(i,j)) + ( E_av(i,j) * E_NCL(i,j) * E_FAD(i,j)) + ( E_NCL(i,j) * E_mv(i,j) * E_FAD(i,j)) + (E_av(i,j) * E_mv(i,j) * E_NCL(i,j))) ;
    end
end
%
%--------------------------------------------------------------------------
%
a = 1   ;
c = 30  ;
disp([E_v1(a);E_v1(c)])
%

% norm stress
figure
plot( E_p/E_p(1) , Sigma(a,:)/Sigma(a,1)   , 'r-' )
hold on
plot( E_p/E_p(1) , Sigma(c,:)/Sigma(a,1)  , 'b-' )
xlabel('Ep')
ylabel('norm stress')
% %

% stress
figure
plot( E_p/E_p(1) , Sigma(a,:)   , 'r-' )
hold on
plot( E_p/E_p(1) , Sigma(c,:)  , 'b-' )
xlabel('Ep')
ylabel('stress')
% %

% strain in nucleus element
figure
plot( E_p/E_p(1) , Epsilon_3(a,:)/Epsilon_3(a,1)   , 'r-' )
% plot( E_p , Epsilon_3(a,:)   , 'b-' )
hold on
plot( E_p/E_p(1) , Epsilon_3(c,:)/Epsilon_3(a,1)   , 'b-' )
% plot( E_p , Epsilon_3(c,:)   , 'r-' )
xlabel('Normalized Ep')
ylabel('Normalized strain in nucleus element')
%

% focal adhesion formation
figure
plot( E_p/E_p(1) , E_FAD(a,:)/E_FAD(a,1)   , 'r-' )
% plot( E_p , E_FAD(a,:)   , 'b-' )
hold on
plot( E_p/E_p(1) , E_FAD(c,:)/E_FAD(a,1)   , 'b-' )
% plot( E_p , E_FAD(c,:)   , 'r-' )
xlabel('Normalized Ep')
ylabel('Normalized focal adhesion formation')
%

% YAP accumulation calculations
% calculate Kin and R
Kin = Kin0 + ((Epsilon_3 ./ Epsilon_3(a,1)).^nn) ./ (1 + alphan * ((Epsilon_3 ./ Epsilon_3(a,1)).^nn));
R = Kin / Kout;

% plot R vs E_p
figure;
plot(E_p / E_p(1), R(a,:), 'b-', 'LineWidth', 1.5);
hold on;
plot(E_p / E_p(1), R(c,:), 'r-', 'LineWidth', 1.5);
xlabel('Normalized Ep');
ylabel('Ratio');
grid on;















