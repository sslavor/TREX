%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   ████████╗██████╗ ███████╗██╗  ██╗
%   ╚══██╔══╝██╔══██╗██╔════╝╚██╗██╔╝
%      ██║   ██████╔╝█████╗   ╚███╔╝ 
%      ██║   ██╔══██╗██╔══╝   ██╔██╗ 
%      ██║   ██║  ██║███████╗██╔╝ ██╗
%      ╚═╝   ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝
%
%   TREX  –  Turbine Rotor EXtractor
%   Supersonic Impulse Blade Profile Generator

%   Version: 1.0
%   Description:
%   Parametric generator for supersonic impulse turbine blade geometry.

%   References:
%   1. Analytical Investigation of supersonic turbomachinery blading , Louis J. Goldman , Vincent J. Scullin
%   2. Design and preliminary optimization of a supersonic turbine for Rotating Detonation Engine , Noraiz MUSHTAQ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUTS 
function TREX_AUTO(nl, Beta_i_UI, M_i_UI, M_L_UI, M_U_UI, T_UI, R_UI, Y_UI)

clc ;
nlines = nl ;

Beta_i = Beta_i_UI ; % deg
Beta_o = -Beta_i ;
M_i = M_i_UI ; % Inlet Mach number
M_o = M_i ; % Outlet Mach number
M_L = M_L_UI ; % Lower surface Mach number (PRESSURE SURFACE)
M_U = M_U_UI ; % Upper surface Mach number (SUCTION SURFACE)

T = T_UI ; % K
R = R_UI ;
Y = Y_UI ;

output = fopen('M-PM.txt','w') ;

% For interpolation of Prandtl-Meyer function and Mach number
fprintf(output,'Ae/At\t\tM\t\tv\n') ;

for i = 1 : 0.0001 : 5
    v_x = PM(i,Y) ;
    fprintf(output,'%f\t%f\n',i,v_x) ;
end

data = readtable('M-PM.txt') ;
M_x = data{:,1} ;
v_x = data{:,2} ;


%% SOLVER
Vcr = sqrt(Y*R*T) ;

M_star_L = M_star(M_L,Y) ;
R_star_L = 1/M_star_L ;

M_star_U = M_star(M_U,Y) ;
R_star_U = 1/M_star_U ;

%fprintf('%f\t%f\n',R_star_L,R_star_U) ;

% LOWER SURFACE (PRESSURE SURFACE)
dtheta_trans_p_i = dtheta_trans(M_i,M_L,Y,true) ;
dtheta_cir_p_i = dtheta_cir(Beta_i,dtheta_trans_p_i,true) ;

dtheta_trans_p_o = dtheta_trans(M_o,M_L,Y,true) ;
dtheta_cir_p_o = dtheta_cir(Beta_o,dtheta_trans_p_o,false) ;

% UPPER SURFACE (SUCTION SURFACE)
dtheta_trans_s_i = dtheta_trans(M_i,M_U,Y,false) ;
dtheta_cir_s_i = dtheta_cir(Beta_i,dtheta_trans_s_i,true) ;

dtheta_trans_s_o = dtheta_trans(M_o,M_U,Y,false) ;
dtheta_cir_s_o = dtheta_cir(Beta_o,dtheta_trans_s_o,false) ;

dv_p = (1/nl) * dtheta_trans_p_i ;
dv_s = (1/nl) * dtheta_trans_s_i ;
v_i = PM(M_i,Y) ;
v_L = PM(M_L,Y) ;
v_U = PM(M_U,Y) ;
kmax = (v_i - v_L)/dv_p ;
jmax = (v_U - v_i)/dv_s ;

% Circular arc points
x_L = -R_star_L * sind(dtheta_cir_p_i) ;
y_L = R_star_L * cosd(dtheta_cir_p_i) ;
x_U = -R_star_U * sind(dtheta_cir_s_i) ;
y_U = R_star_U * cosd(dtheta_cir_s_i) ;

% k -> index for lower surface of the blade
% j -> index for upper surface of the blade
% i -> blade inlet
% o -> blade outlet

gamma = Y ;
syms xx yy RR

%% LOWER SURFACE (PRESSURE SURFACE OR CONCAVE SURFACE)

v_p = [] ;
v_W_p = [] ;
v_I_p = [] ;
theta_p = [] ;
theta_W_p = [] ;
theta_I_p = [] ;
M_p = [] ;
M_W_p = [] ;
M_I_p = [] ;
u_p = [] ;
u_W_p = [] ;
u_I_p = [] ;
X_p = [] ;
X_W_p = [] ;
X_I_p = [] ;
Y_p = [] ;
Y_W_p = [] ;
Y_I_p = [] ;
x_W_p = [] ;
y_W_p = [] ;

dv = dv_p ;
a = 1 ;
b = 1 ;
c = 1 ;
for k = kmax+1:-1:1
    p_or_s = true ;
    % STARTING POINT
    if k == kmax+1
        v_p(a) = PM(M_L,gamma) ;
        theta_p(a) = 0 ;   % Flow angle w.r.t x-axis
        M_p(a) = M_L ;
        u_p(a) = ux(M_p(a)) ; % Corresponding Mach angle
        x = -R_star_L * sind(theta_p(a)) ;
        y = R_star_L * cosd(theta_p(a)) ;
        X_p(a) = x*cosd(dtheta_cir_p_i) - y*sind(dtheta_cir_p_i) ;
        Y_p(a) = x*sind(dtheta_cir_p_i) + y*cosd(dtheta_cir_p_i) ;

        theta_W_p(c) = theta_p(a) ;
        x_W_p(c) = x ;
        y_W_p(c) = y ;
        X_W_p(c) = X_p(a) ;
        Y_W_p(c) = Y_p(a) ;
        
        c = c + 1 ;
        a = a + 1 ;

    else
        % j = 1 -> On the Major-Vortex Expansion Characteristic
        % j = 2 -> On the Wall
        for j = 1:1:2

            % INNER POINTS
            if j == 1
                theta_p(a) = v_i - v_L - (k - 1)*dv ;
                v_p(a) = theta_p(a) + v_p(a-1) - theta_p(a-1) ;
                M_p(a) = Mx(v_x,M_x,v_p(a)) ;
                u_p(a) = ux(M_p(a)) ; % Corresponding Mach angle
                
                R_star = solve_R(v_i, gamma, dv , k ,p_or_s) ;
                x1 = -R_star * sind(theta_p(a)) ;
                y1 = R_star * cosd(theta_p(a)) ;
                X_p(a) = x1*cosd(dtheta_cir_p_i) - y1*sind(dtheta_cir_p_i) ;
                Y_p(a) = x1*sind(dtheta_cir_p_i) + y1*cosd(dtheta_cir_p_i) ;
                
                theta_I_p(b) = theta_p(a) ;
                v_I_p(b) = v_p(a) ;
                M_I_p(b) = M_p(a) ;
                u_I_p(b) = u_p(a) ;
                X_I_p(b) = X_p(a) ;
                Y_I_p(b) = Y_p(a) ;

                b = b + 1 ;
                a = a + 1;

            % WALL POINTS
            else
                theta_p(a) = theta_p(a-1) ;
                v_p(a) = v_p(a-1) ; 
                M_p(a) = Mx(v_x,M_x,v_p(a)) ;
                u_p(a) = ux(M_p(a)) ; % Corresponding Mach angle
                
                m_minus_1 = m_minus(theta_p(a-1),theta_p(a),u_p(a-1),u_p(a)) ;
                tm_minus_1 = t(m_minus_1) ;
                m_plus_1 = theta_W_p(c-1) ;
                tm_plus_1 = t(m_plus_1) ;
                
                eqn1 = yy - y1 == tm_minus_1*(xx-x1) ;
                eqn2 = yy - y_W_p(c-1) == tm_plus_1*(xx-x_W_p(c-1)) ;
                s = solve([eqn1,eqn2],[xx,yy]) ;
                x2 = double(s.xx) ;%s.xx ;
                y2 = double(s.yy) ;%s.yy ;
                X_p(a) = x2*cosd(dtheta_cir_p_i) - y2*sind(dtheta_cir_p_i) ;
                Y_p(a) = x2*sind(dtheta_cir_p_i) + y2*cosd(dtheta_cir_p_i) ;

                theta_W_p(c) = theta_p(a) ;
                v_W_p(c) = v_p(a) ;
                M_W_p(c) = M_p(a) ;
                u_W_p(c) = u_p(a) ;
                x_W_p(c) = x2 ;
                y_W_p(c) = y2 ;
                X_W_p(c) = X_p(a) ;
                Y_W_p(c) = Y_p(a) ;

                c = c + 1 ;
                a = a + 1;

            end
        end
    end
end

%% UPPER SURFACE (SUCTION SURFACE OR CONVEX SURFACE)

v_s = [] ;
v_W_s = [] ;
v_I_s = [] ;
theta_s = [] ;
theta_W_s = [] ;
theta_I_s = [] ;
M_s = [] ;
M_W_s = [] ;
M_I_s = [] ;
u_s = [] ;
u_W_s = [] ;
u_I_s = [] ;
X_s = [] ;
X_W_s = [] ;
X_I_s = [] ;
Y_s = [] ;
Y_W_s = [] ;
Y_I_s = [] ;
x_W_s = [] ;
y_W_s = [] ;

dv = dv_s ;
a = 1 ;
b = 1 ;
c = 1 ;
for j = jmax+1:-1:1
    p_or_s = false ;
    % STARTING POINT
    if j == jmax+1
        v_s(a) = PM(M_U,gamma) ;
        theta_s(a) = 0 ;   % Flow angle w.r.t x-axis
        M_s(a) = M_U ;
        u_s(a) = ux(M_s(a)) ; % Corresponding Mach angle
        x = -R_star_U * sind(theta_s(a)) ;
        y = R_star_U * cosd(theta_s(a)) ;
        X_s(a) = x*cosd(dtheta_cir_s_i) - y*sind(dtheta_cir_s_i) ;
        Y_s(a) = x*sind(dtheta_cir_s_i) + y*cosd(dtheta_cir_s_i) ;

        theta_W_s(c) = theta_s(a) ;
        x_W_s(c) = x ;
        y_W_s(c) = y ;
        X_W_s(c) = X_s(a) ;
        Y_W_s(c) = Y_s(a) ;
        
        c = c + 1 ;
        a = a + 1 ;

    else
        % i = 1 -> On the Major-Vortex Compression Characteristic
        % i = 2 -> On the Wall
        for i = 1:1:2

            % INNER POINTS
            if i == 1
                theta_s(a) = v_U - v_i - (j - 1)*dv ;
                v_s(a) = v_s(a-1) + theta_s(a-1) - theta_s(a) ;
                M_s(a) = Mx(v_x,M_x,v_s(a)) ;
                u_s(a) = ux(M_s(a)) ; % Corresponding Mach angle
                
                R_star = solve_R(v_i, gamma, dv , j ,p_or_s) ;
                x1 = -R_star * sind(theta_s(a)) ;
                y1 = R_star * cosd(theta_s(a)) ;
                X_s(a) = x1*cosd(dtheta_cir_s_i) - y1*sind(dtheta_cir_s_i) ;
                Y_s(a) = x1*sind(dtheta_cir_s_i) + y1*cosd(dtheta_cir_s_i) ;
                
                theta_I_s(b) = theta_s(a) ;
                v_I_s(b) = v_s(a) ;
                M_I_s(b) = M_s(a) ;
                u_I_s(b) = u_s(a) ;
                X_I_s(b) = X_s(a) ;
                Y_I_s(b) = Y_s(a) ;

                b = b + 1 ;
                a = a + 1;

            % WALL POINTS
            else
                theta_s(a) = theta_s(a-1) ;
                v_s(a) = v_s(a-1) ; 
                M_s(a) = Mx(v_x,M_x,v_s(a)) ;
                u_s(a) = ux(M_s(a)) ; % Corresponding Mach angle
                
                m_plus_1_d = m_plus(theta_s(a-1),theta_s(a),u_s(a-1),u_s(a)) ;
                tm_plus_1_d = t(m_plus_1_d) ;
                m_plus_1 = theta_W_s(c-1) ;
                tm_plus_1 = t(m_plus_1) ;
                
                eqn1 = yy - y1 == tm_plus_1_d*(xx-x1) ;
                eqn2 = yy - y_W_s(c-1) == tm_plus_1*(xx-x_W_s(c-1)) ;
                s = solve([eqn1,eqn2],[xx,yy]) ;
                x2 = double(s.xx) ;%s.xx ;
                y2 = double(s.yy) ;%s.yy ;
                X_s(a) = x2*cosd(dtheta_cir_s_i) - y2*sind(dtheta_cir_s_i) ;
                Y_s(a) = x2*sind(dtheta_cir_s_i) + y2*cosd(dtheta_cir_s_i) ;

                theta_W_s(c) = theta_s(a) ;
                v_W_s(c) = v_s(a) ;
                M_W_s(c) = M_s(a) ;
                u_W_s(c) = u_s(a) ;
                x_W_s(c) = x2 ;
                y_W_s(c) = y2 ;
                X_W_s(c) = X_s(a) ;
                Y_W_s(c) = Y_s(a) ;

                c = c + 1 ;
                a = a + 1;

            end
        end
    end
    % Straight line
    if j == 1
        X1 = X_p(end) ;
        Y1 = Y_p(end) ;
        
        %m = t(theta_W_s(end)) ;
        m = t(Beta_i) ;
        X2 = X_s(end) ;
        Y2 = Y1 + m*(X2 - X1) ;
        Y_s_end = Y_s(end) ;
        d_Y_s = Y2 - Y_s_end ;

    end
end

%% PLOTTING

%[G,G_new] = G_star(M_i,Y,v_L,v_U,R_star_L,R_star_U,2*Beta_i,M_star_L,M_star_U) ;

% Circular Arc
figure(1) ;
hold on 
[X_p_cir_1,Y_p_cir_1] = plotArc(R_star_L,dtheta_cir_p_i,0,nl) ;
[X_s_cir_1,Y_s_cir_1] = plotArc(R_star_U ,dtheta_cir_s_i,d_Y_s,nl) ;
[X_p_cir_2,Y_p_cir_2] = plotArc(R_star_L,dtheta_cir_p_i,-d_Y_s,nl) ;
[X_s_cir_2,Y_s_cir_2] = plotArc(R_star_U ,dtheta_cir_s_i,0,nl) ;
title(nlines)

% Points
%{
plot(X_p, Y_p, ".", "Color", "magenta", "MarkerSize", 15);
plot(-X_p, Y_p, ".", "Color", "magenta", "MarkerSize", 15);
plot(X_s, Y_s + d_Y_s, ".", "Color", "green", "MarkerSize", 15);
plot(-X_s, Y_s + d_Y_s, ".", "Color", "green", "MarkerSize", 15);
%}

% Wall
for i = 1:1:nl
    line([X_W_p(i),X_W_p(i+1)],[Y_W_p(i),Y_W_p(i+1)],Linewidth=2,color='red') ;
    line([-X_W_p(i),-X_W_p(i+1)],[Y_W_p(i),Y_W_p(i+1)],Linewidth=2,color='red') ;

    line([X_W_p(i),X_W_p(i+1)],[Y_W_p(i) - d_Y_s,Y_W_p(i+1) - d_Y_s],Linewidth=2,color='red') ;
    line([-X_W_p(i),-X_W_p(i+1)],[Y_W_p(i) - d_Y_s,Y_W_p(i+1) - d_Y_s],Linewidth=2,color='red') ;


    line([X_W_s(i),X_W_s(i+1)],[Y_W_s(i) + d_Y_s,Y_W_s(i+1) + d_Y_s],Linewidth=2,color='red') ;
    line([-X_W_s(i),-X_W_s(i+1)],[Y_W_s(i) + d_Y_s,Y_W_s(i+1) + d_Y_s],Linewidth=2,color='red') ;

    line([X_W_s(i),X_W_s(i+1)],[Y_W_s(i),Y_W_s(i+1)],Linewidth=2,color='red') ;
    line([-X_W_s(i),-X_W_s(i+1)],[Y_W_s(i),Y_W_s(i+1)],Linewidth=2,color='red') ;
end
%{
% Characteristic lines
for i = 1:1:nl+1
    if i == 1
        line([X_W_p(i),X_I_p(i)],[Y_W_p(i),Y_I_p(i)],Linewidth=1,color='black') ;
        line([-X_W_p(i),-X_I_p(i)],[Y_W_p(i),Y_I_p(i)],Linewidth=1,color='black') ;

        line([X_W_p(i),X_I_p(i)],[Y_W_p(i) - d_Y_s,Y_I_p(i) - d_Y_s],Linewidth=1,color='black') ;
        line([-X_W_p(i),-X_I_p(i)],[Y_W_p(i) - d_Y_s,Y_I_p(i) - d_Y_s],Linewidth=1,color='black') ;


        line([X_W_s(i),X_I_s(i)],[Y_W_s(i) + d_Y_s,Y_I_s(i) + d_Y_s],Linewidth=1,color='black') ;
        line([-X_W_s(i),-X_I_s(i)],[Y_W_s(i) + d_Y_s,Y_I_s(i) + d_Y_s],Linewidth=1,color='black') ;

        line([X_W_s(i),X_I_s(i)],[Y_W_s(i),Y_I_s(i)],Linewidth=1,color='black') ;
        line([-X_W_s(i),-X_I_s(i)],[Y_W_s(i),Y_I_s(i)],Linewidth=1,color='black') ;
    else

        % Wall to Internal Points
        line([X_W_p(i),X_I_p(i-1)],[Y_W_p(i),Y_I_p(i-1)],Linewidth=1,color='black') ;
        line([-X_W_p(i),-X_I_p(i-1)],[Y_W_p(i),Y_I_p(i-1)],Linewidth=1,color='black') ;

        line([X_W_p(i),X_I_p(i-1)],[Y_W_p(i) - d_Y_s,Y_I_p(i-1) - d_Y_s],Linewidth=1,color='black') ;
        line([-X_W_p(i),-X_I_p(i-1)],[Y_W_p(i) - d_Y_s,Y_I_p(i-1) - d_Y_s],Linewidth=1,color='black') ;


        line([X_W_s(i),X_I_s(i-1)],[Y_W_s(i) + d_Y_s,Y_I_s(i-1) + d_Y_s],Linewidth=1,color='black') ;
        line([-X_W_s(i),-X_I_s(i-1)],[Y_W_s(i) + d_Y_s,Y_I_s(i-1) + d_Y_s],Linewidth=1,color='black') ;

        line([X_W_s(i),X_I_s(i-1)],[Y_W_s(i),Y_I_s(i-1)],Linewidth=1,color='black') ;
        line([-X_W_s(i),-X_I_s(i-1)],[Y_W_s(i),Y_I_s(i-1)],Linewidth=1,color='black') ;
    end

    % Last Characteristic line
    if i < nl+1 && i > 1
        line([X_I_p(i),X_I_p(i-1)],[Y_I_p(i),Y_I_p(i-1)],Linewidth=1,color='black') ;
        line([-X_I_p(i),-X_I_p(i-1)],[Y_I_p(i),Y_I_p(i-1)],Linewidth=1,color='black') ;

        line([X_I_p(i),X_I_p(i-1)],[Y_I_p(i) - d_Y_s,Y_I_p(i-1) - d_Y_s],Linewidth=1,color='black') ;
        line([-X_I_p(i),-X_I_p(i-1)],[Y_I_p(i) - d_Y_s,Y_I_p(i-1) - d_Y_s],Linewidth=1,color='black') ;


        line([X_I_s(i),X_I_s(i-1)],[Y_I_s(i) + d_Y_s,Y_I_s(i-1) + d_Y_s],Linewidth=1,color='black') ;
        line([-X_I_s(i),-X_I_s(i-1)],[Y_I_s(i) + d_Y_s,Y_I_s(i-1) + d_Y_s],Linewidth=1,color='black') ;

        line([X_I_s(i),X_I_s(i-1)],[Y_I_s(i),Y_I_s(i-1)],Linewidth=1,color='black') ;
        line([-X_I_s(i),-X_I_s(i-1)],[Y_I_s(i),Y_I_s(i-1)],Linewidth=1,color='black') ;
    end
end
%}

% Straight line from Upper arc to Lower arc
line([X_p(end),X_s(end)],[Y_p(end),Y_s(end) + d_Y_s],Linewidth=2,color='black') ;
line([-X_p(end),-X_s(end)],[Y_p(end),Y_s(end) + d_Y_s],Linewidth=2,color='black') ;

line([X_p(end),X_s(end)],[Y_p(end) - d_Y_s,Y_s(end)],Linewidth=2,color='black') ;
line([-X_p(end),-X_s(end)],[Y_p(end) - d_Y_s,Y_s(end)],Linewidth=2,color='black') ;

%% ONE CONTOUR

figure(2) ;
hold on 
[X_p_cir,Y_p_cir] = plotArc(R_star_L,dtheta_cir_p_i,0,nl) ;
[X_s_cir,Y_s_cir] = plotArc(R_star_U ,dtheta_cir_s_i,d_Y_s,nl) ;
title(nlines)

% Points
%{
plot(X_p, Y_p, ".", "Color", "magenta", "MarkerSize", 15);
plot(-X_p, Y_p, ".", "Color", "magenta", "MarkerSize", 15);
plot(X_s, Y_s + d_Y_s, ".", "Color", "green", "MarkerSize", 15);
plot(-X_s, Y_s + d_Y_s, ".", "Color", "green", "MarkerSize", 15);
%}
% Wall
for i = 1:1:nl
    line([X_W_p(i),X_W_p(i+1)],[Y_W_p(i),Y_W_p(i+1)],Linewidth=2,color='red') ;
    line([-X_W_p(i),-X_W_p(i+1)],[Y_W_p(i),Y_W_p(i+1)],Linewidth=2,color='red') ;

    line([X_W_s(i),X_W_s(i+1)],[Y_W_s(i) + d_Y_s,Y_W_s(i+1) + d_Y_s],Linewidth=2,color='red') ;
    line([-X_W_s(i),-X_W_s(i+1)],[Y_W_s(i) + d_Y_s,Y_W_s(i+1) + d_Y_s],Linewidth=2,color='red') ;
end

% Characteristic lines
for i = 1:1:nl+1
    if i == 1
        line([X_W_p(i),X_I_p(i)],[Y_W_p(i),Y_I_p(i)],Linewidth=1,color='black') ;
        line([-X_W_p(i),-X_I_p(i)],[Y_W_p(i),Y_I_p(i)],Linewidth=1,color='black') ;

        line([X_W_s(i),X_I_s(i)],[Y_W_s(i) + d_Y_s,Y_I_s(i) + d_Y_s],Linewidth=1,color='black') ;
        line([-X_W_s(i),-X_I_s(i)],[Y_W_s(i) + d_Y_s,Y_I_s(i) + d_Y_s],Linewidth=1,color='black') ;
    else

        % Wall to Internal Points
        line([X_W_p(i),X_I_p(i-1)],[Y_W_p(i),Y_I_p(i-1)],Linewidth=1,color='black') ;
        line([-X_W_p(i),-X_I_p(i-1)],[Y_W_p(i),Y_I_p(i-1)],Linewidth=1,color='black') ;

        line([X_W_s(i),X_I_s(i-1)],[Y_W_s(i) + d_Y_s,Y_I_s(i-1) + d_Y_s],Linewidth=1,color='black') ;
        line([-X_W_s(i),-X_I_s(i-1)],[Y_W_s(i) + d_Y_s,Y_I_s(i-1) + d_Y_s],Linewidth=1,color='black') ;
    end

    % Last Characteristic line
    if i < nl+1 && i > 1
        line([X_I_p(i),X_I_p(i-1)],[Y_I_p(i),Y_I_p(i-1)],Linewidth=1,color='black') ;
        line([-X_I_p(i),-X_I_p(i-1)],[Y_I_p(i),Y_I_p(i-1)],Linewidth=1,color='black') ;

        line([X_I_s(i),X_I_s(i-1)],[Y_I_s(i) + d_Y_s,Y_I_s(i-1) + d_Y_s],Linewidth=1,color='black') ;
        line([-X_I_s(i),-X_I_s(i-1)],[Y_I_s(i) + d_Y_s,Y_I_s(i-1) + d_Y_s],Linewidth=1,color='black') ;
    end
end

% Straight line from Upper arc to Lower arc
line([X_p(end),X_s(end)],[Y_p(end),Y_s(end) + d_Y_s],Linewidth=2,color='black') ;
line([-X_p(end),-X_s(end)],[Y_p(end),Y_s(end) + d_Y_s],Linewidth=2,color='black') ;
%{
% Wall
for i = 1:1:nl
    line([X_W_p(i),X_W_p(i+1)],[Y_W_p(i),Y_W_p(i+1)],Linewidth=2,color='red') ;
    line([-X_W_p(i),-X_W_p(i+1)],[Y_W_p(i),Y_W_p(i+1)],Linewidth=2,color='red') ;

    line([X_W_s(i),X_W_s(i+1)],[Y_W_s(i) + d_Y_s,Y_W_s(i+1) + d_Y_s],Linewidth=2,color='red') ;
    line([-X_W_s(i),-X_W_s(i+1)],[Y_W_s(i) + d_Y_s,Y_W_s(i+1) + d_Y_s],Linewidth=2,color='red') ;
end

% Characteristic lines
for i = 1:1:nl+1
    if i == 1
        line([X_W_p(i),X_I_p(i)],[Y_W_p(i),Y_I_p(i)],Linewidth=1,color='black') ;
        line([-X_W_p(i),-X_I_p(i)],[Y_W_p(i),Y_I_p(i)],Linewidth=1,color='black') ;

        line([X_W_s(i),X_I_s(i)],[Y_W_s(i) + d_Y_s,Y_I_s(i) + d_Y_s],Linewidth=1,color='black') ;
        line([-X_W_s(i),-X_I_s(i)],[Y_W_s(i) + d_Y_s,Y_I_s(i) + d_Y_s],Linewidth=1,color='black') ;
    else

        % Wall to Internal Points
        line([X_W_p(i),X_I_p(i-1)],[Y_W_p(i),Y_I_p(i-1)],Linewidth=1,color='black') ;
        line([-X_W_p(i),-X_I_p(i-1)],[Y_W_p(i),Y_I_p(i-1)],Linewidth=1,color='black') ;

        line([X_W_s(i),X_I_s(i-1)],[Y_W_s(i) + d_Y_s,Y_I_s(i-1) + d_Y_s],Linewidth=1,color='black') ;
        line([-X_W_s(i),-X_I_s(i-1)],[Y_W_s(i) + d_Y_s,Y_I_s(i-1) + d_Y_s],Linewidth=1,color='black') ;
    end

    % Last Characteristic line
    if i < nl+1 && i > 1
        line([X_I_p(i),X_I_p(i-1)],[Y_I_p(i),Y_I_p(i-1)],Linewidth=1,color='black') ;
        line([-X_I_p(i),-X_I_p(i-1)],[Y_I_p(i),Y_I_p(i-1)],Linewidth=1,color='black') ;

        line([X_I_s(i),X_I_s(i-1)],[Y_I_s(i) + d_Y_s,Y_I_s(i-1) + d_Y_s],Linewidth=1,color='black') ;
        line([-X_I_s(i),-X_I_s(i-1)],[Y_I_s(i) + d_Y_s,Y_I_s(i-1) + d_Y_s],Linewidth=1,color='black') ;
    end
end

% Straight line from Upper arc to Lower arc
line([X_p(end),X_s(end)],[Y_p(end),Y_s(end) + d_Y_s],Linewidth=2,color='black') ;
line([-X_p(end),-X_s(end)],[Y_p(end),Y_s(end) + d_Y_s],Linewidth=2,color='black') ;
%}

hold off

%% WRITING TO AN EXCEL SHEET

X_pressure_1 = [flip(X_W_p(2:end-1))+(-X_W_p(end)) , X_p_cir_1(2:end-1)+(-X_W_p(end)) , -(X_W_p(2:end-1))+(-X_W_p(end))] ;
Y_pressure_1 = [flip(Y_W_p(2:end-1))+(-Y_W_p(end)) , Y_p_cir_1(2:end-1)+(-Y_W_p(end)) , (Y_W_p(2:end-1))+(-Y_W_p(end))] ;

X_pressure_2 = [flip(X_W_p(2:end-1))+(-X_W_p(end)) , X_p_cir_2(2:end-1)+(-X_W_p(end)) , -(X_W_p(2:end-1))+(-X_W_p(end))] ;
%Y_pressure_2 = [flip(Y_W_p(2:end-1))+(-Y_W_p(end))+(-d_Y_s) , Y_p_cir_2(2:end-1)+(-Y_W_p(end))+(-d_Y_s) , (Y_W_p(2:end-1))+(-Y_W_p(end))+(-d_Y_s)] ;
Y_pressure_2 = [flip(Y_W_p(2:end-1))+(-Y_W_p(end))+(-d_Y_s) , Y_p_cir_2(2:end-1)+(-Y_W_p(end)) , (Y_W_p(2:end-1))+(-Y_W_p(end))+(-d_Y_s)] ;

X_suction_1 = [flip(X_W_s(2:end-1)+(-X_W_p(end)))        , (X_s_cir_1(2:end-1))+(-X_W_p(end)) , -(X_W_s(2:end))+(-X_W_p(end)) ] ;
Y_suction_1 = [flip(Y_W_s(2:end-1)+d_Y_s+(-Y_W_p(end)))  , (Y_s_cir_1(2:end-1))+(-Y_W_p(end)) , (Y_W_s(2:end))+d_Y_s+(-Y_W_p(end))] ;

X_suction_2 = [flip(X_W_s(2:end-1)+(-X_W_p(end)))        , (X_s_cir_2(2:end-1))+(-X_W_p(end)) , -(X_W_s(2:end))+(-X_W_p(end)) ] ;
Y_suction_2 = [flip(Y_W_s(2:end-1)+(-Y_W_p(end)))  , (Y_s_cir_2(2:end-1))+(-Y_W_p(end)) , (Y_W_s(2:end))+(-Y_W_p(end))] ;

X_straight_inlet_1 = 0 : ((X_W_s(end))+(-X_W_p(end)))/1e2 : (X_W_s(end))+(-X_W_p(end)) ;
Y_straight_inlet_1 = 0 : ((Y_W_s(end))+d_Y_s+(-Y_W_p(end)))/1e2 : (Y_W_s(end))+d_Y_s+(-Y_W_p(end)) ;

X_straight_inlet_2 = 0 : ((X_W_s(end))+(-X_W_p(end)))/1e2 : (X_W_s(end))+(-X_W_p(end)) ;
Y_straight_inlet_2 = Y_straight_inlet_1 - d_Y_s;

X_straight_outlet_1 = -(X_W_p(end))+(-X_W_p(end)) : -((X_W_s(end))+(-X_W_p(end)))/1e2 : -(X_W_s(end))+(-X_W_p(end)) ; %((-(X_W_s(end))+(-X_W_p(end))) - (-(X_W_p(end))+(-X_W_p(end))))/1e1
Y_straight_outlet_1 = Y_straight_inlet_1;

X_straight_outlet_2 = -(X_W_p(end))+(-X_W_p(end)) : -((X_W_s(end))+(-X_W_p(end)))/1e2 : -(X_W_s(end))+(-X_W_p(end)) ; %((-(X_W_s(end))+(-X_W_p(end))) - (-(X_W_p(end))+(-X_W_p(end))))/1e1
Y_straight_outlet_2 = Y_straight_inlet_1 - d_Y_s;

X_full_1 = [flip(X_suction_1) , flip(X_straight_inlet_1) , X_pressure_1 , X_straight_outlet_1] ;
Y_full_1 = [flip(Y_suction_1) , flip(Y_straight_inlet_1) , Y_pressure_1 , Y_straight_outlet_1] ;
Z_full_1 = zeros(1,length(X_full_1)) ;
first_row = ones(1,length(X_full_1)) ;
second_row = 1:1:length(X_full_1) ;

X_full_2 = [flip(X_suction_2) , flip(X_straight_inlet_2) , X_pressure_2 , X_straight_outlet_2] ;
Y_full_2 = [flip(Y_suction_2) , flip(Y_straight_inlet_2) , Y_pressure_2 , Y_straight_outlet_2] ;
Z_full_2 = zeros(1,length(X_full_2)) ;
first_row = ones(1,length(X_full_2)) ;
second_row = 1:1:length(X_full_2) ;
%{
scaler = 1 ; % Scales the X and Y coordinates to how much ever you want
X_excel_1 = scaler*[flip(X_W_p(2:end))+(-X_W_p(end)) , X_p_cir_1(2:end-1)+(-X_W_p(end)) , -(X_W_p(1:end-1))+(-X_W_p(end)) , -(X_W_p(end))+(-X_W_p(end)) : 0.01 : -flip(X_W_s(end))+(-X_W_p(end)) , -flip(X_W_s)+(-X_W_p(end))             , flip(X_s_cir_1(2:end-1))+(-X_W_p(end)) , X_W_s+(-X_W_p(end))       , X_p(end)+(-X_W_p(end))] ;
Y_excel_1 = scaler*[flip(Y_W_p(2:end))+(-Y_W_p(end)) , Y_p_cir_1(2:end-1)+(-Y_W_p(end)) , (Y_W_p(1:end-1))+(-Y_W_p(end)) , 0 : 0.01 : flip(Y_W_s(1))+d_Y_s+(-Y_W_p(end))                    , flip(Y_W_s(2:end))+d_Y_s+(-Y_W_p(end)) , flip(Y_s_cir_1(2:end-1))+(-Y_W_p(end)) , Y_W_s+d_Y_s+(-Y_W_p(end)) , 0] ;
Z_excel_1 = zeros(1,length(X_excel_1)) ;

X_excel_2 = X_excel_1 ;
Y_excel_2 = Y_excel_1 - d_Y_s ;
Z_excel_2 = zeros(1,length(X_excel_2)) ;
%}

X_pressure_surf_1 = [flip(X_W_p)+(-X_W_p(end)) , X_p_cir_1(2:end-1)+(-X_W_p(end)) , -(X_W_p)+(-X_W_p(end))] ;
Y_pressure_surf_1 = [flip(Y_W_p)+(-Y_W_p(end)) , Y_p_cir_1(2:end-1)+(-Y_W_p(end)) , (Y_W_p)+(-Y_W_p(end))] ;
Z_pressure_surf_1 = zeros(1,length(X_pressure_surf_1)) ;

X_pressure_surf_2 = [flip(X_W_p)+(-X_W_p(end)) , X_p_cir_2(2:end-1)+(-X_W_p(end)) , -(X_W_p)+(-X_W_p(end))] ;
Y_pressure_surf_2 = [flip(Y_W_p)+(-Y_W_p(end))-d_Y_s , Y_p_cir_2(2:end-1)+(-Y_W_p(end)) , (Y_W_p)+(-Y_W_p(end))-d_Y_s] ;
Z_pressure_surf_2 = zeros(1,length(X_pressure_surf_2)) ;

X_suction_surf_1 = [flip(X_W_s+(-X_W_p(end)))        , (X_s_cir_1(2:end-1))+(-X_W_p(end)) , -(X_W_s)+(-X_W_p(end)) ] ;
Y_suction_surf_1 = [flip(Y_W_s)+d_Y_s+(-Y_W_p(end))  , (Y_s_cir_1(2:end-1))+(-Y_W_p(end)) , (Y_W_s)+d_Y_s+(-Y_W_p(end))] ;
Z_suction_surf_1 = zeros(1,length(X_suction_surf_1)) ;

X_suction_surf_2 = [flip(X_W_s+(-X_W_p(end)))        , (X_s_cir_2(2:end-1))+(-X_W_p(end)) , -(X_W_s)+(-X_W_p(end)) ] ;
Y_suction_surf_2 = [flip(Y_W_s)+(-Y_W_p(end))  , (Y_s_cir_2(2:end-1))+(-Y_W_p(end)) , (Y_W_s)+(-Y_W_p(end))] ;
Z_suction_surf_2 = zeros(1,length(X_suction_surf_2)) ;

figure(3);
hold on
plot(X_pressure_surf_1,Y_pressure_surf_1)
plot(X_suction_surf_1,Y_suction_surf_1)
plot(X_pressure_surf_2,Y_pressure_surf_2)
plot(X_suction_surf_2,Y_suction_surf_2)

disp(d_Y_s)

% Combine arrays column-wise
%M = [X_excel_1(:) Y_excel_1(:) Z_excel_1(:)];
M_1 = [X_full_1(:) Y_full_1(:) Z_full_1(:)];
M_2 = [X_full_2(:) Y_full_2(:) Z_full_2(:)];
M_pressure_1 = [X_pressure_surf_1(:) Y_pressure_surf_1(:) Z_pressure_surf_1(:)];
M_pressure_2 = [X_pressure_surf_2(:) Y_pressure_surf_2(:) Z_pressure_surf_2(:)];
M_suction_1 = [X_suction_surf_1(:) Y_suction_surf_1(:) Z_suction_surf_1(:)];
M_suction_2 = [X_suction_surf_2(:) Y_suction_surf_2(:) Z_suction_surf_2(:)];
M_ansys_1 = [first_row(:) second_row(:) X_full_1(:) Y_full_1(:) Z_full_1(:)] ;
M_ansys_2 = [first_row(:) second_row(:) X_full_2(:) Y_full_2(:) Z_full_2(:)] ;

%{
% Write to Excel and notepad
writematrix(M, 'rotor_cords.xlsx');
writematrix(M_pressure,'rotor_p_coords.txt')
writematrix(M_suction,'rotor_s_coords.txt')
writematrix(M, 'rotor_cords.txt');
writematrix(M_ansys, 'rotor_cords_ansys.txt');
%}

writematrix(M_1,          fullfile('profiles','rotor_cords_1.xlsx'));
writematrix(M_2,          fullfile('profiles','rotor_cords_2.xlsx'));
writematrix(M_pressure_1, fullfile('profiles','rotor_p_coords_1.txt'));
writematrix(M_pressure_2, fullfile('profiles','rotor_p_coords_2.txt'));
writematrix(M_suction_1,  fullfile('profiles','rotor_s_coords_1.txt'));
writematrix(M_suction_2,  fullfile('profiles','rotor_s_coords_2.txt'));
writematrix(M_1,          fullfile('profiles','rotor_cords_1.txt'));
writematrix(M_2,          fullfile('profiles','rotor_cords_2.txt'));
%writematrix(M_ansys,    fullfile('profiles','rotor_cords_ansys_1.txt'));

% Open file for writing
fileID1 = fopen(fullfile('profiles','M_ansys_1.txt'),'w');
fileID2 = fopen(fullfile('profiles','M_ansys_2.txt'),'w');

% Define format: integers for first two columns, 6 decimals for the rest
fprintf(fileID1, '%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n', M_ansys_1');
fprintf(fileID2, '%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n', M_ansys_2');

% Close the file
fclose(fileID1);
fclose(fileID2);

end
%% FUNCTIONS

function [x,y] = plotArc(r, theta_c, G, nl)
    % plotArc(r, theta_c)
    % Plots a circular arc centered at (0,0)
    % r       = radius of the arc
    % theta_c = central angle subtended by the arc (in radians)

    % Define arc angles (symmetric about x-axis)
    %theta = linspace(-theta_c/2, theta_c/2, 1000);
    
    %theta = linspace(-theta_c, theta_c, 51);
    theta = linspace(-theta_c, theta_c, (nl*20 + 1));

    % Coordinates of the arc
    x = r * sind(theta);
    y = r * cosd(theta) + G;

    % Plot arc
    plot(x, y, 'b-', 'LineWidth', 2);
    axis equal; grid on;

    % Mark start and end of arc
    %plot(x(1), y(1), 'ro', 'MarkerSize', 8, 'LineWidth', 1);
    %plot(x(end), y(end), 'go', 'MarkerSize', 8, 'LineWidth', 1);
end

function R_solve = solve_R(v_i, Y, dv, k,p_or_s)
    if p_or_s == true   % PRESSURE SURFACE
        Target = f_R_star_with_dv(v_i,Y,dv,k) ;

    else                % SUCTION SURFACE
        Target = f_R_star_with_dv_s(v_i,Y,dv,k) ;
    end
    R_star = sqrt((Y+1)/(Y-1)) ;
    R_lower = 0 ;
    error = 1 ;
    step = 1e-6 ;
    while error > step
        f_R = f_R_star_with_R_star(Y,R_star) ;
        error = abs(f_R - Target) ;
        %fprintf("Target:%.3f\tf_R:%.3f\terror:%.3f\tR_star:%.3f\n",Target,f_R,error,R_star)
        if error < 10*step
            R_solve = R_star ;
            break
        end
        R_star = R_star - step ;
        if R_star < R_lower
            R_solve = nan ;
            break
        end
    end
end

function f_R = f_R_star_with_R_star(Y,R_star)
    term1 = sqrt((Y+1)/(Y-1)) * asin( ((Y-1)/R_star^2) - Y) ;
    term2 = asin( ((Y+1)*R_star^2) - Y) ;
    f_R = term1 + term2 ;
end

function f_R = f_R_star_with_dv(v_i,Y,dv,k)
    v_i = v_i * (pi/180) ; % degrees to radians
    dv = dv * (pi/180) ; % degrees to radians
    term1 = 2*v_i  ;
    term2 = - (pi/2)*(sqrt((Y+1)/(Y-1))-1) ;
    term3 = - 2*(k-1)*dv ;
    f_R = term1 + term2 + term3 ;
end

function f_R = f_R_star_with_dv_s(v_i,Y,dv,k)
    v_i = v_i * (pi/180) ; % degrees to radians
    dv = dv * (pi/180) ; % degrees to radians
    term1 = 2*v_i  ;
    term2 = - (pi/2)*(sqrt((Y+1)/(Y-1))-1) ;
    term3 = 2*(k-1)*dv ;
    f_R = term1 + term2 + term3 ;
end

function d_theta = dtheta_cir(x,y,inlet_condition)
% x -> Beta
% y -> dtheta_trans
% inlet_condition == true , it is the inlet condition (Intlet of rotor)
% else , it is the outlet condition (Outlet of rotor)
    if inlet_condition == true
        d_theta = x - y;
    else
        d_theta = x + y;
    end
end

function d_theta = dtheta_trans(x,y,Y,p_or_s_condition)
% x -> inlet or outlet Mach numbers
% y -> pressure or suction surface Mach numbers
% p_or_s_condition == true , it is a pressure surface (lower or concave)
% else , it is a suction surface (upper or convex)
    if p_or_s_condition == true
        d_theta = PM(x,Y) - PM(y,Y);
    else
        d_theta = PM(y,Y) - PM(x,Y);
    end
end

function m_star = M_star(M,Y)
    m_star = ( ( ((Y+1)/2) * M^2) / ( (1 + ((Y-1)/2) * M^2)) )^(1/2) ;
end

function v = PM(M,Y)
    v = sqrt( (Y+1)/(Y-1) ) * atand( sqrt( ((Y-1)/(Y+1))*(M^2 - 1) ) ) - atand( sqrt(M^2 - 1) ) ;
end

% Mach number and Prandtl-Meyer function
function M2 = Mx(v1,M1,v2)
    M2 = interp1(v1,M1,v2) ;
end

% Mach angle
function u = ux(M)
    u = asind(1/M) ;
end

% K- characteristic line
function z = K_minus(theta1,v1)
    z = theta1 + v1 ;
end

% K+ characteristic line
function z = K_plus(theta1,v1)
    z = theta1 - v1 ;
end

% Angle of C- characteristic slope
function s = m_minus(Ox,Oy,ux,uy)
    s = ( (Ox-ux) + (Oy-uy) ) / 2 ;
end

% Angle of C+ characteristic slope
function s = m_plus(Ox,Oy,ux,uy)
    s = ( (Ox+ux) + (Oy+uy) ) / 2 ;
end

% Slope
function s = t(u)
    s = tand(u) ;
end

% Area ratio
function s = E(M,Y)
    a = (Y+1)/2 ;
    b = (Y-1)/2 ;
    s = (1/M) * ((1/a) * (1 + (b*(M^2))))^(a/(2*b));
end

% Q factor
function Q = Q_factor(v_U,v_L)
% QFACTOR  Computes refined NACA-type Q factor
% as a function of convex (vu) and concave (vl) property angles (degrees).
%
%   Q = Qfactor(vu, vl)
%
%   Inputs:
%       vu - convex surface property angle (degrees)
%       vl - concave surface property angle (degrees)
%   Output:
%       Q  - Q-factor (dimensionless)
%
%   Model:
%       Q = 1 - (A + B*vl) * (1 - exp(-C*vu))

    % Empirical coefficients (fitted to NACA chart)
    A = 0.20;     % base decay amplitude
    B = 0.0085;   % dependence on vl
    C = 0.045;    % exponential saturation rate

    % Compute Q-factor
    Q = 1 - (A + B .* v_L) .* (1 - exp(-C .* v_U));

    % Ensure Q stays within chart bounds
    Q(Q < 0.6) = 0.6;
end

% G* nondimensional blade spacing
function [G,G_new] = G_star(M_i,Y,v_L,v_U,R_star_L,R_star_U,Turning_angle,M_star_L,M_star_U)
    Ae_Astar = E(M_i,Y) ;
    Q = Q_factor(v_U,v_L) ;
    Q_n = Q_new(M_star_L,M_star_U) ;
    G = Ae_Astar * Q * (R_star_L - R_star_U) / cosd(Turning_angle/2) ;
    G_new = Ae_Astar * Q_n * (R_star_L - R_star_U) / cosd(Turning_angle/2) ;
end

function Q = Q_new(M_star_L,M_star_U)
    integrand_2 = Q_integrand(M_star_U) - Q_integrand(M_star_L) ;
    Q = M_star_L*M_star_U/(M_star_U-M_star_L)*integrand_2 ;
end

function integrand = Q_integrand(M_star)
    integrand = ( ((6-M_star^2)^(5/2)/280) + ((6-M_star^2)^(3/2)/28) - ((6-M_star^2)^(1/2)/1.555) - 1.571*log((2.45 + ( (6-M_star^2)^(1/2)))/M_star) ) ;
end




