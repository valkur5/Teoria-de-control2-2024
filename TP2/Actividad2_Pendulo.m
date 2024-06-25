##TP2 Caso 3: Pendulo invertido
pkg load control signal io symbolic;
clc; clear all; close all;
%% variables utiles
Ts=1e-3;
i=1; ki=1;
h = Ts/2.0;
t = 0:h:15;
p_max = floor(15/Ts);
%Parámetros
m  = 0.1;
m_ = m*10;
F  = 0.1;
l  = 1.6;
g  = 9.8;
M  = 1.5;

% Matrices de estados
% RECORRIDO DE 0 A 10: considera m
A=[0       1               0        0;
      0    -F/M            -m*g/M      0;
      0       0               0        1;
      0   -F/(l*M)     -g*(m+M)/(l*M)  0];

% RECORRIDO DE 10 A 0: considera m_
A_=[0       1             0         0;
       0    -F/M         -m_*g/M       0;
       0       0             0         1;
       0   -F/(l*M)  -g*(m_+M)/(l*M)   0];

B=[   0   ;
      1/M  ;
       0   ;
   1/(l*M)];

C = [1 0 0 0;
      0 0 1 0];

D=0;
% Sistema discreto
sysDisc = c2d(ss(A,B,C,D), Ts, 'zoh');
Ad     = sysDisc.a;
Bd     = sysDisc.b;
Cd     = sysDisc.c;
Dd     = sysDisc.d;

sysDisc_ = c2d(ss(A_,B,C,D), Ts, 'zoh');
Ad_     = sysDisc_.a;
Bd_     = sysDisc_.b;
Cd_     = sysDisc_.c;
Dd_     = sysDisc_.d;

% Sistema ampliado
Aa = [Ad , zeros(4,1) ; -Cd(1,:)*Ad, eye(1)];
Ba = [Bd; -Cd(1,:)*Bd];

Aa_ = [Ad_ , zeros(4,1) ; -Cd_(1,:)*Ad_, eye(1)];
Ba_ = [Bd_; -Cd_(1,:)*Bd_];

%%Calculamos el LQR para cada caso
Q1 = diag([1 1 5 5 .000003]); R1=1;
Q2 = diag([400 200 100 100 .000001]); R2=20;

[Klqr1, ~, ~] = dlqr(Aa, Ba, Q1, R1);
[Klqr2, ~, ~] = dlqr(Aa_, Ba_, Q2, R2);

K  = Klqr1(1:4);
KI = -Klqr1(5);

K_  = Klqr2(1:4);
KI_ = -Klqr2(5);

%%Calculamos ahora el observador

% para m
Ao = Ad';
Bo = Cd';
Co = Bd';
Qo = diag([100 5000 50000 10]);
Ro = diag([.01 .01]);
Ko = (dlqr(Ao,Bo,Qo,Ro))';

% para m_
Ao_ = Ad_';
Bo_ = Cd_';
Co_ = Bd_';

Qo_ = diag([1000 5000 9000 10]);
Ro_ = diag([.001 .001]);
Ko_ = (dlqr(Ao_,Bo_,Qo,Ro))';

%%Condiciones de simulación
d(1)     = 0;
d_p(1)   = 0;
phi(1)   = pi;
phi_p(1) = 0;
phi_pp(1) = 0;

posRef = 10;            % referencia de posicion a donde se quiere desplazar el carro

X= [d(1) d_p(1) phi(1) phi_p(1)]';
xop = [0 0 pi 0]'; %Punto de operacion

u = [];
flag  = 0;

ei = 0;
x_hat = [0 0 pi 0]';

for ki=1:p_max

    % Salida de dos componentes
    Ys   = Cd*X;             % salida del sistema
    Y_obs = Cd*(x_hat+xop);   % salida del observador

    ei= ei+posRef-Ys(1);

    %Ley de control
    u1(ki) = -K*(X - xop) + KI*ei;          % sin observador
%     u1(ki)  = -K*(x_hat - xop) + KI*ei;     % con observador

    % Zona Muerta
%     if(abs(u1(ki)) < 0.5)
%         u1(ki) = 0;
%     else
%         u1(ki) = sign(u1(ki))*(abs(u1(ki)) - 0.5);
%     end


    % Integraciones de Euler por paso de simulaci�n
    for kii=1:Ts/h

        u(i) = u1(ki);

        % C�lculo por sistema no lineal
        d_pp       = (1/(M+m))*(u(i) - m*l*phi_pp*cos(phi(i)) + m*l*phi_p(i)^2*sin(phi(i)) - F*d_p(i));
        phi_pp     = (1/l)*(g*sin(phi(i)) - d_pp*cos(phi(i)));
        d_p(i+1)   = d_p(i) + h*d_pp;
        d(i+1)     = d(i) + h*d_p(i);
        phi_p(i+1) = phi_p(i) + h*phi_pp;
        phi(i+1)   = phi(i) + h*phi_p(i);

        if(d(i) >= 9.99)
            if(flag == 0)
                posRef  = 0;
                m  = m_;
                flag = 1;
                K  = K_;
                Ki = KI_;
                Ko = Ko_;
                Ad    = Ad_;
                Bd    = Bd_;
            end
        end
        i=i+1;
    end

    % Actualizaci�n de los estados
    % Estados del sistema
    X     = [d(i-1) d_p(i-1) phi(i-1) phi_p(i-1)]';
    % Estados estimados por el observador
    x_hat = Ad*(x_hat-xop) + Bd*u1(ki) + Ko*(Ys - Y_obs) + xop;
end

u(i) = u1(ki);
figure(1);
subplot(3,2,1); grid on; hold on;
plot(t,phi_p,'LineWidth',1.5);grid on; title('Velocidad angular \phi_p');

subplot(3,2,2); grid on; hold on;
plot(t,phi,'LineWidth',1.5); title('angulo \phi');xlabel('Tiempo');

subplot(3,2,3); grid on; hold on;
plot(t,d,'LineWidth',1.5);title('Posicion grua \delta');xlabel('Tiempo');

subplot(3,2,4); grid on; hold on;
plot(t,d_p,'LineWidth',1.5);title('Velocidad de grua \delta_p');

subplot(3,1,3); grid on; hold on;
plot(t,u,'LineWidth',1.5);title('Accion de control u');xlabel('Tiempo en Seg.');

figure(2);
subplot(2,1,1);grid on; hold on;
plot(phi,phi_p,'LineWidth',1.5);
title('angulo vs Velocidad angular');
xlabel('angulo');ylabel('Velocidad angular');

subplot(2,1,2);grid on; hold on;
plot(d,d_p,'LineWidth',1.5);
title('Distancia vs velocidad');
xlabel('Distancia');ylabel('Velocidad');

