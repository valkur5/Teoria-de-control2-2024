##TP2 Caso 2: Avión
pkg load control signal io symbolic;
clc; clear all; close all;

%Definición de variables útiles
Ts=1e-2; %Sale de dividir la parte real de los polos por 30
h=3e-3;
tiempo=(150/h);
t=0:h:(tiempo*h);
i=1;

%Variables
a=0.07;
w=9;
b=5;
c=150;

%Matrices de spacio de estados en tiempo continuo
A=[-a a 0 0;
    0 0 1 0;
    w^2 -w^2 0 0;
    c 0 0 0]

B=[ 0;
    0;
    b*w^2;
    0]

C=[ 0 0 0 1;
    0 1 0 0];
D=0


%Calculamos el controlador K utilizando Ackerman:
Aa = [A  zeros(4,1) ; -C(1,:) 0];
Ba = [B; 0];

poles = [-15+15i -15-15i -0.0427+0.0492i -0.0427-0.0492i -0.0775];
K_ack = acker(Aa, Ba, poles);
Ka = K_ack(1:4);
KI = -K_ack(end);

%%Cálculo del observador, Se utiliza el método del sistema dual
C_Obs=B';
A_Obs=A';
B_Obs=C';

Qo = 100*diag([1 100 1 100]);
Ro = diag([100000000 100000000]);
Ko=dlqr(A_Obs,B_Obs,Qo,Ro)'; %K del observador


altura(1)=500; %Condición inicial
tita_d(1)=0;tita_p_d(1)=0;
u(1)=1; X_P=[0 0 0 0]';Xang=[0;0;0;500]; ref=100;
X=[0 0 0 altura(1)]';X_p=0;alpha_obs(1)=0;
psi_p=0;psi=0; phi(1)=0;alpha(1)=0;phi_p(1)=0;

for i=1:tiempo

  %Medimos las salidas
  Y_obs = C*(Xang);
  Ys = C*X;

  %Calculo del error
  psi_p=ref-Ys(1);
  psi=psi+psi_p*h;

  %Acción de control
  u(i)=-Ka*X+KI*psi;

  %Zona muerta
  if abs(u(i))<.1
        u(i)=0;
    else
        u(i)=sign(u(i))*(abs(u(i))-.1);
    end
  %%Modelo no lineal
  alpha_p    = a*(phi(i) - alpha(i));
  phi_pp     = -w^2*(phi(i) - alpha(i) -b*u(i));
  altura_p   = c*alpha(i);
  alpha(i+1) = alpha(i) + h*alpha_p;
  phi_p(i+1) = phi_p(i) + h*phi_pp;
  phi(i+1)   = phi(i) + h*phi_p(i);
  altura(i+1)= altura(i) + h*altura_p;

  % actualizacu�n de estados
  X = [alpha(i+1) phi(i+1) phi_p(i+1) altura(i+1)]';

  %Observador
  Xang_p=A*Xang+B*u(i)+Ko*[Ys-Y_obs];
  Xang=Xang+Xang_p*h;
  alpha_obs(i+1)=Xang(1);
end

u(i+1)=-Ka*X+KI*psi;
figure 1;
subplot(3,2,1);plot(t,alpha);grid on; title("Alfa"); hold on, plot(t,alpha_obs,"-.");
legend("sistema real","observador");
subplot(3,2,2);plot(t,phi);grid on;title("tita")
subplot(3,2,3);plot(t,phi_p);grid on;title("tita punto")
subplot(3,2,4);plot(t,u);grid on; title("acción de control")
subplot(3,1,3);plot(t,altura);grid on;title("altura h")


