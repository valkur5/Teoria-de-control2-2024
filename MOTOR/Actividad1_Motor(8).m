##TP1, inciso 8
pkg load control signal io symbolic;
clc; clear all; close all;

%Definición de variables útiles
h=1e-7;
tiempo=(.1/h);
t=0:h:(tiempo*h);
i=1;
retardo_TL=0.186300000000042;

##Los componentes calculados en el inciso ítem 5 anterior son
RA = 28.131
J = 1.8826e-09
Km = 0.060530
Ki = 0.014336
LAA = 1.5349e-03
Bm = 0

%Las matrices son:
A=[-RA/LAA -Km/LAA 0;Ki/J -Bm/J 0;0 1 0]
B=[1/LAA 0;0 -1/J ;0 0]
C_t=[0 0 1]
D=0

#Constantes del PID
Kp=0.1; Ki=0.01;Kd=5;

A1=((2*Kp*h)+(Ki*(h^2))+(2*Kd))/(2*h);
B1=(-2*Kp*h+Ki*(h^2)-4*Kd)/(2*h);
C1=Kd/h;

e_tita(1)=0;

%%creamos las entradas
titaref=(pi/2)*square(2*pi*(t+retardo_TL)/.1); %%1 radian
TL=0;
Va=heaviside(t-0.005)*12;
ia(1)=0; wr(1)=0; tita(1)=0;
X=[ia(1) ; wr(1) ; tita(1)];
u=[Va(1); TL(1)];
acc(1)=0;

%%Creamos el observador
A_dual=A';
B_dual=C_t';
C_dual=B';
Polos_Obs=[-.5e2;-5e2+0.5i;-5e2-0.5i]; % Polos del observador
alfa=poly(Polos_Obs*20);
M=[B_dual A_dual*B_dual (A_dual^2)*B_dual]; %Matriz de controlabilidad dual
autoval=eig(A_dual);
c_ai=poly(autoval);
W=[c_ai(3) c_ai(2) 1; c_ai(2) 1 0; 1 0 0];
T=M*W;
Ko=(fliplr(alfa(2:4)-c_ai(2:4))*inv(T))'; %Ganancia de nuestro observador

X_hat=[0; 0; 0];Yo(2)=0; ia_ob(1)=0;


tic
for(i=2:1:(tiempo+1))
  k=i+2;
  ia(i)=X(1);wr(i)=X(2);tita(i)=X(3);
  if(titaref(i)==(pi/2))
    TL(i)=1.15e-3;
  else
    TL(i)=0;
  end
  e_tita(k)=titaref(i)-tita(i);
  u=[u(1)+A1*e_tita(k)+B1*e_tita(k-1)+C1*e_tita(k-2); TL(i)];
  acc(i)=u(1);
  X_p=A*X+B*u;
  X=X+h*X_p;

  %%ciclo del observador
  X_hat_p=A*X_hat+B*u+Ko*(ia(i-1)-X_hat(1));
  X_hat=X_hat+X_hat_p*h;
  ia_ob(i)=X_hat(1);
end
toc

#Se nos pide medir la corriente con el observador, es por eso que ahora, superponemos las gráficas de la corriente simulada y la corriente medida por el observador
plot(t,ia), hold on, grid on, plot(t,ia_ob,"-."), legend("corriente real","corriente del observador");
