##TP2 Caso 1: Motor
pkg load control signal io symbolic;
clc; clear all; close all;

%Definición de variables útiles
h=10e-7;
tiempo=(.3/h);
t=0:h:(tiempo*h);
i=1;
retardo_TL=0.186300000000042;

##Los componentes calculados en el TP1 son
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
Kp=10; Ki=0.5;Kd=0.001;

A1=((2*Kp*h)+(Ki*(h^2))+(2*Kd))/(2*h);
B1=(-2*Kp*h+Ki*(h^2)-4*Kd)/(2*h);
C1=Kd/h;

e_tita(1)=0;

#K ampliado
Q=diag([1 1e-3 1 1e5]); %%Limita la velocidad de respuesta pero más que nada la intensidad de los efectos del cambio de control.
%%El último valor representa el error, los demás representan las variables de estado en ese orden
R=1e-3; %%Limita la velocidad de respuesta del controlador, refiriéndose a la velocidad a la que X_hat se aproxima a X
A_ampliado=[A zeros(3,1);-C_t 0]
B_ampliado=[B(:,1);0]
Ka=lqr(A_ampliado,B_ampliado,Q,R);

%%creamos las entradas
titaref=(pi/2)*square(2*pi*(t+retardo_TL)/.1); %%varía entre pi/2 y -pi/2
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
Polos_Obs=[-1000;-10000+10i;-10000-10i]; % Polos del observador
alfa=poly(Polos_Obs);
M=[B_dual A_dual*B_dual (A_dual^2)*B_dual]; %Matriz de controlabilidad dual
autoval=eig(A_dual);
c_ai=poly(autoval);
W=[c_ai(3) c_ai(2) 1; c_ai(2) 1 0; 1 0 0];
T=M*W;
Ko=(fliplr(alfa(2:4)-c_ai(2:4))*inv(T))'; %Ganancia de nuestro observador

X_hat=[0; 0; 0];Yo(2)=0; ia_ob(1)=0;
psi=0;

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
  u=[-Ka*[X_hat;psi]; TL(i)];
  if (abs(u(1))<1)
    u(1)=0;
  else
    u(1)=sign(u(1))*(abs(u(1))-1);
  end
  acc(i)=u(1);
  X_p=A*X+B*u;
  X=X+h*X_p;
  psi=psi+h*e_tita(k);

  %%ciclo del observador

  X_hat_p=A*X_hat+B*u+Ko*(ia(i-1)-X_hat(1));
  X_hat=X_hat+X_hat_p*h;
  ia_ob(i)=X_hat(1);
end
toc

#Se nos pide medir la corriente con el observador, es por eso que ahora, superponemos las gráficas de la corriente simulada y la corriente medida por el observador
figure 1;
plot(t,ia), hold on, grid on, plot(t,ia_ob,"-."), title("Corriente de armadura"),legend("corriente real","corriente del observador");

figure 2;
subplot(4,1,1),plot(t,wr),hold on, grid on, title("Velocidad angular");
subplot(4,1,2),plot(t,ia),grid on, title ("Corriente de armadura");
subplot(4,1,3),plot(t,tita),grid on, title("tita");
subplot(4,1,4),plot(t,TL),grid on; title("carga");

figure 3;
plot(t,acc), hold on, grid on, title("Acción de control (V_A)");

