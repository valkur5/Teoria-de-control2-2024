##TP1, inciso 6
pkg load control signal io symbolic;
clc; clear all; close all;

%Definición de variables útiles
h=10e-7;
tiempo=(.1/h);
t=0:h:(tiempo*h);
i=1;
retardo_TL=0.186300000000042;

##Los componentes calculados en el inciso ítem anterior son
RA = 28.131
J = 1.8826e-09
Km = 0.060530
Ki = 0.014336
LAA = 1.5349e-03
Bm = 0

%Las matrices son:
A=[-RA/LAA -Km/LAA 0;Ki/J -Bm/J 0;0 1 0]
B=[1/LAA 0;0 -1/J ;0 0]
C_t=[0 1 0]
D=0

#Constantes del PID
Kp=10; Ki=0.5;Kd=0.001;

%Las adaptamos para trabajar en tiempo discreto
A1=((2*Kp*h)+(Ki*(h^2))+(2*Kd))/(2*h);
B1=(-2*Kp*h+Ki*(h^2)-4*Kd)/(2*h);
C1=Kd/h;

e_tita(1)=0;

%%creamos las entradas
titaref=1; %%1 radian
TL=(1.03e-03/2)*square(2*pi*(t+retardo_TL)/.1)+(1.03e-3/2);
Va=heaviside(t-0.005)*12;
ia(1)=0; wr(1)=0; tita(1)=0;
X=[ia(1) ; wr(1) ; tita(1)];
u=[Va(1); TL(1)];
acc(1)=0;

tic
for(i=1:1:(tiempo+1))
  k=i+2;
  ia(i)=X(1);wr(i)=X(2);tita(i)=X(3);
  e_tita(k)=titaref-tita(i);

  u=[u(1)+A1*e_tita(k)+B1*e_tita(k-1)+C1*e_tita(k-2); 0*TL(i)];
  acc(i)=u(1);
  X_p=A*X+B*u;
  X=X+h*X_p;
end
toc

subplot(4,1,1),plot(t,wr),hold on, grid on, title("Velocidad angular");
subplot(4,1,2),plot(t,ia),grid on, title ("Corriente de armadura");
subplot(4,1,3),plot(t,tita),grid on, title("tita");
subplot(4,1,4),plot(t,TL),grid on; title("carga");
