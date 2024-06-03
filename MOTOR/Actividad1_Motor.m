##TP1 inciso 4
pkg load control signal symbolic;
clc; clear all; close all;

%Definición de variables útiles
h=1e-7;
tiempo=(.3/h);
t=0:h:(tiempo*h);
i=1;

%Definición de variables
LAA=366e-6;
J=5e-9;
RA=55.6;
Bm=0;
Ki=6.49e-3;
Km=6.53e-3;

%Entrada
TL=heaviside(t-0.1)*1.3925e-03;%El valor maximo es 1.3926e-03, sale de hacer ia_max*ki, despejandolo de la ecuación de wr_p cuando wr=0.
Va=heaviside(t-0.01)*12;


ia(1)=0; wr(1)=0; tita(1)=0;

%Definiendo las matrices
A=[-RA/LAA -Km/LAA 0;Ki/J -Bm/J 0;0 1 0]
B=[1/LAA 0;0 -1/J ;0 0]
C_t=[0 1 0]
D=0
X=[ia(1) ; wr(1) ; tita(1)];
u=[Va(1); TL(1)];

tic
for(i=1:1:(tiempo+1))

  ia(i)=X(1);wr(i)=X(2);tita(i)=X(3); u=[Va(i); TL(i)];

  X_p=A*X+B*u;
  X=X+h*X_p;

end
toc

sprintf("La corriente máxima es %f [A], y la carga máxima %f [Nm]", max(ia),max(TL))

subplot(4,1,1);plot(t,ia); title("Corriente");grid on;
subplot(4,1,2);plot(t,wr); title("Velocidad Angular");grid on;
subplot(4,1,3);plot(t,tita); title("Angulo ");grid on;
subplot(4,1,4);plot(t,TL);title("Carga");grid on;
