##TP2 Caso 2: Avión
pkg load control signal io symbolic;
clc; clear all; close all;

%Definición de variables útiles
Ts=1e-4; %Sale de dividir la parte real de los polos por 30
h=3e-3;
tiempo=(100/h);
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

C=[ 0 0 1 0;
    0 0 0 1]
D=0
X=[0 0 0 500]'
%%Matrices de espacio de estados en tiempo discreto
Ad=e^(A*Ts)
Bd=[0;0;0;0];
d_tao=Ts/100; tao=0;
for hh=1:100 %%Hacemos una integral y suponemos 100 muestras, por eso este número
d_Bd=e^(A*tao)*B*d_tao;
Bd=Bd+d_Bd;
tao=tao+d_tao;
end
Bd %Imprimimos B discreto
##ss_discreto=c2d(ss(A,B,C,D),Ts); %esto nos discretiza nuestras matrices de estados solamente con las matrices y el tiempo de muestreo

%Asignamos las matrices correspondientes
##Ad=ss_discreto.a
##Bd=ss_discreto.b

alfa_d(1)=0; tita_d(1)=0;tita_p_d(1)=0; altura_d(1)=0;
u(1)=1; X_P=[0 0 0 0]'; td=0;
for i=1:100
  td=[td i];
  u(i)=1;
  X_p=Ad*X+Bd*u(i);
  alfa_d(i+1)=X_p(1);
  tita_d(i+1)=X_p(2);
  tita_p_d(i+1)=X_p(3);
  altura_d(i+1)=X_p(4);
end
X_p=0;X=0;
for i=1:(tiempo+1)
  X_p=A*X+B*u(1);
  X=X+X_p*h;
  alfa(i)=X(1);
  tita(i)=X(2);
  tita_p(i)=X(3);
  altura(i)=X(4);
end

figure 1;
subplot(4,1,1);plot(td,alfa_d);grid on; title("Alfa")
subplot(4,1,2);plot(td,tita_d);grid on;title("tita")
subplot(4,1,3);plot(td,tita_p_d);grid on;title("tita punto")
subplot(4,1,4);plot(td,altura_d);grid on;title("altura h")

figure 2;
subplot(4,1,1);plot(t,alfa);grid on; title("Alfa")
subplot(4,1,2);plot(t,tita);grid on;title("tita")
subplot(4,1,3);plot(t,tita_p);grid on;title("tita punto")
subplot(4,1,4);plot(t,altura);grid on;title("altura h")
