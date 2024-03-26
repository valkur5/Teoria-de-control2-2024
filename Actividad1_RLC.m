pkg load control signal;
clc; clear all; close all;

%Variables de utilidad
h=0.0001;
tiempo=(0.1/h);
t=0:h:(tiempo*h);
i=1;

%componentes
R=47;
L=1e-3;
C=100e-6;
Ve=12;
I(1)=0;%X1
V_c(1)=0;%X2

%matrices del espacio de estados (ss)
A=[-R/L -1/L; 1/C 0]
B=[1/L; 0]
Ct=[R 0]
D=[0]
X=[0 ; 0];
u(1)=0;
Sistema=ss(A,B,Ct,D);

while(i<(tiempo+1))
  u(i)=Ve;
  I(i)=X(1);V_c(i)=X(2);
  X_P=A*X+B*u(i);%X punto
  X=X+h*X_P;%Esto es el cálculo de la integral como sumatoria
  i=i+1;

end

I(i)=X(1);V_c(i)=X(2);

subplot(2,1,1) ;plot(t,(I/max(I))); title("corriente");
subplot(2,1,2) ;plot(t,V_c/max(V_c)); title("Tension capacitor");

