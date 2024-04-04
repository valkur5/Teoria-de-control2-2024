pkg load control signal;
clc; clear all; close all;

%Variables de utilidad
h=1e-6; %Debe ser más chico que la frecuencia más rápida del sistema
tiempo=(0.2/h);%2e5
t=0:h:(tiempo*h);
i=1;

%componentes
R=47;
L=1e-3;
C=100e-6;
Ve=0;
I(1)=0;%X1
V_c(1)=0;%X2

%matrices del espacio de estados (ss)
A=[-R/L -1/L; 1/C 0]
B=[1/L; 0]
Ct=[R 0]
D=[0]
X=[I ; V_c];
u(1)=0;
[num,den]=ss2tf(A,B,Ct,D);
Y(1)=0;

while(i<(tiempo+1))
  u(i)=Ve;
  I(i)=X(1);V_c(i)=X(2);
  if( mod(i, 10000) == 0) %Cambia el valor de la tensión de alimentación cada 10 mil tics, que es 1ms
    if (Ve==12)
      Ve=-12;
    else
      Ve=12;
    end
  end
  X_P=A*X+B*u(i);%X punto
  X=X+h*X_P;%Esto es el cálculo de la integral como sumatoria

  Y(i)=R*I(i);
  i=i+1;

end

figure 1;
subplot(3,1,1) ;plot(t,I); title("corriente"); grid on;
subplot(3,1,2) ;plot(t,V_c); title("Tension capacitor");grid on;
subplot(3,1,3) ;plot(t,u); title("Tension de entrada");grid on;
figure 2;
plot(t,Y);title("tension de la resistencia"); grid on;
