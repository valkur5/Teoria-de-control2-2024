%%TP1 Inciso 2 y 3
pkg load control signal io;
clc; clear all; close all;


%%Variables de utilidad
h=0.0001;%%Se extrae de analizar el archivo del profesor
tiempo=(0.2/h);
t=0:h:(tiempo*h);
i=1;

%%Datos del archivo exel
data = transpose(csvread("Curvas_Medidas_RLC.csv"));%se transponen los datos porque en el archivo los datos están en columnas, y para separarlos en variables necesitamos esos datos en filas
t_d=data(1,:);
I_d=data(2,:);
V_c_d=data(3,:);
Ve_d=data(4,:);

%% Determinación de componentes
Amplitud=12;
t_inic=0.0007; %se determina este tiempo como el tiempo que tarda la señal en llegar al 60% aprox de su valor máximo, esto se obtiene de analizar las gráficas otorgadas por el profesor
retardo=0.01;%Retardo que hay antes de empezar la simulación

sys_G_ang= Chen(Amplitud, t_d, V_c_d, t_inic,retardo, 12)

[ychen ,tchen,ent]=lsim(sys_G_ang, Ve_d,t_d, [0,0]);

figure 1
plot(t_d,V_c_d); title("Tension del capacitor");grid on;hold on;plot(tchen,ychen,"--"); legend("curva real","curva aproximada");

%%Calculo de la simulación de nuestros valores

corriente=diff(ychen)/h;
C=max(I_d)/max(corriente)
R=sys_G_ang.den{1}(2)/C
L=sys_G_ang.den{1}(1)/C

%Determinacion de matrices de estado
I(1)=0;
V_c(1)=0;
Ve=0;

A=[-R/L -1/L; 1/C 0]
B=[1/L; 0]
Ct=[0 1]
D=[0]
X=[0 ; 0];
u(1)=0;

%Calculo con el método de euler:
while(i<(tiempo+2))
  u(i)=Ve;
  I(i)=X(1);V_c(i)=X(2);
  if( i==100 ) %Genera un retardo de 100 tics, que son aprox. 10ms, la entrada es 0 hasta entonces
  Ve=12;
end
if( mod(i, 500) == 0) %Cambia el valor de la entrada cada 500 tics luego de eso, que son aprox 50 ms
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

%%Impresión de valores


figure 2;
subplot(3,1,1);plot(t_d,I_d); title("Corriente");grid on;hold on;plot(t,I,'--');legend("corriente real","corriente aproximada");
subplot(3,1,2);plot(t_d,V_c_d); title("Tension del capacitor");grid on;hold on;plot(t,V_c,"--");legend("tensión real","tensión aproximada");
subplot(3,1,3);plot(t_d,Ve_d); title("Voltaje de entrada");grid on;hold on;plot(t,u,"--");



