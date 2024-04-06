pkg load control signal io;
clc; clear all; close all;
%%Variables de utilidad
h=0.0001;
tiempo=(0.2/h);
t=0:h:(tiempo*h);
i=1;

%%Datos del archivo exel
data = transpose(csvread("Curvas_Medidas_RLC.csv"));
t_d=data(1,:);
I_d=data(2,:);
V_c_d=data(3,:);
Ve_d=data(4,:);

%% Determinación de componentes
Amplitud=12;
t_inic=0.01;

[val lugar]=min(abs(t_inic-t_d));
y_t1=V_c_d(lugar);
t_t1=t_d(lugar);

[val lugar]=min(abs(2*t_inic-t_d));
y_2t1=V_c_d(lugar);
t_2t1=t_d(lugar);

[val lugar]=min(abs(3*t_inic-t_d));
y_3t1=V_c_d(lugar);
t_3t1=t_d(lugar);

K=12/Amplitud;
k1=(1/Amplitud)*y_t1/K-1;
k2=(1/Amplitud)*y_2t1/K-1;
k3=(1/Amplitud)*y_3t1/K-1;

b=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;

a1=(k1*k2+k3-sqrt(b))/(2*(k1^2+k2));
a2=(k1*k2+k3+sqrt(b))/(2*(k1^2+k2));
beta=(k1+a2)/(a1-a2);

T1_ang=-t_t1/log(a1);
T2_ang=-t_t1/log(a2);
T3_ang=beta*(T1_ang-T2_ang)+T1_ang;
T1(1)=T1_ang;
T2(1)=T2_ang;
T3(1)=T3_ang;
T3_ang=sum(T3/length(T3));
T2_ang=sum(T2/length(T2));
T1_ang=sum(T1/length(T1));

sys_G_ang=tf(K*[T3_ang 1],conv([T1_ang 1],[T2_ang 1]))

%%Calculo de la simulación de nuestros valores
I(1)=0;V_c(1)=0;Ve=0;

##A=[-R/L -1/L; 1/C 0]
##B=[1/L; 0]
##Ct=[0 1]
##D=[0]
##X=[0 ; 0];
##u(1)=0;


%%Impresión de valores

#subplot(3,1,1);plot(t_d,I_d); title("Corriente");grid on;hold on;
#subplot(3,1,2);
plot(t_d,V_c_d); title("Tension del capacitor");grid on;hold on;step(Amplitud*sys_G_ang,
"r",0.2);
#subplot(3,1,3);plot(t_d,Ve_d); title("Voltaje de entrada");grid on;hold on;
