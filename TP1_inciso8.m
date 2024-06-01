##TP1, inciso 8
pkg load control signal io symbolic;
clc; clear all; close all;

%Definición de variables útiles
h=1e-7;
tiempo=(.1/h);
t=0:h:(tiempo*h);
i=1;

%Definición de Constantes

data=transpose(csvread("Curvas_Medidas_Motor.csv"));
t_d=data(1,:);
wr_d=data(2,:);
ia_d=data(3,:);
Va_d=data(4,:);
TL_d=data(5,:);

%cálculo con chen para va
t_inic_va=0.000089;
retardo_va=0.0351;
Amplitud_va=12;

sys_G_ang1= Chen(Amplitud_va, t_d, wr_d, t_inic_va,retardo_va, 198.248802246817);
sys_G_ang1.num{1}(1)=0
[ychen ,tchen,ent]=lsim(sys_G_ang1, Va_d ,t_d, [0,0]);

%Cálculo de chen para TL
t_inic_TL=0.00036;
retardo_TL=0.186300000000042;
Amplitud_TL=1.03e-03;

sys_G_ang2= Chen_TL(Amplitud_TL, t_d, wr_d, t_inic_TL, retardo_TL, 33.39);
sys_G_ang2.den{1}=sys_G_ang1.den{1}
[y2chen ,t2chen,ent2]=lsim(sys_G_ang2, TL_d,t_d, [0,0]);

##Una vez teniendo las gráficas y las funciones de transferencia, podemos determinar los componentes como
RA=12/max(ia_d)
J=(sys_G_ang2.den{1}(2))/(sys_G_ang2.num{1}(2))
Km=1/sys_G_ang1.num{1}(1)
Ki=RA/(Km*sys_G_ang2.num{1}(2))
LAA=sys_G_ang2.num{1}(1)*Ki*Km
Bm=0

%Las matrices son:
A=[-RA/LAA -Km/LAA 0;Ki/J -Bm/J 0;0 1 0]
B=[1/LAA 0;0 -1/J ;0 0]
C_t=[0 1 0]
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
Ko=(fliplr(alfa(2:4)-c_ai(2:4))*inv(T))';

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
  ia_ob=X_hat(1);
end
toc

figure 2;
subplot(4,1,1),plot(t,wr),hold on, grid on, title("Velocidad angular");
subplot(4,1,2),plot(t,ia),grid on, title ("Corriente de armadura");
subplot(4,1,3),plot(t,tita),grid on, title("tita");
subplot(4,1,4),plot(t,TL),grid on; title("carga");

##
##figure 1;
##plot(t_d,wr_d),hold on, plot(t_d,ychen-y2chen), title("comparacion del sistema real con el aproximado"),legend("real","aproximado con chen"); grid on;

##figure 2;
##subplot(4,1,1);plot(t_d,wr_d);grid on; title("velocidad angular"); hold on;
##plot(t_d,ychen-y2chen);
##subplot(4,1,2); plot(t_d,ia_d);grid on; title("Corriente de armadura");
##subplot(4,1,3); plot(t_d,Va_d);grid on; title("Entrada de Va");
##subplot(4,1,4); plot(t_d,TL_d);grid on; title("Entrada de TL");
