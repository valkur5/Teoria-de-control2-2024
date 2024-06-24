##Ejemplo aureliano
clc;clear all;
w=3;a=0.05;b=5;c=100;
%TIEMPOS: de muestreo, de simulación, de Euler y su integración
Ts=1e-2;KMAX=80000;Veces_Euler=100;dt=Ts/Veces_Euler;t_d=(0:KMAX)*Ts;
TamanioFuente=12;
Ref1=0.1;Ref2=500;
%Condiciones iniciales
alfa(1)=0; fi_p(1)=0; fi(1)=0; h(1)=0;
% Expresion matricial
Mat_Ac=[-a a 0 0;0 0 1 0; w^2 -w^2 0 0; c 0 0 0];
Mat_Bc=[0; 0; b*w^2; 0];
Mat_C=[0 0 1 0;0 0 0 1]; I=eye(4);
H=[0;0;0;0]; d_tao=Ts/100;tao=0;
for hh=1:100
dH=expm(Mat_Ac*tao)*Mat_Bc*d_tao;
H=H+dH;
tao=tao+d_tao;
end
Mat_B=H
Mat_A=expm(Mat_Ac*Ts)
Mat_Aa=[Mat_A,zeros(4,2);-Mat_C*Mat_A, eye(2)];
Mat_Ba=[Mat_B;-Mat_C*Mat_B];
Mat_Ma=[Mat_Ba Mat_Aa*Mat_Ba Mat_Aa^2*Mat_Ba Mat_Aa^3*Mat_Ba Mat_Aa^4*Mat_Ba];%Matriz Controlabilidad
rango=rank(Mat_Ma);
Qc=diag([1 1e4 1e4 1e0 1e-4 1e-7]); R=1; %Ts=0.01;
%Contrucción del Hamiltoniano para el cálculo del controlador
H=inv([eye(6) Mat_Ba*inv(R)*Mat_Ba'; zeros(6) Mat_Aa'])*[Mat_Aa zeros(6);-Qc eye(6)];
[V,D]=eig(H);MX1X2=[];
for ii=1:12
if abs(D(ii,ii))<1
MX1X2=[MX1X2 V(:,ii)];
end
end
MX1=MX1X2(1:6,1:6); MX2=MX1X2(7:12,1:6);
Pc=real(MX2*inv(MX1)); % [K1,P,E]=dlqr(Mat_A,Mat_B,Q,R); %En Octave
Ka=inv(R+Mat_Ba'*Pc*Mat_Ba)*Mat_Ba'*Pc*Mat_Aa;
K=Ka(1:4);KI=-Ka(5);
aut_controlador=abs(eig(Mat_Aa-Mat_Ba*Ka))
%____________________________OBSERVADOR______________
%Cálculo del Observador de estados
Mat_Adual=Mat_A';Mat_Bdual=Mat_C';Mat_Cdual=Mat_B';
Mat_Qobs=[Mat_C;Mat_C*Mat_A;Mat_C*Mat_A^2;Mat_C*Mat_A^3];
rango_matriz_obs=rank(Mat_Qobs);
Qobs=diag([1e-2 1e-1 1e2 1e1]);Ro=diag([1e1 1e-1]);
%Contrucción del Hamiltoniano para el cálculo del Observador
Ho=inv([eye(4) Mat_Bdual*inv(Ro)*Mat_Bdual'; zeros(4)
Mat_Adual'])*[Mat_Adual zeros(4);-Qobs eye(4)];
[Vo,Do]=eig(Ho);MX1X2=[];
for ii=1:8
if abs(Do(ii,ii))<1
MX1X2=[MX1X2 Vo(:,ii)];
end
end
MX1o=MX1X2(1:4,:); MX2o=MX1X2(5:8,:);
Po=real(MX2o*inv(MX1o)); % [K1o,Po,Eo]=dlqr(Mat_Adual,Mat_Bdual,Qobs,Ro);
Kobs=(inv(Ro+Mat_Bdual'*Po*Mat_Bdual)*Mat_Bdual'*Po*Mat_Adual)';
p_observador=abs(eig(Mat_A-Kobs*Mat_C)); %Verifica polos de observabilidad
t=0;X=[alfa(1);fi(1);fi_p(1);h(1)];
ve1(1)=0;ve2(1)=0;
u_k(1)=0; Xang=[0;0;0;0]; color='r';for ki=2:KMAX
t=[t ki*Ts];
ys=Mat_C*X; %Acá DEBE medirse y.
ve1(ki)=ve1(ki-1)+Ref1-ys(1);
ve2(ki)=ve2(ki-1)+Ref2-ys(2);
u(ki-1)=-Ka*[Xang;ve1(ki);ve2(ki)]; %Con Observador
if u>1
u=1;
end
X=Mat_A*X+Mat_B*u(ki-1);
alfal(ki)=X(1);
fil(ki)=X(2);
fi_pl(ki)=X(3);
hl(ki)=X(4);
Xang=Mat_A*Xang+Mat_B*u(ki-1)+Kobs*(ys-Mat_C*Xang);%Acá se usa y
if hl(ki)>499.99
Ref2=-500;
end
end
u(ki)=u(ki-1);
figure(1);hold on;
subplot(3,2,1);plot(t,hl,color);grid
on;title('Altura','FontSize',TamanioFuente);hold on;
subplot(3,2,2);plot(t,fi_pl,color);grid on;
title('$\dot{\phi_t}$','Interpreter','latex','FontSize',TamanioFuente);hol
d on;
subplot(3,2,3); plot(t,alfal,color);grid
on;title('\alpha_t','FontSize',TamanioFuente);hold on;
subplot(3,2,4);plot(t,fil,color);grid
on;title('\phi_t','FontSize',TamanioFuente);hold on;
subplot(3,1,3);plot(t,u,color);grid on;title('Acción de control','FontSize',TamanioFuente);xlabel('Tiempo en Seg.','FontSize',TamanioFuente);hold on;
%Verificación de la solución con el modelo en tiempo continuo.
T=t(ki);x=[alfa(1);fi(1);fi_p(1);h(1)];
Veces_Euler;u=[];i=1;Ref2=500;
u_k(1)=0;xang=[0;0;0;0];colorc='b';ve1(1)=0;ve2(1)=0;
for ki=2:KMAX
ys=Mat_C*x; %Acá se mide la salida.
ve1(ki)=ve1(ki-1)+Ref1-ys(1);
ve2(ki)=ve2(ki-1)+Ref2-ys(2);
u1(ki)=-Ka*[xang;ve1(ki);ve2(ki)]; %Con Observador
if u1(ki)>1
u1(ki)=1;
end
if abs(u1(ki))<.95
u1(ki)=0;
else
u1(ki)=(abs(u1(ki))-.95)*sign(u1(ki));
end
for kii=1:Veces_Euler
u(i)=u1(ki);
% Ecuaciones diferenciales
alfa_p = a*(fi(i)-alfa(i));
fi_pp = -(w^2)*(fi(i)-alfa(i)-b*u(i));
h_p = c*alfa(i);
% Integraciones por Euler
alfa(i+1) = alfa(i)+alfa_p*dt;
h(i+1) = h(i)+h_p*dt;
fi_p(i+1) = fi_p(i)+fi_pp*dt;
fi(i+1) = fi(i)+fi_p(i)*dt;
i=i+1;
end
if h(i)>499.99
Ref2=-500;end
x=[alfa(i-1); fi(i-1); fi_p(i-1); h(i-1)];
xang=Mat_A*xang+Mat_B*u1(ki)+Kobs*(ys-Mat_C*xang);
end
u(i)=u1(ki);t=(1:i)*dt;figure(1);hold on;
subplot(3,2,1);plot(t,h,colorc);grid on;
subplot(3,2,2);plot(t,fi_p,colorc);grid on;
subplot(3,2,3);plot(t,alfa,colorc);grid on;
subplot(3,2,4);plot(t,fi,colorc);grid on;
subplot(3,1,3);plot(t,u,colorc);grid on;
