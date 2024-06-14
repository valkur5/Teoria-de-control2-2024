## Copyright (C) 2024 valkur
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.


## Author: valkur <valkur@dialen>
## Created: 2024-04-06

function tf_aprox = Chen_RLC (Amplitud, t, y, t_inic, retardo, last)

[val lugar]=min(abs((t_inic+retardo)-t));
y_t1=y(lugar);
t_t1=t(lugar);

[val lugar]=min(abs(2*t_inic+retardo-t));
y_2t1=y(lugar);
t_2t1=t(lugar);

[val lugar]=min(abs(3*t_inic+retardo-t));
y_3t1=y(lugar);
t_3t1=t(lugar);

K=last/Amplitud;
k1=(1/Amplitud)*y_t1/K-1;
k2=(1/Amplitud)*y_2t1/K-1;
k3=(1/Amplitud)*y_3t1/K-1;

b=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;

a1=(k1*k2+k3-sqrt(b))/(2*(k1^2+k2))
a2=(k1*k2+k3+sqrt(b))/(2*(k1^2+k2))

beta=(k1+a2)/(a1-a2) %%Beta para RLC

T1_ang=real(-t_inic/log(a1));
T2_ang=real(-t_inic/log(a2));
T3_ang=real(beta*(T1_ang-T2_ang)+T1_ang);
T1(1)=T1_ang;
T2(1)=T2_ang;
T3(1)=T3_ang;
T3_ang=sum(T3/length(T3));
T2_ang=sum(T2/length(T2));
T1_ang=sum(T1/length(T1));

plot([t_t1 t_2t1 t_3t1], [y_t1 y_2t1 y_3t1],'o'), hold on
plot(t,y)
##[t_t1 t_2t1 t_3t1]
##[y_t1 y_2t1 y_3t1]

tf_aprox=tf(K*[T3_ang 1],conv([T1_ang 1],[T2_ang 1]));


endfunction
