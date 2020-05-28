% Modelo SIRD con estacionalidad para Lima provincia
% @ktakahashig
% 27/05/2020
%
clear

dt=1; % Paso de tiempo (dias)
Nt=365; % Numero de pasos de tiempo
I0=1000; % Numero inicial de infectados
gamma=0.2; % Tasa de remocion (1/dia)
R0=2.3; % Numero de reproduccion base
pais='peru'; % Pais de matriz de contactos
edades=[0:5:75]'; ncat=length(edades); % Intervalos de edades
dmu=0.26e-2; % Amplitud de ciclo anual de mu
peso=[1 0.5 0.5 0.5]; % Pesos para las matrices de contacto de
                      % 1) hogar, 2) trabajo, 3) escuela, 4) otros
% ~ CFR Italia 
% Adaptado de https://jamanetwork.com/journals/jama/article-abstract/2763667
%     0    5    10   15   20   25   30 35   40    45
ifr=[0.02;0.02;0.03;0.03;0.04;0.04;0.2;0.3; 0.35 ; 0.5 ;
     0.7; 1.3; 3; 4.5; 10 ;20]/100;
%     50  55   60  65  70  75
%
% Matrices de contacto (POLYMOD, https://doi.org/10.1371/journal.pcbi.1005697)
% (Numero de contactos/dia)
home=load(['MU_home.' pais '.txt']); % En el ambito del hogar
work=load(['MU_work.' pais '.txt']); % En el ambito laboral
school=load(['MU_school.' pais '.txt']); % En la escuela, universidad, etc
other=load(['MU_other_locations.' pais '.txt']); % Otros (?)
%
all=peso(1)*home+peso(2)*work+peso(3)*school+peso(4)*other;


% Distribucion de edades (INEI, censo 2017)
% https://censos2017.inei.gob.pe/redatam/
% Pasar de cuentas por año a los intervalos usados
% y Fallecidos confirmados COVID (CDC) en el dpto de Lima
dat=load('Nedades_Lima_prov.txt'); dumedad=[0:115]'; N=[];
dat2=load('fallecidos_minsa_limadpto.txt'); fallecidos=[];
dat3=load('casos_minsa_limaprov.txt'); casos=[];
for n=1:ncat
   if (n<ncat); 
         ii=find(dumedad>=edades(n)&dumedad<edades(n+1));
         jj=length(find(dat2>=edades(n)&dat2<edades(n+1)));
         kk=length(find(dat3>=edades(n)&dat3<edades(n+1)));
   else; ii=find(dumedad>=edades(n)); 
         jj=length(find(dat2>=edades(n)));
         kk=length(find(dat3>=edades(n)));
   end
   N=[N; sum(dat(ii))];
   fallecidos=[fallecidos; jj];
   casos=[casos; kk];
end
pop=sum(N); % Poblacion total



% Matriz de proxima generacion (Next generation matrix) para derivar mu
[Ni Nj]=meshgrid(N/pop*100,N/pop*100);
alltmp=home+work+school+other; % Todos los contactos para el calculo de mu
C=(alltmp.*Ni./Nj); dum=abs(eig(C)); 
mu=gamma*R0/dum(1); % Trasmisividad x contacto 

clear dat dat2 dat3 alltmp C dum ii jj kk n


% Condiciones iniciales
I=zeros(ncat,1); I=I0*N/pop; R=zeros(ncat,1); 
S=N-I; D=zeros(ncat,1);

% Correr el modelo 
C=(all.*Ni./Nj); dum=abs(eig(all)); 
mu2=mu; Rt=mu*dum(1)/gamma; Rt2=mu2*dum(1)/gamma;
%
X=[S;I;R;D]; X2=[S;I;R;D];
%
for i=1:Nt
   X=[X rk4(X(:,end),dt,mu*all',gamma,ifr,N)]; 
   mu2=mu+dmu*(1-cos(2*pi*i/365.24));
   X2=[X2 rk4(X2(:,end),dt,mu2*all',gamma,ifr,N)]; 
   % Calcular R
   C=(all.*Ni./Nj); dum=abs(eig(all)); 
   Rt=[Rt; mu*dum(1)/gamma];
   Rt2=[Rt2; mu2*dum(1)/gamma];
end

% dividir variables
S=X(1:ncat,:); I=X(ncat+1:2*ncat,:); R=X(2*ncat+1:3*ncat,:); D=X(3*ncat+1:4*ncat,:);
S2=X2(1:ncat,:); I2=X2(ncat+1:2*ncat,:); R2=X2(2*ncat+1:3*ncat,:); D2=X2(3*ncat+1:4*ncat,:);



% FIGURA

figure(1,'papersize',[10 14]),clf
%
subplot(2,2,1)
h=plot(0:Nt,sum(I+R+D)/1e6,0:Nt,sum(I2+R2+D2)/1e6,'r--');
set(h,'linewidth',3)
set(gca,'fontsize',14)
title('a\) Infectados acumulados (I+R+D)','fontsize',16)
xlabel('Dias','fontsize',16)
ylabel('Millones de personas','fontsize',16)
%
subplot(2,2,3)
h=plot(0:Nt,sum(D)/1e3,0:Nt,sum(D2)/1e3,'r--');
set(h,'linewidth',3)
set(gca,'fontsize',14)
title('b\) Fallecidos (D)','fontsize',16)
xlabel('Dias','fontsize',16)
ylabel('Miles de personas','fontsize',16)
%
subplot(2,2,2)
h=plot(0:Nt,sum(I)'/1e3,0:Nt,sum(I2)'/1e3,'r--');
set(h,'linewidth',3)
set(gca,'fontsize',14)
title('c\) Infectados activos (I)','fontsize',16)
ylim([0 600])
legend('Sin estacionalidad','Con estacionalidad','location','northeast')
xlabel('Dias','fontsize',16)
ylabel('Miles de personas','fontsize',16)
%
subplot(2,2,4)
h=plot(1:Nt,diff(sum(D)'),1:Nt,diff(sum(D2)'),'r--');
set(h,'linewidth',3)
set(gca,'fontsize',14)
ylim([0 600])
title('d\) Fallecidos por dia (dD/dt)','fontsize',16)
xlabel('Dias','fontsize',16)
ylabel('Personas/dia','fontsize',16)
%
papersize = get (gcf, "papersize"); border = 0.25;
set (gcf, "paperposition", 
     [border, border, (papersize - 2*border)],...
     'paperorientation','landscape')
print('-dpdf','sird_seas.pdf')


figure(2,'papersize',[10 7]),clf
subplot(2,1,1)
h=plot(0:Nt,100*sum(D)./sum(N-S),0:Nt,100*sum(D2)./sum(N-S2),'r--');
set(h,'linewidth',3)
ylim([0 0.6])
set(gca,'fontsize',14)
title('a\) Ratio fallecidos acum / infectados acum (%)','fontsize',16)
xlabel('Dias','fontsize',16)
ylabel('%','fontsize',16)
%
subplot(2,1,2)
h=plot(0:Nt,sum(D(find(edades>=60),:))./sum(D)*100,0:Nt,sum(D2(find(edades>=60),:))./sum(D2)*100,'r--');
set(h,'linewidth',3)
set(gca,'fontsize',14)
set(h,'linewidth',3)
title('b\) Ratio fallecidos acum mayor o igual a 60 años sobre el total (%)','fontsize',16)
xlabel('Dias','fontsize',16)
ylabel('%','fontsize',16)
ylim([40 90])
%
papersize = get (gcf, "papersize"); border = 0.25;
set (gcf, "paperposition",
     [border, border, (papersize - 2*border)],...
     'paperorientation','portrait')
print('-dpdf','death_ratios.pdf')

figure(3,'papersize',[7 7]),clf
h=plot(0:Nt,Rt,0:Nt,Rt2,'r--',[0 Nt],R0*[1 1],'k:');
set(h,'linewidth',3)
ylim([0 3])
set(gca,'fontsize',14)
title('Numero de reproduccion en tiempo real Rt estimado','fontsize',16)
xlabel('Dias','fontsize',16)
legend('Rt sin estacionalidad','Rt con estacionalidad','Numero de reprod. base Ro','location','southwest')
%
papersize = get (gcf, "papersize"); border = 0.25;
set (gcf, "paperposition",
     [border, border, (papersize - 2*border)],...
     'paperorientation','landscape')
print('-dpdf','Rt.pdf')


xlabs=strvcat('0-4','5-9','10-14','15-19','20-24','25-29','30-34',...
 '35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-');


figure(4,'papersize',[11 6]),clf
bar(2:5:77,ifr,'linewidth',3)
set(gca,'XTick',2:5:77,'XTickLabel',xlabs,'fontsize',12)
ylim([0 0.25])
papersize = get (gcf, "papersize"); border = 0.25;
set (gcf, "paperposition",
     [border, border, (papersize - 2*border)],...
     'paperorientation','landscape')
title('Letalidad referencial de las infecciones de COVID-19 usada en el modelo','fontsize',16)
xlabel('Edad','fontsize',14)
print('-dpdf',['ifr.pdf'])

figure(5,'papersize',[11 6]),clf
bar(2:5:77,N/pop*100,'linewidth',3)
set(gca,'XTick',2:5:77,'XTickLabel',xlabs,'fontsize',12)
papersize = get (gcf, "papersize"); border = 0.25;
set (gcf, "paperposition",
     [border, border, (papersize - 2*border)],...
     'paperorientation','landscape')
title('Estructura de la poblacion de la provincia de Lima (Censo 2017)','fontsize',16)
xlabel('Edad','fontsize',14)
ylabel('%','fontsize',14)
print('-dpdf',['pobl_edad.pdf'])




dia=90;
figure(6,'papersize',[20 12]),clf
subplot(2,2,1)
bar(edades+2.5,(R(:,dia)+I(:,dia))/1e3,'stacked')
xlabel('Edad','fontsize',14)
ylabel('Miles de personas','fontsize',14)
set(gca,'XTick',2:5:77,'XTickLabel',xlabs,'fontsize',11)
title(['a\) Infectados acum. simulados (sin estac.) para el dia ' num2str(dia) ...
  ' (total = ' num2str(round(sum(N-S(:,dia)))) ')'],'fontsize',16)
%
subplot(2,2,2)
bar(edades+2.5,D(:,dia))
set(gca,'XTick',2:5:77,'XTickLabel',xlabs,'fontsize',11)
xlabel('Edad','fontsize',14)
ylabel('Personas','fontsize',14)
title(['b\) Fallecidos simulados (sin estac.) para el dia ' num2str(dia) ...
    ' (total = ' num2str(round(sum(D(:,dia)))) ')'],'fontsize',16)
%
subplot(2,2,3)
bar(edades+2.5,casos)
set(gca,'XTick',2:5:77,'XTickLabel',xlabs,'fontsize',11)
title(['c\) Casos confirmados COVID-19 (MINSA) al 25/05 en prov Lima (total = ' ...
     num2str(sum(casos)) ')'],'fontsize',16)
xlabel('Edad','fontsize',14)
ylabel('Personas','fontsize',14)
%
subplot(2,2,4)
bar(edades+2.5,fallecidos)
set(gca,'XTick',2:5:77,'XTickLabel',xlabs,'fontsize',11)
title(['d\) Fallecidos confirmados COVID-19 (MINSA) al 25/05 en dpto Lima (total = ' ...
     num2str(sum(fallecidos)) ')'],'fontsize',16)
xlabel('Edad','fontsize',14)
ylabel('Personas','fontsize',14)
%
papersize = get (gcf, "papersize"); border = 0.25;
set (gcf, "paperposition",
     [border, border, (papersize - 2*border)],...
     'paperorientation','landscape')
print('-dpdf',['fall_edad.pdf'])

