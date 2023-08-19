% %% Recorte vientos
% lat=ncread('vwnd.sig995.mon.mean.nc','lat');
% lon=ncread('vwnd.sig995.mon.mean.nc','lon');
% u_wind=ncread('uwnd.sig995.mon.mean.nc','uwnd');
% v_wind=ncread('vwnd.sig995.mon.mean.nc','vwnd');
% aux_lat=find(lat>=-30 & lat<=0);
% aux_lon=find(lon>=180 & lon<=280);
% lat_wind=lat(aux_lat);
% lon_wind=lon(aux_lon);
% u_wind=u_wind(aux_lon,aux_lat,1309:end);
% v_wind=v_wind(aux_lon,aux_lat,1309:end);
% %% Recorte calor latente
% lat=ncread('lhtfl.mon.mean.nc','lat');
% lon=ncread('lhtfl.mon.mean.nc','lon');
% lh=ncread('lhtfl.mon.mean.nc','lhtfl');
% aux_lat=find(lat>=-29.5234 & lat<=-0.9524);
% aux_lon=find(lon>=180 & lon<=279.3750);
% lat_lh=lat(aux_lat);
% lon_lh=lon(aux_lon);
% lh=lh(aux_lon,aux_lat,1309:end);
% %% C'alculo magnitud del viento para la region y fechas pedidas
% for i=1:51
%     for j=1:16
%         mag_wind(i,j,:)=sqrt(u_wind(i,j,:).^2+v_wind(i,j,:).^2);
%     end
% end
% 
% save('variables.mat','mag_wind','lh','lon_lh','lat_lh','lon_wind','lat_wind')

%% Ejercicio 1:
load("variables.mat")
% promedio climatologico de la magnitud del viento
for i=1:51
    for j=1:16
        mean_wind(i,j)=mean(mag_wind(i,j,:));
    end
end
%promedio climatologico del flujo de calor latente:

for i=1:54
    for j=1:15
        mean_lh(i,j)=mean(lh(i,j,:));
    end
end
h=tiledlayout(2,1)
nexttile
pcolor(lon_wind-360,lat_wind,mean_wind')
title('Magnitud del viento','FontSize',16)
ax1=gca;
ax1.FontSize=12
shading flat
c1=colorbar
c1.Label.String='m/s'
c1.Label.FontSize=18
colormap("turbo")

nexttile

pcolor(lon_lh-360,lat_lh,mean_lh')
shading flat
title('Flujo de calor latente','FontSize',16)
c2=colorbar
c2.Label.String='W/m^2'
c2.Label.FontSize=18
colormap("turbo")
xlabel(h,'Longitud','fontsize',20)
ylabel(h,'Latitud','fontsize',20)
title(h,'Promedio climatol\''ogico desde 1980-2012','interpreter','Latex','Fontsize',20)

ax2=gca;
ax2.FontSize=12;
%% Ejercicio 2:
load("variables.mat")
XXw=length(lon_wind);
YYw=length(lat_wind);
%llevar la matriz espacial lon x lat a una columna (lon*lat,1)
%% MATRIZ DE COVARIANZA DE MAGNITUD DEL VIENTO
for i=1:396
    reorden(:,i)=reshape(squeeze(mag_wind(:,:,i)),XXw*YYw,1);
end
% Cada fila es una serie de tiempo en un punto del mapa. Por lo que se
% tienen 49*14 = M series de tiempo de largo 396 = N

% ciclo anual, se calcula la media y desviaci'on est'andar de cada mes del todos los años, y
% luego se resta a los datos originales.
for i=1:12
    med(:,i)=mean(reorden(:,i:12:end),2);
    est(:,i)=std(reorden(:,i:12:end),0,2);
end

% anomalia estandarizada
c=0;
for i=1:33 % anhos entre 1980 y 2012
    for j=1:12
        c=c+1;
        W(:,c)=(reorden(:,c)-med(:,j))./est(:,j);
    end
end
%% MATRIZ DE COVARIANZE FLUJO CALOR LATENTE
XXl=length(lon_lh);
YYl=length(lat_lh);
%llevar la matriz espacial lon x lat a una columna (lon*lat,1)
reorden=[];
for i=1:396
    reorden(:,i)=reshape(squeeze(lh(:,:,i)),XXl*YYl,1);
end
% Cada fila es una serie de tiempo en un punto del mapa. Por lo que se
% tienen 52*15 = M series de tiempo de largo 396 = N

% ciclo anual, se calcula la media y desviaci'on est'andar de cada mes del todos los años, y
% luego se resta a los datos originales.
med=[];
est=[];
for i=1:12 
    med(:,i)=mean(reorden(:,i:12:end),2);
    est(:,i)=std(reorden(:,i:12:end),0,2);
end

% anomalia estandarizada
c=0;
for i=1:33 % años entre 1980 y 2012
    for j=1:12
        c=c+1;
        LA(:,c)=(reorden(:,c)-med(:,j))./est(:,j);
    end
end

F=[LA;W]; % MATRIZ DE COVARIANZA COMBINADA

% [L,A,E,error]=EOF(F); % entra como (M,N)... (espacio,tiempo)... TRUCO
[L,E,A,error]=EOF(F'); 
% defino vector fecha
c=0;
for i=1980:2012
    for j=1:12
        c=c+1;
        fecha(c)=datenum(i,j,1);
    end
end

modos=6
N=396;
figure
plot(L(1:modos),'o','markersize',10,MarkerFaceColor='b')
hold on
plot(L(1:modos)+L(1:modos)*sqrt(2/N),'+r','LineWidth',3,'markersize',10)
plot(L(1:modos)-L(1:modos)*sqrt(2/N),'+r','LineWidth',3,'Markersize',10)
grid on
title('Diagrama esquem\''atico de los primeros 6 valores propios ','Interpreter','latex','FontSize',18)
legend('Valor Propio','Error del valor propio')
%Se concluye que solo 3 modos no est'an mezclados con otros.

%% Ejercicio 3:
EC=load('EC.txt');
E_enso=EC(1206:1601,3);
C_enso=EC(1206:1601,4);
for modo=1:3
% modo=1;
Mw=length(W(:,1));
Mla=length(LA(:,1));

for i=1:Mw
    corrcoef(W(i,:)',A(:,modo));
    rw(i)=ans(1,2);
end
for i=1:Mla
    corrcoef(LA(i,:)',A(:,modo));
    rla(i)=ans(1,2);
end

figure(modo)
subplot(311)
plot(fecha,A(:,modo)/std(A(:,modo))),datetick
axis tight
title(['Componente principal normalizada ' num2str(modo) ' :' num2str(L(modo)/sum(L)*100) '%'],'FontSize',18)
h=line([datenum(1980,1,1) datenum(2012,12,31)],[0 0])
set(h,'color','k')
subplot(312)

contourf(lon_wind-360,lat_wind,reshape(rw,51,16)'),colorbar
caxis([-0.8 0.8])
title('Campo de correlaci\''on del viento con la componente principal','Interpreter','latex',FontSize=18)
subplot(313)
contourf(lon_lh-360,lat_lh,reshape(rla,54,15)'),colorbar
caxis([-0.8 0.8])
title('Campo de correlaci\''on del flujo de calor latente con la componente principal','Interpreter','latex',FontSize=18)

figure(modo+3)
h2=tiledlayout(1,2)
nexttile
plot(fecha,A(:,modo),LineWidth=1.5),datetick
coef=polyfit(fecha,A(:,modo),1);
hold on
plot(fecha,polyval(coef,fecha),"Color",'red',LineWidth=1.5)
axis tight
legend('Componente principal','Tendencia lineal','Fontsize',15)
grid on
title(['Tendencia lineal de la componente principal modo ' num2str(modo)],'Fontsize',18)
xlabel('Años','Fontsize',15)

%montecarlo
for m=1:1000
    princ_rem=remuestreo(A(:,modo));
    coef_mc=polyfit(fecha,princ_rem,1);
    pendiente(m)=coef_mc(1)*365;
end
nexttile
histogram(pendiente)
hold on
xline([prctile(pendiente,97.5) prctile(pendiente,2.5)],'LineWidth',3,'Color','k','LineStyle','--')
xline(coef(1)*365,'LineWidth',3,'Color','r','LineStyle','--')
grid on
legend('','','','Pendiente observada','Fontsize',15)
xlabel('Tasa de cambio/año componente principal','Fontsize',15)
title(['Significancia de la tendencia con 95% de confianza del modo ' num2str(modo)],'Fontsize',18)
end

%% 3 c)
for modo=1:3
    A_m=detrend(A(:,modo));
    figure(modo)
    h=tiledlayout(2,2);
    nexttile
    scatter(A_m,C_enso)
    grid on
    hold on
    plot(A_m,polyval(polyfit(A_m,C_enso,1),A_m),'Color','r')
    correC=corr(C_enso,A_m);
    ylabel('C ENSO','FontSize',20)
    text(1,1,['\rho=' num2str(correC)],'Fontsize',18)
    title(['Diagrama de dispersión modo ' num2str(modo)],'Fontsize',20)
%% ----------------
    nexttile
    for m=1:1000
        princ_rem=remuestreo(A_m);
        C_enso_rem=remuestreo(C_enso);
        corre_mc(m)=corr(princ_rem',C_enso_rem');
    end
    histogram(corre_mc)
    hold on
    xline([prctile(corre_mc,97.5) prctile(corre_mc,2.5)],'LineWidth',3,'Color','k','LineStyle','--')
    xline(correC,'LineWidth',3,'Color','r','LineStyle','--')
    grid on
    legend('','','','\rho observado','Fontsize',15)
    title(['Significancia de \rho'],'Fontsize',20)

%% ----------------
    nexttile
    scatter(A_m,E_enso)
    hold on
    plot(A_m,polyval(polyfit(A_m,E_enso,1),A_m),'Color','r')
    ylabel('E ENSO','FontSize',20)
    correE=corr(E_enso,A_m);
    text(1,1,['\rho=' num2str(correE)],'Fontsize',18)
    xlabel('Componente principal','FontSize',18)
    grid on

%% -----------------
    nexttile
    for m=1:1000
        princ_rem=remuestreo(A_m); %remuestreo
        E_enso_rem=remuestreo(E_enso); %remuestreo
        corre_mc(m)=corr(princ_rem',E_enso_rem');
    end
    histogram(corre_mc)
    hold on
    xline([prctile(corre_mc,97.5) prctile(corre_mc,2.5)],'LineWidth',3,'Color','k','LineStyle','--')
    xline(correE,'LineWidth',3,'Color','r','LineStyle','--')
    grid on
    legend('','','','\rho observado','Fontsize',15)
    xlabel('Coeficiente de correlación','FontSize',18)
end

figure()
geoscatter(-20,-150)
geolimits([-30 0],[-180 -80])
geobasemap('landcover')
hold on
geoplot([0 -30],[-180,-180],'r-','LineWidth',3)
geoplot([-30 -30],[-180, -80],'r-','LineWidth',3)
geoplot([-30 0],[-80,-80],'r-','LineWidth',3)
geoplot([0 0],[-80,-180],'r-','LineWidth',3)


