clear;close all;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%0 = marsh
%1 = natural, fixed channel
%2 = man-made ditch
%3 = isolated pond
%4 = connected pond
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
rng(1)

%initialization
load A
load h
dx=1;
[N,M]=size(A);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Time paramters
dt=1;% years

%General parmeters
rmud=2600;Wmud=0.25;emud=rmud*Wmud;
rorg=1200;Worg=0.06;eorg=rorg*Worg;
forgSUSP=0.14;

slr=1/1000;
amp=3.6/2;
springneap=0.1*amp;
T=12.5/365/24; %in 1/yr

%Inorganic sediment
Co=30/1000;% reference seidment concentration in the flood water [g/l] = [kg/m3]
beta=1/20;%1/m decay of sediment with distance
facseddist=0.3;%strom factor background SSC 

%Creep
Ko=0.1;%m2/yr%soil diffusion with vegetation0.1;%
hch=0.;%channel depth, incluences creep rate
hditcho=1.1;%elevation of ditches. INFLUENCES BOTH CREEP AND DRAINAGE POTENTIAL

%Marsh parameters and organic accretion
Dmin=0;Dmax=amp;
Orgmax=6/1000;

%pond paramteters
Kseed=0.0004;%[#of ponds/m2/yr] seeding
Er=0.015;%pond expansion rate [m/yr]  . The value is 0.02/0.03!!!
dpondo=amp*0.25;% [m] initial depth of new pond compared to existing marsh
zlimfornewpondformation=Dmin;%do not /seed a new pond if there is not a marsh already there!!!
Dbasepond=Dmin;%% [m a.m.s.l.]base layer below which new ponds to do scour
Rateponddeep=3/1000;%active pond deepening 
pondtoglitutto=0;%when you do ponding, do you take the elevation from everything or only from the plant material? [1] is everything. Does not change elevation!!

%Buffer distance from channel for pond drainage
RAD=20;%[m] offset of channel driange for pond drinage

%ditches; elevation loss due to water table lowering due to dithes
DOfac=3/1000;% mm/yr%how much to reduce because of the oxidation. if zero there is no oxidation
DOx=20; %[m]

elevationdynamic=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%Domain boundary
A(1,:)=-1;A(end,:)=-1;A(:,1)=-1;A(:,end)=-1;
bnd=find(A==-1);%find the indexes of all boundary cells
A(bnd)=1;%make the boundary as fake channels that do not evelove BUT ALSO DO NOT DRAIN (sse Ddrn AA=0 later)

%initialize the ditch depth and keep track of it
hditch=A*0+hditcho;; %keep track of the ditch elevation

%Channels and creeks
h(A==1)=hch;
%Ditches
h(A==2)=hditch(A==2);
%Isolated ponds
h(A==3)=0.8;

% %reset the initial ponds
A(A==3)=0;
h=h*0+amp*0.9;

%Initilize the book-keeping
zINRo=0.5;%intial thinkess of mud
zORGo=1;%intial thinkess of plant organic
zORG=A*0+zORGo;
zINR=A*0+zINRo;
Ttot_acc=A*0;
Torg_acc=A*0;
Ttot_cr=A*0;
Torg_cr=A*0;
Ttot_pond=A*0;
Torg_pond=A*0;
Ttot_ditchwl=A*0;
Torg_ditchwl=A*0;
SLRcum=A*0;
hstore=h;
hinitial=h-zORG-zINR;%-zMUDORG;
ditchmaterail=A*0;%keep track when you scour ditches

%Create the matrix with the distance from channels and from ditches
%find distance from connected channels (NOT a ditch, they do not convey sediments!!!)
AA=A;AA(:,1)=0;AA(:,end)=0;AA(1,:)=0;AA(end,:)=0;
[Ddrn,idx]=bwdist(AA==1);Ddrn=double(dx*Ddrn);%distance from natural channels
[Ddrn2,idx2]=bwdist(AA==2);Ddrn2=double(dx*Ddrn2);%distance from ditches
clear AA


%OPTIONAL:load the steady state configuration
%This overwrite all previous geometry
load Asteady;load hsteady;h=hsteady;A=Asteady;


hinitial=h-zORG-zINR;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1]);tic
%assume that simulation sstart at year 1850
timeditch=80;%ditches dug in 1930
tmax=160;%simulation ends at year 2010
time=[1:tmax]*dt;
plotint=1;%only for visualization!!!
tpondanlaysis=plotint;%%does not affect simulations
k=0;SS=0;

for t=1:tmax;tic
    if t>50;slr=2.9/1000;end%slr accelreates at year 1900
    
%OPTIONAL: Make the ditches%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if t==timeditch;
    DI=A*0;
    load ditch2;DI(ditch2)=1;load ditch3;DI(ditch3)=1;  
    
    %scour the ditches
    A((A==0 | A==3) & DI==1)=2;%make ditches only where there is an existing marsh (A==0)
    ditchmaterail=A*0;ditchmaterail(A==2)=h(A==2)-hditch(A==2);
    h(A==2)=hditch(A==2);

    %update the drianage mask to allow for ditches to drain ponds!!!
    AA=A;AA(:,1)=0;AA(:,end)=0;AA(1,:)=0;AA(end,:)=0;
    [Ddrn2,idx2]=bwdist(AA==2);Ddrn2=double(dx*Ddrn2);clear AA
end


%Pond reconnection
G=(A==0 & h<Dmin);%the drowned marsh (elevation below the vegetation limit) that is not a pond
CC = bwconncomp(A==3,4); %remember that the marsh is A==0  %CC = bwconncomp(A==3,4);
for i=1:CC.NumObjects;%find the ponds that interesect a channel (A==1)
xxx=CC.PixelIdxList{i};
       %1)drainage due to natrual channels
        a=find(Ddrn(xxx)<=RAD);  %is there any part of the pond in the drianage zone?
        if length(a)>0   
        A(xxx)=0;%A(xxx(A(xxx)==3))=0;%drain the whole pond!    
        end
       %2)drainage due to ditches
        a=find(Ddrn2(xxx)<=RAD);  %is there any part of the pond in the drianage zone?
        if length(a)>0  
        A(xxx)=0;%A(xxx(A(xxx)==3))=0;%drain the whole pond!
        end
       %3)drainage of ponds due to intersecting a drowned marsh, which is assumed to be connected
        pxxx=xxx(h(xxx)<Dmin);%only select the subset of isolated ponds that are below the vegetation limit
        if length(pxxx)>0
            %check if any of the neighouring ponds is a drowned marsh (as identified above)
             a=find(  G(max(1,min([pxxx+1 pxxx-1 pxxx+N pxxx-N],N*M)))==1 );%|   G(circshift(CC.PixelIdxList{i},[0 -1]))==1 |  G(circshift(CC.PixelIdxList{i},[1 0]))==1 |  G(circshift(CC.PixelIdxList{i},[-1 0]))==1  );  %is there any part of the pond in the drianage zone?
             if length(a)>0 
             A(pxxx)=0;%A(xxx(A(xxx)==3))=0;%drain the whole pond!
             end
        end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ELEVATION DYNAMICS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if elevationdynamic==1;
    
%Inorganic depositon
C=Co*(facseddist+(1-facseddist)*exp(-beta*Ddrn));%C=Co*(facseddist+(1-facseddist)*exp(-beta*min(Ddrn,Ddrn2)));
dm=max(0,amp-h);%water depth at high tide % dmM=ampM-h; % dmM2=ampM2-h; % hydroperiod=(max(0,min(2*amp,dm))+max(0,min(2*ampM,dmM))+max(0,min(2*ampM2,dmM2)))/3;
hydroperiod=max(0,min(2*amp,dm)).*min(1,dm/(2*amp*springneap));%the second part modifyes the hyudrperiod to account for spring neap modulation
dep=C.*hydroperiod/T/(rmud*Wmud);
dep(A==1)=0;%no deposition in channels BUT YES ON DITCHES
dep(A==3)=0;%no deposition in isolated ponds

%Organic depositions
Bpeak=4*(Dmax-h).*(h-Dmin)/(Dmax-Dmin)^2;Bpeak(Bpeak<0)=0;
%The connected ponds cannot grow vegetation
Bpeak(A==3 | A==1 | A==2)=0;
org=Bpeak*Orgmax;

P=A;
P(bnd)=4;
P(P==1 | P==2)=0;%channels and ditches
PP=[P(:,2:end) P(:,end)]; FR=(h-[h(:,2:end) h(:,end)]).*(1-(P==0 & PP==3)).*(1-(P==3 & PP==0)).*(P~=4 & PP~=4);%do not creep at the pond boundary!
PP=[P(:,1) P(:,1:end-1)]; FL=([h(:,1) h(:,1:end-1)]-h).*(1-(P==0 & PP==3)).*(1-(P==3 & PP==0)).*(P~=4 & PP~=4);%do not creep at the pond boundary!
PP=[P(2:end,:); P(end,:)];FD=(h-[h(2:end,:); h(end,:)]).*(1-(P==0 & PP==3)).*(1-(P==3 & PP==0)).*(P~=4 & PP~=4);%do not creep at the pond boundary!
PP=[P(1,:); P(1:end-1,:)];FU=([h(1,:); h(1:end-1,:)]-h).*(1-(P==0 & PP==3)).*(1-(P==3 & PP==0)).*(P~=4 & PP~=4);%do not creep at the pond boundary!
cr=(FL-FR+FU-FD)/dx^2*Ko.*(P~=3);%explicit


%%%%%%%%%%%%%%%
    minimo=min(0,min(zINR,zORG));
    zORGi=zORG-minimo+0.00001;
    zINRi=zINR-minimo+0.00001;
    ratio=zORGi./(zORGi+zINRi);
    ratio=A*0+0.5;
%%%%%%%%%%%%%%



%factor to reduce organic accretion with water table drop
 DO=Ddrn2;
 ditchdriangeeffect=DOfac*exp(-DO/DOx); %DOcorrection=Ddrn2<DOx;%max(0,1-Ddrn2/DOx); %ditchdriangeeffect(Ddrn<20)=0;
 ditchdriangeeffect(A==2 | A==1)=0; %hte ditch do not do a lowering of themselves by the water level effect!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main evolution of the elevation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if pondtoglitutto==1;Fporg=ratio;else;Fporg=A*0+1;end

%pond deepening
ponddeepening=(A==3)*Rateponddeep.*(h>Dmin);   

htemp=h;
h=h+(org+dep-slr+cr-ponddeepening-ditchdriangeeffect)*dt;

%Impose fixed depth in the channels and ditches
h(A==1)=hch;%Channels
hpreditch=h;%store ditch elevation before maintanance
h(A==2)=hditch(A==2);%Ditche maintanance
zORG(A==2)=zORG(A==2)-(hpreditch(A==2)-h(A==2));  %ASSUMPTION #1 THAT DITCH MAINTANANCE REMOVED ONLY FROM zORG!!!
end  %END OF ELEVATION DYNAMICS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






CHarea=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Store the values MAPS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if elevationdynamic==1

    SLRcum=SLRcum+(A==0 | A==3 | A==2)*dt*slr;
    
    %update z: elevation + creep       
    zORG=zORG+(org+(A==0 | A==3 | A==2).*cr.*ratio)*dt; %ASSUMPTION #3 THAT CREEP REMOVED FROM BOTH zORG and zINR
    zINR=zINR+(dep+(A==0 | A==3 | A==2).*cr.*(1-ratio))*dt;   %ASSUMPTION #3 THAT CREEP REMOVED FROM BOTH zORG and zINR
    %water table effect
    zORG=zORG-ditchdriangeeffect*dt.*Fporg;%ASSUMPTION #2 THAT DITCH MAINTANANCE REMOVED ONLY FROM zORG!!! 
    %update z: ponding
    zORG=zORG-ponddeepening*dt.*Fporg;
    zINR=zINR-ponddeepening*dt.*(1-Fporg);   
    
    %accretion
    Ttot_acc=Ttot_acc+(org+dep)*dt;                           
    Torg_acc=Torg_acc+(org*eorg+dep*emud*forgSUSP)*dt;
    %creep
    Ttot_cr=Ttot_cr+(A==0 | A==3 | A==2).*cr.*dt;
    Torg_cr=Torg_cr+(A==0 | A==3 | A==2).*cr.*(ratio*eorg +(1-ratio)*emud*forgSUSP)*dt;    
    %active pond deepening.It all consume the organic material  
    Ttot_pond=Ttot_pond+ponddeepening.*dt;
    Torg_pond=Torg_pond+ponddeepening.*(eorg.*Fporg +emud*forgSUSP*(1-Fporg))*dt;
    %ditchwatertable effect
    Ttot_ditchwl=Ttot_ditchwl+ditchdriangeeffect.*dt;
    Torg_ditchwl=Torg_ditchwl+ditchdriangeeffect.*eorg*dt;%ASSUMPTION #2 THAT DITCH MAINTANANCE REMOVED ONLY FROM zORG!!!
    

celle=(A==3 | A==0 | A==2) & Ddrn>=CHarea;    
%Store the values time series
marsharea=sum((A(:)==3 | A(:)==0 | A(:)==2) & Ddrn(:)>=CHarea);
%marsharea=sum(A(:)==3 | A(:)==0 )
Dorg(t)=sum(org(celle))/marsharea;%only make the mass balance on the marsh and ponds
Dinr(t)=sum(dep(celle))/marsharea;%only make the mass balance on the marsh and ponds
SLR(t)=slr;%only make the mass balance on the marsh and ponds

Dcreeplost(t)=+sum(cr(A==1 & Ddrn>=CHarea))/marsharea;%the creep that goes into the channels  %Dcreeplost(t)=-sum(cr(A==0 | A==3))/marsharea;
%ASSUMPTION #3 THAT CREEP REMOVED FROM BOTH zORG and zINR
DcreeplostORG(t)=(sum(cr(A==1 & Ddrn>=CHarea).*ratio(A==1 & Ddrn>=CHarea))*eorg     +sum(cr(A==1 & Ddrn>=CHarea).*(1-ratio(A==1 & Ddrn>=CHarea)))*emud*forgSUSP)      /marsharea;%the creep that goes into the channels
                %how much goes in the ditches and needs to be removed to keep them at the same height
                Dcreeplostditch(t)=+sum(cr(A==2 & Ddrn>=CHarea))/marsharea  +sum(dep(A==2 & Ddrn>=CHarea))/marsharea  -sum(A(:)==2  & Ddrn(:)>=CHarea)*dt*slr/marsharea;%the creep that goes into the channels
                
EPdeepen(t)=sum(ponddeepening(A==3 & Ddrn>=CHarea))/marsharea;
EPdeepenORG(t)=sum(ponddeepening(A==3 & Ddrn>=CHarea).*Fporg(A==3 & Ddrn>=CHarea)*eorg)/marsharea  +sum(ponddeepening(A==3 & Ddrn>=CHarea).*((1-Fporg(A==3 & Ddrn>=CHarea))*emud*forgSUSP))/marsharea;
Morg(t)=sum(org(celle))*eorg/marsharea   +sum(dep(celle)*emud*forgSUSP)/marsharea;
Dditchdriangeeffect(t)=sum(ditchdriangeeffect((A==0 | A==3) & Ddrn>=CHarea))/marsharea;
DditchdriangeeffectORG(t)=sum(ditchdriangeeffect((A==0 | A==3) & Ddrn>=CHarea))*eorg/marsharea;

%keep track of much material rmoved due to maintance
ditchremoval(t)=sum(hpreditch(A==2 & Ddrn>=CHarea)-h(A==2 & Ddrn>=CHarea))/marsharea;
ditchremovalORG(t)=sum(hpreditch(A==2 & Ddrn>=CHarea)-h(A==2 & Ddrn>=CHarea))*eorg/marsharea;%ASSUMPTION #1 THAT DITCH MAINTANANCE REMOVED ONLY FROM zORG!!!


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%New ponds formation
a=find(A==0 & Ddrn>RAD & Ddrn2>RAD);%find stuff that can become a pond!& D>RAD & D>RAD
r=rand(size(A));r=r(a);
seed=a(r<Kseed*dt*dx^2);  %.*(max(0,D(a))/100)
A(seed)=3;%transform to pond
if elevationdynamic==1;
dpond=min(dpondo,max(0,h(seed)-Dbasepond));%how much you can deepen!
EPseed(t)=(sum(dpond))/(dt)/marsharea;%Store elevation lost due to pond formation
EPseedORG(t)=sum(dpond.*(Fporg(seed)*eorg + (1-Fporg(seed))*emud.*forgSUSP))/(dt)/marsharea ;%Store elevation lost due to pond formation
    %Update the material
    zORG(seed)=zORG(seed)-dpond.*Fporg(seed);
    zINR(seed)=zINR(seed)-dpond.*(1-Fporg(seed));    
    %Keep track of organic material pond loss
    Torg_pond(seed)=Torg_pond(seed)+dpond.*Fporg(seed)*eorg;
    Ttot_pond(seed)=Ttot_pond(seed)+dpond;
    h(seed)=h(seed)-dpond;
end

%Isolated pond expansion
%find a cell that is a marsh and has at least one pond nearby. record cell indices
I=(circshift(A,[1 0])==3) + (circshift(A,[-1 0])==3) + (circshift(A,[0 1])==3) + (circshift(A,[0 -1])==3);
a=find(h>Dmin & A==0 & I>0 ); %find the places where the pond can expand to. It need to chew up a marsh (h>Dmin!!!)
r=rand(size(A));r=r(a);
trn=a(r<dt/dx*Er*I(a)); %.*(max(0,D(a))/100)
A(trn)=3;%transform to pond
if elevationdynamic==1;
dpond=min(dpondo,max(0,h(trn)-Dbasepond));%how much you can deepen!
EPexp(t)=(sum(dpond))/(dt)/marsharea;%Store elevation lost due to pond formation
EPexpORG(t)=sum(dpond.*(Fporg(trn)*eorg + (1-Fporg(trn))*emud.*forgSUSP))/(dt)/marsharea ;%Store elevation lost due to pond formation
    %Update the materail
    zORG(trn)=zORG(trn)-dpond.*Fporg(trn);
    zINR(trn)=zINR(trn)-dpond.*(1-Fporg(trn)); 
    %Keep track of pond loss
    Torg_pond(trn)=Torg_pond(trn)+dpond.*Fporg(trn)*eorg;
    Ttot_pond(trn)=Ttot_pond(trn)+dpond;
    h(trn)=h(trn)-dpond;
end


if elevationdynamic==1;
%Keep track of elevation change
deltah(t)=(sum(hstore(A==0 | A==3)-h(A==0 | A==3)))/(dt);%per unit of time
pozze(t)=sum(Bpeak(:)==0 & A(:)==3)/marsharea;
mudflats(t)=sum(Bpeak(:)==0 & A(:)==0)/marsharea;

totalorganic(t)=sum(zORG(celle)*eorg )  + sum(zINR(celle)*emud*forgSUSP ) ;
totalelevation(t)=sum(h(celle) );
end


%%%PLOTTING DURING THE CYCLE!!!
if mod(t,plotint)==0
   
subplot(1,2,1);
imagesc(h);caxis([1 amp]);axis equal
subplot(1,2,2);
G=A;
G(A==1)=0;%channels
G(A==2)=1.4;%ditches
G(A==0)=1.9;%marsh
G(A==3)=0.8;%ponds
G(A==0 & h<Dmin)=2.5;%the drowned marsh (elevation below the vegetation limit) that is not a pond
G(A==3 & h<Dmin)=3.5;%A pond with elevation below the vegeation limit
imagesc(G);axis equal;caxis([0 4]);colormap('jet');
title(strcat(num2str(0*1850+t*dt),' years   ',num2str(slr*1000)))
pause(0.001)

end

toc
end;%end of time cycle







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plots at the end of the simulation

figure
subplot(2,1,1)
diffel=(totalelevation-[NaN totalelevation(1:end-1)])/marsharea/dt*1000;
plot(time,(Dorg+Dinr)*1000,time,(Dorg+Dinr-Dcreeplost-EPexp-EPseed-EPdeepen-Dditchdriangeeffect-ditchremoval)*1000,time,Dorg*1000,time,Dinr*1000,time,-Dcreeplost*1000,time,-(EPexp+EPseed+EPdeepen)*1000,time,SLR*1000,'k');%,time,Dcreeplost)
hold on;plot(time,-Dditchdriangeeffect*1000,'-.m')
hold on;plot(time,-Dcreeplostditch*1000,'-.g')
hold on;plot(time,-ditchremoval*1000,'--k')
hold on;plot(time,diffel+SLR*1000,'--r')%--*k
ylim([-5 8])
legend('Gross accretion','Net accretion','Org','Inr','bank creep','pond loss','SLR','ditch subsidence','material removed from ditches')
ylabel('Spatially averaged vertical rates [mm/yr]')
xlim([time(1) time(end)])

subplot(2,1,2)
diffORG=(totalorganic-[NaN totalorganic(1:end-1)])/marsharea/dt*1000;
plot(time,(Morg)*1000,time,(Morg-DcreeplostORG-EPexpORG-EPseedORG-EPdeepenORG-DditchdriangeeffectORG-ditchremovalORG)*1000,time,-DcreeplostORG*1000,time,-(EPexpORG+EPseedORG+EPdeepenORG)*1000)
hold on;plot(time,diffORG,'o')%,time,SLR*1000,'k');%,time,Dcreeplost)
hold on
plot(time,-DditchdriangeeffectORG*1000,'-.m')
ylim([-0 200]);
ylim([-700 700])
legend('Gross Org accretion (Org)','Net Org accretion','bank creep Org','pond loss Org','NetAccretionOrganic')
xlabel('Time [years]')
ylabel('Spatially averaged vertical rates [g/m^2/yr]')


%Organic ratio by VOLUME
% ratio=zORG./(zORG+zINR);ratio(A==1)=NaN;
% ratioI=org./(dep+org);ratio(A==1)=NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Depth average values
LOI=(max(0,zORG)*eorg+zINR*forgSUSP*emud)./(max(0,zORG)*eorg+zINR*emud);%atio(isnan(ratioI))=1;ratio(A==1)=NaN;
bulkdensity=(max(0,zORG)*eorg+zINR*emud)./(max(0,zORG)+zINR);

%Instantaneous values
LOII=(org*eorg+dep*forgSUSP*emud)./(org*eorg+dep*emud);%atio(isnan(ratioI))=1;ratio(A==1)=NaN;
bulkdensityI=(org*eorg+dep*emud)./(org+dep);

LOII(A==3)=NaN;
bulkdensityI(A==3)=NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,2,1);
IM=imagesc(LOI);axis equal
colormap('jet');caxis([0.1 0.5]);colorbar
subplot(2,2,2);
IM=imagesc(bulkdensity);axis equal
colormap('jet');colorbar;%IM=imagesc(bulkdensity.*ratio);colorbar;colormap('jet');colorbar;
caxis([100 650])

figure
subplot(2,1,1);
IM=imagesc(LOII);axis equal;set(IM,'alphadata',~(A==3));%colorbar%colormap('jet');
colormap('jet');caxis([0.1 0.4]);colorbar;
subplot(2,1,2)
IM=imagesc(bulkdensityI);axis equal;set(IM,'alphadata',~(A==3));%colorbar%colormap('jet');colorbar;colormap('jet');%IM=imagesc(bulkdensityI.*ratioI);colorbar;colormap('jet');
colormap('jet');caxis([100 500]);colorbar


figure;imagesc(dep*1000);axis equal;caxis([0 20])
figure;imagesc(org*1000);axis equal;caxis([0 6])

figure;imagesc(dep./(dep+org));colormap('jet');caxis([0 1])

figure
imagesc((dep+org-ditchdriangeeffect)*1000);colormap('jet');caxis([0 10]);%colorbar


mean(h((A==0 | A==3) & Ddrn>CHarea))
median(h((A==0 | A==3) & Ddrn>CHarea))
s=(h((A==0 | A==3) & Ddrn>CHarea));prctile(s,90)
figure;hist(s,100);xlim([-0.25 1.6])



