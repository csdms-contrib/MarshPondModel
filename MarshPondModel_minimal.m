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
hstore=h;

%Create the matrix with the distance from channels and from ditches
%find distance from connected channels (NOT a ditch, they do not convey sediments!!!)
AA=A;AA(:,1)=0;AA(:,end)=0;AA(1,:)=0;AA(end,:)=0;
[Ddrn,idx]=bwdist(AA==1);Ddrn=double(dx*Ddrn);%distance from natural channels
[Ddrn2,idx2]=bwdist(AA==2);Ddrn2=double(dx*Ddrn2);%distance from ditches
clear AA


%OPTIONAL:load the steady state configuration
%This overwrite all previous geometry
load Asteady;load hsteady;h=hsteady;A=Asteady;










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1]);
%assume that simulation starts at year 1850
timeditch=80;%ditches dug in 1930
tmax=160;%simulation ends at year 2010
time=[1:tmax]*dt;
plotint=1;%only for visualization!!!
tpondanlaysis=plotint;%%does not affect simulations
k=0;SS=0;

for t=1:tmax;
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


%factor to reduce organic accretion with water table drop
 DO=Ddrn2;
 ditchdriangeeffect=DOfac*exp(-DO/DOx); %DOcorrection=Ddrn2<DOx;%max(0,1-Ddrn2/DOx); %ditchdriangeeffect(Ddrn<20)=0;
 ditchdriangeeffect(A==2 | A==1)=0; %hte ditch do not do a lowering of themselves by the water level effect!

%Pond deepening
ponddeepening=(A==3)*Rateponddeep.*(h>Dmin);   
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main evolution of the elevation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Total elevation change
h=h+(org+dep-slr+cr-ponddeepening-ditchdriangeeffect)*dt;

%Impose fixed depth in the channels and ditches
h(A==1)=hch;%Channels
hpreditch=h;%store ditch elevation before maintanance
h(A==2)=hditch(A==2);%Ditche maintanance
end  %END OF ELEVATION DYNAMICS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%New ponds formation
a=find(A==0 & Ddrn>RAD & Ddrn2>RAD);%find stuff that can become a pond!& D>RAD & D>RAD
r=rand(size(A));r=r(a);
seed=a(r<Kseed*dt*dx^2);  %.*(max(0,D(a))/100)
A(seed)=3;%transform to pond
if elevationdynamic==1;
dpond=min(dpondo,max(0,h(seed)-Dbasepond));%how much you can deepen!
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
    h(trn)=h(trn)-dpond;
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
end%end of plotting


end;%end of time cycle







