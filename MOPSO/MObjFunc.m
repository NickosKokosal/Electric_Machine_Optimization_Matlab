function  opti = MObjFunc(X)
% clear all; % close all;
% enter here the path corresponding to the installation directory of FEMM
addpath C:\femm42\mfiles
savepath; openfemm(1); newdocument(0);

deg        = pi/180;   % conversion factor deg rad
% Geometric data (in mm)
length     = 1085;     % active axial length 
airgap     = 1.6;      % Air gap
Rbor       = 400;      % Radius bore
beta       = 0.66;     % pole-pitch ratio
betad      = X(1);     % tooth to slot ratio %betad=0.42;
ThickyokeR = X(2);     % rotor yoke thickness %ThickyokeR=20.4;
ThickyokeS = X(3);     % stator yoke thickness %ThickyokeS=20.4;
Thickmagnet= 1.6;      % magnet thickness
Dslot      = 29.9;     % Slot depth
np         = 8;        % number of pole pairs
Pbec       = 3;        % thickness of the nozzles
ferm       = 0.5;      % slot closing coef
nd         = 2*np;     % ratio of machine to field of study
% Current
n          = 80;       % number of conductors
Nt         = n/(2*np); % number of conductors per slot 
npas       = 30;       % number of rotation steps in half a period set to 30
dtpas      = 180/(npas*np);
tta0       = dtpas/2; 
Al         = 60000;    % Al
I          = 0;        % Current
Nn         = 753;      % Speed
kr         = 0.6;      % kr

% calculation of the main radius
RintCR     = Rbor-airgap-Thickmagnet-ThickyokeR;
RextCR     = Rbor-airgap-Thickmagnet; 
RextA      = Rbor-airgap;
Rag1       = RextA+0.9*airgap/2;
Rag2       = RextA+1.1*airgap/2;
Rag0       = RextA+airgap/2;
Rfe        = Rbor+Pbec+Dslot;
Rde        = Rbor+Pbec;
RextS      = Rbor+Dslot+ThickyokeS+Pbec;
Rme        =(Rbor+Pbec+Rfe)/2;
RmCR       =(RintCR+RextCR)/2;
StepTooth  = 180/(3*np);

% ddefinition of materials and circuits
mi_probdef(0,'millimeters','planar', 1.e-8,length,30);
mi_addcircprop('A+',sqrt(2)*I*cos(-5*np*dtpas-0*pi/3),1);
mi_addcircprop('B+',sqrt(2)*I*cos(-5*np*dtpas-2*pi/3),1);
mi_addcircprop('C+',sqrt(2)*I*cos(-5*np*dtpas-4*pi/3),1);
mi_addcircprop('A-',-sqrt(2)*I*cos(-5*np*dtpas-0*pi/3),1);
mi_addcircprop('B-',-sqrt(2)*I*cos(-5*np*dtpas-2*pi/3),1);
mi_addcircprop('C-',-sqrt(2)*I*cos(-5*np*dtpas-4*pi/3),1);
mi_getmaterial('Air');
mi_addmaterial('Coil',1,1,0,0,58,0,0,1,3,0,0,1,7.77);
mi_getmaterial('Pure Iron');
mi_getmaterial('NdFeB 32 MGOe');
%definition of boundary conditions
mi_addboundprop('A=0', 0, 0, 0, 0, 0, 0, 0, 0, 0);
mi_addboundprop('Ap1', 0, 0, 0, 0, 0, 0, 0, 0, 5);
mi_addboundprop('Ap2', 0, 0, 0, 0, 0, 0, 0, 0, 5);
mi_addboundprop('Ap3', 0, 0, 0, 0, 0, 0, 0, 0, 5);
mi_addboundprop('Ap4', 0, 0, 0, 0, 0, 0, 0, 0, 5);
mi_addboundprop('Ap5', 0, 0, 0, 0, 0, 0, 0, 0, 5);
mi_addboundprop('Ap6', 0, 0, 0, 0, 0, 0, 0, 0, 5);
% rotor hThickmagnetd design
mi_addnode(RintCR, 0);
mi_addnode(RintCR*cos(2*pi/(2*np)),RintCR*sin(2*pi/(2*np)));
mi_addarc(RintCR, 0, RintCR*cos(2*pi/(2*np)),RintCR*sin(2*pi/(2*np)),360/(2*np),1);
mi_addnode(RextCR, 0);
mi_addnode(RextCR*cos(((1-beta)/2)*2*pi/(2*np)), RextCR*sin(((1-beta)/2)*2*pi/(2*np)));
mi_addnode(RextCR*cos(2*pi/(2*np)-((1-beta)/2)*2*pi/(2*np)), RextCR*sin(2*pi/(2*np)-((1-beta)/2)*2*pi/(2*np)));
mi_addnode(RextCR*cos(2*pi/(2*np)),RextCR*sin(2*pi/(2*np)));
mi_addarc(RextCR, 0, RextCR*cos(((1-beta)/2)*2*pi/(2*np)), RextCR*sin(((1-beta)/2)*2*pi/(2*np)),((1-beta)/2)*2*180/(2*np) ,1);
xint1=RextCR*cos(((1-beta)/2)*2*pi/(2*np));
yint1=RextCR*sin(((1-beta)/2)*2*pi/(2*np)); 
xint2=RextCR*cos(2*pi/(2*np)-((1-beta)/2)*2*pi/(2*np));
yint2=RextCR*sin(2*pi/(2*np)-((1-beta)/2)*2*pi/(2*np));
mi_addarc(xint1,yint1,xint2,yint2,beta*360/(2*np),1);
xint3=RintCR*cos(2*pi/(2*np));
yint3=RintCR*sin(2*pi/(2*np));
xo3=RextCR*cos(2*pi/(2*np));
yo3=RextCR*sin(2*pi/(2*np));
mi_addarc(xint2,yint2,xo3,yo3,((1-beta)/2)*2*180/(2*np),1);
mi_addsegment(RintCR, 0, RextCR, 0);
mi_addsegment(xint3,yint3,xo3,yo3);

% design of magnets 
xr31=RextA*cos(((1-beta)/2)*2*pi/(2*np));
yr31=RextA*sin(((1-beta)/2)*2*pi/(2*np));
xr32=RextA*cos(2*pi/(2*np)-((1-beta)/2)*2*pi/(2*np));
yr32=RextA*sin(2*pi/(2*np)-((1-beta)/2)*2*pi/(2*np));
mi_addnode(xr31,yr31);
mi_addnode(xr32,yr32);
mi_addarc(xr31,yr31,xr32,yr32,beta*360/(2*np),1);
mi_addsegment(xint1,yint1,xr31,yr31);
mi_addsegment(xint2,yint2,xr32,yr32);

% drawing of the air rotor part
mi_addnode(Rag1,0);
mi_addsegment(RextCR,0,Rag1,0); 
mi_addnode(Rag1*cos(2*pi/(2*np)),Rag1*sin(2*pi/(2*np)));
mi_addsegment(xo3,yo3,Rag1*cos(2*pi/(2*np)),Rag1*sin(2*pi/(2*np)));

%assignment of rotor geometry to group 2 and CL rotor
%nodes
mi_selectnode(RintCR, 0);
mi_selectnode(RintCR*cos(2*pi/(2*np)),RintCR*sin(2*pi/(2*np)));
mi_selectnode(RextCR, 0);
mi_selectnode(RextCR*cos(((1-beta)/2)*2*pi/(2*np)), RextCR*sin(((1-beta)/2)*2*pi/(2*np)));
mi_selectnode(RextCR*cos(2*pi/(2*np)-((1-beta)/2)*2*pi/(2*np)), RextCR*sin(2*pi/(2*np)-((1-beta)/2)*2*pi/(2*np)));
mi_selectnode(RextCR*cos(2*pi/(2*np)),RextCR*sin(2*pi/(2*np)));
mi_selectnode(xr31,yr31);
mi_selectnode(xr32,yr32);
mi_setnodeprop('<None>',2);
mi_clearselected;
%arcs 
%interieur arc
mi_selectarcsegment(RintCR*cos(2*pi/(4*np)),RintCR*sin(2*pi/(4*np)));
mi_setarcsegmentprop(1,'A=0',0,2);
mi_clearselected;
%other arcs
mi_selectarcsegment(RextCR*cos(2*pi/(4*np)),RextCR*sin(2*pi/(4*np)));
mi_selectarcsegment(RextCR*cos(0.5*(0.5-beta/2)*2*pi/(2*np)),RextCR*sin(0.5*(0.5-beta/2)*2*pi/(2*np)));
mi_selectarcsegment(RextCR*cos(0.5*(1.5+beta/2)*2*pi/(2*np)),RextCR*sin(0.5*(1.5+beta/2)*2*pi/(2*np)));
mi_selectarcsegment(RextA*cos(2*pi/(4*np)),RextA*sin(2*pi/(4*np)));
mi_setarcsegmentprop(1,'<None>',0,2);
mi_clearselected;
%segments
mi_selectsegment(0.5*(xint1+xr31),0.5*(yr31+yint1));
mi_selectsegment(0.5*(xint2+xr32),0.5*(yr32+yint2));
mi_setsegmentprop('<None>',0,1,0,2);
mi_clearselected;
%external segments 
mi_selectsegment(0.5*(RintCR+RextCR),0);
mi_selectsegment(0.5*(RextCR+RintCR)*cos(2*pi/(2*np)),0.5*(RextCR+RintCR)*sin(2*pi/(2*np)));
mi_setsegmentprop('Ap1',0,1,0,2);
mi_clearselected;

mi_selectsegment(0.5*(Rag1+RextCR),0);
mi_selectsegment(0.5*(Rag1+RextCR)*cos(2*pi/(2*np)),0.5*(Rag1+RextCR)*sin(2*pi/(2*np)));
mi_setsegmentprop('Ap2',0,1,0,2);
mi_clearselected;


%material assignment rotor
%yoke
RmCr=(RintCR+RextCR)/2;
mi_addblocklabel(RmCr*cos(pi/(2*np)),RmCr*sin(pi/(2*np)));
mi_selectlabel(RmCr*cos(pi/(2*np)),RmCr*sin(pi/(2*np)));
mi_setblockprop('Pure Iron', 1, 0, '<None>', 0, 2, 1);
mi_clearselected;
%magnet
RmA=(RextA+RextCR)/2;
mi_addblocklabel(RmA*cos(pi/(2*np)),RmA*sin(pi/(2*np)));
mi_selectlabel(RmA*cos(pi/(2*np)),RmA*sin(pi/(2*np)));
mi_setblockprop('NdFeB 32 MGOe', 1, 0, '<None>', 180/(2*np), 2, 1);
mi_clearselected;
%air gap 
Rma=(RextA+Rag1)/2;
mi_addblocklabel(Rma*cos(pi/(2*np)),Rma*sin(pi/(2*np)));
mi_selectlabel(Rma*cos(pi/(2*np)),Rma*sin(pi/(2*np)));
mi_setblockprop('Air', 1, 0, '<None>', 0, 2, 1);
mi_clearselected;

% drawing of the air stator part
mi_addnode(Rag2,0);
mi_addnode(Rbor,0);
mi_addsegment(Rag2,0,Rbor,0); 
mi_addnode(Rag2*cos(2*pi/(2*np)),Rag2*sin(2*pi/(2*np)));
mi_addnode(Rbor*cos(2*pi/(2*np)),Rbor*sin(2*pi/(2*np)));
mi_addsegment(Rag2*cos(2*pi/(2*np)),Rag2*sin(2*pi/(2*np)),Rbor*cos(2*pi/(2*np)),Rbor*sin(2*pi/(2*np)) );
%stator design 
%slot zone design 
mi_addnode(Rde*cos(0.5*betad*StepTooth*deg),Rde*sin(0.5*betad*StepTooth*deg));
mi_addnode(Rfe*cos(0.5*betad*StepTooth*deg),Rfe*sin(0.5*betad*StepTooth*deg));
mi_addnode(Rde*cos((1-0.5*betad)*StepTooth*deg),Rde*sin((1-0.5*betad)*StepTooth*deg));
mi_addnode(Rfe*cos((1-0.5*betad)*StepTooth*deg),Rfe*sin((1-0.5*betad)*StepTooth*deg));
mi_addnode(Rde*cos(0.5*StepTooth*deg),Rde*sin(0.5*StepTooth*deg));
mi_addnode(Rfe*cos(0.5*StepTooth*deg),Rfe*sin(0.5*StepTooth*deg));

angno=0.5*betad*StepTooth*deg+0.5*ferm*(1-betad)*StepTooth*deg;
mi_addnode(Rde*cos(angno),Rde*sin(angno));
mi_addnode(Rbor*cos(angno),Rbor*sin(angno));
angno1=(1-0.5*betad)*StepTooth*deg-0.5*ferm*(1-betad)*StepTooth*deg;
mi_addnode(Rde*cos(angno1),Rde*sin(angno1));
mi_addnode(Rbor*cos(angno1),Rbor*sin(angno1));

mi_addsegment(Rde*cos((1-0.5*betad)*StepTooth*deg),Rde*sin((1-0.5*betad)*StepTooth*deg),Rfe*cos((1-0.5*betad)*StepTooth*deg),Rfe*sin((1-0.5*betad)*StepTooth*deg));
mi_addsegment(Rde*cos(0.5*betad*StepTooth*deg),Rde*sin(0.5*betad*StepTooth*deg),Rfe*cos(0.5*betad*StepTooth*deg),Rfe*sin(0.5*betad*StepTooth*deg));
mi_addsegment(Rde*cos(angno),Rde*sin(angno), Rbor*cos(angno),Rbor*sin(angno));
mi_addsegment(Rde*cos(angno1),Rde*sin(angno1), Rbor*cos(angno1),Rbor*sin(angno1));
mi_addsegment(Rde*cos(0.5*StepTooth*deg),Rde*sin(0.5*StepTooth*deg),Rfe*cos(0.5*StepTooth*deg),Rfe*sin(0.5*StepTooth*deg));
mi_zoomnatural();

mi_addarc(Rde*cos(angno),Rde*sin(angno),Rde*cos(0.5*StepTooth*deg),Rde*sin(0.5*StepTooth*deg),(1-ferm)*(1-betad)*StepTooth,1);
mi_addarc(Rde*cos(0.5*StepTooth*deg),Rde*sin(0.5*StepTooth*deg),Rde*cos(angno1),Rde*sin(angno1),(1-ferm)*(1-betad)*StepTooth,1);
mi_addarc(Rfe*cos(0.5*betad*StepTooth*deg),Rfe*sin(0.5*betad*StepTooth*deg), Rfe*cos(0.5*StepTooth*deg),Rfe*sin(0.5*StepTooth*deg), betad*StepTooth,1);
mi_addarc(Rfe*cos(0.5*StepTooth*deg),Rfe*sin(0.5*StepTooth*deg), Rfe*cos((1-0.5*betad)*StepTooth*deg),Rfe*sin((1-0.5*betad)*StepTooth*deg), betad*StepTooth,1);
mi_addarc(Rde*cos(0.5*betad*StepTooth*deg),Rde*sin(0.5*betad*StepTooth*deg),Rde*cos(angno),Rde*sin(angno),0.5*ferm*(1-betad)*StepTooth,1);
mi_addarc(Rde*cos(angno1),Rde*sin(angno1),Rde*cos((1-0.5*betad)*StepTooth*deg),Rde*sin((1-0.5*betad)*StepTooth*deg),0.5*ferm*(1-betad)*StepTooth,1);
mi_addarc(Rbor,0,Rbor*cos(angno),Rbor*sin(angno),angno/deg,1);


%slot selection and rotation

mi_selectnode(Rde*cos(0.5*betad*StepTooth*deg),Rde*sin(0.5*betad*StepTooth*deg));
mi_selectnode(Rfe*cos(0.5*betad*StepTooth*deg),Rfe*sin(0.5*betad*StepTooth*deg));
mi_selectnode(Rde*cos((1-0.5*betad)*StepTooth*deg),Rde*sin((1-0.5*betad)*StepTooth*deg));
mi_selectnode(Rfe*cos((1-0.5*betad)*StepTooth*deg),Rfe*sin((1-0.5*betad)*StepTooth*deg));
mi_selectnode(Rde*cos(angno),Rde*sin(angno));
mi_selectnode(Rbor*cos(angno),Rbor*sin(angno));
mi_selectnode(Rde*cos(angno1),Rde*sin(angno1));
mi_selectnode(Rbor*cos(angno1),Rbor*sin(angno1));
mi_selectnode(Rbor*cos(pi/(3*np)),Rbor*sin(pi/(3*np)));
mi_selectnode(Rbor,0);
mi_selectnode(Rde*cos(0.5*StepTooth*deg),Rde*sin(0.5*StepTooth*deg));
mi_selectnode(Rfe*cos(0.5*StepTooth*deg),Rfe*sin(0.5*StepTooth*deg));
mi_setnodeprop('<None>',3);
mi_copyrotate2(0,0,StepTooth,2,0)
mi_clearselected;

mi_selectsegment(Rme*cos((1-0.5*betad)*StepTooth*deg),Rme*sin((1-0.5*betad)*StepTooth*deg));
mi_selectsegment(Rme*cos(0.5*betad*StepTooth*deg),Rme*sin(0.5*betad*StepTooth*deg));
mi_selectsegment(Rme*cos(0.5*StepTooth*deg),Rme*sin(0.5*StepTooth*deg));
mi_selectsegment(0.5*(Rde+Rbor)*cos(angno),0.5*(Rde+Rbor)*sin(angno));
mi_selectsegment(0.5*(Rde+Rbor)*cos(angno1),0.5*(Rde+Rbor)*sin(angno1));
mi_setsegmentprop('<None>',0,1,0,3);
mi_copyrotate2(0,0,StepTooth,2,1);
mi_clearselected;

mi_selectarcsegment(Rfe*cos((1.5-0.5*betad)*StepTooth*deg/2),Rfe*sin((1.5-0.5*betad)*StepTooth*deg/2));
mi_selectarcsegment(Rfe*cos((1+betad)*StepTooth*deg/4), Rfe*sin((1+betad)*StepTooth*deg/4));
mi_selectarcsegment(Rbor*cos(angno/2),Rbor*sin(angno/2));
mi_selectarcsegment(Rde*cos((StepTooth*(1-0.5*betad))*deg/2+angno1/2),Rde*sin((StepTooth*(1-0.5*betad))*deg/2+angno1/2));
mi_selectarcsegment(Rde*cos((0.5*StepTooth)*deg/2+angno1/2),Rde*sin((0.5*StepTooth)*deg/2+angno1/2));
mi_selectarcsegment(Rde*cos((0.5*StepTooth)*deg/2+angno/2),Rde*sin((0.5*StepTooth)*deg/2+angno/2));
mi_selectarcsegment(Rde*cos((0.5*StepTooth*betad)*deg/2+angno/2),Rde*sin((0.5*StepTooth*betad)*deg/2+angno/2));
mi_setarcsegmentprop(1,'<None>',0,3);
mi_copyrotate2(0,0,StepTooth,2,3);
mi_clearselected;

mi_selectarcsegment(Rbor*cos(angno/2),Rbor*sin(angno/2));
mi_setarcsegmentprop(1,'<None>',0,3);
mi_copyrotate2(0,0,angno1/deg,1,3);
mi_clearselected;

mi_selectarcsegment(Rbor*cos(0.5*(angno1+StepTooth*deg)),Rbor*sin(0.5*(angno1+StepTooth*deg)));
mi_setarcsegmentprop(1,'<None>',0,3);
mi_copyrotate2(0,0,StepTooth,2,3);
mi_clearselected;

%stator yoke design
mi_addnode(RextS,0); 
mi_addnode(RextS*cos(2*pi/(2*np)),RextS*sin(2*pi/(2*np)));
mi_addarc(RextS,0,RextS*cos(2*pi/(2*np)),RextS*sin(2*pi/(2*np)),360/(2*np),1);
mi_addsegment(Rbor,0,RextS,0);
mi_addsegment(Rbor*cos(2*pi/(2*np)),Rbor*sin(2*pi/(2*np)),RextS*cos(2*pi/(2*np)),RextS*sin(2*pi/(2*np)));

%Little square 1 for calculation of flux density
mi_addnode(0.999*Rme*cos(0.99*2*StepTooth*deg), 0.999*Rme*sin(0.99*2*StepTooth*deg));
mi_addnode(1.001*Rme*cos(0.99*2*StepTooth*deg), 1.001*Rme*sin(0.99*2*StepTooth*deg));
mi_addnode(0.999*Rme*cos(1.01*2*StepTooth*deg), 0.999*Rme*sin(1.01*2*StepTooth*deg));
mi_addnode(1.001*Rme*cos(1.01*2*StepTooth*deg), 1.001*Rme*sin(1.01*2*StepTooth*deg));
mi_addsegment(0.999*Rme*cos(0.99*2*StepTooth*deg), 0.999*Rme*sin(0.99*2*StepTooth*deg),1.001*Rme*cos(0.99*2*StepTooth*deg), 1.001*Rme*sin(0.99*2*StepTooth*deg));
mi_addsegment(1.001*Rme*cos(0.99*2*StepTooth*deg), 1.001*Rme*sin(0.99*2*StepTooth*deg),1.001*Rme*cos(1.01*2*StepTooth*deg), 1.001*Rme*sin(1.01*2*StepTooth*deg));
mi_addsegment(1.001*Rme*cos(1.01*2*StepTooth*deg), 1.001*Rme*sin(1.01*2*StepTooth*deg),0.999*Rme*cos(1.01*2*StepTooth*deg), 0.999*Rme*sin(1.01*2*StepTooth*deg));
mi_addsegment(0.999*Rme*cos(0.99*2*StepTooth*deg), 0.999*Rme*sin(0.99*2*StepTooth*deg),0.999*Rme*cos(1.01*2*StepTooth*deg), 0.999*Rme*sin(1.01*2*StepTooth*deg));

mi_addblocklabel(Rme*cos(2*StepTooth*deg),Rme*sin(2*StepTooth*deg));
mi_selectlabel(Rme*cos(2*StepTooth*deg),Rme*sin(2*StepTooth*deg));
mi_setblockprop('Pure Iron', 1, 0, '<None>', 0, 1, 1);
mi_clearselected;

%Litttle square 2 for caclulation of flux density
RmCs=(Rfe+RextS)/2;
mi_addnode(0.999*RmCs*cos(0.99*2.5*StepTooth*deg), 0.999*RmCs*sin(0.99*2.5*StepTooth*deg));
mi_addnode(1.001*RmCs*cos(0.99*2.5*StepTooth*deg), 1.001*RmCs*sin(0.99*2.5*StepTooth*deg));
mi_addnode(0.999*RmCs*cos(1.01*2.5*StepTooth*deg), 0.999*RmCs*sin(1.01*2.5*StepTooth*deg));
mi_addnode(1.001*RmCs*cos(1.01*2.5*StepTooth*deg), 1.001*RmCs*sin(1.01*2.5*StepTooth*deg));
mi_addsegment(0.999*RmCs*cos(0.99*2.5*StepTooth*deg), 0.999*RmCs*sin(0.99*2.5*StepTooth*deg),1.001*RmCs*cos(0.99*2.5*StepTooth*deg), 1.001*RmCs*sin(0.99*2.5*StepTooth*deg));
mi_addsegment(1.001*RmCs*cos(0.99*2.5*StepTooth*deg), 1.001*RmCs*sin(0.99*2.5*StepTooth*deg),1.001*RmCs*cos(1.01*2.5*StepTooth*deg), 1.001*RmCs*sin(1.01*2.5*StepTooth*deg));
mi_addsegment(1.001*RmCs*cos(1.01*2.5*StepTooth*deg), 1.001*RmCs*sin(1.01*2.5*StepTooth*deg),0.999*RmCs*cos(1.01*2.5*StepTooth*deg), 0.999*RmCs*sin(1.01*2.5*StepTooth*deg));
mi_addsegment(0.999*RmCs*cos(1.01*2.5*StepTooth*deg), 0.999*RmCs*sin(1.01*2.5*StepTooth*deg),0.999*RmCs*cos(0.99*2.5*StepTooth*deg), 0.999*RmCs*sin(0.99*2.5*StepTooth*deg));

mi_addblocklabel(RmCs*cos(2.5*StepTooth*deg),RmCs*sin(2.5*StepTooth*deg));
mi_selectlabel(RmCs*cos(2.5*StepTooth*deg),RmCs*sin(2.5*StepTooth*deg));
mi_setblockprop('Pure Iron', 1, 0, '<None>', 0, 1, 1);
mi_clearselected;


%stator material assignment
%yoke
RmCs=(Rfe+RextS)/2;
mi_addblocklabel(RmCs*cos(pi/(2*np)),RmCs*sin(pi/(2*np)));
mi_selectlabel(RmCs*cos(pi/(2*np)),RmCs*sin(pi/(2*np)));
mi_setblockprop('Pure Iron', 1, 0, '<None>', 0, 1, 1);
mi_clearselected;
% Current

Lup=repmat(['B+';'A-';'C-';'B-';'A+';'C+'],np,1);
Ldown=repmat(['C-';'B-';'A+';'C+';'B+';'A-'],np,1);
%Upper slots
for i=1:3
    mi_addblocklabel(Rme*cos(StepTooth*(1.5-0.5*betad+2*(i-1))*deg/2),Rme*sin(StepTooth*(1.5-0.5*betad+2*(i-1))*deg/2));
    mi_selectlabel(Rme*cos(StepTooth*(1.5-0.5*betad+2*(i-1))*deg/2),Rme*sin(StepTooth*(1.5-0.5*betad+2*(i-1))*deg/2));
    mi_setblockprop('Coil',1,0,Lup(i,:),0,3,Nt);
    mi_clearselected;
end
%Downer slots
for j=1:3
    mi_addblocklabel(Rme*cos(StepTooth*(1+betad+4*(j-1))*deg/4),Rme*sin(StepTooth*(1+betad+4*(j-1))*deg/4));
    mi_selectlabel(Rme*cos(StepTooth*(1+betad+4*(j-1))*deg/4),Rme*sin(StepTooth*(1+betad+4*(j-1))*deg/4));
    mi_setblockprop('Coil',1,0,Ldown(j,:),0,3,Nt);
    mi_clearselected;
end

%stator CL assignment 
mi_selectsegment(0.5*(Rag2+Rbor),0);
mi_selectsegment(0.5*(Rag2+Rbor)*cos(2*pi/(2*np)),0.5*(Rag2+Rbor)*sin(2*pi/(2*np)));
mi_setsegmentprop('Ap4',0,1,0,1);
mi_clearselected;
mi_selectsegment(0.5*(Rbor+RextS),0);
mi_selectsegment(0.5*(Rbor+RextS)*cos(2*pi/(2*np)),0.5*(Rbor+RextS)*sin(2*pi/(2*np)));
mi_setsegmentprop('Ap5',0,1,0,1);
mi_clearselected;
%interior arc
mi_selectarcsegment(RextS*cos(2*pi/(4*np)),RextS*sin(2*pi/(4*np)));
mi_setarcsegmentprop(1,'A=0',0,1);
mi_clearselected;
%geometry backup
mi_saveas('geometry.fem');
mi_close;

SurfaceE=Rme*StepTooth*deg*(1-betad)*Dslot*1e-6/2;
Mmagnet=2*StepTooth*deg*RmA*Thickmagnet*length*7500*10^(-9);
Mcu=0.6*6*SurfaceE*length*8960*10^(-3)*(1-0.42)/(1-betad);

for i=0:npas-1
%for i=0:5
    
    tta=i*dtpas+tta0;
    ang(i+1)=tta;
    opendocument('geometry.fem');
    %No changing the currents because I=0.
    mi_selectgroup(2);
    mi_moverotate(0,0,tta);
    mi_clearselected;
    mi_addarc(Rag1,0,Rag2*cos(tta*deg),Rag2*sin(tta*deg),tta,1);
    mi_addarc(Rag1*cos(pi/np),Rag1*sin(pi/np),Rag2*cos(tta*deg+pi/np),Rag2*sin(tta*deg+pi/np),tta,1);
    mi_selectarcsegment(Rag0*cos(tta*deg/2),Rag0*sin(tta*deg/2));
    mi_selectarcsegment(Rag0*cos(pi/np+tta*deg/2),Rag0*sin(pi/np+tta*deg/2));
    mi_setarcsegmentprop(1,'Ap3',0,1);
    mi_clearselected;
    mi_zoomnatural();
    mi_saveas('temp1.fem');
    mi_analyse;
    mi_loadsolution;
    
    
    %trigger torque storage
    mo_selectblock(RmA*cos(0.5*pi/np+tta*deg),RmA*sin(0.5*pi/np+tta*deg));
    mo_selectblock(RmCR*cos(0.5*pi/np+tta*deg),RmCR*sin(0.5*pi/np+tta*deg));
    det(i+1)=nd*mo_blockintegral(22);
    mo_clearblock;

    %flows Phase 1
    mo_selectblock(Rme*cos(0.5*StepTooth*(1+betad)*deg/2),Rme*sin(0.5*StepTooth*(1+betad)*deg/2));
    a1=nd*Nt*mo_blockintegral(1)/SurfaceE;
    mo_clearblock;
    mo_selectblock(Rme*cos(StepTooth*(1.5-0.5*betad)*deg/2+2*StepTooth*deg),Rme*sin(StepTooth*(1.5-0.5*betad)*deg/2+2*StepTooth*deg));
    a6=nd*Nt*mo_blockintegral(1)/SurfaceE;
    mo_clearblock;
    Flux1(i+1)=a1+a6;
    %flow Phase 2
    mo_selectblock(Rme*cos(StepTooth*(1.5-0.5*betad)*deg/2),Rme*sin(StepTooth*(1.5-0.5*betad)*deg/2));
    a2=nd*Nt*mo_blockintegral(1)/SurfaceE;
    mo_clearblock;
    mo_selectblock(Rme*cos(0.5*StepTooth*(1+betad)*deg/2+StepTooth*deg),Rme*sin(0.5*StepTooth*(1+betad)*deg/2+StepTooth*deg));
    a3=nd*Nt*mo_blockintegral(1)/SurfaceE;
    Flux2(i+1)=a2+a3;
    mo_clearblock;
    %flow Phase 3
    mo_selectblock(Rme*cos(StepTooth*(1.5-0.5*betad)*deg/2+StepTooth*deg),Rme*sin(StepTooth*(1.5-0.5*betad)*deg/2+StepTooth*deg));
    a4=nd*Nt*mo_blockintegral(1)/SurfaceE;
    mo_clearblock;
    mo_selectblock(Rme*cos(0.5*StepTooth*(1+betad)*deg/2+2*StepTooth*deg),Rme*sin(0.5*StepTooth*(1+betad)*deg/2+2*StepTooth*deg)); 
    a5=nd*Nt*mo_blockintegral(1)/SurfaceE;
    Flux3(i+1)=a4+a5;
    mo_clearblock;

    %induction in the stator yoke
    Sb=0.00008*StepTooth*deg*(Rme^2)*10^(-6);
    mo_selectblock(Rme*cos(2*StepTooth*deg),Rme*sin(2*StepTooth*deg));
    Indx(i+1)=mo_blockintegral(8)/Sb;
    Indy(i+1)=mo_blockintegral(9)/Sb;
    Ind(i+1)=sqrt(Indx(i+1)*Indx(i+1)+Indy(i+1)*Indy(i+1));
    mo_clearblock;
    
    Sb2=0.0001*StepTooth*deg*(RmCs^2)*10^(-6);
    mo_selectblock(RmCs*cos(2.5*StepTooth*deg),RmCs*sin(2.5*StepTooth*deg));
    Indx2(i+1)=mo_blockintegral(8)/Sb2;
    Indy2(i+1)=mo_blockintegral(9)/Sb2;
    Ind2(i+1)=sqrt(Indx2(i+1)*Indx2(i+1)+Indy2(i+1)*Indy2(i+1));
    mo_clearblock;
    mi_close;
    mo_close;
end
for j=0:npas-1
    ang(npas+j+1)=(npas+j)*dtpas+tta0;
    Flux1(npas+j+1)=-Flux1(j+1);
    Flux2(npas+j+1)=-Flux2(j+1);
    Flux3(npas+j+1)=-Flux3(j+1);
    Ind(npas+j+1)=Ind(j+1);
    Ind2(npas+j+1)=Ind2(j+1);
    det(npas+j+1)=det(j+1);
end

for j=1:2*npas-2
    fem1(j+1)=(Flux1(j+2)-Flux1(j))/(2*dtpas*deg);
    fem2(j+1)=(Flux2(j+2)-Flux2(j))/(2*dtpas*deg);
    fem3(j+1)=(Flux3(j+2)-Flux3(j))/(2*dtpas*deg);
end
fem1(1)= (Flux1(2)-Flux1(2*npas))/(2*dtpas*deg);   
fem1(2*npas)=(Flux1(1)-Flux1(2*npas-1))/(2*dtpas*deg);
fem2(1)= (Flux2(2)-Flux2(2*npas))/(2*dtpas*deg)   ;
fem2(2*npas)=(Flux2(1)-Flux2(2*npas-1))/(2*dtpas*deg);
fem3(1)= (Flux3(2)-Flux3(2*npas))/(2*dtpas*deg);  
fem3(2*npas)=(Flux3(1)-Flux3(2*npas-1))/(2*dtpas*deg);

%closefemm;


% calculation at position 0 NeutRborized magnet and unit current in phase
% 1
% opendocument('geometrieMAP.fem');
% mi_modifycircprop('A',1,1);
% mi_modifycircprop('B',1,0);
% mi_modifycircprop('C',1,0);
% 
% mi_selectlabel(RmA*cos(pi/(2*np)),RmA*sin(pi/(2*np)));
% mi_setblockprop('Air', 1, 0, '<None>', 0, 2, 1);
% mi_clearselected;
% 
% mi_addsegment(Rag1,0,Rag2,0);
% mi_addsegment(Rag1*cos(pi/np),Rag1*sin(pi/np),Rag2*cos(pi/np),Rag2*sin(pi/np));
% mi_selectsegment(Rag0,0);
% mi_selectsegment(Rag0*cos(pi/np),Rag0*sin(pi/np));
% mi_setarcsegmentprop(1,'Ap3',0,1);
% mi_setsegmentprop('Ap3', 0, 1, 0, 1)
% mi_clearselected;
% mi_zoomnatural();
% mi_saveas('temp2.fem');
% mi_analyse;
% mi_loadsolution;
% mo_selectblock(Rme*cos(0.5*StepTooth*deg),Rme*sin(0.5*StepTooth*deg));
% Lxx=nd*Nt*mo_blockintegral(1)/SurfaceE;
% mo_clearblock;
% mo_selectblock(Rme*cos(2.5*StepTooth*deg),Rme*sin(2.5*StepTooth*deg));
% Mxx=nd*Nt*mo_blockintegral(1)/SurfaceE;
% mo_showdensityplot(1,0,0,2,'mag')

%Plot of the field map at position 0 


opendocument('geometry.fem');
I=Al*pi*2*Rbor*0.001/(6*np*n);
mi_modifycircprop('A+',1,sqrt(2)*I*cos(-5*np*dtpas-0*pi/3));
mi_modifycircprop('B+',1,sqrt(2)*I*cos(-5*np*dtpas-2*pi/3));
mi_modifycircprop('C+',1,sqrt(2)*I*cos(-5*np*dtpas-4*pi/3));
mi_modifycircprop('A-',1,-sqrt(2)*I*cos(-5*np*dtpas-0*pi/3));
mi_modifycircprop('B-',1,-sqrt(2)*I*cos(-5*np*dtpas-2*pi/3));
mi_modifycircprop('C-',1,-sqrt(2)*I*cos(-5*np*dtpas-4*pi/3));

mi_selectlabel(RmA*cos(pi/(2*np)),RmA*sin(pi/(2*np)));
mi_setblockprop('NdFeB 32 MGOe', 1, 0, '<None>', 180/(2*np), 2, 1);
mi_clearselected;

mi_addsegment(Rag1,0,Rag2,0);
mi_addsegment(Rag1*cos(pi/np),Rag1*sin(pi/np),Rag2*cos(pi/np),Rag2*sin(pi/np));
mi_selectsegment(Rag0,0);
mi_selectsegment(Rag0*cos(pi/np),Rag0*sin(pi/np));
mi_setarcsegmentprop(1,'Ap3',0,1);
mi_setsegmentprop('Ap3', 0, 1, 0, 1)
mi_clearselected;

%Little square for calculation of flux density in the rotor
Sb3=0.0001*StepTooth*deg*(RmCR^2)*10^(-6);
mi_addnode(0.999*RmCR*cos(0.99*2.5*StepTooth*deg), 0.999*RmCR*sin(0.99*2.5*StepTooth*deg));
mi_addnode(1.001*RmCR*cos(0.99*2.5*StepTooth*deg), 1.001*RmCR*sin(0.99*2.5*StepTooth*deg));
mi_addnode(0.999*RmCR*cos(1.01*2.5*StepTooth*deg), 0.999*RmCR*sin(1.01*2.5*StepTooth*deg));
mi_addnode(1.001*RmCR*cos(1.01*2.5*StepTooth*deg), 1.001*RmCR*sin(1.01*2.5*StepTooth*deg));
mi_addsegment(0.999*RmCR*cos(0.99*2.5*StepTooth*deg), 0.999*RmCR*sin(0.99*2.5*StepTooth*deg),1.001*RmCR*cos(0.99*2.5*StepTooth*deg), 1.001*RmCR*sin(0.99*2.5*StepTooth*deg));
mi_addsegment(1.001*RmCR*cos(0.99*2.5*StepTooth*deg), 1.001*RmCR*sin(0.99*2.5*StepTooth*deg),1.001*RmCR*cos(1.01*2.5*StepTooth*deg), 1.001*RmCR*sin(1.01*2.5*StepTooth*deg));
mi_addsegment(1.001*RmCR*cos(1.01*2.5*StepTooth*deg), 1.001*RmCR*sin(1.01*2.5*StepTooth*deg),0.999*RmCR*cos(1.01*2.5*StepTooth*deg), 0.999*RmCR*sin(1.01*2.5*StepTooth*deg));
mi_addsegment(0.999*RmCR*cos(1.01*2.5*StepTooth*deg), 0.999*RmCR*sin(1.01*2.5*StepTooth*deg),0.999*RmCR*cos(0.99*2.5*StepTooth*deg), 0.999*RmCR*sin(0.99*2.5*StepTooth*deg));

mi_addblocklabel(RmCR*cos(2.5*StepTooth*deg),RmCR*sin(2.5*StepTooth*deg));
mi_selectlabel(RmCR*cos(2.5*StepTooth*deg),RmCR*sin(2.5*StepTooth*deg));
mi_setblockprop('Pure Iron', 1, 0, '<None>', 0, 1, 1);
mi_clearselected;


mi_zoomnatural();
mi_saveas('temp3.fem');
mi_analyse;
mi_loadsolution;

mo_selectblock(Rme*cos(2*StepTooth*deg),Rme*sin(2*StepTooth*deg));
Inx=mo_blockintegral(8)/Sb;
Iny=mo_blockintegral(9)/Sb;
Indtooth=sqrt(Inx^2+Iny^2);
mo_clearblock;

mo_selectblock(RmCs*cos(2.5*StepTooth*deg),RmCs*sin(2.5*StepTooth*deg));
Inx2=mo_blockintegral(8)/Sb2;
Iny2=mo_blockintegral(9)/Sb2;
IndyokeS=sqrt(Inx2^2+Iny2^2);
mo_clearblock;

mo_selectblock(RmCR*cos(2.5*StepTooth*deg),RmCR*sin(2.5*StepTooth*deg));
Inx3=mo_blockintegral(8)/Sb3;
Iny3=mo_blockintegral(9)/Sb3;
IndyokeR=sqrt(Inx3^2+Iny3^2);
mo_clearblock;

%mo_showdensityplot(1,0,0,2,'mag');



%Flow chart and EMF
% figure(1);
% plot(ang,Flux1,'b+-.');
% hold on
% plot(ang,Flux2,'k');
% hold on
% plot(ang,Flux3,'c+');
% legend('Flow in the phases');
% figure(2)
% plot(ang,fem1);
%hold on
%plot(ang,fem2,'k');
%hold on
%plot(ang,fem3,'c+');
%legend('FEM at Omega_m =1rad/s');
%figure(3)
%plot(ang,det);
%legend('Detent torque');
% figure(4)
% plot(ang,Ind);
% legend('|Flux density Tooth|');
% figure(5)
% plot(ang,Ind2);
% legend('|Flux density Yoke|');
% disp(['Natural inductance of a phase Lxx = ', num2str(Lxx)]);
% disp(['Mutual inductance Mxy = ', num2str(Mxx)]);
% 
% disp(['Flux density tooth (t=0) = ', num2str(Indtooth)]);
% disp(['Flux density yoke (t=0) = ', num2str(Indyoke)]);

fem1max= max(abs(fftshift(fft(fem1))))/npas;
E=fem1max*(2*pi*Nn/60)*2*np/sqrt(2);
P=3*E*I;
%Test efficiency
Spp=1; %Nombres d'encoches par p?le et par phase
RhoCu=1.725*(10^(-8)); %R?sistivit? du Cuivre
miron=7860; %Masse volumique du fer
Po=5.637; %Pertes volumique du fer ? 50Hz, 1,5T
Mtooth=miron*pi*0.001*length*(0.001*Dslot*betad*0.001*(2*Rbor+2*Pbec+Dslot)+0.001*(Pbec*0.001*(1-betad)*(1-ferm)*0.001*(2*Rbor+Pbec))); %Masse dent
MyokeS=miron*0.001*ThickyokeR*pi*0.001*(2*Rbor+2*Dslot+ThickyokeR+2*Pbec)*0.001*length;
MyokeR=miron*0.001*ThickyokeS*pi*0.001*(2*Rbor-2*airgap-2*Thickmagnet)*0.001*length;
%Masse cul
Pftooth=3*Po*(((Nn/(50*60))^1.5)*((Indtooth/1.5)^2.2))*Mtooth; %Iron Losses tooth
PfyokeS=1.5*Po*(((Nn/(50*60))^1.5)*((IndyokeS/1.5)^2.2))*MyokeS; %Iron Losses Yoke Stator
PfyokeR=1.5*Po*(((Nn/(50*60))^1.5)*((IndyokeR/1.5)^2.2))*MyokeR; %Iron Losses Yoke Rotor
Piron=Pftooth+PfyokeS+PfyokeR; %Totale Iron losses
%Joules losses
Sc=kr*Dslot*0.001*(1-0.42)*pi*2*Rbor*0.001/(6*n*np*Spp);
Ra=RhoCu*6*n*np*Spp*length*0.001/Sc;
Pj=3*Ra*(I)^2; 
Losses=Pj+Piron;
eff=P/(P+Losses);
%M=1.8*10^3/(Mtooth+Myoke+Mcu+Mmagnet);
effnew=(exp(eff)/0.06306)-40.995784;
M=541/(Mtooth+MyokeS+MyokeR+Mcu+Mmagnet);
Mnew=(sqrt(M)/0.426349)- 1.3454961;
%opti=1/(0.6*effnew+0.4*Mnew);
inveff=1/eff;
invM=1/M;
opti=[inveff invM];
end 