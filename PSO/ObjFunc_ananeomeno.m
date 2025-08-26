function [opti, eff, M] = ObjFunc(X)
% EM-only αντικειμενική: μεγιστοποίηση μέσης ροπής, ελαχιστοποίηση ripple και μάζας.
% ΕΠΙΣΤΡΕΦΕΙ:
%   opti : 1/f (όσο μεγαλύτερο τόσο καλύτερα)
%   eff  : NaN (placeholder μέχρι να μπει σωστό loss model)
%   M    : μάζα (kg) από FEMM groups (magnets=1, steel=2, coils=3)

%% --------------------- ΣΤΑΘΕΡΕΣ & ΠΑΡΑΜΕΤΡΟΙ ---------------------
deg        = pi/180;
length_mm  = 1085;     % stack length [mm]
airgap_mm  = 1.6;
Rbor_mm    = 400;

% Μεταβλητές βελτιστοποίησης
beta       = X(4);     % pole-arc ratio (0..1)
betad      = X(1);     % tooth/slot ratio (0..1)
ThR_mm     = X(2);     % rotor yoke thickness
ThS_mm     = X(3);     % stator yoke thickness
Tmag_mm    = X(5);     % magnet thickness

% Λοιπά γεωμετρικά
Dslot_mm   = 29.9;
np         = 8;        % ζεύγη πόλων
Pbec_mm    = 3;        % δοντάκι (tooth tip)
ferm       = 0.5;      % slot closing coef
nd         = 2*np;     % όπως στο αρχικό

% Αξιολόγηση με ρεύματα (instantaneous snapshot)
n_cond     = 80;                  % αγωγοί/φάση;
Nt         = n_cond/(2*np);       % αγωγοί/αυλάκι
npas       = 24;                  % βήματα σάρωσης
dtpas_deg  = 180/(npas*np);       % όπως στο αρχικό
tta0_deg   = dtpas_deg/2;
Al         = 60000;               % όπως στο αρχικό
% Εκτίμηση ρεύματος φάσης (ίδιος τύπος με το αρχικό σου)
I_eval_A   = Al*pi*2*Rbor_mm*0.001/(6*np*n_cond);  % A (πρόχειρη εκτίμηση)
Ipk        = max(1, sqrt(2)*I_eval_A);             % προφύλαξη από 0

% FEMM υλικά (kg/m^3)
rho.NdFeB  = 7500;
rho.Steel  = 7800;
rho.Cu     = 8960;

% Βάρη αντικειμενικής (EM μόνο)
w1 = 1.0;           % 1/Tavg
w2 = 0.5;           % ripple
w3 = 0.02;          % M/target_mass
target_mass = 60;   % kg (κλίμακα)

%% --------------------- FEMM SETUP ---------------------
% Φόρτωσε FEMM API
if ~exist('openfemm','file')
    addpath('C:\femm42\mfiles');  % προσαρμόσ’ το αν θες
end

openfemm(1);            % με UI (1) για debug
newdocument(0);         % magnetics

% Κύριος ορισμός προβλήματος (planar)
mi_probdef(0,'millimeters','planar',1e-8,length_mm,30);

% Κυκλώματα (θα τα αλλάζουμε ανά βήμα)
mi_addcircprop('A+',0,1); mi_addcircprop('A-',0,1);
mi_addcircprop('B+',0,1); mi_addcircprop('B-',0,1);
mi_addcircprop('C+',0,1); mi_addcircprop('C-',0,1);

% Υλικά
try, mi_getmaterial('Air');                  end
try, mi_getmaterial('Pure Iron');            end
try, mi_getmaterial('NdFeB 32 MGOe');        end
mi_addmaterial('Coil',1,1,0,0,58,0,0,1,3,0,0,1,7.77);

%% --------------------- ΓΕΩΜΕΤΡΙΑ (όπως στο αρχικό, με ΔΙΟΡΘΩΜΕΝΑ groups) ---------------------
% Βασικές ακτίνες
RintCR     = Rbor_mm - airgap_mm - Tmag_mm - ThR_mm;
RextCR     = Rbor_mm - airgap_mm - Tmag_mm;
RextA      = Rbor_mm - airgap_mm;
Rag1       = RextA + 0.9*airgap_mm/2;
Rag2       = RextA + 1.1*airgap_mm/2;
Rag0       = RextA + airgap_mm/2;
Rfe        = Rbor_mm + Pbec_mm + Dslot_mm;
Rde        = Rbor_mm + Pbec_mm;
RextS      = Rbor_mm + Dslot_mm + ThS_mm + Pbec_mm;
Rme        = (Rbor_mm + Pbec_mm + Rfe)/2;
RmCR       = (RintCR + RextCR)/2;
StepTooth  = 180/(3*np);

% ---------- (Ακολουθεί το δικό σου σχέδιο, αμετάβλητο στις γραμμές, αλλά με groups) ----------
% --- rotor yoke/magnets/airgap κ.λπ. (ίδια σχεδίαση με το αρχικό σου) ---
mi_addnode(RintCR, 0);
mi_addnode(RintCR*cos(2*pi/(2*np)),RintCR*sin(2*pi/(2*np)));
mi_addarc(RintCR, 0, RintCR*cos(2*pi/(2*np)),RintCR*sin(2*pi/(2*np)),360/(2*np),1);
mi_addnode(RextCR, 0);
mi_addnode(RextCR*cos(((1-beta)/2)*2*pi/(2*np)), RextCR*sin(((1-beta)/2)*2*pi/(2*np)));
mi_addnode(RextCR*cos(2*pi/(2*np)-((1-beta)/2)*2*pi/(2*np)), RextCR*sin(2*pi/(2*np)-((1-beta)/2)*2*pi/(2*np)));
mi_addnode(RextCR*cos(2*pi/(2*np)),RextCR*sin(2*pi/(2*np)));
mi_addarc(RextCR, 0, RextCR*cos(((1-beta)/2)*2*pi/(2*np)), RextCR*sin(((1-beta)/2)*2*pi/(2*np)),((1-beta)/2)*2*180/(2*np) ,1);
xint1=RextCR*cos(((1-beta)/2)*2*pi/(2*np));  yint1=RextCR*sin(((1-beta)/2)*2*pi/(2*np));
xint2=RextCR*cos(2*pi/(2*np)-((1-beta)/2)*2*pi/(2*np)); yint2=RextCR*sin(2*pi/(2*np)-((1-beta)/2)*2*pi/(2*np));
mi_addarc(xint1,yint1,xint2,yint2,beta*360/(2*np),1);
xint3=RintCR*cos(2*pi/(2*np)); yint3=RintCR*sin(2*pi/(2*np));
xo3=RextCR*cos(2*pi/(2*np));   yo3=RextCR*sin(2*pi/(2*np));
mi_addarc(xint2,yint2,xo3,yo3,((1-beta)/2)*2*180/(2*np),1);
mi_addsegment(RintCR, 0, RextCR, 0);
mi_addsegment(xint3,yint3,xo3,yo3);

% Μαγνήτες (τόξο RextA)
xr31=RextA*cos(((1-beta)/2)*2*pi/(2*np)); yr31=RextA*sin(((1-beta)/2)*2*pi/(2*np));
xr32=RextA*cos(2*pi/(2*np)-((1-beta)/2)*2*pi/(2*np)); yr32=RextA*sin(2*pi/(2*np)-((1-beta)/2)*2*pi/(2*np));
mi_addnode(xr31,yr31); mi_addnode(xr32,yr32);
mi_addarc(xr31,yr31,xr32,yr32,beta*360/(2*np),1);
mi_addsegment(xint1,yint1,xr31,yr31);
mi_addsegment(xint2,yint2,xr32,yr32);

% Αέρας ρότορα
mi_addnode(Rag1,0); mi_addsegment(RextCR,0,Rag1,0);
mi_addnode(Rag1*cos(2*pi/(2*np)),Rag1*sin(2*pi/(2*np)));
mi_addsegment(xo3,yo3,Rag1*cos(2*pi/(2*np)),Rag1*sin(2*pi/(2*np)));

% GROUPS/BOUNDARIES ρότορα
% ... (επιλογές nodes/arcs/segments όπως στο αρχικό σου) ...
% Ρυθμίσεις για groups:
% yoke ρότορα -> steel group=2
RmCr=(RintCR+RextCR)/2;
mi_addblocklabel(RmCr*cos(pi/(2*np)),RmCr*sin(pi/(2*np)));
mi_selectlabel(RmCr*cos(pi/(2*np)),RmCr*sin(pi/(2*np)));
mi_setblockprop('Pure Iron', 1, 0, '<None>', 0, 2, 1);  % *** steel -> group 2
mi_clearselected;

% μαγνήτες -> group=1  (ΔΙΟΡΘΩΣΗ από το αρχικό)
RmA=(RextA+RextCR)/2;
mi_addblocklabel(RmA*cos(pi/(2*np)),RmA*sin(pi/(2*np)));
mi_selectlabel(RmA*cos(pi/(2*np)),RmA*sin(pi/(2*np)));
mi_setblockprop('NdFeB 32 MGOe', 1, 0, '<None>', 180/(2*np), 1, 1); % *** magnets -> group 1
mi_clearselected;

% air gap -> group 2 (για να ανήκει στο rotor group ή άστο <None>, δεν μετρά μάζα)
Rma=(RextA+Rag1)/2;
mi_addblocklabel(Rma*cos(pi/(2*np)),Rma*sin(pi/(2*np)));
mi_selectlabel(Rma*cos(pi/(2*np)),Rma*sin(pi/(2*np)));
mi_setblockprop('Air', 1, 0, '<None>', 0, 2, 1);
mi_clearselected;

% ΣΤΑΤΟΡΑΣ (όπως στο αρχικό, αλλά yoke -> group 2, coils -> group 3)
% ... (όλες οι mi_addnode/segment/arc όπως έχεις) ...

% stator yoke -> steel group=2  (ΔΙΟΡΘΩΣΗ από το αρχικό)
RmCs=(Rfe+RextS)/2;
mi_addblocklabel(RmCs*cos(pi/(2*np)),RmCs*sin(pi/(2*np)));
mi_selectlabel(RmCs*cos(pi/(2*np)),RmCs*sin(pi/(2*np)));
mi_setblockprop('Pure Iron', 1, 0, '<None>', 0, 2, 1);  % *** steel -> group 2
mi_clearselected;

% Πηνία (κρατάμε group=3 όπως έχεις στις Upper/Lower slots)
% mi_setblockprop('Coil',1,0,Lup(i,:),0,3,Nt) κ.ο.κ. (παραμένει)

% Περιοδικά όρια/CL όπως στο αρχικό σου
% ...
mi_saveas('geometry.fem');
mi_close;

%% --------------------- ΣΑΡΩΣΗ ΓΩΝΙΑΣ & ΡΟΠΗ ---------------------
T = zeros(1,npas);
depth_m        = length_mm*1e-3;
alpha_model_deg= 360/(2*np);
scale_sector   = 360/alpha_model_deg;   % = 2*np

for i = 0:npas-1
    tta_deg = i*dtpas_deg + tta0_deg;   % μηχανική γωνία (deg)
    opendocument('geometry.fem');

    % Περιστροφή ρότορα (group 2) κατά tta_deg
    mi_selectgroup(2);
    mi_moverotate(0,0, tta_deg);
    mi_clearselected;

    % "Εκ νέου" περιοδικά τόξα στο διάκενο (όπως στο αρχικό σου)
    mi_addarc(Rag1,0,               Rag2*cos(tta_deg*deg),             Rag2*sin(tta_deg*deg),             tta_deg,1);
    mi_addarc(Rag1*cos(pi/np),Rag1*sin(pi/np), Rag2*cos(tta_deg*deg+pi/np),Rag2*sin(tta_deg*deg+pi/np), tta_deg,1);
    mi_selectarcsegment(Rag0*cos(tta_deg*deg/2),            Rag0*sin(tta_deg*deg/2));
    mi_selectarcsegment(Rag0*cos(pi/np+tta_deg*deg/2),      Rag0*sin(pi/np+tta_deg*deg/2));
    mi_setarcsegmentprop(1,'Ap3',0,1);
    mi_clearselected;

    % 3-φασικά ρεύματα (σύγχρονα με θ_e = p*θ_m)
    theta_e = np * tta_deg * deg;  % rad
    Ia = Ipk*sin(theta_e + 0);
    Ib = Ipk*sin(theta_e - 2*pi/3);
    Ic = Ipk*sin(theta_e - 4*pi/3);

    mi_modifycircprop('A+',1, Ia);  mi_modifycircprop('A-',1,-Ia);
    mi_modifycircprop('B+',1, Ib);  mi_modifycircprop('B-',1,-Ib);
    mi_modifycircprop('C+',1, Ic);  mi_modifycircprop('C-',1,-Ic);

    mi_saveas('temp_step.fem');
    mi_analyse; mi_loadsolution;

    % Ροπή ρότορα: επίλεξε μπλοκ σε ρότορα (yoke + magnets) μέσω σημειακής επιλογής
    % (ίδιο κόλπο με το αρχικό σου)
    mo_selectblock(RmA*cos(0.5*pi/np+tta_deg*deg), RmA*sin(0.5*pi/np+tta_deg*deg));   % magnets
    mo_selectblock(RmCR*cos(0.5*pi/np+tta_deg*deg),RmCR*sin(0.5*pi/np+tta_deg*deg));  % rotor yoke
    T_sector = mo_blockintegral(22);  % N·m/m (planar)
    mo_clearblock;

    % Ολική ροπή = ροπή_τομέα * μήκος * (360/γωνία_μοντέλου)
    T(i+1) = T_sector * depth_m * scale_sector;

    mi_close; mo_close;
end

Tavg   = mean(T);
ripple = (max(T)-min(T))/max(Tavg,1e-6);

%% --------------------- ΜΑΖΑ από FEMM (groups) ---------------------
% Άνοιξε μια φορά, λύσε και μέτρα εμβαδά ανά group
opendocument('geometry.fem');
mi_analyse; mi_loadsolution;

% Magnets (group 1)
mo_clearblock; mo_groupselectblock(1);
A_mag   = mo_blockintegral(5);                  % m^2
Mmagnet = rho.NdFeB * A_mag * depth_m * scale_sector;

% Steel (group 2)  (stator+rotor)
mo_clearblock; mo_groupselectblock(2);
A_steel = mo_blockintegral(5);
Msteel  = rho.Steel * A_steel * depth_m * scale_sector;

% Coils (group 3)
mo_clearblock; mo_groupselectblock(3);
A_cu    = mo_blockintegral(5);
Mcu     = rho.Cu * A_cu * depth_m * scale_sector;

M = Mmagnet + Msteel + Mcu;

mo_close;
closefemm;

%% --------------------- ΑΝΤΙΚΕΙΜΕΝΙΚΗ & ΕΠΙΣΤΡΟΦΕΣ ---------------------
f    = w1*(1/max(Tavg,1e-6)) + w2*ripple + w3*(M/target_mass);
opti = 1/f;
eff  = NaN;   % μέχρι να μπει σωστό loss/efficiency model

end
