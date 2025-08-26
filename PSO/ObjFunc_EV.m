function [f, Tavg, ripple] = ObjFunc(X)
% ΜΟΝΟ ΗΛΕΚΤΡΟΜΑΓΝΗΤΙΚΟΙ ΣΤΟΧΟΙ: μεγιστοποίηση μέσης ροπής, ελαχιστοποίηση ripple,
% ελαχιστοποίηση απωλειών σιδήρου & μάζας μαγνήτη (μέσω κόστους).
%
% ΕΠΙΣΤΡΕΦΕΙ:
%   f      : scalar objective για PSO (minimize)
%   Tavg   : μέση ροπή (Nm)
%   ripple : (max(T)-min(T))/Tavg   (αδιάστατο)

% --------- ΣΤΑΘΕΡΕΣ ΜΟΝΤΕΛΟΥ ---------
deg    = pi/180;
np     = 8;          % ζεύγη πόλων/2? -> εδώ έχεις np=8: άρα p=8/2? Στον δικό σου κώδικα nd=2*np
p      = np;         % ηλεκτρικοί πόλοι/2 = ζεύγη πόλων (ταιριάζει με θ_e = p*θ_m)
npas   = 18;         % βήματα σάρωσης (κράτα μικρό για ταχύτητα, 12–24)
I_rms  = 200;        % Ρεύμα φάσης RMS για αξιολόγηση (μόνο ηλεκτρομαγνητικό)
I_pk   = sqrt(2)*I_rms;

% ---------- FEMM SETUP ----------
if ~exist('openfemm','file')
    addpath('C:\femm42\mfiles'); % προσαρμόσ’ το
end
openfemm(1);
newdocument(0);
length = 1085;               % mm (έχεις ήδη)
mi_probdef(0,'millimeters','planar',1e-8,length,30);

% Κυκλώματα (μηδενικά αρχικά, θα τα αλλάζουμε στη σάρωση)
mi_addcircprop('A+', 0, 1);
mi_addcircprop('B+', 0, 1);
mi_addcircprop('C+', 0, 1);
mi_addcircprop('A-', 0, 1);
mi_addcircprop('B-', 0, 1);
mi_addcircprop('C-', 0, 1);

% Υλικά (safe)
try, mi_getmaterial('Air'); end
try, mi_getmaterial('Pure Iron'); end
try, mi_getmaterial('NdFeB 32 MGOe'); end
mi_addmaterial('Coil',1,1,0,0,58,0,0,1,3,0,0,1,7.77);

% ---------- ΓΕΩΜΕΤΡΙΑ ----------
% Χρησιμοποίησε την υπάρχουσα γεωμετρία σου (ό,τι έχεις πιο πάνω),
% ΑΛΛΑ βγάλε τα “I=0” και τις περίπλοκες καμπύλες για ταχύτητα.
% Ιδανικά βάλε τη δική σου build_συνάρτηση. Για συντομία, καλούμε helper:
build_geom_minimal(X, np);   % Δες helper στο τέλος (πατάει στα δικά σου R*, betad, beta, κ.λπ.)

% Fallback Air label για σιγουριά
mi_addblocklabel(0,0); mi_selectlabel(0,0);
mi_setblockprop('Air',1,0,'',0,0,0); mi_clearselected;

% SAVE πριν από mesh/solve
mi_saveas(fullfile(tempdir,'model_tmp.fem'));

% ---------- ΣΑΡΩΣΗ ΡΟΠΤΙΑΣ ----------
T = zeros(1, npas);
for k = 1:npas
    theta_m = (k-1) * (pi/npas);  % μηχανική γωνία (rad) σε μισό περίοδο για παράδειγμα
    theta_e = p * theta_m;        % ηλεκτρική γωνία
    
    % 3-φασικά ρεύματα (instantaneous)
    Ia = I_pk*sin(theta_e + 0); 
    Ib = I_pk*sin(theta_e - 2*pi/3);
    Ic = I_pk*sin(theta_e - 4*pi/3);
    
    % Εφάρμοσε ρεύματα στα circuits (A+/A- κτλ)
    mi_modifycircprop('A+',1, Ia);  mi_modifycircprop('A-',1,-Ia);
    mi_modifycircprop('B+',1, Ib);  mi_modifycircprop('B-',1,-Ib);
    mi_modifycircprop('C+',1, Ic);  mi_modifycircprop('C-',1,-Ic);

    % Περιστροφή ρότορα κατά theta_m (χρησιμοποίησε group=2 που έχεις για rotor)
    if k==1
        % τίποτα
    else
        mi_selectgroup(2);
        mi_moverotate(0,0, (theta_m - ((k-2)*(pi/npas)))/deg ); % σε μοίρες
        mi_clearselected;
    end

    % Λύση
    mi_saveas(fullfile(tempdir,'model_tmp.fem'));
    mi_createmesh; mi_analyze(1); pause(0.02); mi_loadsolution;

    % Ροπή (Nm): planar -> Nm/m * depth(m)
    depth_m = length*1e-3;
    try
        mo_groupselectblock(2);                 % rotor group
        T_2D = mo_blockintegral(22);            % Nm/m
        mo_clearblock;
        T(k) = T_2D*depth_m;
    catch
        T(k) = 0;
    end

    mo_close;  % κλείσε post για να μπορείς να μετακινήσεις/αναλύσεις ξανά
end

% Μετρικές ροπής
Tavg   = mean(T);
if Tavg<=0, Tavg = max(Tavg, 1e-6); end
ripple = (max(T)-min(T))/Tavg;

% ---------- Proxy απωλειών σιδήρου & μάζα μαγνήτη ----------
% (ελαφρύς υπολογισμός: δειγματοληψία |B| σε 2-3 σημεία ανά δόντι/ζυγό)
% εδώ δείχνω μόνο απλό proxy = (B_tooth^2 + B_yoke^2)
Pfe_proxy = 0; Mmagnet = 0;
try
    opendocument(fullfile(tempdir,'model_tmp.fem'));
    mi_analyze(1); mi_loadsolution;
    % Δείγμα σημείων (προσαρμόσ’ το στα δικά σου)
    +% π.χ. ένα σημείο tooth & ένα yoke:
    Rtooth = 0.7; Ryoke = 0.85; Rrot = 0.6;  % ΣΧΕΤΙΚΑ προς RextS (π.χ. normalize, εδώ ενδεικτικά)
    [Bx1,By1] = mo_getb(Rtooth*100, 0);      % ΠΡΟΣΑΡΜΟΣΕ συντεταγμένες στο δικό σου μοντέλο!
    [Bx2,By2] = mo_getb(Ryoke*100,  0);
    B1 = hypot(Bx1,By1); B2 = hypot(Bx2,By2);
    Pfe_proxy = B1^2 + B2^2;                 % απλό proxy, όσο μικρότερο τόσο καλύτερο
    % Μάζα μαγνήτη (αν έχεις group=magnet=1, μπορείς από area)
    mo_clearblock; mo_groupselectblock(1);
    A_mag_mm2 = mo_blockintegral(5); mo_clearblock;
    Mmagnet = 7500 * (A_mag_mm2*1e-6) * (length*1e-3); % kg ~ ρ*V
    mo_close;
catch
    Pfe_proxy = 10; Mmagnet = 0.1;  % ασφαλή defaults
end

closefemm;

% ---------- SCALAR OBJECTIVE (μόνο Η/Μ) ----------
% Μικρότερο f = καλύτερο:
w1 = 1.0;   % βάρος για 1/Tavg
w2 = 0.5;   % βάρος για ripple
w3 = 0.05;  % βάρος για Pfe_proxy
w4 = 0.01;  % βάρος για Mmagnet
f = w1*(1/max(Tavg,1e-6)) + w2*ripple + w3*Pfe_proxy + w4*Mmagnet;

end

% ===== helper γεωμετρίας (πολύ απλοποιημένη, κράτα τη δική σου καλύτερα) =====
function build_geom_minimal(X, np)
% X = [beta_d, ThR, ThS, beta, Tmag] κατά το δικό σου σχήμα ορίων.
% Εδώ βάλε την "δική σου" γεωμετρία από το μεγάλο script σου.
% ΣΗΜΑΝΤΙΚΟ: Βάλε σωστά group=2 για ROTOR, group=1 για MAGNETS, group=3 για COILS, group=1/2 για steel κ.λπ.
% Για συντομία παραλείπω λεπτομέρειες — χρησιμοποίησε το ήδη δοκιμασμένο τμήμα σου
% (απλώς φρόντισε να μην ορίζεις I=0 εδώ, τα ρεύματα τα αλλάζουμε στη σάρωση).
%
% Τip: κράτα όλα τα mi_addblocklabel(...); mi_setblockprop(..., group, turns)
% ώστε να παίζει mo_groupselectblock(2) για torque και group=1 για μαγνήτες.
end
