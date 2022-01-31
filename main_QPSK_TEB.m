%----------------------------------------------------------------------------------------
%         TP TS113:   Télécoms 2020-2021
%         Simulation d'une QPSK en bande de base en présence de bruit
%         Par:        Kadiri Asma - Abrouk Nada  (T1 - G1)
%----------------------------------------------------------------------------------------



clear all;
close all;
clc;

%% -------------   Initialisation des parametres   ----------------------------
%------------------------------------------------------------------------------

fe=4000000;                                               % Frequence d'echantillonage
Te=1/fe;                                                % Periode d'echantillonage
M=4;                                                    % Nombre de symboles dans la modulation
N=log2(M);                                              % Nombre de bits par symbole
Ds=1000000;                                                % Debit symboles (symbole/s)
Ts=1/Ds;                                                % Periode entre symboles
fse=Ts * fe;                                            % Rythme d'echantillonage
nfft=512;                                               % Nombre d'echantillons pour la fft
Ns=5000;                                                % Nombre de symboles par paquet
axe_temps=linspace(0,10*Ts -Te,100);                    % Axe de temps
% g_t=fonction_g_t(axe_temps,Ns,Ts);                      % Filtre de mise en forme (porte)
% g_t_a=fliplr(g_t);                                      % Filtre adapté à g_t
%Eg=sum((abs(g_t).^2));                                  % Energie du filtre de mise en forme
roll_off=0.35;                                           % Roll_off du filtre
span=8;                                                 % Nombre de symboles (car retard de groupe = 4Ts))
sps=fse;                                                 % Nombre d'echantillons par symbole
g_t_2=rcosdesign(roll_off,span,sps,'sqrt');             % Filtre de mise en forme (racine de cosinus surélevé)
g_t_a_2=fliplr(g_t_2);                                  % Filtre adapté à g_t_2
Eg_2=sum((abs(g_t_2).^2));                              % Energie du filtre de mise en forme 2
sigA2=1;                                                % variance théorique des symboles (Car variables aleatoires gaussiennes de Variance 1 par randi())
eb_n0_db= 0:0.5:10;                                     % Liste des Eb/N0 en db
eb_n0 = 10.^(eb_n0_db/10);                              % Liste des Eb/N0
%sigma2 = sigA2 * Eg ./(N * eb_n0);                      % Variance de bruit complexe en bande de base
sigma2_2 = sigA2 * Eg_2 ./(N * eb_n0);                  % Variance de bruit complexe en bande de base (avec energie du 2eme filtre)
TEB = zeros(size(eb_n0));                               % Tableau des TEB (resultats)
TEB_2 = zeros(size(eb_n0));                             % Tableau des TEB (resultats du 2eme filtre)
Pb= qfunc(sqrt(2*eb_n0));                               % Tableau des probabilités d'erreurs théorique
window=ones(nfft,1);                                    % Fenêtre rectangulaire pour le calcul de la DSP


% %% -----------------------------------------------------------------------------
% %                     Debut de la simulation pour le filtre de mise en
% %                     forme porte
% % ------------------------------------------------------------------------------
% 
% 
% for i=1:length(eb_n0)                                   % Début de la boucle principale de calcul
%     
%     %Initiation des compteurs
%     error_cnt=0;                                        % Compteur d'erreurs
%     bit_cnt = 0;                                        % Compteur de bits transmis
%     
%     while error_cnt <100                                % On boucle tant qu'on a pas atteint 100 erreurs, TEB probable si il est obtenu avec 100 erreurs min
%        
%         %% ----------------    Emetteur  --------------------------------------------------
%         
%         s_b_TX=randi([0,3],5000,1);                         % Generation d'une sequence de bit uniformement repartie (transformé en base 10)
%         s_s=pskmod(s_b_TX,M,pi/M,'gray');                   % Generation des symboles
%         s_s_sur_ech=upsample(s_s,fse);                      % Sur-echantillonage 
%         s_l=conv(s_s_sur_ech,g_t);                          % Filtrage par le filtre de mise en forme
%         
%         %% -----------------    Canal    --------------------------------------------------
%         nl = sqrt(sigma2(i)/2) * (randn(size(s_l)) +1j*randn(size(s_l)));           % Bruit blanc gaussien centré (BBGC)
%        
%         
%         %% -----------------  Recepteur  -------------------------------------------------
%         
%         s_l_b=s_l + nl;                                                             % Ajout du BBGC
%         r_l=conv(g_t_a,s_l_b);                                                      % Filrage par le filtre adapté
%         r_l_ech=r_l(length(g_t):fse :Ns * 10);                                      % Sous_echantillonnage du signal r_l 
%         s_b_RX=pskdemod(r_l_ech,M,pi/M,'gray');                                     % Demodulation des symboles: sequence de bits reçus
%         
%         s_b2_RX=de2bi(s_b_RX);                                                      % Conversion decimale/binaire de la sequence reçue                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
%         s_b2_TX=de2bi(s_b_TX);                                                      % Conversion decimale/binaire de la sequence émise
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
%         
%         %% calcul de proba d'erreur binaire et incrementation des compteurs
%         
%         matrice_erreur=sum(abs(s_b2_TX - s_b2_RX));                                 % Matrice des différences entre les deux sequences TX et RX
%         pb=(matrice_erreur(1) + matrice_erreur(2));                                 % Probabilité d'erreur binaire
%         
%         error_cnt=error_cnt + pb;
%         bit_cnt = bit_cnt + 2*Ns;
%        
%        
%     end
%     
%     TEB(i)=error_cnt/bit_cnt;                                                       % Remplissage du tableau contenant les TEB
%     
% end



%% -----------------------------------------------------------------------------
%                     Debut de la simulation pour le filtre de mise en
%                     forme racine de cosinus surelevé
% ------------------------------------------------------------------------------


for i=1:length(eb_n0)                                   % Début de la boucle principale de calcul
    
    %Initiation des compteurs
    error_cnt=0;                                        % Compteur d'erreurs
    bit_cnt = 0;                                        % Compteur de bits transmis
    
    while error_cnt <100                                % On boucle tant qu'on a pas atteint 100 erreurs, TEB probable si il est obtenu avec 100 erreurs min
       
        %% ----------------    Emetteur  --------------------------------------------------
        
        s_b_TX=randi([0,3],1,5000);                         % Generation d'une sequence de bit uniformement repartie (transformé en base 10)
        s_s=pskmod(s_b_TX,M,pi/M,'gray');                   % Generation des symboles
        %s_s=(association_bit_symbole(s_b_TX))';
        s_s_sur_ech=upsample(s_s,fse);                      % Sur-echantillonage 
        s_l_2=conv(s_s_sur_ech,g_t_2);                      % Filtrage par le filtre de mise en forme
        
        %% -----------------    Canal    --------------------------------------------------
        nl = sqrt(sigma2_2(i)/2) * (randn(size(s_l_2)) +1j*randn(size(s_l_2)));           % Bruit blanc gaussien centré (BBGC)
        %h = sinc((-12:12)-2) + sinc((-12:12)-5) ;   % chaque sinc presente un trajet
        h=1;
        %h = h.*hann(length(h)).';
        y=conv(h,s_l_2);
        %nl = randn(1,length(y));
        %nl = 0;
        
        %% -----------------  Recepteur  -------------------------------------------------
        
        s_l_b=y + nl;    % Ajout du BBGC
        r_l=conv(g_t_a_2,s_l_b);                                                    % Filrage par le filtre adapté
        r_l_ech=r_l(length(g_t_2):fse :Ns * fse+length(g_t_2)-1);                    % Sous_echantillonnage du signal r_l 
        %r_l_ech=r_l(fse:fse :length(s_l_b)); 
        s_b_RX=pskdemod(r_l_ech,M,pi/M,'gray');                                     % Demodulation des symboles: sequence de bits reçus
        %s_b_RX = association_symbole_bit(r_l_ech(1:5000));
        
        s_b2_RX=de2bi(s_b_RX);                                                      % Conversion decimale/binaire de la sequence reçue                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
        s_b2_TX=de2bi(s_b_TX);                                                      % Conversion decimale/binaire de la sequence émise
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
        
        %% calcul de proba d'erreur binaire et incrementation des compteurs
        
        matrice_erreur=sum(abs(s_b2_TX - s_b2_RX));                                 % Matrice des différences entre les deux sequences TX et RX
        pb=(matrice_erreur(1) + matrice_erreur(2));                                 % Probabilité d'erreur binaire
        
        error_cnt=error_cnt + pb;
        bit_cnt = bit_cnt + 2*Ns;
       
       
    end
    
    TEB_2(i)=error_cnt/bit_cnt;                                                     % Remplissage du tableau contenant les TEB
    
end

%% calcul de la DSP
% 
% [DSP,f] = pwelch(s_l,window,0,nfft,'centered',fe);                                  % Densite Spectrale de Puissance experimentale filtre porte
[DSP_2,f_2] = pwelch(r_l,window,0,nfft,'centered');                            % Densite Spectrale de Puissance experimentale filtre racine cos surélevé

%% Efficacité spectrale

%ef = 

%% -------------------------------  Figures  ------------------------------------
%--------------------------------------------------------------------------------

figure,
semilogy(eb_n0_db,TEB_2,'k','linewidth',2);                                           % Tracé du TEB experimentale en fonction de Eb/N0

% hold on;
% semilogy(eb_n0_db,TEB_2,'r','linewidth',2);                                         % Tracé du TEB experimentale avec filtre cos surélevé en fonction de Eb/N0

hold on;
semilogy(eb_n0_db,Pb,'g','linewidth',2);                                            % Tracé du TEB theorique en fonction de Eb/N0

legend('TEB empirique filtre racine cos surélevé','TEB théorique');
title('TEB en fonction du rapport erreur sur bruit EB/N0');
xlabel('rapport signal / bruit (dB)');
ylabel('TEB');
grid();
% 
% figure(2);
% semilogy(f,pow2db(DSP),'r','linewidth',2);                                          % Tracé de la Densité Spectrale de Puissance experimentale filtre porte
% 
% hold on;
% semilogy(f_2,pow2db(DSP_2),'b','linewidth',2);                                      % Tracé de la Densité Spectrale de Puissance theorique filtre racine cos surélevé
% 
% title('Représentation de la DSP théorique et Empirique');
% xlabel('Frequence (Hz)');
% ylabel('DSP (dB/Hz)');
% legend('DSP Empirique filtre porte', 'DSP Empirique filtre racine cos surélevé');
% grid();



figure,
plot(abs(conv(h,conj(h(end:-1:1)))))
figure,
stem(h)
figure,
semilogy(f_2,DSP_2);  
title('DSP Empirique filtre racine cos surélevé')
fvtool(h)
