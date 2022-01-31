clear;
close all;
clc;

%% Parametres
Fse = 4;
N_fft=512;
M = 4;
n_b = log2(M);
Nombre_symbole=5000;
Nombre_bits=Nombre_symbole*n_b;
S_b=randi([0,1],1,Nombre_bits);
alpha = 0.35;
span = 8;
p_n = 25;

%% Emetteur
g = rcosdesign(alpha, span, Fse, 'sqrt');
S_s=Emetteur(S_b,g,Fse);
%S_s = pskmod(S_b, M, pi/M, 'gray');

%% Canal
h = sinc((-12:12)-1.7) + sinc((-12:12)-5.7);    %délai entre les deux chemins augmenté : BC est étroite
h_2 = sinc((-12:12)-1.3) + sinc((-12:12)-5.7); 
h_2 = h_2.*hann(length(h_2)).';
h = h.*hann(length(h)).';                   %etalement difference entre les peaks 
fvtool(h)
fvtool(h_2)
%s=upsample(S_s,Fse);
s=S_s;
S_l=conv(s,g);

%% echantillonnage pour un res plus rigoureux
%% Echantillonnage
% S_l=upsample(S_l,Fse);

%% Filtre racine cos surélevé
S = conv(S_l, h);
S_2 = conv(S_l, h_2);

%% Reception
%S = S ;
S_2 = S_2 + randn(1, length(S_2));

r_l = conv(S,g);
r_l_ech = r_l(Fse:Fse:length(S));

r_l_2 = conv(S_2,g);
r_l_ech_2 = r_l_2(Fse:Fse:length(S_2));

[Pxx, frequences] = pwelch(r_l,N_fft,0,N_fft,'centered');
[Pxx_2, frequences_2] = pwelch(r_l_2,N_fft,0,N_fft,'centered');

%% Figures
% figure, plot(g)
% axis([0 100 0 0.6])
% 
% figure, plot(abs(S_l))
% axis([0 100 0 0.6])

% figure, plot(abs(S))
% axis([0 100 0 0.4])

% plot(abs(r_l_ech))
% axis([0 100 0 1.5])

figure, plot(abs(conv(h, conj(h(end:-1:1)))))


figure
subplot(2,1,1);
plot(frequences,20*log(Pxx));
xlabel('frequence (Hz)');
ylabel('magnitude (dcb)');
title('DSP théorique du filtre racine cosinus surélevé sans bruit');

subplot(2,1,2); 
plot(frequences_2,20*log(Pxx_2));
xlabel('frequence (Hz)');
ylabel('magnitude (dcb)');
title('DSP théorique du filtre racine cosinus surélevé avec bruit');
