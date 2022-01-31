function bits_recus=recepteur(y_l,g_a,Fse)
fe=1e4;
M = 4;
n_b = log2(M);
periode_symbole=1e-3;
Nombre_symbole=5000;
Fs=(periode_symbole*fe);
Nombre_bits=Nombre_symbole*n_b;

r_l = conv(y_l,g_a);
%échantillonnage au rythme Ts
%r_l_ech = upsample(r_l_cor,1/(Fse*periode_symbole));
r_l_ech = r_l(Fse:Fs:length(y_l));
%figure(1)
%plot(temps,r_l_ech)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bloc de décision
A_n_estime=zeros(1,Nombre_symbole);

for(i=1:Nombre_symbole)
    if real(r_l_ech(i))>0 && imag(r_l_ech(i))>0
        A_n_estime(i)=(1/sqrt(2))*(1+j);
    elseif real(r_l_ech(i))>0 && imag(r_l_ech(i))<0
        A_n_estime(i)=(1/sqrt(2))*(1-j);
    elseif real(r_l_ech(i))<0 && imag(r_l_ech(i))>0
        A_n_estime(i)=(1/sqrt(2))*(-1+j);
    elseif real(r_l_ech(i))<0 && imag(r_l_ech(i))<0
        A_n_estime(i)=-(1/sqrt(2))*(1+j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bloc association bits symboles
bits_recus=zeros(1,Nombre_bits);

a=1;
b=1;
while(b<Nombre_bits)
    if A_n_estime(a)==(1/sqrt(2))*(1+j)
        bits_recus(b)=0;
        bits_recus(b+1)=0;
    elseif A_n_estime(a)==(1/sqrt(2))*(-1+j)
        bits_recus(b)=0;
        bits_recus(b+1)=1;  
    elseif A_n_estime(a)==(1/sqrt(2))*(1-j)
        bits_recus(b)=1;
        bits_recus(b+1)=0;
    elseif A_n_estime(a)==-(1/sqrt(2))*(1+j)
        bits_recus(b)=1;
        bits_recus(b+1)=1;
    end
    b=b+2;
    a=a+1;
end
%retourne bits recus