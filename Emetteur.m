function S_l=Emetteur(S_b,g,Fse)

M = 4;
n_b = log2(M);
Nombre_symbole=5000;
Nombre_bits=Nombre_symbole*n_b;
S_s=zeros(1,Nombre_symbole);

i=1;
k=1;
while(i<=Nombre_bits) 
    if S_b(i)==0 && S_b(i+1)==0
        S_s(k)=(1/sqrt(2))*(1+1j);
    elseif S_b(i)==0 && S_b(i+1)==1
        S_s(k)=(1/sqrt(2))*(-1+1j);   
    elseif S_b(i)==1 && S_b(i+1)==0
        S_s(k)=(1/sqrt(2))*(1-1j);
    elseif S_b(i)==1 && S_b(i+1)==1
        S_s(k)=-(1/sqrt(2))*(1+1j);
    end
    i=i+2;
    k=k+1;
end
%retourne S_s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%filtre de mise en forme
s=upsample(S_s,Fse);
S_l=conv(s,g);

end