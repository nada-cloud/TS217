function compteur_erreur=Erreur(bits_recus,S_b)
M = 4;
n_b = log2(M);
Nombre_symbole=5000;
Nombre_bits=Nombre_symbole*n_b;

compteur_erreur=0;
for i=1:Nombre_bits
    if S_b(i) ~= bits_recus(i)
        compteur_erreur = compteur_erreur+1;
    end
end