#include <iostream>
#include <cmath>
#include "Vecteur.h"
#include "Matrice.h"
#include "Tenseur.h"


using namespace std;

//
//CLASSE TENSEUR :
//Classe qui represente un tenseur d'ordre aribtraire de reels au format float
//


//retourne phi_inv(i) = (i1,...,id) sous fome de vecteur /!\ ce sont des float
Vecteur phi_inv(int i, int ordre, int *d){
    int ind[ordre];
    int f[ordre];
    Vecteur ind_vect(ordre);

    f[ordre-1] = i;
    ind[ordre-1] = f[ordre - 1] % d[ordre - 1];
    if (ind[ordre-1] == 0){ind[ordre-1] = d[ordre - 1];}

    for (int t = ordre - 2; t >= 0; t--){
        f[t] = ((f[t+1] - ind[t+1])/(d[t+1])) + 1;
        ind[t] = f[t] % d[t];
        if (ind[t] == 0){ind[t] = d[t];}
    }

    for (int i = 0; i < ordre; i++){
        ind_vect[i] = ind[i];
    }
    return ind_vect;
};

int prod_n_i(int t, int ordre, int *d){
    int prod_n = 1;
    for (int j = ordre - 1; j >= t; j--){
        prod_n = prod_n * d[j];
    }
    return prod_n;
};

//retourne ph(i_1,..i_d) sous forme d'int
int phi(int *ind, int ordre, int *d){
    int i = ind[ordre - 1];
    for (int t = ordre - 2; t >= 0; t--){
        i = i + prod_n_i(t+1, ordre, d)*(ind[t] - 1);
    }
    return i;
};

Tenseur:: Tenseur(){};

//Constructeur prenant en arguments un tableau d’entiers ainsi que sa taille, et construisant un tenseur d’ordre et de dimensions correspondantes dont tous les coefficients seront initialises a 0
Tenseur:: Tenseur(int *tab, int size){
    d = new int[size];
    ordre = size;
    nbelts = 1;
    for (int i = 0; i < size; i++){
        nbelts = nbelts * tab[i];
        d[i] = tab[i];
    }
    vect_tens = Vecteur(nbelts);
    for (int i = 0; i < nbelts; i++){
        vect_tens.tab[i] = 0;
    }
};
void Tenseur :: affiche_vect(){
    vect_tens.affiche();
};

//Constructeur prenant en arguments un tableau d’entiers, sa taille d, ainsi qu’un objet de type Vecteur et construisant un tenseur d’ordre d dont la version vectorisee sera initialisee avec l’objet Vecteur
Tenseur:: Tenseur(int *tab, int size, Vecteur Vect){
    d = new int[size];
    ordre = size;
    nbelts = 1;
    for (int i = 0; i < ordre; i++){
        nbelts = nbelts * tab[i];
        d[i] = tab[i];
    }
    vect_tens = Vecteur(nbelts);
    for (int i = 0; i < nbelts; i++){
        vect_tens[i] = Vect.tab[i];
    }
};

//Constructeur prenant en arguments un tableau d’entiers, sa taille d, un entier k entre 1 et d et un objet de classe Matrice representant une matrice A : ce constructeur construit un tenseur d’ordre d dont la version vectorisee est initialisee avec les coefficients de A de sorte que A represente le k-ieme mode du tenseur
Tenseur:: Tenseur(int *tab, int size, int k, Matrice Mat){
    d = new int[size];
    ordre = size;
    nbelts = 1;
    for (int i = 0; i < ordre; i++){
        d[i] = tab[i];
        nbelts = nbelts * d[i];
    }
    vect_tens = Vecteur(nbelts);

    int nbelts_moins_k = nbelts/d[k-1];

    int d_moins_k_2[ordre - 1];
    for (int i=0; i < k; i++){
        d_moins_k_2[i] = d[i];
    }
    for (int i = k - 1; i < ordre - 1; i++){
        d_moins_k_2[i] = d[i+1];
    }

    for (int i = 0; i < d[k-1] ; i ++){
        for (int j = 0; j < nbelts_moins_k; j++){
            int i_k = i;
            int j_k = j;

            Vecteur I_moins_ik;
            int I[ordre];
            I_moins_ik = Vecteur(ordre - 1);

            I_moins_ik = phi_inv(j_k+1, ordre-1 ,d_moins_k_2);

            for (int l = 0; l < k-1; l++){
                I[l] = I_moins_ik[l];
            }
            for (int l = k; l < ordre; l++){
                I[l] = I_moins_ik[l-1];
            }
            I[k - 1] = i_k;

            vect_tens[phi(I,ordre,d)+1] = Mat[j][i];
        }
    }
};

//Destructeur
Tenseur:: ~Tenseur(){
    delete [] d;
};

//Surcharge de l’operateur [ ] pour qu’il permette d’acceder aux coefficients de la version vectorisee du tenseur
float & Tenseur:: operator[](int i){
    return vect_tens[i];
};

//Constructeur de recopie
Tenseur:: Tenseur(const Tenseur & Tens){
    ordre = Tens.ordre;
    nbelts = Tens.nbelts;
    d = new int[ordre];
    for (int k = 0; k < ordre; k ++){
        d[k] = Tens.d[k];
    }
    for (int i=0; i < nbelts; i++){
        vect_tens[i] = Tens.vect_tens.tab[i];
    }
};

//Surcharge de l'opérateur d'affectation
Tenseur & Tenseur:: operator = (const Tenseur &Tens){
    if (this !=&Tens){
        //delete [] d;
        ordre = Tens.ordre;
        nbelts = Tens.nbelts;
        d = new int[ordre];
        for (int k = 0; k < ordre; k ++){
            d[k] = Tens.d[k];
        }
        for (int i=0; i < nbelts; i++){
            vect_tens[i] = Tens.vect_tens.tab[i];
        }
    }
    return *this;
};

//Surcharge des operateurs + et - pour qu’ils permettent d’additionner ou de soustraire deux tenseurs de memes dimensions (et de meme ordre)
Tenseur Tenseur :: operator+ (Tenseur &Tens2){
    Tenseur Tens (Tens2.d,Tens2.ordre);
    for (int i =0; i< Tens2.nbelts; i++){
        Tens.vect_tens[i] = vect_tens[i] + Tens2[i];
    }
    return(Tens);
};

Tenseur Tenseur :: operator- (Tenseur &Tens2){
    Tenseur Tens (Tens2.d,Tens2.ordre);
    for (int i =0; i< Tens2.nbelts; i++){
        Tens.vect_tens[i] = vect_tens[i] - Tens2[i];
    }
    return(Tens);
};
//Fonction membre mode prenant en argument un entier k entre 1 et d et renvoyant un objet de classe Matrice representant le k-i eme mode du tenseur appelant

Matrice Tenseur:: mode (int k){
    int nbelts_moins_k = nbelts/d[k-1];
    Matrice Mode_k (d[k-1],nbelts_moins_k);
    Vecteur I_moins_ik;
    I_moins_ik = Vecteur(ordre - 1);

    for (int t = 0; t < nbelts + 1; t++){
        Vecteur I_t = phi_inv(t, ordre, d);
        int I_t_moins_ik[ordre - 1];
        for (int i=0; i < k; i++){
            I_t_moins_ik[i] = I_t[i];
        }
        for (int i = k - 1; i < ordre - 1; i++){
            I_t_moins_ik[i] = I_t[i+1];
        }

        int d_moins_k_2[ordre - 1];
        for (int i=0; i < k; i++){
            d_moins_k_2[i] = d[i];
        }
        for (int i = k - 1; i < ordre - 1; i++){
            d_moins_k_2[i] = d[i+1];
        }

        int line = I_t[k-1];
        int col = phi(I_t_moins_ik,ordre-1,d_moins_k_2);

        Mode_k[col-1][line-1] = vect_tens[t-1];
    }
    return Mode_k;
};
