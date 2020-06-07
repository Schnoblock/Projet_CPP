#include <iostream>
#include <cmath>
#include "Vecteur.h"
#include "Matrice.h"
#include "Tenseur.h"
#include "TenseurSVD.h"


using namespace std;


TenseurSVD :: TenseurSVD(){};

TenseurSVD :: TenseurSVD(int *tab, int size, Matrice *TablMat):Tenseur(tab, size){
    TableauMatrice = new Matrice[size];
    for (int i = 0; i < size; i ++){
        TableauMatrice[i] = TablMat[i];
    }
};

TenseurSVD :: TenseurSVD(int *tab, int size, Vecteur Vect, Matrice *TablMat):Tenseur(tab, size, Vect){
    TableauMatrice = new Matrice[size];
    for (int i = 0; i < size; i ++){
        TableauMatrice[i] = TablMat[i];
    }
};

TenseurSVD :: TenseurSVD(int *tab, int size, int k, Matrice Mat, Matrice *TablMat):Tenseur(tab, size, k, Mat){
    TableauMatrice = new Matrice[size];
    for (int i = 0; i < size; i ++){
        TableauMatrice[i] = TablMat[i];
    }
};

TenseurSVD :: ~TenseurSVD(){
    delete [] TableauMatrice;
};

TenseurSVD :: TenseurSVD(const TenseurSVD & Tens):Tenseur(Tens){
    TableauMatrice = new Matrice[ordre];
    for (int i = 0; i < ordre; i ++){
        TableauMatrice[i] = Tens.TableauMatrice[i];
    }
};

TenseurSVD & TenseurSVD :: operator = (const TenseurSVD & Tens){
    if (this !=&Tens){
        ordre = Tens.ordre;
        nbelts = Tens.nbelts;
        d = new int[ordre];
        for (int k = 0; k < ordre; k ++){
            d[k] = Tens.d[k];
        }
        for (int i=0; i < nbelts; i++){
            vect_tens[i] = Tens.vect_tens.tab[i];
        }
        TableauMatrice = new Matrice[Tens.ordre];
        for (int i = 0; i < ordre; i ++){
            TableauMatrice[i] = Tens.TableauMatrice[i];
        }
    }
    return *this;
};

Tenseur pmod(Tenseur Tens, Matrice M, int k);

//fonction membre tenseurtotal sera egalement implementee : celle-ci renvoit un objet de la classe Tenseur contenant le produit du coeur de tenseur avec les diferents facteurs (c'est-a-dire le tenseur en forme non factorisee)
Tenseur TenseurSVD :: tenseurtotal(){
    Tenseur S(d, ordre, vect_tens);
    for (int k = 0; k < ordre; k++){
        Matrice U_k (TableauMatrice[k]);
        Tenseur prod_k(pmod(S,U_k,k = k+1));
        S = prod_k;
    }
    return S;
};
