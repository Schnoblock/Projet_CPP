#include <iostream>
#include <cmath>
#include "Vecteur.h"
#include "Matrice.h"
#include "Tenseur.h"

using namespace std;

void Vecteur::affiche(){
    std::cout << "dim = " << dim << std::endl;
    std::cout << "tab = ";
    std::cout << "[";
    for (int i = 0; i < dim - 1; i++) {
        std::cout << tab[i] << " ";
    }
    std::cout << tab[dim-1];
    std::cout << "]";
    std::cout << std::endl;
};

//Constructeur
Vecteur::Vecteur(){};

//Constructeur avec un argument entier qui cree un vecteur representant
//un tableau de reels au format float de la taille correspondant a l’argument,
//ou les coefficients seront initalises a 0
Vecteur::Vecteur(int n):dim(n), tab(new float[n]){
    for (int i = 0; i < dim; i++){
        tab[i] = 0;
    }
};

//Constructeur avec deux arguments, respectivement de type float * et int,
//qui cree un vecteur representant un tableau de reels au format float de
//la taille correspondant au second argument, ou les coefficients du tableau
//seront ceux du tableau passe en premier argument
Vecteur::Vecteur(float *val, int n):dim(n), tab(new float[n]){
    for (int i = 0; i < dim; i++){
        tab[i] = val[i];
    }
};

//Destructeur
Vecteur:: ~Vecteur() {
    delete[] tab;
};

//Constructeur de recopie
Vecteur::Vecteur(const Vecteur & v){
    dim = v.dim ;
    tab = new float[dim];
    for (int i = 0; i < dim; i++){tab[i] = v.tab[i];}
};

//Surdefinition de l’operateur d’affectation
Vecteur &Vecteur :: operator = (const Vecteur & v){
    if (this !=&v){
        //delete [] tab;
        dim = v.dim ;
        tab = new float [dim];
        for( int i =0;i< dim ;i++){tab [i] = v.tab [i];}
    }
    return *this ;
};
//Surdefinition des operateurs + et - pour effectuer la somme et la soustraction de deux vecteurs de meme taille
Vecteur Vecteur :: operator+ (Vecteur &v2){
    Vecteur v (v2.dim);
    for (int i =0; i<dim; i++){
        v.tab[i] = tab[i] + v2.tab[i];
    }
    return(v);
};

Vecteur Vecteur :: operator- (Vecteur &v2){
    Vecteur v (v2.dim);
    for (int i =0; i<dim; i++){
        v.tab[i] = tab[i] - v2.tab[i];
    }
    return(v);
};

//Surdefinition de l’operateur [ ] pour pouvoir acceder aux elements du
//tableau de facon a pouvoir les modifier
float & Vecteur:: operator[](int i){
    return tab[i];
};

//Fonction subvec prenant en argument deux entiers i ≤ j et renvoyant un
//nouvel objet de type Vecteur dont le tableau contiendra les composantes
//i, i + 1, . . . , j du tableau du vecteur appelant
Vecteur Vecteur:: subvec(int i, int j){
    int newdim;
    float *newtab;
    newdim = j - i + 1;
    Vecteur subVecteur(newdim);
    for (int k = 0; k < newdim; k++){
        subVecteur[k] = tab[i + k];
    }
    return subVecteur;
};
