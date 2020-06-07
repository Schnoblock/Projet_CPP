#include <iostream>
#include <cmath>
#include "Vecteur.h"
#include "Matrice.h"

using namespace std;

//Fonction affiche sans argument ni type de retour, qui affiche les dimensions
//et les coefficients de la matrice
void Matrice:: affiche(){
    std::cout << "dims = " << dims[0] << "," << dims[1] << std::endl;
    std::cout << "mat = " << std::endl;
    for (int i = 0; i < dims[0] ; i++){
        std::cout << "[ ";
        for (int j = 0; j < dims[1] ; j++){
            std::cout << mat[j][i] << " ";
        }
        std::cout << "]" << std::endl;
    }
};

//
Matrice:: Matrice(){};

int Matrice :: get_dims(int x){
    if (x == 0){return dims[0];}
    if (x == 1){return dims[1];}
};

//Constructeur prenant en arguments deux entiers correspondant au nombre
//de lignes et de colonnes de la matrice a representer et construisant un
//tableau rempli de vecteurs nuls
Matrice:: Matrice(int rows, int cols){
    dims[0] = rows;
    dims[1] = cols;
    mat = new Vecteur[cols];
    for (int j = 0; j < cols; j++){
        mat[j] = Vecteur(rows);
    }
};


//Constructeur prenant en arguments un objet de type Vecteur et construisant
//une matrice carree diagonale dont les elements diagonaux seront donnees
//par les composantes de ce vecteur
Matrice:: Matrice(Vecteur v2){
    dims[0] = v2.dim;
    dims[1] = v2.dim;
    mat = new Vecteur[v2.dim];
    for (int j = 0; j < v2.dim; j++){
        mat[j] = Vecteur(v2.dim);
    }
    for (int i=0; i < v2.dim; i++){
        mat[i][i] = v2.tab[i];
    }
};

//Constructeur prenant en arguments un tableau (de type Vecteur *) ainsi
//que sa taille et construisant une matrice dont les colonnes seront donnees
//par les vecteurs du tableau
Matrice:: Matrice(Vecteur *mat2, int rows, int cols){
    dims[0] = rows;
    dims[1] = cols;
    for(int j = 0; j < cols; j++){
        mat[j] = mat2[j];
    }
};

//Destructeur
Matrice:: ~Matrice(){delete [] mat;};

//Constructeur de recopie
Matrice:: Matrice(const Matrice & M){
    dims[0] = M.dims[0];
    dims[1] = M.dims[1];
    mat = new Vecteur[dims[1]];
    for (int j = 0 ; j < dims[1]; j ++){
        mat[j] = Vecteur(dims[0]);
        for (int i = 0 ; i < dims[0]; i ++){
            mat[j][i] = M.mat[j][i];
        }
    }
};

//Surdefinition de l’operateur d’affectation
Matrice & Matrice:: operator = (const Matrice &M){
    if (this !=&M){
        //delete [] mat;
        dims[0] = M.dims[0];
        dims[1] = M.dims[1];
        mat = new Vecteur [dims[1]];
    for( int i = 0; i < dims[1] ; i++){mat [i] = M.mat [i];}
    }
    return *this ;
};

//Surdefinition de l’operateur [ ] afin de pouvoir acceder a chacun des
//objets de type Vecteur du tableau mat
Vecteur & Matrice:: operator[] (int j){
    return mat[j];
};

//Surdefinition des operateurs +, - et * pour coder respectivement la somme,
//la difference et le produit de deux matrices de tailles compatibles
Matrice Matrice::operator+ (Matrice & M2){
    Matrice M (M2.dims[0], M2.dims[1]);
    for (int i = 0; i<dims[0]; i++){
        for (int j = 0; j < dims[1]; j++){
            M[j][i] = mat[j][i] + M2[j][i];
        }
    }
    return M;
};


Matrice Matrice::operator- (Matrice &M2){
    Matrice M (M2.dims[0], M2.dims[1]);
    for (int i = 0; i<dims[0]; i++){
        for (int j = 0; j < dims[1]; j++){
            M[j][i] = mat[j][i] - M2[j][i];
        }
    }
    return M;
};


Matrice Matrice::operator* (Matrice &M2){
    Matrice M (dims[0], M2.dims[1]);
    for (int i = 0; i<dims[0]; i++){
        Vecteur line_i (dims[1]);
        for (int j = 0; j < dims[1]; j++){
            line_i[j] = mat[j][i];
        }
        for (int j = 0; j < M2.dims[1]; j++){
            M[j][i] = dot(line_i,M2[j]);
        }
    }
    return M;
};


Matrice Matrice::operator* (float k){
    Matrice Mprod (dims[0],dims[1]);
    for (int i = 0; i < dims[0]; i++){
        for(int j=0; j < dims[1]; j++){
            Mprod[j][i] = mat[j][i] * k;
        }
    }
    return Mprod;
};

//Fonction mvprod prenant en argument un objet de type Vecteur dont la
//dimension correspond au nombre de colonnes de la matrice appelante, et
//retournant un Vecteur correspondant au produit de la matrice appelante
//avec le vecteur en argument
Vecteur Matrice:: mvprod(Vecteur V){
    Vecteur Vprod (dims[1]);
    for (int i = 0; i < dims[0]; i++){
        Vecteur line_i (dims[1]);
        for (int j = 0; j < dims[0]; j++){
            line_i[j] = mat[j][i];
        }
        Vprod[i] = dot(line_i,V);
    }
    return Vprod;
};

//Fonction transpose permettant de construire un nouvel objet de type Matrice
//correspondant a la transposee de la matrice appelante (on rappelle que si
//B = AT, alors Bi,j = Aj,i pour tout couple (i, j) indexant un element de B)
Matrice Matrice:: transpose(){
    Matrice M (dims[0], dims[1]);
    for (int i = 0; i<dims[0]; i++){
        for (int j = 0; j < dims[1]; j++){
            M[j][i] = mat[i][j];
        }
    }
    return M;
};

//Fonction submat prenant en argument quatre entiers il ≤ jl et ic ≤ jc,
//et renvoyant un objet de type Matrice dont le tableau sera forme par les
//lignes il a jl et les colonnes ic a jc de la matrice appelante
Matrice Matrice:: submat(int il, int jl, int ic, int jc){
    int newdims[2];
    newdims[0] = jl - il + 1;
    newdims[1] = jc - ic + 1;
    Matrice M (newdims[0], newdims[1]);
    for (int c = 0; c < newdims[1]; c++){
        Vecteur col_c (newdims[0]);
        for (int l = 0; l < newdims[0]; l++){
            col_c = mat[ic + c].subvec(il, jl);
        }
        M[c] = col_c;
    }
    return M;
};
