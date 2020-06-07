#ifndef MATRICE_H_INCLUDED
#define MATRICE_H_INCLUDED

//
//CLASSE MATRICE :
//Classe qui repressente un tesnseur d'ordre 2
//

class Matrice{

    Vecteur *mat;

    int dims[2];

    public:
        //Fonction affiche sans argument ni type de retour, qui affiche les dimensions
        //et les coefficients de la matrice
        void affiche();

        int get_dims(int x);

        //
        Matrice();

        //Constructeur prenant en arguments deux entiers correspondant au nombre
        //de lignes et de colonnes de la matrice a representer et construisant un
        //tableau rempli de vecteurs nuls
        Matrice(int rows, int cols);


        //GOOD
        //Constructeur prenant en arguments un objet de type Vecteur et construisant
        //une matrice carree diagonale dont les elements diagonaux seront donnees
        //par les composantes de ce vecteur
        Matrice(Vecteur v2);

        //GOOD
        //Constructeur prenant en arguments un tableau (de type Vecteur *) ainsi
        //que sa taille et construisant une matrice dont les colonnes seront donnees
        //par les vecteurs du tableau
        Matrice(Vecteur *mat2, int rows, int cols);

        //GOOD
        //Destructeur
        virtual ~Matrice();

        //GOOD
        //Constructeur de recopie
        Matrice(const Matrice & M);

        //JAMAIS TESTE
        //Surdefinition de l’operateur d’affectation
        Matrice & operator = (const Matrice &M);

        //GOOD
        //Surdefinition de l’operateur [ ] afin de pouvoir acceder a chacun des
        //objets de type Vecteur du tableau mat
        Vecteur & operator[] (int j);

        //GOOD
        //Surdefinition des operateurs +, - et * pour coder respectivement la somme,
        //la difference et le produit de deux matrices de tailles compatibles
        Matrice operator+ (Matrice & M2);

        //GOOD

        Matrice operator- (Matrice &M2);

        // GOOD

        Matrice operator* (Matrice &M2);


        // ??

        Matrice operator* (float k);

        // GOOD
        //Fonction mvprod prenant en argument un objet de type Vecteur dont la
        //dimension correspond au nombre de colonnes de la matrice appelante, et
        //retournant un Vecteur correspondant au produit de la matrice appelante
        //avec le vecteur en argument
        Vecteur mvprod(Vecteur V);

        //GOOD
        //Fonction transpose permettant de construire un nouvel objet de type Matrice
        //correspondant a la transposee de la matrice appelante (on rappelle que si
        //B = AT, alors Bi,j = Aj,i pour tout couple (i, j) indexant un element de B)
        Matrice transpose();

        //GOOD
        //Fonction submat prenant en argument quatre entiers il ≤ jl et ic ≤ jc,
        //et renvoyant un objet de type Matrice dont le tableau sera forme par les
        //lignes il a jl et les colonnes ic a jc de la matrice appelante
        Matrice submat(int il, int jl, int ic, int jc);


        //Fonctions amies :

        //Fonction norm
        friend float norm(Matrice M);

        //Surdefinition de l’operateur *
        //friend Matrice operator* (const Matrice &M2, float k);

        //Fonction outer
        friend Matrice outer(Vecteur v1, Vecteur v2);

        friend class Vecteur;

        friend class Tenseur;

        friend Matrice Id(int n);

        friend Matrice reductridiag(Matrice &D);

        friend Matrice qrsym(Matrice &A, Matrice &Q);

        friend bool test_diag(Matrice D);

        friend void p_q(Matrice T,int &p,int &q);

        friend bool test_nulle(Matrice M);

        friend Matrice qrpivot(Matrice A, Matrice &Q);

        friend void svd(Matrice A, Matrice *TableauMatrice);


};

float norm(Matrice M);

//renvoie l'outer de v1 et v2 : v1v2T
Matrice outer(Vecteur v1, Vecteur v2);

//Renvoie l'indicatrice de taille n
Matrice Id(int n);

#endif // MATRICE_H_INCLUDED
