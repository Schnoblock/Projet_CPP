#ifndef VECTEUR_H_INCLUDED
#define VECTEUR_H_INCLUDED

//
//CLASSE VECTEUR :
//Classe qui represente un tenseur d'ordre 1 de reels au format float
//


class Vecteur{
    //Membres donnees

    int dim;
    float * tab;


    public:
        //Fonction affiche sans argument ni valeur de retour qui affiche les valeurs
        //de dim ainsi que les elements de tab
        void affiche();

        //Constructeur
        Vecteur();

        //Constructeur avec un argument entier qui cree un vecteur representant
        //un tableau de reels au format float de la taille correspondant a l’argument,
        //ou les coefficients seront initalises a 0
        Vecteur(int n);

        //Constructeur avec deux arguments, respectivement de type float * et int,
        //qui cree un vecteur representant un tableau de reels au format float de
        //la taille correspondant au second argument, ou les coefficients du tableau
        //seront ceux du tableau passe en premier argument
        Vecteur(float *val, int n);

        //Destructeur
        virtual ~Vecteur();

        //Constructeur de recopie
        Vecteur(const Vecteur & v);

        //Surdefinition de l’operateur d’affectation
        Vecteur & operator = (const Vecteur & v);

        //Surdefinition des operateurs + et - pour effectuer la somme et la soustraction de deux vecteurs de meme taille
        Vecteur operator+ (Vecteur &v2);

        Vecteur operator- (Vecteur &v2);

        //Surdefinition de l’operateur [ ] pour pouvoir acceder aux elements du
        //tableau de facon a pouvoir les modifier
        float & operator[](int i);

        //Fonction subvec prenant en argument deux entiers i ≤ j et renvoyant un
        //nouvel objet de type Vecteur dont le tableau contiendra les composantes
        //i, i + 1, . . . , j du tableau du vecteur appelant
        Vecteur subvec(int i, int j);

        //Fonctions amies :

        //Fonction sum (non demandée mais utile pour le calcul du produit scalaire)
        friend float sum(Vecteur v);

        //Fonction dot
        friend float dot(Vecteur v1, Vecteur v2);

        //Surcharge de l'operateur *
        friend Vecteur operator* (const Vecteur &v, float k);
        friend Vecteur operator* (float k, const Vecteur &v);

        //Fonction norm
        friend float norm(Vecteur v);

        friend Vecteur householder(Vecteur x, float &beta);

        friend class Matrice;

        friend float max(Vecteur v);

        friend class Tenseur;

        friend class TenseurSVD;
};

#endif // VECTEUR_H_INCLUDED
