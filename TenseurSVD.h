#ifndef TENSEURSVD_H_INCLUDED
#define TENSEURSVD_H_INCLUDED

//
//CLASSE TENSEUR :
//Classe qui represente un tenseur d'ordre aribtraire de reels au format float
//


class TenseurSVD : public Tenseur{
    //Memberes donnees
    protected :

        Matrice *TableauMatrice;


    //Fonctions membres
    public :

        TenseurSVD();

        TenseurSVD(int *tab, int size, Matrice *TablMat);

        TenseurSVD(int *tab, int size, Vecteur Vect,Matrice *TableauMatrice);

        TenseurSVD(int *tab, int size, int k, Matrice Mat, Matrice *TableauMatrice);

        virtual ~TenseurSVD();

        TenseurSVD(const TenseurSVD & Tens);

        TenseurSVD & operator = (const TenseurSVD &Tens);

        Tenseur tenseurtotal();

        friend class Matrice;

        friend class Vecteur;

        friend int indice_i_T(int i, Tenseur T);

        friend TenseurSVD hosvd(Tenseur T);

};

#endif // TENSEUR_H_INCLUDED
