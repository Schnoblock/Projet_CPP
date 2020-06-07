#ifndef TENSEUR_H_INCLUDED
#define TENSEUR_H_INCLUDED

//
//CLASSE TENSEUR :
//Classe qui represente un tenseur d'ordre aribtraire de reels au format float
//


class Tenseur{
    protected :
        //Memberes donnees

        //ordre du tenseur
        int ordre;

        //tableau contenant les valeurs des d dimensions du tenseur
        int *d;

        //Produit des dimensions du tenseur

        int nbelts;

    //Version vectorisee du tenseur, objet de type Vecteur representant un tableau de taille nbelts contenant les elements du tenseur dans l’ordre defini en (4.1) par la fonction Phi;

    Vecteur vect_tens;

    //Fonctions membres
    public :

        void affiche_vect();

        Tenseur();

        //Constructeur prenant en arguments un tableau d’entiers ainsi que sa taille, et construisant un tenseur d’ordre et de dimensions correspondantes dont tous les coefficients seront initialises a 0
        Tenseur(int *tab, int size);

        //Constructeur prenant en arguments un tableau d’entiers, sa taille d, ainsi qu’un objet de type Vecteur et construisant un tenseur d’ordre d dont la version vectorisee sera initialisee avec l’objet Vecteur
        Tenseur(int *tab, int size, Vecteur Vect);

        //Constructeur prenant en arguments un tableau d’entiers, sa taille d, un entier k entre 1 et d et un objet de classe Matrice representant une matrice A : ce constructeur construit un tenseur d’ordre d dont la version vectorisee est initialisee avec les coefficients de A de sorte que A represente le k-ieme mode du tenseur
        Tenseur(int *tab, int size, int k, Matrice Mat);

        //Destructeur
        virtual ~Tenseur();

        //Constructeur de recopie
        Tenseur(const Tenseur & Tens);

        //Surcharge de l'opérateur d'affectation
        Tenseur & operator = (const Tenseur &Tens);

        //Surcharge de l’operateur [ ] pour qu’il permette d’acceder aux coefficients de la version vectorisee du tenseur
        float & operator[](int i);

        //Surcharge des operateurs + et - pour qu’ils permettent d’additionner ou de soustraire deux tenseurs de memes dimensions (et de meme ordre)
        Tenseur operator- (Tenseur &Tens2);
        Tenseur operator+ (Tenseur &Tens2);
        //Fonction membre mode prenant en argument un entier k entre 1 et d et renvoyant un objet de classe Matrice representant le k-i eme mode du tenseur appelant
        Matrice mode (int k);

        friend class Matrice;

        friend class Vecteur;

        friend int indice_i_T(int i, Tenseur T);

        friend Tenseur pmod(Tenseur Tens, Matrice M, int k);

};

Vecteur phi_inv(int i, int ordre, int *d);

int prod_n_i(int t, int ordre, int *d);

//retourne ph(i_1,..i_d) sous forme d'int
int phi(int *ind, int ordre, int *d);

int phi_v(Vecteur ind, int ordre, int *d);

#endif // TENSEUR_H_INCLUDED
