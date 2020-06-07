#include <iostream>
#include <cmath>

#include "Vecteur.h"
#include "Matrice.h"
#include "Tenseur.h"
#include "TenseurSVD.h"

using namespace std;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
};

//FONCTIONS VECTEUR

//Fonction qui somme toutes les valeurs d'un vecteur
float sum(Vecteur v){
    float somme = 0;
    for (int i=0; i < v.dim; i++){
        somme = somme + v[i];
    }
    return somme;
};


//Fonction dot qui prend en arguments deux objets de type Vecteur correspondant a deux vecteurs de meme taille et renvoit leur produit scalaire
float dot(Vecteur v1, Vecteur v2){
    float scalaire;
    Vecteur vprod (v2.dim);
    for (int i = 0; i < v2.dim; i++){
        vprod[i] = v1[i] * v2[i];
    }
    scalaire = sum(vprod);
    return scalaire;
};

//Surcharge de l’operateur * pour deux arguments de type float et Vecteur, qui renvoie un objet de type Vecteur dont le tableau correspondra a la multiplication du tableau du second argument par le reel en premier argument
Vecteur operator* (const Vecteur &v, float k){
    Vecteur vprodk (v.dim);
    for (int i = 0; i < v.dim; i++){
        vprodk[i] = v.tab[i] * k;
    }
    return vprodk;
};

Vecteur operator* (float k, const Vecteur &v){
    Vecteur vprodk (v.dim);
    for (int i = 0; i < v.dim; i++){
        vprodk[i] = v.tab[i] * k;
    }
    return vprodk;
};

//Fonction norm qui prend en argument un objet de type Vecteur et renvoit sa norme euclidienne, ou norme l2
float norm(Vecteur v){
    float normsquare;
    float norm;
    normsquare = dot(v,v);
    norm = sqrt(normsquare);
    return norm;
};


//FONCTIONS MATRICE

//Renvoie la norme de Frobenuis d'une matrice M (calculée comme la racine de la trace de Mt * M)
float norm(Matrice M){
    Matrice M2 (M.dims[0],M.dims[0]);
    M2 = M.transpose() * M;
    float doublesum=0;
    for (int i = 0; i < M.dims[0]; i++){
            doublesum = doublesum + M2[i][i];
    }
    return sqrt(doublesum);
};

//renvoie l'outer de v1 et v2 : v1v2T
Matrice outer(Vecteur v1, Vecteur v2){
    Matrice Mv1 (v1);
    Matrice Mv2 (v2);
    Matrice Mouter (Mv1.dims[0], Mv2.dims[1]);
    for (int i = 0; i < Mouter.dims[0]; i++){
        for(int j=0; j < Mouter.dims[1]; j++){
            Mouter[j][i] = v1[i] * v2[j];
        }
    }
    return Mouter;
};

//Renvoie l'indicatrice de taille n
Matrice Id(int n){
    Vecteur I(n);
    for (int i = 0; i < n; i++){
        I[i] = 1;
    }
    Matrice Id(I);
    return Id;
};

//FONCTIONS SVC

//Max d'un vecteur
float max(Vecteur v){
    int n = v.dim;
    float max = v[0];
    for (int i = 1; i < n-1; i++){
        if (v[i] > max){
            max = v[i];
        }
    }
    return max;
};

//Test si la matrice est diagonale
bool test_diag(Matrice D){
    int n = D.dims[0];
    bool diag = false;
    int p = 0;
    for (int line = 1; line < n-1; line++){
        for (int col = 0; col < line ; col++ ){
            if (abs(D[col][line]) >= 10e-6){
                p = p + 1;
            }
        }
    }
    for (int line = 0; line < n-2; line++){
        for (int col = line + 1; col < n - 1 ; col++ ){
            if (abs(D[col][line]) >= 10e-6){
                p = p + 1;
            }
        }
    }
    if (p == 0){
        diag = true;
    }
    return diag;
};

//Test si une matrice M est nulle
bool test_nulle(Matrice M){
    bool nulle = false;
    int rows = M.dims[0];
    int cols = M.dims[1];
    int p = 0;
    for (int i = 0; i < rows - 1; i++){
        for (int j = 0; j < cols - 1; j++){
            if (abs(M[j][i]) >= 10e-6){
                p = p + 1;
            }
        }
    }
    if (p == 0){
        nulle = true;
    }
    return nulle;
};

//p et q dans la matrice T telle qu'elle s'écrive
//[T_1  0   0 ]
//[ 0  T_2  0 ] T_1 pxp T_3 qxq
//[ 0   0  T_3]

void p_q(Matrice T,int &p,int &q){
    int n = T.dims[0];
    p = 0;
    q = 0;
    Matrice T_1(T.submat(0,p,0,p));
    Matrice T_3(T.submat(n-1-q,n-1,n-1-q,n-1));
    while (test_diag(T_1)){
        p = p + 1;
    }
    while (test_diag(T_3) && n - 1 - q > p){
        q = q + 1;
    }
};

//Fonction cppidcolgivens prenant quatre flottants en param`etres et modifiant les valeurs des deux derniers (correspondant a c et s dans l’algorithme 1)
void cppidcolgivens(float x, float z, float &c, float &s){
    float tau;
    if (z == 0){
        c = 1;
        s = 0;
    }
    else{
        if(abs(z) > abs(x)){
            tau = -x/z;
            s = 1 / (sqrt(1 + tau * tau));
            c = s * tau;
        }
        else{
            tau = -z/x;
            c = 1 / (sqrt(1 + tau * tau));
            s = c * tau;
        }
    }
};


//Fonction householder prenant en argument un objet de type Vecteur (x) ainsi qu’un flottant passe par reference (beta) et renvoyant un objet de type Vecteur (v) tout en modifiant la valeur du flottant basée sur l'algorithme 2)
Vecteur householder(Vecteur x, float &beta){
    int n = x.dim;
    float x0 = x[0];
    float sigma;
    float mu;

    Vecteur x2n (n - 1);
    Vecteur v (n);

    x2n = x.subvec(1, n-1);
    sigma = dot(x2n,x2n);

    v[0] = 1;
    for (int i = 0; i < n - 1 ; i++){
        v[i+1] = x2n[i];
    }
    if(sigma == 0 && x0 >= 0){
        beta = 0;
    }
    else{
        if(sigma == 0 && x0 < 0){
            beta = 2;
        }
        else{
            mu = sqrt(x[0] * x0 + sigma);
            if(x[0] <= 0){
                v[0] = x[0] - mu;
            }
            else{
                v[0] = - (sigma)/(x[0] + mu);
            }
            beta = (2 * v[0] * v[0])/(sigma * v[0] * v[0]);
            v = v * (1/v[0]);
        }
    }
    return v;
};

//Fonction reductridiag basee sur l’algorithme 3, qui a en parametre un matrice modifiable et en sortie un objet de type Matrice
Matrice reductridiag(Matrice &D){
    int l = D.dims[0];
    Matrice Dsquared (D.transpose() * D);
    float d = (1/2)*(D[l-2][l-2] - D[l-1][l-1]);
    float inter = d + sgn(d) * sqrt(d * d + Dsquared[l-2][l-1]);
    float mu = D[l-1][l-1] - Dsquared[l-2][l-1] / inter;
    float x = D[0][0] - mu;
    float z = D[0][1];
    Matrice Z (Id(l));
    for (int k = 0; k < l-2; k++){
        float c,s;
        cppidcolgivens(x,z,c,s);
        for (int j = 0; j < l-1; j++){
            float tau_1 = D[k][j];
            float tau_2 = D[k+1][j];
            D[k][j] = c * tau_1 - s * tau_2;
            D[k+1][j] = s * tau_1 + c * tau_2;
            tau_1 = Z[k][j];
            tau_2 = Z[k +1][j];
            Z[k][j] = c * tau_1 - s * tau_2;
            Z[k+1][j] = s * tau_1 + c * tau_2;
        }
        for (int j = 0; j < l-1; j++){
            float tau_1 = D[j][k];
            float tau_2 = D[j][k+1];
            D[j][k] = c * tau_1 - s * tau_2;
            D[j][k+1] = s * tau_1 + c * tau_2;
        }
        if (k < l-1){
            x = D[k][k+1];
            z = D[k][k+2];
        }
    }
    return Z;
};

//Factorisation QR symetrique decrite par l’algorithme 4 via une fonction qrsym. Cette fonction prend deux objets de type Matrice modifiables en arguments repr´esentant la matrice a factoriser et le facteur Q a calculer. Elle fait appel a la foction reductdiag.
Matrice qrsym(Matrice &A, Matrice &Q){
    int n = A.dims[0];
    Q = Id(n);

    //Tridiag de la Matrice

    //2
    for (int k = 0; k < n - 2; k++){

        //3
        float beta;
        Matrice Asub (A.submat(k+1,n-1,k,k));
        Vecteur Asubvect (Asub[0]);
        Vecteur v (householder(Asubvect , beta));

        //4
        Matrice AA (A.submat(k+1, n-1, k+1, n-1));
        Vecteur p ((AA * beta).mvprod(v));

        //5
        Vecteur prod = v * (beta / 2) * dot(p,v);
        Vecteur w (p - prod);

        //6
        A[k][k+1] = norm(Asubvect);

        //7
        A[k+1][k] = A[k][k + 1];

        //8
        Matrice out_v_w = outer(v,w);
        Matrice out_w_v = outer(w,v);
        Matrice sum_out = out_v_w + out_w_v;
        AA = AA - sum_out;

        //9
        Matrice QQ (Q.submat(k+1, n-1, k+1, n-1));
        Matrice out_v_v (outer(v,v));
        Matrice vv_QQ (out_v_v * QQ);
        Matrice prod_out(vv_QQ * beta);
        QQ = QQ - prod_out;
    }

    //11
    Matrice T (n,n);

    //12
    for (int j = n-1; j >=0; j--){

        //13
        T[j][j] = A[j][j];

        //14
        if (j > 0){
            T[j][j-1] = A[j-1][j];
            T[j-1][j] = T[j][j-1];
        }
    }

    //Diagonalisation de T et mise à jour de Q

    //16
    while(test_diag(T) == 0){
        for (int i = 0; i < n-1;i++){
            if(abs(T[i+1][i]) + abs(T[i][i+1]) <= 10e-9 * (abs(T[i][i]) - abs(T[i+1][i+1]))){
                T[i+1][i] = 0;
                T[i][i+1] = 0;
            }
        }

        //22
        int p = 0, q = 0;
        p_q(T, p , q);

        //23
        if (p+1 + q+1 < n){

            //24
            Matrice T_2 (T.submat(p+1,n-q-2,p+1,n-q-2));
            Matrice Z (reductridiag(T_2));

            //25
            Matrice T_hat;
            T_hat = T;
            Matrice prod_Q(Id(n));
            for (int i = p+1; i < n-q-1; i++){
                for (int j = p+1; j < n-q-1; j++){
                    T_hat[j][j] = T_2[j - p - 1][i - p - 1];
                    prod_Q[j][i] = Z[j - p - 1][i - p - 1];
                }
            }
            Q = Q * prod_Q;
            Matrice T_hat_T (T_hat.transpose());
            T = (T_hat + T_hat_T) * 0.5;
            Matrice T_hat_plus_T((T_hat + T_hat_T) * 0.5);
        }
    }
    return Q;
};

//Factorisation QR non symetrique avec pivotage (algorithme 5). Cette fonction a deux objets de type Matrice (la matrice a factoriser ainsi que son facteur Q) en argument, qu’elle devra modifier. Elle renvoit egalement une matrice correspondant aux differents pivots de colonnes effectues par l’algorithme
Matrice qrpivot(Matrice A, Matrice &Q){
    int n = A.dims[1];
    int m = A.dims[0];
    Vecteur c (n);

    //1
    Matrice PIV (Id(n));

    //2
    for (int j = 0; n-1; j++){

        //3
        Vecteur A_j (A[j]);
        c[j] = dot(A_j,A_j);
    }

    //5
    int r = -1;
    //6
    float tau = max(c);
    //7
    while (tau > 0 && r < n){
        //8
        r = r + 1;
        //9
        //Defnir k comme le plus petit entier tel que r <= k <= n et ck = tau
        int k = r;
        while(k < n && c[k] != tau){k += 1;}

        //10
        Vecteur perm_col_A (A[k]);
        A[k] = A[r];
        A[r] = perm_col_A;

        //11
        float perm_vect (c[k]);
        c[k] = c[r];
        c[r] = perm_vect;

        //12
        Vecteur perm_col_PIV (PIV[k]);
        PIV[k] = PIV[r];
        PIV[r] = perm_col_PIV;

        //13
        Vecteur A_col_r (A[r]);
        Vecteur A_r_r (A_col_r.subvec(r, m - 1));
        float beta;
        Vecteur v (householder(A_r_r, beta));

        //14
        Matrice A_r_m_r_n (A.submat(r, m-1, r, n-1));
        Matrice vvT (outer(v,v));
        Matrice betavvT (vvT * beta);
        Matrice Result_1(betavvT * A_r_m_r_n);
        Matrice Result_2(A_r_m_r_n - Result_1);
        for (int i = r; i < n; i++){
            for (int j = r; j < m; j++){
                A[j][i] = A[j][i] - Result_2[j-r][i-r];
            }
        }

        //15
        Vecteur subv (v.subvec(1,m-r));
        for (int i = r+1; i < m; i++){
            A[r][i] = subv[i - r - 1];
        }

        //16
        Matrice Asquared (A.transpose() * A);
        for (int i = r+1; i < n; i ++){
            //17
            c[i] = c[i] - Asquared[i][r];
            //18
        }

        //19
        if (r < n-1){
            //20
            tau = max(c.subvec(r+1,n-1));
        }
        //21
        else{
            //22
            tau = 0;
        }
    }

    //Calcul de Q


    //25
    Q = Id(m);

    //26
    Vecteur v_2 (m);

    //27
    for (int j = n-1; j >= 0; j--){

        //28

        Vecteur new_vjm(m-j+1);
        new_vjm = A[j].subvec(j, m-1);
        new_vjm[0] = 1;
        for(int i = j; i < m; i++){
            v_2[i] = new_vjm[i-j];
        }

        //29

        float beta_2 = 1 + norm(A[j].subvec(j,m-1)) * norm(A[j].subvec(j,m-1));
        beta_2 = 1 / beta_2;

        //30

        Matrice Q_j_m (Q.submat(j, m-1, j, m-1));
        Matrice vvT_2 (outer(v_2,v_2));
        Matrice betavvT (vvT_2 * beta_2);
        Matrice Result_1(betavvT * Q_j_m);
        Matrice Result_2(Q_j_m - Result_1);
        for (int i = j; i < m; i++){
            for (int j_2 = j; j_2 < m; j_2++){
                Q[j_2][i] = Q[j_2][i] - Result_2[j_2-j][i-j];
            }
        }
    }

    //31
    return PIV;
};

//Fonction SVD, decomposition en valeurs singulieres. Prend en arguments un objet de type Matrice representant la matrice d’origine, ainsi qu’un tableau d’objets de type Matrice (de taille 3) qui representera les Matrices de la decomposition U, Sigma  et V dans cet ordre.
void svd(Matrice A, Matrice *TableauMatrice){
    int n = A.dims[1];
    int m = A.dims[0];
    Matrice U = TableauMatrice[0];
    Matrice Sigma = TableauMatrice[1];
    Matrice V = TableauMatrice[2];

    if (m >= n){
        Matrice Q_1(n,n);
        Matrice ATA (A.transpose() * A);
        qrsym(ATA, Q_1);
        Matrice AQ_1 (A * Q_1);
        Matrice Q_2(m,m);
        Matrice PIV = qrpivot(AQ_1, Q_2);
        Matrice R = Q_2.transpose() * AQ_1 * PIV;
        for (int j = 0; j < n; j++){
            if (R[j][j] < 0){
                Q_2[j] = - 1 *Q_2[j];
            }
        }
        R = Q_2.transpose() * AQ_1 * PIV;
        U = Q_2;
        Sigma = R;
        V = Q_1 * PIV;
    }
    else{
        Matrice Q_1(m,m);
        Matrice AT = A.transpose();
        Matrice AAT (A * AT);
        qrsym(AAT, Q_1);
        Matrice Q_2(n,n);
        Matrice ATQ_1 (AT * Q_1);
        Matrice PIV = qrpivot(ATQ_1, Q_2);
        Matrice R = Q_2.transpose() * ATQ_1 * PIV;
        for (int i = 0; i < m; i++){
            if (R[i][i] < 0){
                Q_2[i] = - 1 *Q_2[i];
            }
        }
        R = Q_2.transpose() * ATQ_1 * PIV;
        V = Q_2;
        Sigma = R.transpose();
        U = Q_1 * PIV;
    }
};










//FONCTIONS TENSEUR


//Fonction pmod prenant en arguments un objet de classe Tenseur representant un tenseur d’ordre d, un objet de classe Matrice et un entier k entre 1 et d, et effectuant le produit k-modal du tenseur par cette matrice.
Tenseur pmod(Tenseur Tens, Matrice M, int k){
    int m_k = M.get_dims(1);
    int n_k = Tens.d[k-1];

    int ordre_pmod = Tens.ordre;

    int d_pmod[ordre_pmod];
    for (int i = 0; i < ordre_pmod; i++){
        d_pmod[i] = Tens.d[i];
    }
    d_pmod[k-1] = m_k;

    Tenseur pmod(d_pmod, Tens.ordre);

    int nbelts_pmod = 1;
    for (int i = 0; i < ordre_pmod; i++){
        nbelts_pmod = nbelts_pmod * d_pmod[i];
    }

    for (int t = 0; t < nbelts_pmod; t ++){
        Vecteur I_z(ordre_pmod);
        I_z = phi_inv(t+1, ordre_pmod, d_pmod);
        float sum = 0;
        int i_k = I_z[k-1];
        for (int j = 0; j < n_k; j++){
            int I_s[Tens.ordre];
            for (int m = 0; m < ordre_pmod; m++){
                I_s[m] = I_z[m];
            }
            I_s[k] = j;
            sum += M[j][i_k]*Tens[phi(I_s, Tens.ordre,Tens.d) - 1];
        }
        pmod[t] = sum;
    }
    return pmod;
};

//HOSVD

TenseurSVD hosvd(Tenseur T){
    Matrice TablMat[T.ordre];
    TenseurSVD TSVD(T.d, T.ordre, T.vect_tens, TablMat));
    for (int k = 0; k < T.ordre; k++){
        Matrice U, Sigma, V;
        Matrice  TablMatrice_SVD[3];
        TablMatrice_SVD[0] = U;
        TablMatrice_SVD[1] = Sigma;
        TablMatrice_SVD[2] = V;

        svd(T.mode(k+1), TablMatrice_SVD);

        TSVD.TableauMatrice[k] = U;
    }
    return TSVD;
};

int main() {
    //Validation Partie 1 : classe Vecteur
    cout << endl << endl << "#####";
    cout << endl << "#####   Validation Partie 1 : classe Vecteur" << endl;
    cout << "#####" << endl;

    //Q1
    cout << endl << "-----------" << "Part I : Question 1" << endl;
    float tableu[3] = {1,1,1};
    float tablev[4] = {3,4,0,0};
    Vecteur u (tableu,3);
    Vecteur v (tablev,4);
    cout << endl << "u = " << endl;
    u.affiche();
    cout << endl << "v = " << endl;
    v.affiche();

    //Q2
    cout << endl << "-----------" << "Part I : Question 2" << endl;
    Vecteur t(u);
    cout << endl << "t = " << endl;
    t.affiche();
    cout << endl << "u = " << endl;
    u.affiche();

    //Q3
    cout << endl << "-----------" << "Part I : Question 3" << endl;
    u[2] = 0;
    cout << endl << "t = " << endl;
    t.affiche();
    cout << endl << "u = " << endl;
    u.affiche();

    //Q4
    cout << endl << "-----------" << "Part I : Question 4" << endl;
    cout << endl << "dot(v,v) : " << dot(v,v) << endl;
    cout << endl << "norm(v) : " << norm(v) << endl;

    //Q5
    cout << endl << "-----------" << "Part I : Question 5" << endl;
    float unsurnorm;
    unsurnorm = 1/norm(v);
    cout << endl << "(1/|v|)*v = " << endl;
    Vecteur produnsurnorm (v);
    produnsurnorm = v * unsurnorm;
    produnsurnorm.affiche();

    //Q6
    cout << endl << "-----------" << "Part I : Question 6" << endl;
    Vecteur w (3);
    w = v.subvec(1, 3);
    cout << endl << "w = " << endl;
    w.affiche();
    cout << endl << "v = " << endl;
    v.affiche();

    //Q7
    cout << endl << "-----------" << "Part I : Question 7" << endl;
    Vecteur somme (w);
    Vecteur soust (w);
    somme = u+w;
    soust = u-w;
    cout << "u + v" << endl;
    somme.affiche();
    cout << "u - v" << endl;
    soust.affiche();

    //Validation Partie 2 : classe Matrice
    cout << endl << endl << "#####";
    cout << endl << "#####   Validation Partie 2 : classe Matrice" << endl;
    cout << "#####" << endl;


    //Q1
    cout << endl << "-----------" << "Part II : Question 1" << endl;
    float c_1[3] = {1,1,0};
    float c_2[3] = {-0.5,2,-1};
    float c_3[3] = {0,-1,1};
    Vecteur A_c_1(c_1,3);
    Vecteur A_c_2(c_2,3);
    Vecteur A_c_3(c_3,3);
    Matrice A(3,3);
    A[0] = A_c_1;
    A[1] = A_c_2;
    A[2] = A_c_3;
    cout << endl << "A = " << endl;
    A.affiche();

    //Q2
    cout << endl << "-----------" << "Part II : Question 2" << endl;
    float c_4[3] = {-2,0};
    float c_5[3] = {3,1};
    Vecteur B_c_1(c_4,2);
    Vecteur B_c_2(c_5,2);
    Matrice B(2,2);
    B[0] = B_c_1;
    B[1] = B_c_2;
    cout << endl << "B = " << endl;
    B.affiche();

    Matrice C (B);
    B[1][0] = 0;
    cout << endl << "B = " << endl;
    B.affiche();
    cout << endl << "C = " << endl;
    C.affiche();

    //Q3
    cout << endl << "-----------" << "Part II : Question 3" << endl;
    Matrice D (A.submat(0,2,0,1));
    cout << endl << "D = " << endl;
    D.affiche();

    //Q4
    cout << endl << "-----------" << "Part II : Question 4" << endl;
    float diag[3] = {3,2,1};
    Vecteur v_2 (diag,3);
    cout << endl << "v = " << endl;
    v_2.affiche();
    Matrice DI (v_2);
    cout << endl << "DI = " << endl;
    DI.affiche();

    //Q5
    cout << endl << "-----------" << "Part II : Question 5" << endl;
    Matrice SumBC (B + C);
    Matrice SousBC (B -C);
    Matrice ProdDC (D * C);
    cout << endl << "B + C = " << endl;
    SumBC.affiche();
    cout << endl << "B - C = " << endl;
    SousBC.affiche();
    cout << endl << "D * C = " << endl;
    ProdDC.affiche();

    //Question 6
    cout << endl << "-----------" << "Part II : Question 6" << endl;
    cout << endl << "C = " << endl;
    C.affiche();
    Matrice CT (C.transpose());
    cout << endl << "CT = " << endl;
    CT.affiche();
    Matrice C2 (C * CT);
    cout << endl << "C * CT = " << endl;
    C2.affiche();
    float Frob = norm(C);
    cout << endl << "norm(C) = " << endl;
    cout << Frob << endl;

    //TESTS
    cout << endl << "-----------" << "Part II : Autres validations personnelles" << endl;

    cout << endl << "DI * v_2 = " << endl;
    Vecteur Mvprod (DI.mvprod(v_2));
    Mvprod.affiche();

    cout << endl << "DI * 2 = " << endl;
    float k = 2;
    Matrice Mat (DI * k);
    Mat.affiche();

    cout << endl << "outer(v_2,v_2) = " << endl;
    Matrice OUT (outer(v_2,v_2));
    OUT.affiche();

    //Validation Partie 3 : SVD
    cout << endl << endl << "#####";
    cout << endl << "#####   Validation Partie 3 : SVD" << endl;
    cout << "#####" << endl;
    cout << endl << "-----------" << "Part III : Question 1" << endl;
    float c, s;
    float x = 1;
    float z = 2;
    cppidcolgivens(x, z, c, s);
    cout << endl << "c = " << c << endl << "s = " << s << endl;

    cout << endl << "-----------" << "Part III : Question 2" << endl;
    float vect_x[2] = {-1,0};
    float sqareroot = 1/sqrt(2);
    float vect_y[2] = {sqareroot,sqareroot};
    float vect_z[1] = {-4};
    Vecteur X(vect_x,2);
    Vecteur Y(vect_y,2);
    Vecteur Z(vect_z,1);
    cout << endl << "X = " << endl;
    X.affiche();
    cout << endl << "Y = " << endl;
    Y.affiche();
    cout << endl << "Z = " << endl;
    Z.affiche();



    cout << endl << " Avec Householder sur : " << endl;

    float beta_X;
    Vecteur v_hous_X (householder(X , beta_X));
    cout << endl << "X = " << endl;
    v_hous_X.affiche();
    cout << endl << "beta_X = " << endl;
    cout << beta_X << endl;


    cout << endl << "Y = " << endl;
    float beta_Y = 0;
    Vecteur v_hous_Y (householder(Y , beta_Y));
    v_hous_Y.affiche();
    cout << endl << "beta_Y = " << endl;
    cout << beta_Y << endl;


    cout << endl << "Z = " << endl;
    float beta_Z = 0;
    Vecteur v_hous_Z (householder(Z , beta_Z));
    v_hous_Z.affiche();
    cout << endl << "beta_Z = " << endl;
    cout << beta_Z << endl;

    cout << endl << "-----------" << "Part III : Question 3" << endl;

    float c_10[2] = {10,-6};
    float c_11[2] = {-6,10};
    Vecteur M_c_1(c_10,2);
    Vecteur M_c_2(c_11,2);
    Matrice M (2,2);
    M[0] = M_c_1;
    M[1] = M_c_2;
    cout << endl << "M = " << endl;
    M.affiche();

    cout << endl << "-----------" << "Part III : Validation Test Diag" << endl;
    bool bool_test = test_diag(A);
    cout << endl << bool_test << endl;
    bool_test = test_diag(DI);
    cout << endl << bool_test << endl;

    cout << endl << "-----------" << "Part III : Question 4" << endl;

    float colo_1[3] = {1,0.5,0.33};
    float colo_2[3] = {0.5,0.33,0.25};
    float colo_3[3] = {0.33,0.25,0.20};

    Vecteur vect_1 (colo_1,3);
    Vecteur vect_2 (colo_2,3);
    Vecteur vect_3 (colo_3,3);

    Matrice GrandM(3,3);
    GrandM[0] = vect_1;
    GrandM[1] = vect_2;
    GrandM[2] = vect_3;

    Matrice Q;

    Matrice QRPIV_M (qrpivot(GrandM, Q));

    cout << endl << "-----------" << "Part III : Question 5" << endl;

    Matrice U, Sigma, V;
    Matrice  TablMatrice[3];
    TablMatrice[0] = U;
    TablMatrice[1] = Sigma;
    TablMatrice[2] = V;

    float colo_1_A[2] = {1,0};
    float colo_2_A[2] = {0,-1};

    Vecteur vect_1_A (colo_1_A,2);
    Vecteur vect_2_A (colo_2_A,2);

    Matrice AA(2,2);
    AA[0] = vect_1_A;
    AA[1] = vect_2_A;

    float colo_1_B[2] = {float(2*sqrt(2)),float(-1*sqrt(2))};
    float colo_2_B[2] = {float(-2*sqrt(2)),float(-1*sqrt(2))};

    Vecteur vect_1_B (colo_1_B,2);
    Vecteur vect_2_B (colo_2_B,2);

    Matrice BB(2,2);
    BB[0] = vect_1_A;
    BB[1] = vect_2_A;

    float colo_1_C[2] = {0.5,float(sqrt(3)/2)};
    float colo_2_C[2] = {float(3*sqrt(3)/2),float(-3/2)};
    float colo_3_C[2] = {0,0};


    Vecteur vect_1_C (colo_1_C,2);
    Vecteur vect_2_C (colo_3_C,2);
    Vecteur vect_3_C (colo_3_C,2);


    Matrice CC(2,3);
    CC[0] = vect_1_C;
    CC[1] = vect_2_C;
    CC[2] = vect_3_C;

    svd(AA, TablMatrice);

    TablMatrice[0].affiche();
    TablMatrice[1].affiche();
    TablMatrice[2].affiche();

    svd(BB, TablMatrice);


    TablMatrice[0].affiche();
    TablMatrice[1].affiche();
    TablMatrice[2].affiche();

    svd(CC, TablMatrice);

    TablMatrice[0].affiche();
    TablMatrice[1].affiche();
    TablMatrice[2].affiche();


    cout << endl << endl << "#####";
    cout << endl << "#####   Validation Partie 4 : Tenseur" << endl;
    cout << "#####" << endl;
    //Q1
    cout << endl << "-----------" << "Part IV : Question 1" << endl;

    int d_T[3];
    d_T[0] = 2;
    d_T[1] = 2;
    d_T[2] = 2;
    Tenseur calT(d_T,3);
    cout << endl << "calT = " << endl;
    calT.affiche_vect();

    //Q2
    cout << endl << "-----------" << "Part IV : Question 2" << endl;

    int nbelt = d_T[0] * d_T[1] * d_T[2] ;
    Tenseur calU(d_T,3);
    for (int i=0; i < nbelt; i++){
        calU[i] = 1;
    }
    cout << endl << "calU = " << endl;
    calU.affiche_vect();

    //Q3
    cout << endl << "-----------" << "Part IV : Question 3" << endl;
    Tenseur calV(calT + calU);
    Tenseur calW(calU - calT);

    cout << endl << "calV = calU + calT = " << endl;
    calV.affiche_vect();
    cout << endl << "calW = calU - calT = " << endl;
    calW.affiche_vect();

    //Q4
    cout << endl << "-----------" << "Part IV : Question 4" << endl;

    int I_z[3] = {2,2,2};

    calU[phi(I_z,3,d_T)-1] = -1;

    cout << endl << "calU = " << endl;
    calU.affiche_vect();

    cout << endl << "calV = " << endl;
    calV.affiche_vect();

    //Q5
    cout << endl << "-----------" << "Part IV : Question 5" << endl;

    cout << endl << "calT^(k) = " << endl;

    calT.mode(1).affiche();

    //Q6
    cout << endl << "-----------" << "Part IV : Question 6" << endl;

    float col_1[2] = {1,4};
    float col_2[2] = {3,0.3333};
    float col_3[2] = {0,1.5};
    float col_4[2] = {-1,2};

    Vecteur vec_1 (col_1,2);
    Vecteur vec_2 (col_2,2);
    Vecteur vec_3 (col_3,2);
    Vecteur vec_4 (col_4,2);

    Matrice Mm(2,4);
    Mm[0] = vec_1;
    Mm[1] = vec_2;
    Mm[2] = vec_3;
    Mm[3] = vec_4;

    cout << endl << "Matrice redifinissant calT = " << endl;

    Mm.affiche();

    Tenseur calT_2 (d_T, 3, 2, Mm);
    calT = calT_2;

    cout << endl << "calT = " << endl;

    calT.affiche_vect();

    //Q7
    cout << endl << "-----------" << "Part IV : Question 7" << endl;

    float col_5[3] = {3,0,0};
    float col_6[3] = {-1,6,-3};

    Vecteur vec_5 (col_5,3);
    Vecteur vec_6 (col_6,3);

    Matrice calA(3,2);
    calA[0] = vec_5;
    calA[1] = vec_6;

    Tenseur S (pmod(calT, calA, 3));

    S.affiche_vect();


    cout << endl << "-----------" << "Brouillon" << endl;

    v_2.affiche();

    int a = v_2[0] + 0.9 ;

    cout << a << endl;

    int du[3] = {3,4,2};
    int rp = prod_n_i(0, 3, du);
    cout << rp << endl;

    int I[3] = {3,4,2};

    int yo = phi(I,3,du);

    cout << "yo = ";
    cout << yo << endl;

    Vecteur pos_yo (phi_inv(yo, 3, du));
    cout << endl << "postition de yo : ";
    pos_yo.affiche();

    int I_2[3];
    for (int i = 0; i < 3; i++){
        I_2[i] = pos_yo[i];
    }

    int yo_2 = phi(I_2,3,du);

    cout << "yo_2 = " << yo_2 << endl;


    Vecteur T(24);
    for (int i = 0; i < 24; i++){
        T[i] = i+1;
    }

    cout << "ordre = " << du[0] * du[1] * du[2] << endl;

    T.affiche();

    Tenseur TENS(du, 3, T);

    TENS.affiche_vect();

    Matrice TENS_1 = TENS.mode(2);

    TENS_1.affiche();

    Tenseur TENS_via_Mat(du, 3, 2, TENS_1);

    TENS_via_Mat.affiche_vect();

    return 0;
}
