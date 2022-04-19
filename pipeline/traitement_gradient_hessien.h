#ifndef GRADIENT_H
#define GRADIENT_H

#include "pipeline_first.h"

/*ON DECIDE DE CE QUE L'ON VA ENVOYER AU TRAITEMENT SUIVANT*/
//Carte de gradient : (utile pour le hessien ensuite)
// - HAUT GAUCHE :  (h_haut*h_haut + w_gauche*w_gauche)
// - HAUT DROITE :  (h_haut*h_haut + w_droite*w_droite)
// - BAS GAUCHE :   (h_bas*h_bas + w_gauche*w_gauche)
// - BAS DROITE :   (h_bas*h_bas + w_droite*w_droite)
#define FLG_GRADIENT_RES_CARTE 3

//4 différents gradients : 
// - HAUT GAUCHE :  Img_gradient_vertical_gauche
// - HAUT DROITE :  Img_gradient_horizontal_haut
// - BAS GAUCHE :   Img_gradient_horizontal_bas
// - BAS DROITE :   Img_gradient_vertical_droite
#define FLG_GRADIENT_RES_GRADIENT 2
//Les deux images
#define FLG_GRADIENT_RES_BOTH 1


/*ON DECIDE DU HESSIEN QUE L'ON VA CALCULER*/
#define FLG_HESSIEN_GAUCHE_BAS 8
#define FLG_HESSIEN_GAUCHE_HAUT 7
#define FLG_HESSIEN_DROITE_BAS 6
#define FLG_HESSIEN_DROITE_HAUT 5
//Tous les hessiens gauches (2)
#define FLG_HESSIEN_GAUCHE 4
//Tous les hessiens droits (2)
#define FLG_HESSIEN_DROITE 3
//Tous les hessiens bas (2)
#define FLG_HESSIEN_BAS 2
//Tous les hessiens hauts (2)
#define FLG_HESSIEN_HAUT 1
//Tous les hessiens (4)
#define FLG_HESSIEN_TOUS 0


//Dire de qu'elle dérivation on prend pour harris
#define FLG_HARRIS_GAUCHE_BAS 8
#define FLG_HARRIS_GAUCHE_HAUT 7
#define FLG_HARRIS_DROITE_BAS 6
#define FLG_HARRIS_DROITE_HAUT 5


/**   ==================== ==================== ==================== ====================         
 *              ==================== ==================== ====================
 *    ==================== ==================== ==================== ====================**/
/**   ====================         
 *     TRAITEMENT GRADIENT        
 *    ====================**/
class Traitement_gradient : public Traitement_image{
    private :
    int flg_grd;
    void gradient(int nH, int nW, OCTET *ImgIn,OCTET *Img_carte_gradient,OCTET *Img_gradient_horizontal_bas,OCTET *Img_gradient_horizontal_haut,OCTET *Img_gradient_vertical_droite,OCTET *Img_gradient_vertical_gauche,OCTET *Img_gradient);
    
    public :
        Traitement_gradient(int flg = FLG_GRADIENT_RES_BOTH);
        virtual ~Traitement_gradient();
        void appliquer();
        virtual void afficher(std::ostream& os);
};

void Traitement_gradient::gradient(int nH, int nW, OCTET *ImgIn,OCTET *Img_carte_gradient,OCTET *Img_gradient_horizontal_bas,OCTET *Img_gradient_horizontal_haut,OCTET *Img_gradient_vertical_droite,OCTET *Img_gradient_vertical_gauche,OCTET *Img_gradient){
    float valeur_max_hg = 0.0f;
    float valeur_max_hd = 0.0f;
    float valeur_max_bg = 0.0f;
    float valeur_max_bd = 0.0f;
   
    for (int i = 1; i < nH-1; i++){
        for (int j = 1; j < nW-1; j++){
            float h_bas = ImgIn[(i+1)*nW+j] - ImgIn[i*nW+j];
            float w_droite = ImgIn[i*nW+j+1] - ImgIn[i*nW+j];
            float h_haut = ImgIn[(i-1)*nW+j] - ImgIn[i*nW+j];
            float w_gauche = ImgIn[i*nW+j-1] - ImgIn[i*nW+j];
            
            Img_gradient_horizontal_bas[i*nW+j] = h_bas;
            Img_gradient_horizontal_haut[i*nW+j] = h_haut;
            Img_gradient_vertical_gauche[i*nW+j] = w_gauche;
            Img_gradient_vertical_droite[i*nW+j] = w_droite;

            int largeur_grad = nW*2;
            Img_gradient[i*largeur_grad+j] = w_gauche;
            Img_gradient[i*largeur_grad+j+nW] = h_haut;
            Img_gradient[(i+nH)*largeur_grad+j] = h_bas;
            Img_gradient[(i+nH)*largeur_grad+j+nW] = w_droite;
            
            float val_hg = sqrt(h_haut*h_haut + w_gauche*w_gauche);
            float val_hd = sqrt(h_haut*h_haut + w_droite*w_droite);
            float val_bg = sqrt(h_bas*h_bas + w_gauche*w_gauche);
            float val_bd = sqrt(h_bas*h_bas + w_droite*w_droite);


            Img_carte_gradient[i*largeur_grad+j] = (val_hg<0)? 0 : (val_hg>255)? 255 : val_hg;
            Img_carte_gradient[i*largeur_grad+j+nW] = (val_hd<0)? 0 : (val_hd>255)? 255 : val_hd;
            Img_carte_gradient[(i+nH)*largeur_grad+j] = (val_bg<0)? 0 : (val_bg>255)? 255 : val_bg;
            Img_carte_gradient[(i+nH)*largeur_grad+j+nW] = (val_bd<0)? 0 : (val_bd>255)? 255 : val_bd;

            if(val_hg>valeur_max_hg){
                valeur_max_hg = val_hg;
            }
            if(val_hd>valeur_max_hd){
                valeur_max_hd = val_hd;
            }
            if(val_bg>valeur_max_bg){
                valeur_max_bg = val_bg;
            }
            if(val_bd>valeur_max_bd){
                valeur_max_bd = val_bd;
            }

        }
    }

    if(AFFICHAGE_DEBUG){
        std::cout << "================================================" << std::endl;
        std::cout << "VALEUR MAXIMALE HAUT GAUCHE : " << valeur_max_hg << std::endl;
        std::cout << "VALEUR MAXIMALE HAUT DROITE : " << valeur_max_hd << std::endl;
        std::cout << "VALEUR MAXIMALE BAS GAUCHE : " << valeur_max_bg << std::endl;
        std::cout << "VALEUR MAXIMALE BAS DROITE : " << valeur_max_bd << std::endl;
    }

}
Traitement_gradient::Traitement_gradient(int flg) : Traitement_image(){
    this->flg_grd = flg;
}
Traitement_gradient::~Traitement_gradient(){}
void Traitement_gradient::appliquer(){
    bool write_img = getTraitementEcriture();

    int flg_carte_suite = (flg_grd == FLG_GRADIENT_RES_BOTH || flg_grd == FLG_GRADIENT_RES_CARTE)? FLG_IMG_SUITE : FLG_IMG_NO_SUITE;
    int flg_gradient_suite = (flg_grd == FLG_GRADIENT_RES_BOTH || flg_grd == FLG_GRADIENT_RES_GRADIENT)? FLG_IMG_SUITE : FLG_IMG_NO_SUITE;

    for (size_t i = 0; i < img_in.size(); i++){
        IMG* img_cur = img_in[i];
        if(AFFICHAGE_DEBUG){
            std::cout << "APPLIQUER[" << i << "] = " << *img_cur << std::endl;
        }

        OCTET *Img_carte_gradient;
        OCTET *Img_gradient_horizontal_bas;
        OCTET *Img_gradient_horizontal_haut;
        OCTET *Img_gradient_vertical_droite;
        OCTET *Img_gradient_vertical_gauche;
        OCTET *Img_gradient;

        int nTaille = img_cur->nH*img_cur->nW;
        allocation_tableau(Img_carte_gradient, OCTET, nTaille*4);
        allocation_tableau(Img_gradient_horizontal_bas, OCTET, nTaille);
        allocation_tableau(Img_gradient_horizontal_haut, OCTET, nTaille);
        allocation_tableau(Img_gradient_vertical_droite, OCTET, nTaille);
        allocation_tableau(Img_gradient_vertical_gauche, OCTET, nTaille);
        allocation_tableau(Img_gradient, OCTET, nTaille*4);

        this->gradient(img_cur->nH,img_cur->nW,img_cur->img_data,Img_carte_gradient,Img_gradient_horizontal_bas,Img_gradient_horizontal_haut,Img_gradient_vertical_droite,Img_gradient_vertical_gauche,Img_gradient);
        
        char* nom_carte_gradient = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_carte_gradient,"%s_carte_gradient", img_cur->nom);    
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH*2, img_cur->nW*2, Img_carte_gradient, write_img,nom_carte_gradient,flg_carte_suite));

        char* nom_gradientHB = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_gradientHB,"%s_gradient_hori_bas", img_cur->nom);    
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_gradient_horizontal_bas, write_img,nom_gradientHB,FLG_IMG_NO_SUITE));

        char* nom_gradientHH = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_gradientHH,"%s_gradient_hori_haut", img_cur->nom);    
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_gradient_horizontal_haut, write_img,nom_gradientHH,FLG_IMG_NO_SUITE));

        char* nom_gradientVD = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_gradientVD,"%s_gradient_verti_droite", img_cur->nom);    
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_gradient_vertical_droite, write_img,nom_gradientVD,FLG_IMG_NO_SUITE));

        char* nom_gradientVG = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_gradientVG,"%s_gradient_verti_gauche", img_cur->nom);    
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_gradient_vertical_gauche, write_img,nom_gradientVG,FLG_IMG_NO_SUITE));

        char* nom_gradient = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_gradient,"%s_gradient_mix", img_cur->nom);    
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH*2, img_cur->nW*2, Img_gradient, write_img,nom_gradient,flg_gradient_suite));
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_gradient::afficher(std::ostream& os){
    os << "TRAITEMENT GRADIENT" << std::endl;
    Traitement_image::afficher(os);
}


/**   ====================         
 *     TRAITEMENT HESSIEN        
 *    ====================**/
class Traitement_hessien : public Traitement_image{
    private :
    int flg_hes;
    void hessien(int nH, int nW, OCTET *Img_deriv_x,OCTET *Img_deriv_y, OCTET *Img_hessien, int vx,int vy);
    void extraire_gradients(int nH, int nW,OCTET *Img_gradient_horizontal_bas,OCTET *Img_gradient_horizontal_haut,OCTET *Img_gradient_vertical_droite,OCTET *Img_gradient_vertical_gauche,OCTET *Img_gradient);

    public :
        Traitement_hessien(int flg = FLG_HESSIEN_TOUS);
        virtual ~Traitement_hessien();
        void appliquer();
        virtual void afficher(std::ostream& os);
};

void Traitement_hessien::hessien(int nH, int nW, OCTET *Img_deriv_x,OCTET *Img_deriv_y, OCTET *Img_hessien, int vx,int vy){
    for (int i = 1; i < nH-1; i++){
        for (int j = 1; j < nW-1; j++){
            float xy = Img_deriv_x[(i+vy)*nW+j] - Img_deriv_x[i*nW+j];
            float yx = Img_deriv_y[i*nW+j+vx] - Img_deriv_y[i*nW+j];
            float xx = Img_deriv_x[i*nW+j+vx] - Img_deriv_x[i*nW+j];
            float yy = Img_deriv_y[(i+vy)*nW+j] - Img_deriv_x[i*nW+j];

            Img_hessien[i*nW*2+j] = xx;
            Img_hessien[i*nW*2+j+nW] = xy;
            Img_hessien[(i+nH)*(nW*2)+j] = yx;
            Img_hessien[(i+nH)*(nW*2)+j+nW] = yy;
        }
    }
}
void Traitement_hessien::extraire_gradients(int nH, int nW,OCTET *Img_gradient_horizontal_bas,OCTET *Img_gradient_horizontal_haut,OCTET *Img_gradient_vertical_droite,OCTET *Img_gradient_vertical_gauche,OCTET *Img_gradient){
    int largeur_grad = nW/2;
    int hauteur_grad = nH/2;

    for (int i = 0; i < hauteur_grad; i++){
        for (int j = 0; j < largeur_grad; j++){
            Img_gradient_horizontal_bas[i*largeur_grad+j] = Img_gradient[(i+hauteur_grad)*nW+j];
            Img_gradient_horizontal_haut[i*largeur_grad+j] = Img_gradient[i*nW+j+largeur_grad];
            Img_gradient_vertical_gauche[i*largeur_grad+j] = Img_gradient[i*nW+j];
            Img_gradient_vertical_droite[i*largeur_grad+j] = Img_gradient[(i+hauteur_grad)*nW+j+largeur_grad];

        }
    }   
}
Traitement_hessien::Traitement_hessien(int flg): Traitement_image(){
    this->flg_hes = flg;
}
Traitement_hessien::~Traitement_hessien(){}
void Traitement_hessien::appliquer(){
    bool write_img = getTraitementEcriture();
    if(AFFICHAGE_DEBUG){
        std::cout << "HESSIEN APPLIQUER" << std::endl;
        std::cout << "WRITE IMG : " << write_img << std::endl;
    }

    for (size_t i = 0; i < this->img_in.size(); i++){
        if(AFFICHAGE_DEBUG){
            std::cout << "IMG IN " << std::endl;
            std::cout << *img_in[i];
        }


        IMG * img_cur = this->img_in[i];
        
        OCTET *Img_gradient_horizontal_bas;
        OCTET *Img_gradient_horizontal_haut;
        OCTET *Img_gradient_vertical_droite;
        OCTET *Img_gradient_vertical_gauche;
        int nTaille = img_cur->nH*img_cur->nW;
        allocation_tableau(Img_gradient_horizontal_bas, OCTET, nTaille);
        allocation_tableau(Img_gradient_horizontal_haut, OCTET, nTaille);
        allocation_tableau(Img_gradient_vertical_droite, OCTET, nTaille);
        allocation_tableau(Img_gradient_vertical_gauche, OCTET, nTaille);

        extraire_gradients(img_cur->nH, img_cur->nW,Img_gradient_horizontal_bas,Img_gradient_horizontal_haut,Img_gradient_vertical_droite,Img_gradient_vertical_gauche,img_cur->img_data);

        OCTET *Img_hessien_gauche_bas;
        OCTET *Img_hessien_gauche_haut;
        OCTET *Img_hessien_droite_bas;
        OCTET *Img_hessien_droite_haut;

        char* nom_hessien_gb = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_hessien_gb,"%s_hessien_gb", img_cur->nom); 
        char* nom_hessien_gh = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_hessien_gh,"%s_hessien_gh", img_cur->nom);    
        char* nom_hessien_db = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_hessien_db,"%s_hessien_db", img_cur->nom);  
        char* nom_hessien_dh = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_hessien_dh,"%s_hessien_dh", img_cur->nom);

        switch (this->flg_hes){
            case FLG_HESSIEN_TOUS:
                if(AFFICHAGE_DEBUG) std::cout << "HESSIEN TOUS" << std::endl;

                allocation_tableau(Img_hessien_gauche_bas, OCTET, nTaille);
                allocation_tableau(Img_hessien_gauche_haut, OCTET, nTaille);
                allocation_tableau(Img_hessien_droite_bas, OCTET, nTaille);
                allocation_tableau(Img_hessien_droite_haut, OCTET, nTaille);

                this->hessien(img_cur->nH/2, img_cur->nW/2, Img_gradient_vertical_gauche,Img_gradient_horizontal_bas, Img_hessien_gauche_bas, 1,1);
                this->hessien(img_cur->nH/2, img_cur->nW/2, Img_gradient_vertical_gauche,Img_gradient_horizontal_haut, Img_hessien_gauche_haut, 1,1);
                this->hessien(img_cur->nH/2, img_cur->nW/2, Img_gradient_vertical_droite,Img_gradient_horizontal_bas, Img_hessien_droite_bas, 1,1);
                this->hessien(img_cur->nH/2, img_cur->nW/2, Img_gradient_vertical_droite,Img_gradient_horizontal_haut, Img_hessien_droite_haut, 1,1);
      
                this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_hessien_gauche_bas, write_img,nom_hessien_gb));
                this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_hessien_gauche_haut, write_img,nom_hessien_gh));
                this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_hessien_droite_bas, write_img,nom_hessien_db));
                this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_hessien_droite_haut, write_img,nom_hessien_dh));
                break;    
            case FLG_HESSIEN_GAUCHE:
                if(AFFICHAGE_DEBUG) std::cout << "HESSIEN GAUCHE" << std::endl;

                allocation_tableau(Img_hessien_gauche_bas, OCTET, nTaille);
                allocation_tableau(Img_hessien_gauche_haut, OCTET, nTaille);
                this->hessien(img_cur->nH/2, img_cur->nW/2, Img_gradient_vertical_gauche,Img_gradient_horizontal_bas, Img_hessien_gauche_bas, 1,1);
                this->hessien(img_cur->nH/2, img_cur->nW/2, Img_gradient_vertical_gauche,Img_gradient_horizontal_haut, Img_hessien_gauche_haut, 1,1);
      
                this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_hessien_gauche_bas, write_img,nom_hessien_gb));
                this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_hessien_gauche_haut, write_img,nom_hessien_gh)); 
                break;
            case FLG_HESSIEN_GAUCHE_BAS:
                if(AFFICHAGE_DEBUG) std::cout << "HESSIEN GAUCHE BAS" << std::endl;

                allocation_tableau(Img_hessien_gauche_bas, OCTET, nTaille);
                this->hessien(img_cur->nH/2, img_cur->nW/2, Img_gradient_vertical_gauche,Img_gradient_horizontal_bas, Img_hessien_gauche_bas, 1,1);  
                this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_hessien_gauche_bas, write_img,nom_hessien_gb));
                break;
            case FLG_HESSIEN_GAUCHE_HAUT:
                if(AFFICHAGE_DEBUG) std::cout << "HESSIEN GAUCHE HAUT" << std::endl;

                allocation_tableau(Img_hessien_gauche_haut, OCTET, nTaille);
                this->hessien(img_cur->nH/2, img_cur->nW/2, Img_gradient_vertical_gauche,Img_gradient_horizontal_haut, Img_hessien_gauche_haut, 1,1);
                this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_hessien_gauche_haut, write_img,nom_hessien_gh)); 
                break;
            case FLG_HESSIEN_DROITE:
                if(AFFICHAGE_DEBUG) std::cout << "HESSIEN DROITE" << std::endl;

                allocation_tableau(Img_hessien_droite_bas, OCTET, nTaille);
                allocation_tableau(Img_hessien_droite_haut, OCTET, nTaille);
                this->hessien(img_cur->nH/2, img_cur->nW/2, Img_gradient_vertical_droite,Img_gradient_horizontal_bas, Img_hessien_droite_bas, 1,1);
                this->hessien(img_cur->nH/2, img_cur->nW/2, Img_gradient_vertical_droite,Img_gradient_horizontal_haut, Img_hessien_droite_haut, 1,1);
      
                this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_hessien_droite_bas, write_img,nom_hessien_db));
                this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_hessien_droite_haut, write_img,nom_hessien_dh));
                break;
            case FLG_HESSIEN_DROITE_HAUT:
                if(AFFICHAGE_DEBUG) std::cout << "HESSIEN DROITE HAUT" << std::endl;

                allocation_tableau(Img_hessien_droite_haut, OCTET, nTaille);
                this->hessien(img_cur->nH/2, img_cur->nW/2, Img_gradient_vertical_droite,Img_gradient_horizontal_haut, Img_hessien_droite_haut, 1,1);
                this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_hessien_droite_haut, write_img,nom_hessien_dh));
                break;
            case FLG_HESSIEN_DROITE_BAS:
                if(AFFICHAGE_DEBUG) std::cout << "HESSIEN DROITE BAS" << std::endl;

                allocation_tableau(Img_hessien_droite_bas, OCTET, nTaille);
                this->hessien(img_cur->nH/2, img_cur->nW/2, Img_gradient_vertical_droite,Img_gradient_horizontal_bas, Img_hessien_droite_bas, 1,1);
                this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_hessien_droite_bas, write_img,nom_hessien_db));
                break;
            case FLG_HESSIEN_BAS:
                if(AFFICHAGE_DEBUG) std::cout << "HESSIEN BAS" << std::endl;

                allocation_tableau(Img_hessien_gauche_bas, OCTET, nTaille);
                allocation_tableau(Img_hessien_droite_bas, OCTET, nTaille);
                
                this->hessien(img_cur->nH/2, img_cur->nW/2, Img_gradient_vertical_gauche,Img_gradient_horizontal_bas, Img_hessien_gauche_bas, 1,1);
                this->hessien(img_cur->nH/2, img_cur->nW/2, Img_gradient_vertical_droite,Img_gradient_horizontal_bas, Img_hessien_droite_bas, 1,1);
                this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_hessien_gauche_bas, write_img,nom_hessien_gb));
                this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_hessien_droite_bas, write_img,nom_hessien_db));
                break;  
            case FLG_HESSIEN_HAUT:
                if(AFFICHAGE_DEBUG) std::cout << "HESSIEN HAUT" << std::endl;

                allocation_tableau(Img_hessien_gauche_haut, OCTET, nTaille);
                allocation_tableau(Img_hessien_droite_haut, OCTET, nTaille);

                this->hessien(img_cur->nH/2, img_cur->nW/2, Img_gradient_vertical_gauche,Img_gradient_horizontal_haut, Img_hessien_gauche_haut, 1,1);
                this->hessien(img_cur->nH/2, img_cur->nW/2, Img_gradient_vertical_droite,Img_gradient_horizontal_haut, Img_hessien_droite_haut, 1,1);
                this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_hessien_gauche_haut, write_img,nom_hessien_gh));
                this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_hessien_droite_haut, write_img,nom_hessien_dh));
                break;
        }
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_hessien::afficher(std::ostream& os){
    os << "TRAITEMENT HESSIEN" << std::endl;
    Traitement_image::afficher(os);
}


/**   ====================         
 *    TRAITEMENT ROSSENFELD        
 *    ====================**/
class Traitement_rossenfeld : public Traitement_image{
    private :
    void rossenfeld(int nH, int nW, OCTET *Img_deriv_x,OCTET *Img_deriv_y, OCTET *Img_deriv_XX,OCTET *Img_deriv_YY,OCTET *Img_deriv_XY,OCTET *Img_deriv_YX, OCTET *Img_out);

    bool effectue = false;
    public :
        Traitement_rossenfeld();
        virtual ~Traitement_rossenfeld();
        void appliquer();
        virtual void afficher(std::ostream& os);
};

void Traitement_rossenfeld::rossenfeld(int nH, int nW, OCTET *Img_deriv_x,OCTET *Img_deriv_y, OCTET *Img_deriv_XX,OCTET *Img_deriv_YY,OCTET *Img_deriv_XY,OCTET *Img_deriv_YX, OCTET *Img_out){
    if(AFFICHAGE_DEBUG) {
        std::cout << "ROSSENFELD TRAITEMENT" << std::endl;
        std::cout << "  ->NH : " << nH << std::endl;
        std::cout << "  ->NW : " << nW << std::endl;
    }

    for (int i = 0; i < nH; i++){
        for (int j = 0; j < nW; j++){
            float deno = (Img_deriv_x[i*nW+j] * Img_deriv_x[i*nW+j] * Img_deriv_y[i*nW+j] * Img_deriv_y[i*nW+j]);
            Img_out[i*nW+j] = (float)((Img_deriv_XX[i*nW+j]*Img_deriv_y[i*nW+j]*Img_deriv_y[i*nW+j]) + (Img_deriv_YY[i*nW+j] * Img_deriv_x[i*nW+j]*Img_deriv_x[i*nW+j]) - (2 * Img_deriv_YX[i*nW+j] * Img_deriv_x[i*nW+j] * Img_deriv_y[i*nW+j])) / deno;
        }
    }
}
Traitement_rossenfeld::Traitement_rossenfeld(): Traitement_image(){}
Traitement_rossenfeld::~Traitement_rossenfeld(){}
void Traitement_rossenfeld::appliquer(){
    bool write_img = getTraitementEcriture();
    if(AFFICHAGE_DEBUG){
        std::cout << "ROSSENFLED APPLIQUER" << std::endl;
        std::cout << "WRITE IMG : " << write_img << std::endl;
    }

    if(this->img_in.size() == 2 && !effectue){
        this->effectue = true;
        IMG * img_gradient = this->img_in[0];
        IMG * img_hessien = this->img_in[1];


        OCTET * Img_deriv_x, *Img_deriv_y,*Img_deriv_XX,*Img_deriv_YY,*Img_deriv_XY,*Img_deriv_YX,*Img_out;
        
        int nTaille = img_gradient->nH/2 * img_gradient->nW/2;

        allocation_tableau(Img_deriv_x, OCTET, nTaille);
        allocation_tableau(Img_deriv_y, OCTET, nTaille);

        allocation_tableau(Img_deriv_XX, OCTET, nTaille);
        allocation_tableau(Img_deriv_YY, OCTET, nTaille);
        allocation_tableau(Img_deriv_XY, OCTET, nTaille);
        allocation_tableau(Img_deriv_YX, OCTET, nTaille);
        allocation_tableau(Img_out, OCTET, nTaille);

        int m_xw = 0;
        int m_xh = 0;

        int m_yh = 0;
        int m_yw = 0;

        for (int i = (m_xh *img_gradient->nH/2); i < ((m_xh+1) *img_gradient->nH/2); i++){
            for (int j = (m_xw*img_gradient->nW/2); j < ((m_xw +1)*img_gradient->nW/2); j++){
                Img_deriv_x[i*img_gradient->nW/2+j] = img_gradient->img_data[i*img_gradient->nW+j];
            }   
        }
        for (int i = (m_yh *img_gradient->nH/2); i < ((m_yh+1) *img_gradient->nH/2); i++){
            for (int j = (m_yw*img_gradient->nW/2); j < ((m_yw +1)*img_gradient->nW/2); j++){
                Img_deriv_y[i*img_gradient->nW/2+j] = img_gradient->img_data[i*img_gradient->nW+j];
            }   
        }

        for (int i = 0; i < img_hessien->nH/2; i++){
            for (int j = 0; j < img_hessien->nW/2; j++){
                Img_deriv_XX[i*img_hessien->nW/2+j] = img_hessien->img_data[i*img_hessien->nW + j];
                Img_deriv_XY[i*img_hessien->nW/2+j] = img_hessien->img_data[i*img_hessien->nW + j + img_hessien->nW/2];
                Img_deriv_YX[i*img_hessien->nW/2+j] = img_hessien->img_data[(i+img_hessien->nH/2)*img_hessien->nW + j];
                Img_deriv_YY[i*img_hessien->nW/2+j] = img_hessien->img_data[(i+img_hessien->nH/2)*img_hessien->nW + j + img_hessien->nW/2];
            }   
        }

        rossenfeld(img_gradient->nH/2, img_gradient->nW/2,Img_deriv_x,Img_deriv_y,Img_deriv_XX,Img_deriv_YY,Img_deriv_XY,Img_deriv_YX,Img_out);
        
        char* nom_rossenfeld = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_rossenfeld,"%s+%s_rossenfeld", img_gradient->nom,img_hessien->nom);  
        this->img_out.push_back(new IMG(this->format_res,img_gradient->nH/2, img_gradient->nW/2, Img_out, write_img,nom_rossenfeld));

    }else{
        this->img_out.clear();
        if(AFFICHAGE_DEBUG){
            std::cout << "ROSSENFLED IMPOSSIBLE D'APPLIQUER : Gradient + Hessien : " << std::endl;
            std::cout << "      -> TAILLE : " << this->img_in.size() << std::endl;
            std::cout << "      -> EFFECTUE : " << this->effectue << std::endl;
        }
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_rossenfeld::afficher(std::ostream& os){
    os << "TRAITEMENT ROSSENFELD" << std::endl;
    Traitement_image::afficher(os);
}


/**   ====================         
 *     TRAITEMENT BEAUDET        
 *    ====================**/
class Traitement_beaudet : public Traitement_image{
    private :
    void beaudet(int nH, int nW, OCTET *Img_deriv_XX,OCTET *Img_deriv_YY,OCTET *Img_deriv_XY,OCTET *Img_deriv_YX, OCTET *Img_out);

    public :
        Traitement_beaudet();
        virtual ~Traitement_beaudet();
        void appliquer();
        virtual void afficher(std::ostream& os);
};

void Traitement_beaudet::beaudet(int nH, int nW, OCTET *Img_deriv_XX,OCTET *Img_deriv_YY,OCTET *Img_deriv_XY,OCTET *Img_deriv_YX, OCTET *Img_out){
    if(AFFICHAGE_DEBUG) {
        std::cout << "BEAUDET TRAITEMENT" << std::endl;
        std::cout << "  ->NH : " << nH << std::endl;
        std::cout << "  ->NW : " << nW << std::endl;
    }

    for (int i = 0; i < nH; i++){
        for (int j = 0; j < nW; j++){
            Img_out[i*nW+j] = (Img_deriv_XX[i*nW+j]*Img_deriv_YY[i*nW+j]) - (Img_deriv_XY[i*nW+j] * Img_deriv_YX[i*nW+j]);
        }
    }
}
Traitement_beaudet::Traitement_beaudet(): Traitement_image(){}
Traitement_beaudet::~Traitement_beaudet(){}
void Traitement_beaudet::appliquer(){
    bool write_img = getTraitementEcriture();
    if(AFFICHAGE_DEBUG){
        std::cout << "BEAUDET APPLIQUER" << std::endl;
        std::cout << "WRITE IMG : " << write_img << std::endl;
    }

    for (int k = 0; k < this->img_in.size(); k++){
        if(AFFICHAGE_DEBUG){
            std::cout << "IMG IN " << std::endl;
            std::cout << *img_in[k];
        }

        IMG * img_cur = this->img_in[k];
        OCTET *Img_deriv_XX,*Img_deriv_YY,*Img_deriv_XY,*Img_deriv_YX,*Img_out;
        
        int nTaille = img_cur->nH/2 * img_cur->nW/2;

        allocation_tableau(Img_deriv_XX, OCTET, nTaille);
        allocation_tableau(Img_deriv_YY, OCTET, nTaille);
        allocation_tableau(Img_deriv_XY, OCTET, nTaille);
        allocation_tableau(Img_deriv_YX, OCTET, nTaille);
        allocation_tableau(Img_out, OCTET, nTaille);


        for (int i = 0; i < img_cur->nH/2; i++){
            for (int j = 0; j < img_cur->nW/2; j++){
                Img_deriv_XX[i*img_cur->nW/2+j] = img_cur->img_data[i*img_cur->nW + j];
                Img_deriv_XY[i*img_cur->nW/2+j] = img_cur->img_data[i*img_cur->nW + j + img_cur->nW/2];
                Img_deriv_YX[i*img_cur->nW/2+j] = img_cur->img_data[(i+img_cur->nH/2)*img_cur->nW + j];
                Img_deriv_YY[i*img_cur->nW/2+j] = img_cur->img_data[(i+img_cur->nH/2)*img_cur->nW + j + img_cur->nW/2];
            }   
        }

        beaudet(img_cur->nH/2, img_cur->nW/2,Img_deriv_XX,Img_deriv_YY,Img_deriv_XY,Img_deriv_YX,Img_out);
        
        char* nom_beaudet = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_beaudet,"%s_beaudet", img_cur->nom);  
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH/2, img_cur->nW/2, Img_out, write_img,nom_beaudet));
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_beaudet::afficher(std::ostream& os){
    os << "TRAITEMENT BEAUDET" << std::endl;
    Traitement_image::afficher(os);
}


/**   ====================         
 *  TRAITEMENT HARRIS-STEPHAN        
 *    ====================**/
class Traitement_HarrisStephan : public Traitement_image{
    private :
    float alpha;
    int flg_harris;
    void harris_stephan(int nH, int nW, OCTET *Img_deriv_X,OCTET *Img_deriv_Y, OCTET *Img_out);

    public :
        Traitement_HarrisStephan(float alpha, int flg = FLG_HARRIS_GAUCHE_BAS);
        virtual ~Traitement_HarrisStephan();
        void appliquer();
        virtual void afficher(std::ostream& os);
};

void Traitement_HarrisStephan::harris_stephan(int nH, int nW, OCTET *Img_deriv_X,OCTET *Img_deriv_Y, OCTET *Img_out){
    if(AFFICHAGE_DEBUG) {
        std::cout << "HARRIS STEPHAN TRAITEMENT" << std::endl;
        std::cout << "  ->NH : " << nH << std::endl;
        std::cout << "  ->NW : " << nW << std::endl;
    }

    for (int i = 0; i < nH; i++){
        for (int j = 0; j < nW; j++){
            int FIxIx = Img_deriv_X[i*nW+j]*Img_deriv_X[i*nW+j];
            int FIyIy = Img_deriv_Y[i*nW+j]*Img_deriv_Y[i*nW+j];
            int FIyIx = Img_deriv_X[i*nW+j]*Img_deriv_Y[i*nW+j];

            Img_out[i*nW+j] = FIxIx * FIyIy - FIyIx*FIyIx - (alpha * (FIxIx + FIyIy));
        }
    }
}
Traitement_HarrisStephan::Traitement_HarrisStephan(float alpha, int flg): Traitement_image(){
    this->alpha = alpha;
    this->flg_harris = flg;
}
Traitement_HarrisStephan::~Traitement_HarrisStephan(){}
void Traitement_HarrisStephan::appliquer(){
    bool write_img = getTraitementEcriture();
     if(AFFICHAGE_DEBUG){
        std::cout << "HARRIS APPLIQUER" << std::endl;
        std::cout << "WRITE IMG : " << write_img << std::endl;
    }

    for (int k = 0; k < this->img_in.size(); k++){
        if(AFFICHAGE_DEBUG){
            std::cout << "IMG IN " << std::endl;
            std::cout << *img_in[k];
        }

        IMG * img_cur = this->img_in[k];
        
        OCTET *Img_deriv_X,*Img_deriv_Y,*Img_out;
        
        int nTaille = img_cur->nH * img_cur->nW;

        allocation_tableau(Img_deriv_X, OCTET, nTaille);
        allocation_tableau(Img_deriv_Y, OCTET, nTaille);
        allocation_tableau(Img_out, OCTET, nTaille);

        for (int i = 1; i < img_cur->nH-1; i++){
            for (int j = 1; j < img_cur->nW-1; j++){
                float h_bas =       img_cur->img_data[(i+1)*img_cur->nW+j] - img_cur->img_data[i*img_cur->nW+j];
                float w_droite =    img_cur->img_data[i*img_cur->nW+j+1] - img_cur->img_data[i*img_cur->nW+j];
                float h_haut =      img_cur->img_data[(i-1)*img_cur->nW+j] - img_cur->img_data[i*img_cur->nW+j];
                float w_gauche =    img_cur->img_data[i*img_cur->nW+j-1] - img_cur->img_data[i*img_cur->nW+j];
                
                switch (this->flg_harris){
                case FLG_HARRIS_GAUCHE_BAS : 
                    Img_deriv_X[i*img_cur->nW+j] = w_gauche;
                    Img_deriv_Y[i*img_cur->nW+j] = h_bas;
                    break;
                case FLG_HARRIS_GAUCHE_HAUT :
                    Img_deriv_X[i*img_cur->nW+j] = w_gauche;
                    Img_deriv_Y[i*img_cur->nW+j] = h_haut;
                    break;
                case FLG_HARRIS_DROITE_BAS :
                    Img_deriv_X[i*img_cur->nW+j] = w_droite;
                    Img_deriv_Y[i*img_cur->nW+j] = h_bas;
                    break;
                case FLG_HARRIS_DROITE_HAUT :
                    Img_deriv_X[i*img_cur->nW+j] = w_droite;
                    Img_deriv_Y[i*img_cur->nW+j] = h_haut;
                    break;
                }
            }
        }

        harris_stephan(img_cur->nH, img_cur->nW,Img_deriv_X,Img_deriv_Y,Img_out);
        
        char* nom_harris = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_harris,"%s_harrisStephan", img_cur->nom);  
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, Img_out, write_img,nom_harris));
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_HarrisStephan::afficher(std::ostream& os){
    os << "TRAITEMENT HARRIS" << std::endl;
    os << "ALPHA : "  << alpha << std::endl;
    Traitement_image::afficher(os);
}

#endif