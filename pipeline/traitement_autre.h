#ifndef TRAITEMENTAUTRE_H
#define TRAITEMENTAUTRE_H

#include "pipeline_first.h"

#define FLG_EROSIONDILATATION_CROIX 0
#define FLG_EROSIONDILATATION_ORDRE_DEUX 1
#define FLG_EROSION 2
#define FLG_DILATATION 3


/**   ==================== ==================== ==================== ====================         
 *              ==================== ==================== ====================
 *    ==================== ==================== ==================== ====================**/
/**   ====================         
 *    TRAITEMENT SOUSTRACTION        
 *    ====================**/

class Traitement_soustraction : public Traitement_image{
    private :
    void soustraction(int nH, int nW, const OCTET* ImgIn,OCTET* Imgres);
    public :
        Traitement_soustraction();
        virtual ~Traitement_soustraction();
        void appliquer();
        virtual void afficher(std::ostream& os);
};
 
void Traitement_soustraction::soustraction(int nH, int nW, const OCTET* ImgIn, OCTET* Imgres){
    for (int i = 0; i < nH; i++){
        for (int j = 0; j < nW/2; j++){
            Imgres[i*nW/2 + j] = ImgIn[i*nW + j] - ImgIn[i*nW + j+nW/2];
        }
    }
}
Traitement_soustraction::Traitement_soustraction(): Traitement_image(){

}
Traitement_soustraction::~Traitement_soustraction(){
    
}
//On va soustraire la première partie de l'image (partie gauche) à la seconde partie de l'image (partie droite)
void Traitement_soustraction::appliquer(){
    bool write_img = getTraitementEcriture();
    if(AFFICHAGE_DEBUG){
        std::cout << "SOUSTRACTION APPLIQUER" << std::endl;
        std::cout << "WRITE IMG : " << write_img << std::endl;
    }
    for (size_t k = 0; k < this->img_in.size(); k++){
        if(AFFICHAGE_DEBUG){
            std::cout << "IMG IN " << std::endl;
            std::cout << *img_in[k];
        }

        IMG * img_cur = this->img_in[k];

        OCTET *img_sous;
        int nTaille = img_cur->nH*(img_cur->nW/2);
        allocation_tableau(img_sous, OCTET, nTaille);
        soustraction(img_cur->nH, img_cur->nW, img_cur->img_data , img_sous);

        char* nom_sous = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_sous,"%s_soustraction", img_cur->nom); 
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW/2, img_sous, write_img,nom_sous));
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_soustraction::afficher(std::ostream& os){
    os << "TRAITEMENT SOUSTRACTION" << std::endl;
    Traitement_image::afficher(os);
}


/**   ====================         
 *    TRAITEMENT UNIFORME        
 *    ====================**/
class Traitement_uniforme : public Traitement_image{
    private :
    int add;
    public :
        Traitement_uniforme(int add);
        virtual ~Traitement_uniforme();
        void appliquer();
        virtual void afficher(std::ostream& os);
};
    
Traitement_uniforme::Traitement_uniforme(int add) : Traitement_image(){
    this->add = add;
}
Traitement_uniforme::~Traitement_uniforme(){
}
void Traitement_uniforme::appliquer(){
    bool write_img = getTraitementEcriture();
    if(AFFICHAGE_DEBUG){
        std::cout << "UNIFORME APPLIQUER" << std::endl;
        std::cout << "WRITE IMG : " << write_img << std::endl;
    }
    for (size_t k = 0; k < this->img_in.size(); k++){
        if(AFFICHAGE_DEBUG){
            std::cout << "IMG IN " << std::endl;
            std::cout << *img_in[k];
        }

        IMG * img_cur = this->img_in[k];
        IMG * img_add = new IMG(img_cur);
        for (size_t i = 0; i < img_cur->nH*img_cur->nW; i++){
            int val = img_cur->img_data[i] + add;
            img_add->img_data[i] = ((val<0)? 0 : ((val>255)? 255 : val));
        }
        char*nom = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom,"%s_uniforme%d",img_add->nom,add);
        img_add->nom = nom;
        this->img_out.push_back(img_add);
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }    
}
void Traitement_uniforme::afficher(std::ostream& os){
    os << "TRAITEMENT UNIFORME : " << add << std::endl;
    Traitement_image::afficher(os);
}


/**   ====================         
 * TRAITEMENT EROSION DILATATION        
 *    ====================**/
class Traitement_erosiondilatation : public Traitement_image{
    private :
        int couleur_fond;
        int flg;
        int ordre;
        int flg_ero_dil;
        const char* nom;

        int voisins_croix(int i, int j,int n, OCTET* imint);
        int voisins_ordredeux(int i, int j,int n, OCTET* imint);
    public :
        Traitement_erosiondilatation(int c,int ordre,int flg_erodil, int flg = FLG_EROSIONDILATATION_CROIX);
        virtual ~Traitement_erosiondilatation();
        void appliquer();
        virtual void afficher(std::ostream& os);
};

int Traitement_erosiondilatation::voisins_croix(int i, int j,int n, OCTET* imint){
    int val = couleur_fond;
    if(this->flg_ero_dil == FLG_DILATATION){
        //Si on est pas en bas
        if((i+1)<n && val<imint[(i+1)*n+j]){
            val = imint[(i+1)*n+j];
        }

        //Si on est pas à droite
        if((j+1)<n && val < imint[i*n+j+1]){
            val = imint[i*n+j+1];
        }

        //Si on est pas a gauche
        if(j!=0 && val <imint[i*n+j-1]){
            val = imint[i*n+j-1];
        }

        //Si on est pas en haut
        if(i!=0 && val < imint[(i-1)*n+j]){
            val = imint[(i-1)*n+j];
        }
    }

    if(this->flg_ero_dil == FLG_EROSION){
        //Si on est pas en bas
        if((i+1)<n && val>imint[(i+1)*n+j]){
            val = imint[(i+1)*n+j];
        }

        //Si on est pas à droite
        if((j+1)<n && val > imint[i*n+j+1]){
            val = imint[i*n+j+1];
        }

        //Si on est pas a gauche
        if(j!=0 && val >imint[i*n+j-1]){
            val = imint[i*n+j-1];
        }

        //Si on est pas en haut
        if(i!=0 && val > imint[(i-1)*n+j]){
            val = imint[(i-1)*n+j];
        }
    }
    return val;
}
int Traitement_erosiondilatation::voisins_ordredeux(int i, int j,int n, OCTET* imint){
    int val = this->couleur_fond;
    if(this->flg_ero_dil == FLG_DILATATION){
        //Si on est pas en bas
        if((i+1)<n && val<imint[(i+1)*n+j]){
            val = imint[(i+1)*n+j];
        }

        //Si on est pas à droite
        if((j+1)<n && val < imint[i*n+j+1]){
            val = imint[i*n+j+1];
        }

        //Si on est pas a gauche
        if(j!=0 && val <imint[i*n+j-1]){
            val = imint[i*n+j-1];
        }

        //Si on est pas en haut
        if(i!=0 && val < imint[(i-1)*n+j]){
            val = imint[(i-1)*n+j];
        }

        //Si on est pas en bas à gauche
        if((i+1)<n && j!=0 && val < imint[(i+1)*n+(j-1)]){
            val = imint[(i+1)*n+(j-1)];
        }

        //Si on est pas en bas a droite
        if((i+1)<n && (j+1)<n && val< imint[(i+1)*n+(j+1)]){
            val = imint[(i+1)*n+(j+1)];
        }

        //Si on est pas en haut a gauche
        if(i!=0 && j!=0 && val <imint[(i-1)*n+(j-1)]){
            val = imint[(i-1)*n+(j-1)];
        }

        //Si on est pas en haut a droite
        if(i!=0 && (j+1)<n && val <imint[(i-1)*n+(j+1)]){
            val  = imint[(i-1)*n+(j+1)];
        }
    }

    if(this->flg_ero_dil == FLG_EROSION){
        //Si on est pas en bas
        if((i+1)<n && val>imint[(i+1)*n+j]){
            val = imint[(i+1)*n+j];
        }

        //Si on est pas à droite
        if((j+1)<n && val > imint[i*n+j+1]){
            val = imint[i*n+j+1];
        }

        //Si on est pas a gauche
        if(j!=0 && val >imint[i*n+j-1]){
            val = imint[i*n+j-1];
        }

        //Si on est pas en haut
        if(i!=0 && val > imint[(i-1)*n+j]){
            val = imint[(i-1)*n+j];
        }

        //Si on est pas en bas à gauche
        if((i+1)<n && j!=0 && val > imint[(i+1)*n+(j-1)]){
            val = imint[(i+1)*n+(j-1)];
        }

        //Si on est pas en bas a droite
        if((i+1)<n && (j+1)<n && val> imint[(i+1)*n+(j+1)]){
            val = imint[(i+1)*n+(j+1)];
        }

        //Si on est pas en haut a gauche
        if(i!=0 && j!=0 && val >imint[(i-1)*n+(j-1)]){
            val = imint[(i-1)*n+(j-1)];
        }

        //Si on est pas en haut a droite
        if(i!=0 && (j+1)<n && val >imint[(i-1)*n+(j+1)]){
            val  = imint[(i-1)*n+(j+1)];
        }
    }
    return val;
}
Traitement_erosiondilatation::Traitement_erosiondilatation(int c,int ordre, int flg_erodil, int flg) : Traitement_image(){
    this->couleur_fond = c;
    this->flg = flg;
    this->ordre = ordre;
    this->flg_ero_dil = flg_erodil;
    this->nom = ((this->flg_ero_dil == FLG_DILATATION)? "dilatation" : (this->flg_ero_dil == FLG_EROSION)? "erosion" : "nieronidila");
}
Traitement_erosiondilatation::~Traitement_erosiondilatation(){
}
void Traitement_erosiondilatation::appliquer(){
    bool write_img = getTraitementEcriture();
    if(AFFICHAGE_DEBUG){
        std::cout << this->nom <<" APPLIQUER" << std::endl;
        std::cout << "WRITE IMG : " << write_img << std::endl;
    }

    for (size_t k = 0; k < this->img_in.size(); k++){
        if(AFFICHAGE_DEBUG){
            std::cout << "IMG IN " << std::endl;
            std::cout << *img_in[k];
        }

        IMG * img_cur = this->img_in[k];
        OCTET *img_erosiondil;
        int nTaille = img_cur->nH*img_cur->nW;
        allocation_tableau(img_erosiondil, OCTET, nTaille);

        OCTET *buffer;
        allocation_tableau(buffer, OCTET, nTaille);
        for (int i = 0; i < nTaille; i++){
            buffer[i] = img_cur->img_data[i]; 
        }

        for (int o = 0; o < ordre; o++){
            for (int i=0; i < img_cur->nH; i++){
                for (int j=0; j < img_cur->nW; j++){
                    img_erosiondil[i*img_cur->nW+j] = ((flg==FLG_EROSIONDILATATION_CROIX)? voisins_croix(i,j,img_cur->nW,buffer) : voisins_ordredeux(i,j,img_cur->nW,buffer));
                }
            }
            for (int i = 0; i < nTaille; i++){
                buffer[i] = img_erosiondil[i];
            }
        }      

        char* nom_erosiondil = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_erosiondil,"%s_%s%d", img_cur->nom,this->nom,this->ordre); 
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, img_erosiondil, write_img,nom_erosiondil));
    }  
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_erosiondilatation::afficher(std::ostream& os){
    os << "TRAITEMENT " << nom << " : " << couleur_fond << std::endl;
    os << "ORDRE : " << ordre << std::endl;
    Traitement_image::afficher(os);
}

#endif