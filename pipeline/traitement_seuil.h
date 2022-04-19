#ifndef SEUIL_H
#define SEUIL_H

#include "pipeline_first.h"

//Seuillage inférieur ou supérieur au seuil
#define FLG_SEUILLAGE_INF 0
#define FLG_SEUILLAGE_SUP 1

#define FLG_HISTERESIS_4V 0
#define FLG_HISTERESIS_8V 1


/**   ==================== ==================== ==================== ====================         
 *              ==================== ==================== ====================
 *    ==================== ==================== ==================== ====================**/
/**   ====================         
 *    TRAITEMENT SEUILLAGE        
 *    ====================**/
class Traitement_seuillage : public Traitement_image{
    private :
    int flg_seuil;
    int val_NDG;

    public :
        Traitement_seuillage(int val_NDG, int flg = FLG_SEUILLAGE_INF);
        virtual ~Traitement_seuillage();
        void appliquer();
        virtual void afficher(std::ostream& os);
};

Traitement_seuillage::Traitement_seuillage(int val_NDG, int flg): Traitement_image(){
    this->flg_seuil = flg;
    this->val_NDG = val_NDG;
}
Traitement_seuillage::~Traitement_seuillage(){}
void Traitement_seuillage::appliquer(){
    bool write_img = getTraitementEcriture();
    if(AFFICHAGE_DEBUG){
        std::cout << "SEUILLAGE APPLIQUER" << std::endl;
        std::cout << "WRITE IMG : " << write_img << std::endl;
    }

    for (size_t i = 0; i < this->img_in.size(); i++){
        if(AFFICHAGE_DEBUG){
            std::cout << "IMG IN " << std::endl;
            std::cout << *img_in[i];
        }

        IMG * img_cur = this->img_in[i];
        OCTET *img_seuillee;
        int nTaille = img_cur->nH*img_cur->nW;
        allocation_tableau(img_seuillee, OCTET, nTaille);

        for (size_t j = 0; j < nTaille; j++){
            img_seuillee[j] = ((this->flg_seuil == FLG_SEUILLAGE_INF)? ((img_cur->img_data[j]<this->val_NDG)? 255 : 0) : ((this->flg_seuil == FLG_SEUILLAGE_SUP)? ((img_cur->img_data[j]>this->val_NDG)? 255 : 0) : img_cur->img_data[j]));
        }

        char* nom_seuil = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_seuil,"%s_seuillee%d", img_cur->nom, this->val_NDG); 
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, img_seuillee, write_img,nom_seuil));
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_seuillage::afficher(std::ostream& os){
    os << "TRAITEMENT SEUILLAGE : " << this->val_NDG << std::endl;
    os << ((this->flg_seuil == FLG_SEUILLAGE_INF)? "INF" : (this->flg_seuil == FLG_SEUILLAGE_SUP)? "SUP" : "NSP") << std::endl;
    Traitement_image::afficher(os);
}


//HYSTERESIS
class Traitement_hysteresis : public Traitement_image{
    private :
    int val_NDG_min;
    int val_NDG_max;
    int couleur_contour;
    int couleur_fond;

    int flg_hist;

    int getvalVoisins8(OCTET* img, int indice, int nH, int nW);
    int getvalVoisins4(OCTET* img, int indice, int nH, int nW);


    public :
        Traitement_hysteresis(int val_NDG_min, int val_NDG_max, int couleur_contour = 255, int flg_hist = FLG_HISTERESIS_8V);
        virtual ~Traitement_hysteresis();
        void appliquer();
        virtual void afficher(std::ostream& os);
};

int Traitement_hysteresis::getvalVoisins8(OCTET* img, int indice, int nH, int nW){
    //Si le px n'est pas à gauche
    if(indice%nW != 0 && img[indice-1] == couleur_contour){
        return couleur_contour;
    }

    //Si le px n'est pas en haut
    if(indice/nW != 0 && img[indice-nW] == couleur_contour){
         return couleur_contour;
    }

    //Si le px n'est pas en bas
    if(indice/nW != (nH-1) && img[indice+nW] == couleur_contour){
        return couleur_contour;
    }

    //Si le px n'est pas à droite
    if((indice+1)%nW != 0 && img[indice+1] == couleur_contour){
        return couleur_contour;
    }

    //Si on est pas en bas à gauche
    if(indice%nW != 0 && indice/nW != (nH-1) && img[indice-1+nW] == couleur_contour){
        return couleur_contour;
    }

    //Si on est pas en bas a droite
    if(indice/nW != (nH-1) && (indice+1)%nW != 0 && img[indice+1+nW] == couleur_contour){
        return couleur_contour;
    }

    //Si on est pas en haut a gauche
    if(indice/nW != 0 && indice%nW != 0 && img[indice-nW-1] == couleur_contour){
        return couleur_contour;
    }

    //Si on est pas en haut a droite
    if(indice/nW != 0 && (indice+1)%nW != 0 && img[indice-nW+1] == couleur_contour){
        return couleur_contour;
    }

    return couleur_fond;
}
int Traitement_hysteresis::getvalVoisins4(OCTET* img, int indice, int nH, int nW){
    //Si le px n'est pas à gauche
    if(indice%nW != 0 && img[indice-1] == couleur_contour){
        return couleur_contour;
    }

    //Si le px n'est pas en haut
    if(indice/nW != 0 && img[indice-nW] == couleur_contour){
         return couleur_contour;
    }

    //Si le px n'est pas en bas
    if(indice/nW != (nH-1) && img[indice+nW] == couleur_contour){
        return couleur_contour;
    }

    //Si le px n'est pas à droite
    if((indice+1)%nW != 0 && img[indice+1] == couleur_contour){
        return couleur_contour;
    }

    return couleur_fond;
}
Traitement_hysteresis::Traitement_hysteresis(int val_NDG_min, int val_NDG_max,int couleur_contour, int flg_hist): Traitement_image(){
   this->val_NDG_min = val_NDG_min;
   this->val_NDG_max = val_NDG_max;
   this->couleur_contour = couleur_contour;
   this->couleur_fond = 255-this->couleur_contour;
   this->flg_hist = flg_hist;
}
Traitement_hysteresis::~Traitement_hysteresis(){}
void Traitement_hysteresis::appliquer(){
    bool write_img = getTraitementEcriture();
    if(AFFICHAGE_DEBUG){
        std::cout << "HYSTERESIS APPLIQUER" << std::endl;
        std::cout << "WRITE IMG : " << write_img << std::endl;
    }

    for (size_t k = 0; k < this->img_in.size(); k++){
        if(AFFICHAGE_DEBUG){
            std::cout << "IMG IN " << std::endl;
            std::cout << *img_in[k];
        }

        IMG * img_cur = this->img_in[k];
        OCTET *img_seuillee;
        int nTaille = img_cur->nH*img_cur->nW;
        allocation_tableau(img_seuillee, OCTET, nTaille);

        for (int i=0; i < nTaille; i++){
            if (img_cur->img_data[i]<val_NDG_min) {
                img_seuillee[i]=couleur_fond;
            }
            if (img_cur->img_data[i]>val_NDG_max) {
                img_seuillee[i]=couleur_contour;
            }
        }
        for (int i=0; i < img_cur->nH; i++){
            for (int j=0; j < img_cur->nW; j++){
                if ( img_cur->img_data[i*img_cur->nW+j] <= this->val_NDG_max && img_cur->img_data[i*img_cur->nW+j] >= this->val_NDG_min){ 
                    img_seuillee[i*img_cur->nW+j] = ((flg_hist==FLG_HISTERESIS_4V)? getvalVoisins4(img_seuillee, i*img_cur->nW+j,img_cur->nH,img_cur->nW) : getvalVoisins8(img_seuillee, i*img_cur->nW+j,img_cur->nH,img_cur->nW));
                }
            }
        }
        
        char* nom_seuil = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_seuil,"%s_hysteresis%dx%d", img_cur->nom, this->val_NDG_min,this->val_NDG_max); 
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, img_seuillee, write_img,nom_seuil));
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_hysteresis::afficher(std::ostream& os){
    os << "TRAITEMENT HYSTERESIS" << std::endl;
    os << "Smin : " << this->val_NDG_min <<  std::endl;
    os << "Smax : " << this->val_NDG_max <<  std::endl;
    os << "couleur fond : " << this->couleur_fond <<  std::endl;
    os << "couleur contour : " << this->couleur_contour <<  std::endl;
    os << "NB VOIS : " << ((this->flg_hist ==FLG_HISTERESIS_4V )? "4": ((this->flg_hist ==FLG_HISTERESIS_8V)? "8" : "NSP")) <<  std::endl;
    Traitement_image::afficher(os);
}


//EXTREMUM : min et max
class Traitement_seuillageExtremums : public Traitement_image{
    private :
    int val_NDG_min;
    int val_NDG_max;

    public :
        Traitement_seuillageExtremums(int val_NDGmin, int val_NDGmax);
        virtual ~Traitement_seuillageExtremums();
        void appliquer();
        virtual void afficher(std::ostream& os);
};

Traitement_seuillageExtremums::Traitement_seuillageExtremums(int val_NDGmin, int val_NDGmax): Traitement_image(){
    this->val_NDG_min = val_NDGmin;
    this->val_NDG_max = val_NDGmax;
}
Traitement_seuillageExtremums::~Traitement_seuillageExtremums(){}
void Traitement_seuillageExtremums::appliquer(){
    bool write_img = getTraitementEcriture();
    if(AFFICHAGE_DEBUG){
        std::cout << "SEUILLAGE APPLIQUER" << std::endl;
        std::cout << "WRITE IMG : " << write_img << std::endl;
    }

    for (size_t i = 0; i < this->img_in.size(); i++){
        if(AFFICHAGE_DEBUG){
            std::cout << "IMG IN " << std::endl;
            std::cout << *img_in[i];
        }

        IMG * img_cur = this->img_in[i];
        OCTET *img_seuillee;
        int nTaille = img_cur->nH*img_cur->nW;
        allocation_tableau(img_seuillee, OCTET, nTaille);

        for (int j = 0; j < nTaille; j++){
            if(img_cur->img_data[j] < val_NDG_min){
                img_seuillee[j] = val_NDG_min;
            }else if(img_cur->img_data[j] > val_NDG_max){
                img_seuillee[j] = val_NDG_max;
            }else{
                img_seuillee[j] = img_cur->img_data[j];
            }
        }

        char* nom_seuil = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_seuil,"%s_seuilleeMin%dMax%d", img_cur->nom, this->val_NDG_min,this->val_NDG_max); 
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, img_seuillee, write_img,nom_seuil));
        
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_seuillageExtremums::afficher(std::ostream& os){
    os << "TRAITEMENT SEUILLAGE EXTREMUMS "  << std::endl;
    os << "Smin : " << this->val_NDG_min <<  std::endl;
    os << "Smax : " << this->val_NDG_max <<  std::endl;
    Traitement_image::afficher(os);
}


//MOYENNE
class Traitement_seuillageMoy : public Traitement_image{
    private :
    int flg_seuil;

    public :
        Traitement_seuillageMoy(int flg = FLG_SEUILLAGE_INF);
        virtual ~Traitement_seuillageMoy();
        void appliquer();
        virtual void afficher(std::ostream& os);
};

Traitement_seuillageMoy::Traitement_seuillageMoy(int flg): Traitement_image(){
    this->flg_seuil = flg;
}
Traitement_seuillageMoy::~Traitement_seuillageMoy(){}
void Traitement_seuillageMoy::appliquer(){
    bool write_img = getTraitementEcriture();
    if(AFFICHAGE_DEBUG){
        std::cout << "SEUILLAGE APPLIQUER" << std::endl;
        std::cout << "WRITE IMG : " << write_img << std::endl;
    }

    for (size_t i = 0; i < this->img_in.size(); i++){
        if(AFFICHAGE_DEBUG){
            std::cout << "IMG IN " << std::endl;
            std::cout << *img_in[i];
        }

        IMG * img_cur = this->img_in[i];
        OCTET *img_seuillee;
        int nTaille = img_cur->nH*img_cur->nW;
        allocation_tableau(img_seuillee, OCTET, nTaille);

        float S = 0;
        for (int j=0; j < nTaille; j++){
            S+= img_cur->img_data[i];
        }
        S/=nTaille;
        if(AFFICHAGE_DEBUG){
            std::cout << "MOYENNE : " << S << std::endl;
        }
        int seuil = S;

        for (size_t j = 0; j < nTaille; j++){
            img_seuillee[j] = ((this->flg_seuil == FLG_SEUILLAGE_INF)? ((img_cur->img_data[j]<seuil)? 255 : 0) : ((this->flg_seuil == FLG_SEUILLAGE_SUP)? ((img_cur->img_data[j]>seuil)? 255 : 0) : img_cur->img_data[j]));
        }

        char* nom_seuil = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_seuil,"%s_seuillee%d", img_cur->nom, seuil); 
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, img_seuillee, write_img,nom_seuil));
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_seuillageMoy::afficher(std::ostream& os){
    os << "TRAITEMENT SEUILLAGE MOYEN : " << std::endl;
    os << ((this->flg_seuil == FLG_SEUILLAGE_INF)? "INF" : (this->flg_seuil == FLG_SEUILLAGE_SUP)? "SUP" : "NSP") << std::endl;
    Traitement_image::afficher(os);
}


//MEDIAN
class Traitement_seuillageMedian : public Traitement_image{
    private :
    int nombre_seuil;
    float* proportions;

    void seuiller(OCTET* ImgIn, int nH, int nW,OCTET* ImgOut);

    public :
        Traitement_seuillageMedian(int nbS);
        virtual ~Traitement_seuillageMedian();
        void appliquer();

        virtual void afficher(std::ostream& os);
        Traitement_seuillageMedian * ajouterProportion(int ind, float p);
};

void Traitement_seuillageMedian::seuiller(OCTET* ImgIn, int nH, int nW,OCTET* ImgOut){
    int nTaille = nH * nW;
    
    float prop = 0.f;
    for(int i = 0; i<nombre_seuil; i++){
        if(proportions[i]<0){
            prop = (1.f-prop) * 1.f/(float)((float) nombre_seuil - (float)i);
            if( AFFICHAGE_DEBUG){
                std::cout << "PROPORTION RESTANTE POUR [" << i << "," << nombre_seuil-1 << "] : " << prop << std::endl;
            }
            break;
        }else{
            prop+=proportions[i];
        }
        if( AFFICHAGE_DEBUG){
            std::cout << "PROPORTION POUR [" << i << "] : " << proportions[i] << std::endl;
        }
    }

    for(int i = 0; i< nombre_seuil; i++){
        if(proportions[i]<0){
            proportions[i] = prop;
        }
        if( AFFICHAGE_DEBUG){
            std::cout << "PROPORTION POUR [" << i << "] : " << proportions[i] << std::endl; 
        }
    }

    
    int* seuil = new int[nombre_seuil]{0};
    seuil[0] = nTaille*proportions[0];
    for(int i = 1; i< nombre_seuil-1; i++){
        seuil[i] = nTaille*proportions[i] + seuil[i-1];
    }

    //HISTOGRAMME
    int* nbElementsLus = new int[256]{};
    for (int i=0; i < nH; i++){
        for (int j=0; j < nW; j++){
            nbElementsLus[ImgIn[i*nW+j]]++;
        }
    }

    for(int i = 0; i<nombre_seuil-1; i++){
        int nb = 0;
        for(int j = 0; j<256;j++){
            nb+=nbElementsLus[j];
            if(nb>=seuil[i]){
                seuil[i] = j;
                break;
            }
        }
    }

    int* valeur_seuil = new int[nombre_seuil]{0};
    valeur_seuil[0] = 0;
    valeur_seuil[nombre_seuil-1] = 255;

    for(int i = 1; i<nombre_seuil-1; i++){
        valeur_seuil[i] = (seuil[i-1]+seuil[i])/2;
    }

    if(AFFICHAGE_DEBUG){
        std::cout << "ON A LES SEUILS SUIVANTS : " << std::endl;
        for(int i = 0; i< nombre_seuil;i++){
            std::cout << "  seuil[" << i+1 << "]" << std::endl; 
            std::cout << "      -> val : " << seuil[i] << std::endl;
            std::cout << "      -> pix : " << valeur_seuil[i] << std::endl; 
        }
    }

    for (int i=0; i < nH; i++){
        for (int j=0; j < nW; j++){
            bool changement = false;
            for(int k = 0; k<nombre_seuil-1; k++){
                if ( ImgIn[i*nW+j] < seuil[k]){
                    ImgOut[i*nW+j]= valeur_seuil[k];
                    changement = true;
                    break;
                }
                if(changement){
                    break;
                }
            }
            if(!changement){
                ImgOut[i*nW+j] = 255;
            }
        }
    }
}
Traitement_seuillageMedian::Traitement_seuillageMedian(int nbS): Traitement_image(){
    this->nombre_seuil = nbS;
    this->proportions = new float[nbS];
    for(int i = 0; i<this->nombre_seuil; i++) this->proportions[i] = -1.f;
}
Traitement_seuillageMedian::~Traitement_seuillageMedian(){
    delete this->proportions;
}
void Traitement_seuillageMedian::appliquer(){
    bool write_img = getTraitementEcriture();
    if(AFFICHAGE_DEBUG){
        std::cout << "SEUILLAGE MED APPLIQUER" << std::endl;
        std::cout << "WRITE IMG : " << write_img << std::endl;
    }

    for (size_t i = 0; i < this->img_in.size(); i++){
        if(AFFICHAGE_DEBUG){
            std::cout << "IMG IN " << std::endl;
            std::cout << *img_in[i];
        }

        IMG * img_cur = this->img_in[i];
        OCTET *img_seuillee;
        int nTaille = img_cur->nH*img_cur->nW;
        allocation_tableau(img_seuillee, OCTET, nTaille);

        seuiller(img_cur->img_data,img_cur->nH,img_cur->nW,img_seuillee);

        char* nom_seuil = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_seuil,"%s_seuilleeMED%d", img_cur->nom,nombre_seuil); 
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, img_seuillee, write_img,nom_seuil));
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_seuillageMedian::afficher(std::ostream& os){
    os << "TRAITEMENT SEUILLAGE MED : " << nombre_seuil << std::endl;
    Traitement_image::afficher(os);
}
Traitement_seuillageMedian * Traitement_seuillageMedian::ajouterProportion(int ind, float p){
    if(ind >= 0 && ind < nombre_seuil) this->proportions[ind] = p;
    return this;
}


//MINIMUM LOCAUX
//NdgSimp : proportion de valeurs à prendre en compte : nombre de plages dans l'histogramme avant détection des minimums locaux
class Traitement_seuillageMinLocaux : public Traitement_image{
    private :
    int NDGSimp;

    void seuiller(OCTET* ImgIn, int nH, int nW,OCTET* ImgOut);

    public :
        Traitement_seuillageMinLocaux(int ndgSimp);
        virtual ~Traitement_seuillageMinLocaux();
        void appliquer();
        virtual void afficher(std::ostream& os);
};

void Traitement_seuillageMinLocaux::seuiller(OCTET* ImgIn, int nH, int nW,OCTET* ImgOut){
    int nTaille = nH * nW;
    int nb_pas = 256/NDGSimp;

    int* histoSeuille = new int[NDGSimp]{};
    //HISTOGRAMME
    for (int i=0; i < nH; i++){
        for (int j=0; j < nW; j++){
            histoSeuille[(int)floor(((double)ImgIn[i*nW+j]/256.0) * NDGSimp)]++;
        }
    }

    //Détermination des minimums locaux
    int * seuil = new int[NDGSimp]{};
    seuil[0] = seuil[NDGSimp-1] = -1;
    for (int i = 1; i < NDGSimp-1; i++){
        if(histoSeuille[i]<histoSeuille[i+1] && histoSeuille[i]<histoSeuille[i-1]){
            seuil[i] = 1;
        }else{
            seuil[i] = -1;
        }
    }

    int nbMinimLoc = 0;
    for (int i = 0; i < NDGSimp; i++){
        if(AFFICHAGE_DEBUG) std::cout << "MINIMUM LOCAL DE " << i*nb_pas << " ? ";
        if(seuil[i]>0){
            if(AFFICHAGE_DEBUG)  std::cout << "OUI" << std::endl;
            nbMinimLoc++;
        }else{
            if(AFFICHAGE_DEBUG) std::cout << "NON" << std::endl;
        }
    }

    int* seuils = new int[nbMinimLoc]{};
    int cpt = 0;
    for (int i = 0; i < NDGSimp; i++){
        if(seuil[i]>0){
            seuils[cpt] = i*nb_pas;
            cpt++; 
        }
    }
    
    if(AFFICHAGE_DEBUG) {
        std::cout << " SEUILS DE L'IMAGE : " << std::endl;
        for(int i = 0; i<nbMinimLoc; i++){
            std::cout << "   " << seuils[i] << std::endl;
        }
    }

    int* valSeuils = new int[nbMinimLoc+2]{};
    valSeuils[0] = 0;
    valSeuils[nbMinimLoc+1] = 255;
    for (int i = 0; i < nbMinimLoc-1; i++){
        int vs=(seuils[i] + seuils[i+1])/2;
        valSeuils[i] = (vs<0)? 0: (vs>255)? 255 : vs;
    }

    if(AFFICHAGE_DEBUG) {
        std::cout << "VALEUR DE x<" << seuils[0] << " : " << 0 <<  std::endl;
        for (int i = 0; i < nbMinimLoc; i++){
            std::cout << "VALEUR DE x c [" << seuils[i] << "," << seuils[i+1] << "] : " << (seuils[i] + seuils[i+1])/2 << std::endl;
        }
        std::cout << "VALEUR DE x>" << seuils[nbMinimLoc-1] << " : " << 255 <<  std::endl;
        std::cout << "NB DE SEUILS : " << nbMinimLoc <<  std::endl;
    }

    for (int i=0; i < nH; i++){
        for (int j=0; j < nW; j++){
            bool changement = false;
            for(int k = 0; k<nbMinimLoc; k++){
                if ( ImgIn[i*nW+j] < seuils[k]){
                    ImgOut[i*nW+j] = valSeuils[k];
                    changement = true;
                    break;
                }
                if(changement){
                    break;
                }
            }
            if(!changement){
                ImgOut[i*nW+j] = 255;
            }
        }
    }
}
Traitement_seuillageMinLocaux::Traitement_seuillageMinLocaux(int ndgSimp): Traitement_image(){
    this->NDGSimp = ndgSimp;
}
Traitement_seuillageMinLocaux::~Traitement_seuillageMinLocaux(){}
void Traitement_seuillageMinLocaux::appliquer(){
    bool write_img = getTraitementEcriture();
    if(AFFICHAGE_DEBUG){
        std::cout << "SEUILLAGE MIN LOCAUX APPLIQUER" << std::endl;
        std::cout << "WRITE IMG : " << write_img << std::endl;
    }

    for (size_t i = 0; i < this->img_in.size(); i++){
        if(AFFICHAGE_DEBUG){
            std::cout << "IMG IN " << std::endl;
            std::cout << *img_in[i];
        }

        IMG * img_cur = this->img_in[i];
        OCTET *img_seuillee;
        int nTaille = img_cur->nH*img_cur->nW;
        allocation_tableau(img_seuillee, OCTET, nTaille);

        seuiller(img_cur->img_data,img_cur->nH,img_cur->nW,img_seuillee);

        char* nom_seuil = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_seuil,"%s_seuilleeMinLocaux%d", img_cur->nom,NDGSimp); 
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, img_seuillee, write_img,nom_seuil));
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_seuillageMinLocaux::afficher(std::ostream& os){
    os << "TRAITEMENT SEUILLAGE MIN LOCAUX : " << NDGSimp << std::endl;
    Traitement_image::afficher(os);
}

#endif