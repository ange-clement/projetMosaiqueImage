#ifndef HISTOGRAMME_H
#define HISTOGRAMME_H

#include "pipeline_first.h"


/**   ==================== ==================== ==================== ====================         
 *              ==================== ==================== ====================
 *    ==================== ==================== ==================== ====================**/
/**   ====================         
 *TRAITEMENT EGALISATION HISTO        
 *    ====================**/
class Traitement_egalisationHistogramme : public Traitement_image{
    private :
    void egalisation(int nW, int nH,const OCTET * ImgIn, OCTET * ImgOut);
    public :
        Traitement_egalisationHistogramme();
        virtual ~Traitement_egalisationHistogramme();
        void appliquer();
        virtual void afficher(std::ostream& os);
};

void Traitement_egalisationHistogramme::egalisation(int nW, int nH, const OCTET * ImgIn, OCTET * ImgOut){
    int nTaille = nW * nH;
    int* nbElementsLus = new int[256]{};
    for (int i=0; i < nH; i++){
        for (int j=0; j < nW; j++){
            nbElementsLus[ImgIn[i*nW+j]]++;
        }
    }
    float* ddp = new float[256]{};
    for(int i = 0; i<256; i++){
        ddp[i] = (float) nbElementsLus[i]/(float) nTaille;
    }
    float* Frep = new float[256]{};
    Frep[0] = ddp[0];
    for(int i = 1; i<256; i++){
        Frep[i] = Frep[i-1] + ddp[i];
    }
    for (int i = 0; i < nTaille; i++){
        int val = 255 * Frep[ImgIn[i]];
        ImgOut[i] = (val<0)? 0 : (val>255)? 255 : val;
    }

    delete nbElementsLus;
    delete ddp;
    delete Frep;
}
Traitement_egalisationHistogramme::Traitement_egalisationHistogramme() : Traitement_image(){}
Traitement_egalisationHistogramme::~Traitement_egalisationHistogramme(){}
void Traitement_egalisationHistogramme::appliquer(){
     bool write_img = getTraitementEcriture();
    if(AFFICHAGE_DEBUG){
        std::cout << "EGALISATION APPLIQUER" << std::endl;
        std::cout << "WRITE IMG : " << write_img << std::endl;
    }

    for (size_t k = 0; k < this->img_in.size(); k++){
        if(AFFICHAGE_DEBUG){
            std::cout << "IMG IN " << std::endl;
            std::cout << *img_in[k];
        }

        IMG * img_cur = this->img_in[k];
        OCTET *img_egal;
        int nTaille = img_cur->nH*img_cur->nW;
        allocation_tableau(img_egal, OCTET, nTaille);

        egalisation(img_cur->nW, img_cur->nH,img_cur->img_data, img_egal);

        char* nom_spe = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_spe,"%s_egalisationHISTO", img_cur->nom); 
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, img_egal, write_img,nom_spe));
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_egalisationHistogramme::afficher(std::ostream& os){
    os << "TRAITEMENT EGALISATION HISTOGRAMME" << std::endl;
    Traitement_image::afficher(os);
}


/**   ====================         
 *TRAITEMENT EXPENSION DYNAMIQUE        
 *    ====================**/
class Traitement_expensionDynamique : public Traitement_image{
    private :
    void expension(int nW, int nH,const OCTET * ImgIn, OCTET * ImgOut);
    public :
        Traitement_expensionDynamique();
        virtual ~Traitement_expensionDynamique();
        void appliquer();
        virtual void afficher(std::ostream& os);
};

void Traitement_expensionDynamique::expension(int nW, int nH, const OCTET * ImgIn, OCTET * ImgOut){
    int nTaille = nW * nH;
    int* nbElementsLus = new int[256]{};
    for (int i=0; i < nH; i++){
        for (int j=0; j < nW; j++){
            nbElementsLus[ImgIn[i*nW+j]]++;
        }
    }

    int min = 256;
    int max = 0;

    for (int i = 0; i < 256; i++){
        if(nbElementsLus[i] != 0 && min>i){
            min = i;
        }
        if(nbElementsLus[i] !=0 && max<i){
            max = i;
        }
    }
    max++;

    float alpha = 255.0f * (float) min/((float) max-min);
    float beta = 255.0f/((float) max-min);

    for (int i=0; i < nTaille; i++){
         int val = alpha + beta *  ImgIn[i];
         ImgOut[i] = (val<0)? 0 : (val>255)? 255 : val;
    }
}
Traitement_expensionDynamique::Traitement_expensionDynamique() : Traitement_image(){}
Traitement_expensionDynamique::~Traitement_expensionDynamique(){}
void Traitement_expensionDynamique::appliquer(){
    bool write_img = getTraitementEcriture();
    if(AFFICHAGE_DEBUG){
        std::cout << "EXPENSION DYNAMIQUE APPLIQUER" << std::endl;
        std::cout << "WRITE IMG : " << write_img << std::endl;
    }

    for (size_t k = 0; k < this->img_in.size(); k++){
        if(AFFICHAGE_DEBUG){
            std::cout << "IMG IN " << std::endl;
            std::cout << *img_in[k];
        }

        IMG * img_cur = this->img_in[k];
        OCTET *img_exp;
        int nTaille = img_cur->nH*img_cur->nW;
        allocation_tableau(img_exp, OCTET, nTaille);

        expension(img_cur->nW, img_cur->nH,img_cur->img_data, img_exp);

        char* nom_spe = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_spe,"%s_expensionDynamique", img_cur->nom); 
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, img_exp, write_img,nom_spe));
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_expensionDynamique::afficher(std::ostream& os){
    os << "TRAITEMENT EXPENSION DYNAMIQUE" << std::endl;
    Traitement_image::afficher(os);
}


/**    ====================         
 *TRAITEMENT SPECIFICATION HISTOGRAMME        
 *     ====================**/
class Traitement_specificationHistogramme : public Traitement_image{
    private :
    IMG * img_specification;
    void specification(int nW, int nH,const OCTET * ImgIn, OCTET * ImgOut);
    public :
        Traitement_specificationHistogramme(IMG * img_specification);
        virtual ~Traitement_specificationHistogramme();
        void appliquer();
        virtual void afficher(std::ostream& os);
};

void Traitement_specificationHistogramme::specification(int nW, int nH,const OCTET * ImgIn, OCTET * ImgOut){
    OCTET *ImgIn_Beg;
    int nTaille = nW * nH;
    allocation_tableau(ImgIn_Beg, OCTET, img_specification->nW*img_specification->nH);

    int* nbElementsLus = new int[256]{};
    for (int i=0; i < img_specification->nH; i++){
        for (int j=0; j < img_specification->nW; j++){
            nbElementsLus[img_specification->img_data[i*img_specification->nW+j]]++;
        }
    }

    float* ddp = new float[256]{};
    for(int i = 0; i<256; i++){
        ddp[i] = (float) nbElementsLus[i]/(float) nTaille;
    }

    float* Frep = new float[256]{};
    Frep[0] = ddp[0];

    for(int i = 1; i<256; i++){
        Frep[i] = Frep[i-1] + ddp[i];
    }

    for (int i = 0; i < img_specification->nH*img_specification->nW; i++){
        int val = 255 * Frep[img_specification->img_data[i]];
        ImgIn_Beg[i] = (val<0)? 0 : (val>255)? 255 : val;
    }


    int* nbElementsLus_R = new int[256]{};
    for (int i=0; i < nH; i++){
        for (int j=0; j < nW; j++){
            nbElementsLus_R[ImgIn[i*nW+j]]++;
        }
    }

    float* ddp_B = new float[256]{};
    for(int i = 0; i<256; i++){
        ddp_B[i] = (float) nbElementsLus_R[i]/(float) nTaille;
    }

    float* Frep_B = new float[256]{};
    Frep_B[0] = ddp_B[0];

    for(int i = 1; i<256; i++){
        Frep_B[i] = Frep_B[i-1] + ddp_B[i];
    }

    for (int i = 0; i < img_specification->nH*img_specification->nW; i++){
        float valeur = ImgIn_Beg[i]/255.0f;
        
        float v = 0.0f;
        int j = 0;
        for (; j < 256 && v<valeur; j++){
            v=Frep_B[j];
        }
        ImgOut[i] = j;
    }
}
Traitement_specificationHistogramme::Traitement_specificationHistogramme(IMG * img_specification) : Traitement_image(){
    this->img_specification = img_specification;
}
Traitement_specificationHistogramme::~Traitement_specificationHistogramme(){
    delete img_specification;
}
void Traitement_specificationHistogramme::appliquer(){
    bool write_img = getTraitementEcriture();
    if(AFFICHAGE_DEBUG){
        std::cout << "SPECIFICATION APPLIQUER" << std::endl;
        std::cout << "WRITE IMG : " << write_img << std::endl;
    }

    for (size_t k = 0; k < this->img_in.size(); k++){
        if(AFFICHAGE_DEBUG){
            std::cout << "IMG IN " << std::endl;
            std::cout << *img_in[k];
        }

        IMG * img_cur = this->img_in[k];
        OCTET *img_spe;
        int nTaille = img_cur->nH*img_cur->nW;
        allocation_tableau(img_spe, OCTET, nTaille);

        specification(img_cur->nW, img_cur->nH,img_cur->img_data, img_spe);

        char* nom_spe = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_spe,"%s_specificationHISTO", img_cur->nom); 
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, img_spe, write_img,nom_spe));
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_specificationHistogramme::afficher(std::ostream& os){
    os << "TRAITEMENT SPECIFICATION HISTOGRAMME" << std::endl;
    os << *img_specification << std::endl;
    Traitement_image::afficher(os);
}


#endif