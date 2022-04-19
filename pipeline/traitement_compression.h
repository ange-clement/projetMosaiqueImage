#ifndef COMPRESSION_H
#define COMPRESSION_H

#include "pipeline_first.h"

#define FLG_ONDELETTE_4RES 0
#define FLG_ONDELETTE_RES 1

/**   ==================== ==================== ==================== ====================         
 *              ==================== ==================== ====================
 *    ==================== ==================== ==================== ====================**/
/**   ====================         
 *   TRAITEMENT DECORELATION        
 *    ====================**/
class Traitement_decorrelation : public Traitement_image{
    private :
        int difference_voisin_gauche_haut(int i, int j,int n, OCTET* imgbase);
        int DPCM(int i, int j,int n, OCTET* imgbase);
        int correction_erreur(int val);
    public :
        Traitement_decorrelation();
        virtual ~Traitement_decorrelation();
        void appliquer();
        virtual void afficher(std::ostream& os);
};

int Traitement_decorrelation::difference_voisin_gauche_haut(int i, int j,int n, OCTET* imgbase){
    return (j!=0)? imgbase[i*n+j] - imgbase[i*n+j-1] : (i!=0)? imgbase[i*n+j] - imgbase[(i-1)*n+j] :  imgbase[i*n+j];
}
int Traitement_decorrelation::DPCM(int i, int j,int n, OCTET* imgbase){
    if(j!=0){
        if(i!=0){
            int v_AB = abs(imgbase[i*n+j-1] - imgbase[(i-1)*n+j-1]);
            int v_BC = abs(imgbase[(i-1)*n+j-1] - imgbase[(i-1)*n+j]);

            return (v_AB<v_BC)? imgbase[i*n+j] - imgbase[(i-1)*n+j] : imgbase[i*n+j] - imgbase[i*n+j-1];
        }else{
            //En haut et pas à gauche
            return  imgbase[i*n+j] - imgbase[i*n+j-1];
        }
    }else{
        //A gauche
        return (i!=0)? imgbase[i*n+j] - imgbase[(i-1)*n+j] :  imgbase[i*n+j];
    }
}
int Traitement_decorrelation::correction_erreur(int val){
    val = (val<-128)? -128 : (val>=128)? 127 : val;
    return val + 128;
}
Traitement_decorrelation::Traitement_decorrelation(): Traitement_image(){}
Traitement_decorrelation::~Traitement_decorrelation(){}
void Traitement_decorrelation::appliquer(){
    bool write_img = getTraitementEcriture();
    if(AFFICHAGE_DEBUG){
        std::cout << "DECORRELATION APPLIQUER" << std::endl;
        std::cout << "WRITE IMG : " << write_img << std::endl;
    }

    for (size_t i = 0; i < this->img_in.size(); i++){
        if(AFFICHAGE_DEBUG){
            std::cout << "IMG IN " << std::endl;
            std::cout << *img_in[i];
        }

        IMG * img_cur = this->img_in[i];
        OCTET *img_decorrelation;
        int nTaille = img_cur->nH*img_cur->nW;
        allocation_tableau(img_decorrelation, OCTET, nTaille);

        for (int i=0; i < img_cur->nH; i++){
            for (int j=0; j < img_cur->nW; j++){
                img_decorrelation[i*img_cur->nW+j] = correction_erreur(DPCM(i,j,img_cur->nW,img_cur->img_data));
            }
        }

        char* nom_decorrelation = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_decorrelation,"%s_decorrelation", img_cur->nom); 
        this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, img_decorrelation, write_img,nom_decorrelation));
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_decorrelation::afficher(std::ostream& os){
    os << "TRAITEMENT DECORRELATION" << std::endl;
    Traitement_image::afficher(os);
}


/**   ====================         
 *   TRAITEMENT ONDELETTE        
 *    ====================**/
class Traitement_ondelette : public Traitement_image{
    private :
        int N;
        int flg_traitement;
        float BF(int i, int j, OCTET* img, int n);
        float MFH(int i, int j, OCTET* img, int n);
        float MFV(int i, int j, OCTET* img, int n);
        float HF(int i, int j, OCTET* img, int n);

        void appliquer_4res(IMG* img);
        void appliquer_res(OCTET* ImgIn, int nH, int nW,OCTET* ImgOut);

    public :
        Traitement_ondelette(int N, int flg = FLG_ONDELETTE_RES);
        virtual ~Traitement_ondelette();
        void appliquer();
        virtual void afficher(std::ostream& os);
};

float Traitement_ondelette::BF(int i, int j, OCTET* img, int n){
    return ((float) (img[i*n+j] + img[(i+1)*n+j] +img[i*n+j+1] +img[(i+1)*n+j+1])) / 4.0f;
}
float Traitement_ondelette::MFH(int i, int j, OCTET* img, int n){
    return ((float) (img[i*n+j] + img[(i+1)*n+j] -img[i*n+j+1] -img[(i+1)*n+j+1])) / 2.0f;
}
float Traitement_ondelette::MFV(int i, int j, OCTET* img, int n){
    return ((float) (img[i*n+j] - img[(i+1)*n+j] +img[i*n+j+1] -img[(i+1)*n+j+1])) / 2.0f;
}
float Traitement_ondelette::HF(int i, int j, OCTET* img, int n){
    return ((float) (img[i*n+j] - img[(i+1)*n+j] -img[i*n+j+1] +img[(i+1)*n+j+1]));
}
void Traitement_ondelette::appliquer_4res(IMG* img){
    OCTET* ImgIn, *ImgOut_BF,*ImgOut_MFH,*ImgOut_MFV,*ImgOut_HF;

    int nH = img->nH;
    int nW = img->nW;
    allocation_tableau(ImgIn, OCTET, nH*nW);

    for (size_t i = 0; i < nH*nW; i++){
        ImgIn[i] = img->img_data[i];
    }

    bool write_img = getTraitementEcriture();

    int nbdec = 0;
    int taille_H = nH;
    int taille_W = nW;

    while (nbdec < N){
        taille_H/=2;
        taille_W/=2;

        //On fait une decomposition par ondelettes : méthode de HAN
        int out_i = 0;
        if(AFFICHAGE_DEBUG){
            std::cout << "ON COMMENCE TRAITEMENT : " << nbdec << std::endl;
            std::cout << "TAILLE H : " << taille_H << std::endl;
            std::cout << "TAILLE W : " << taille_W << std::endl;
        }

        allocation_tableau(ImgOut_BF, OCTET, taille_W*taille_H);
        allocation_tableau(ImgOut_MFH, OCTET, taille_W*taille_H);
        allocation_tableau(ImgOut_MFV, OCTET, taille_W*taille_H);
        allocation_tableau(ImgOut_HF, OCTET, taille_W*taille_H);


        for (int i=0; i < taille_H*2; i+=2){
            int out_j = 0;
            for (int j=0; j < taille_W*2; j+=2){
                int bf = BF(i,j,ImgIn, taille_W*2);
                int mfh = MFH(i,j,ImgIn, taille_W*2) + 128;
                int mfv = MFV(i,j,ImgIn, taille_W*2) + 128;
                int hf = HF(i,j,ImgIn, taille_W*2) + 128;
                

                ImgOut_BF[out_i*taille_W+out_j] = (bf<0)? 0 : (bf>255)? 255 : bf;
                ImgOut_MFH[out_i*taille_W+out_j] = (mfh<0)? 0 : (mfh>255)? 255 : mfh;
                ImgOut_MFV[out_i*taille_W+out_j] = (mfv<0)? 0 : (mfv>255)? 255 : mfv;
                ImgOut_HF[out_i*taille_W+out_j] = (hf<0)? 0 : (hf>255)? 255 : hf;
                out_j++;
            }
            out_i++;
        }

        char* BF_char = new char[TAILLE_MAX_NAME_IMG];
        sprintf(BF_char, "%s_OndeletteBF%d",img->nom, N);
        this->img_out.push_back(new IMG(this->format_res,taille_H, taille_W, ImgOut_BF, write_img,BF_char));

        char* MFH_char = new char[TAILLE_MAX_NAME_IMG];
        sprintf(MFH_char, "%s_OndeletteMFH%d",img->nom, N);
        this->img_out.push_back(new IMG(this->format_res,taille_H, taille_W, ImgOut_MFH, write_img,MFH_char));

        char* MFV_char = new char[TAILLE_MAX_NAME_IMG];
        sprintf(MFV_char, "%s_OndeletteMFV%d",img->nom ,N);
        this->img_out.push_back(new IMG(this->format_res,taille_H, taille_W, ImgOut_MFV, write_img,MFV_char));

        char* HF_char = new char[TAILLE_MAX_NAME_IMG];
        sprintf(HF_char, "%s_OndeletteHF%d",img->nom, N);
        this->img_out.push_back(new IMG(this->format_res,taille_H, taille_W, ImgOut_HF, write_img,HF_char));

        nbdec++;

        for (int i = 0; i < taille_W*taille_H; i++){
            ImgIn[i] = ImgOut_BF[i];
        }
    }
}
void Traitement_ondelette::appliquer_res(OCTET* ImgIn, int nH, int nW,OCTET* ImgOut){
    int nbdec = 0;
    int taille_H = nH;
    int taille_W = nW;
    int nTaille = nH * nW;

    OCTET* buffer;
    allocation_tableau(buffer, OCTET, nH*nW);
    for (size_t i = 0; i < nH*nW; i++){
        buffer[i] = ImgIn[i];
    }

    while (nbdec < N){
        taille_H/=2;
        taille_W/=2;

        //On fait une decomposition par ondelettes : méthode de HAN
        int out_i = 0;
        if(AFFICHAGE_DEBUG){
            std::cout << "ON COMMENCE TRAITEMENT : " << nbdec << std::endl;
            std::cout << "TAILLE H : " << taille_H << std::endl;
            std::cout << "TAILLE W : " << taille_W << std::endl;
        }

        for (int i=0; i < taille_H*2; i+=2){
            int out_j = 0;
            for (int j=0; j < taille_W*2; j+=2){
                int bf = BF(i,j,buffer,nW);
                int mfh = MFH(i,j,buffer,nW) + 128;
                int mfv = MFV(i,j,buffer,nW) + 128;
                int hf = HF(i,j,buffer,nW) + 128;
                

                ImgOut[out_i*nW+out_j] = (bf<0)? 0 : (bf>255)? 255 : bf;
                ImgOut[(out_i+taille_H)*nW+out_j] = (mfh<0)? 0 : (mfh>255)? 255 : mfh;
                ImgOut[out_i*nW+out_j + taille_W] = (mfv<0)? 0 : (mfv>255)? 255 : mfv;
                ImgOut[(out_i+taille_H)*nW+out_j+taille_W] = (hf<0)? 0 : (hf>255)? 255 : hf;
                out_j++;
            }
            out_i++;
        }
        nbdec++;
        for (int i = 0; i < nTaille; i++){
            buffer[i] = ImgOut[i];
        }
    }
}
Traitement_ondelette::Traitement_ondelette(int N, int flg): Traitement_image(){
    this->N = N;
    this->flg_traitement = flg;
}
Traitement_ondelette::~Traitement_ondelette(){}
void Traitement_ondelette::appliquer(){
bool write_img = getTraitementEcriture();
    if(AFFICHAGE_DEBUG){
        std::cout << "ONDELETTE APPLIQUER" << std::endl;
        std::cout << "WRITE IMG : " << write_img << std::endl;
    }

    for (size_t i = 0; i < this->img_in.size(); i++){
        if(AFFICHAGE_DEBUG){
            std::cout << "IMG IN " << std::endl;
            std::cout << *img_in[i];
        }

        IMG * img_cur = this->img_in[i];
        
        if(this->flg_traitement == FLG_ONDELETTE_4RES){
            if(AFFICHAGE_DEBUG) std::cout << "ONDELETTE 4RES" << std::endl;
            appliquer_4res(img_cur);
        }else if(this->flg_traitement == FLG_ONDELETTE_RES){
            if(AFFICHAGE_DEBUG) std::cout << "ONDELETTE RES" << std::endl;
            OCTET *img_ondelette;
            int nTaille = img_cur->nH*img_cur->nW;
            allocation_tableau(img_ondelette, OCTET, nTaille);

            appliquer_res(img_cur->img_data,img_cur->nH,img_cur->nW,img_ondelette);

            char* nom_ondelette = new char[TAILLE_MAX_NAME_IMG];
            sprintf(nom_ondelette,"%s_Ondelette%d", img_cur->nom,N); 
            this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW, img_ondelette, write_img,nom_ondelette));
        }else{
            if(AFFICHAGE_DEBUG){
                std::cout << "AUCUN RESULTAT CAR FLG INCONNU" << std::endl;
            }
        }
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_ondelette::afficher(std::ostream& os){
    os << "TRAITEMENT ONDELETTE : " << N << std::endl;
    os << "RES : " << ((flg_traitement == FLG_ONDELETTE_4RES)? "4" : (flg_traitement == FLG_ONDELETTE_RES)? "1" : "NSP")<< std::endl;
    Traitement_image::afficher(os);
}


#endif