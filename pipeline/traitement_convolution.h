#ifndef TRAITEMENTCONVO_H
#define TRAITEMENTCONVO_H

#include "pipeline_first.h"

typedef struct masque masque;

#define TAILLE_MAX_NOM_MASK 100

/**   ==================== ==================== ==================== ====================         
 *              ==================== ==================== ====================
 *    ==================== ==================== ==================== ====================**/
/**   ====================         
 *      TRAITEMENT CONVO        
 *    ====================**/
struct masque{
    float* masque;
    int nb_x;
    int nb_y;
    const char* name;

    float gaussian(float x, float mu,float sigma){
        const float a = (x-mu)/sigma;
        return std::exp(-.5f * a * a); 
    }

    void onDestroy(){
        delete[] masque;
        delete[] name;
    }

    struct::masque gaussian(int k, float sigma){
        onDestroy();
        int taille = k*2+1;
        nb_x = taille;
        nb_y = taille;

        masque = new float[taille*taille];
        float sum = 0.0f;
        for (int i = 0; i < taille; i++){
            for (int j = 0; j < taille; j++){
                float elem = gaussian(i,k,sigma) * gaussian(j,k,sigma);
                masque[i* taille+j] = elem;
                sum+=elem;
            }
        }
        for (int i = 0; i < taille; i++){
            for (int j = 0; j < taille; j++){
                masque[i* taille+j]/=sum;
            }
        }
        
        char* nom = new char[TAILLE_MAX_NOM_MASK];
        sprintf(nom, "Gaussian%d|%f", k, sigma);
        this->name = nom;

        return *this;
    }

    struct::masque gaussian(int k){
        onDestroy();
        int taille = k*2+1;
        float sigma = (float) k /2.0f;
        nb_x = taille;
        nb_y = taille;

        masque = new float[taille*taille];
        float sum = 0.0f;
        for (int i = 0; i < taille; i++){
            for (int j = 0; j < taille; j++){
                float elem = gaussian(i,k,sigma) * gaussian(j,k,sigma);
                masque[i* taille+j] = elem;
                sum+=elem;
            }
        }
        for (int i = 0; i < taille; i++){
            for (int j = 0; j < taille; j++){
                masque[i* taille+j]/=sum;
            }
        }

        char* nom = new char[TAILLE_MAX_NOM_MASK];
        sprintf(nom, "Gaussian%d|%f", k, sigma);
        this->name = nom;

        return *this;
    }

    struct::masque moyenneur(int taille){
        onDestroy();
        masque = new float[taille*taille];
        for (size_t i = 0; i < taille*taille; i++){
            masque[i] = 1.0f/(float)(taille*taille);
        }
        nb_x = taille;
        nb_y = taille;

        char* nom = new char[TAILLE_MAX_NOM_MASK];
        sprintf(nom, "Moyenneur%d", taille);
        this->name = nom;

        return *this;
    }

    struct::masque binomial_newton_ordre(){
        onDestroy();
        masque = new float[25];
        masque[0] = 1.0f/256.0f;
        masque[1] = 4.0f/256.0f;
        masque[2] = 6.0f/256.0f;
        masque[3] = 4.0f/256.0f;
        masque[4] = 1.0f/256.0f;

        masque[5] = 4.0f/256.0f;
        masque[6] = 16.0f/256.0f;
        masque[7] = 24.0f/256.0f;
        masque[8] = 16.0f/256.0f;
        masque[9] = 4.0f/256.0f;

        masque[10] = 6.0f/256.0f;
        masque[11] = 24.0f/256.0f;
        masque[12] = 36.0f/256.0f;
        masque[13] = 24.0f/256.0f;
        masque[14] = 6.0f/256.0f;

        masque[15] = 4.0f/256.0f;
        masque[16] = 16.0f/256.0f;
        masque[17] = 24.0f/256.0f;
        masque[18] = 16.0f/256.0f;
        masque[19] = 4.0f/256.0f;

        masque[20] = 1.0f/256.0f;
        masque[21] = 4.0f/256.0f;
        masque[22] = 6.0f/256.0f;
        masque[23] = 4.0f/256.0f;
        masque[24] = 1.0f/256.0f;

        nb_x = 5;
        nb_y = 5;
        
        char* nom = new char[TAILLE_MAX_NOM_MASK];
        sprintf(nom, "BinoNewton%dx%d", nb_x, nb_y);
        this->name = nom;

        return *this;
    }

    struct::masque sobel_prewit_y(bool reverse, float val){
        onDestroy();
        masque = new float[9];
        if(reverse){
            masque[0] = 1.0f;
            masque[1] = val;
            masque[2] = 1.0f;
        }else{
            masque[0] = -1.0f;
            masque[1] = -val;
            masque[2] = -1.0f;
        }

        masque[3] = 0;
        masque[4] = 0;
        masque[5] = 0;
        if(reverse){
            masque[6] = -1.0f;
            masque[7] = -val;
            masque[8] = -1.0f;
        }else{
            masque[6] = 1.0f;
            masque[7] = val;
            masque[8] = 1.0f;
        }

        nb_x = 3;
        nb_y= 3;

        char* nom = new char[TAILLE_MAX_NOM_MASK];
        sprintf(nom, "SobelPrewitY%f%s", val, ((reverse)? "REVERSED" : "NOTREVERSED"));
        this->name = nom;

        return *this;
    }

    struct::masque sobel_prewit_x(bool reverse, float val){
        onDestroy();
        masque = new float[9];
        if(reverse){
            masque[0] = 1.0f;
            masque[1] = 0;
            masque[2] = -1.0f;

            masque[3] = val;
            masque[4] = 0;
            masque[5] = -val;
    
            masque[6] = 1.0f;
            masque[7] = 0;
            masque[8] = -1.0f;
        }else{
            masque[0] = -1.0f;
            masque[1] = 0;
            masque[2] = 1.0f;

            masque[3] = -val;
            masque[4] = 0;
            masque[5] = val;
    
            masque[6] = -1.0f;
            masque[7] = 0;
            masque[8] = 1.0f;
        }

        nb_x = 3;
        nb_y = 3;

        char* nom = new char[TAILLE_MAX_NOM_MASK];
        sprintf(nom, "SobelPrewitX%f%s", val, ((reverse)? "REVERSED" : "NOTREVERSED"));
        this->name = nom;

        return *this;
    }

    struct::masque laplacien_4connexes(){
        onDestroy();
        masque = new float[9];
        masque[0] = 0.0f;
        masque[1] = 1.0f;
        masque[2] = 0.0f;

        masque[3] = 1.0f;
        masque[4] = -4.0f;
        masque[5] = 1.0f;
    
        masque[6] = 0.0f;
        masque[7] = 1.0f;
        masque[8] = 0.0f;

        nb_x = 3;
        nb_y= 3;

        this->name = "Laplacien4CO";

        return *this;
    }

    struct::masque laplacien_8connexes(){
        onDestroy();
        masque = new float[9];
        masque[0] = 1.0f;
        masque[1] = 1.0f;
        masque[2] = 1.0f;

        masque[3] = 1.0f;
        masque[4] = -8.0f;
        masque[5] = 1.0f;
    
        masque[6] = 1.0f;
        masque[7] = 1.0f;
        masque[8] = 1.0f;

        nb_x = 3;
        nb_y= 3;

        this->name = "Laplacien8CO";

        return *this;
    }

    struct::masque kirsh(){
        onDestroy();
        masque = new float[9];
        masque[0] = 5.0f;
        masque[1] = 5.0f;
        masque[2] = 5.0f;

        masque[3] = -3.0f;
        masque[4] = 0.0f;
        masque[5] = -3.0f;
    
        masque[6] = -3.0f;
        masque[7] = -3.0f;
        masque[8] = -3.0f;

        nb_x = 3;
        nb_y= 3;

        this->name = "Kirsh3x3";

        return *this;
    }

    struct::masque robinson(){
        onDestroy();
        masque = new float[9];
        masque[0] = 1.0f;
        masque[1] = 1.0f;
        masque[2] = 1.0f;

        masque[3] = 1.0f;
        masque[4] = -2.0f;
        masque[5] = 1.0f;
    
        masque[6] = -1.0f;
        masque[7] = -1.0f;
        masque[8] = -1.0f;

        nb_x = 3;
        nb_y= 3;

        this->name = "Robinson3x3";

        return *this;
    }

    struct::masque robinson_laplacien(){
        onDestroy();
        masque = new float[9];
        masque[0] = 1.0f;
        masque[1] = -2.0f;
        masque[2] = 1.0f;

        masque[3] = -2.0f;
        masque[4] = 4.0f;
        masque[5] = -2.0f;
    
        masque[6] = 1.0f;
        masque[7] = -2.0f;
        masque[8] = 1.0f;

        nb_x = 3;
        nb_y= 3;

        this->name = "RobinsonLaplacien3x3";

        return *this;
    }

    void afficher(std::ostream& os){
        os << "*** MASQUE ***" << std::endl;
        for (size_t i = 0; i < nb_y; i++){
            for (size_t j = 0; j < nb_x; j++){
                os << "["<<masque[j+i*nb_x]<<"]";
            }
            os << std::endl;
        }
    }
};

//Traitement de l'image normale (pas de traitement quand on atteind les bords)
#define FLG_NORMAL 0
#define FLG_MIROR 1

//Pour avoir image non convoluée et image convoluée
#define FLG_MIXING 2
#define FLG_NO_MIXING 3

class Traitement_convo : public Traitement_image{
    private :
        masque masque_convo;
        int flag = FLG_NORMAL;
        int flg_mix = FLG_NO_MIXING;

        void convolution(int nH, int nW,int taille_x, int taille_y, OCTET *ImgIn,OCTET *ImgRes,float* filtre){
            if(AFFICHAGE_DEBUG) {
                std::cout << "TRAITEMENT CONVOLUTION" << std::endl;
                std::cout << "HEIGHT : " << nH << std::endl;
                std::cout << "WIDTH : " << nW << std::endl;
                std::cout << "Taille X : " << taille_x << std::endl;
                std::cout << "Taille Y : " << taille_y << std::endl;
            }

            int ecart_x = (taille_x-1)/2;
            int ecart_y = (taille_y-1)/2;

            for (int i = ecart_y; i < nH-(taille_y-ecart_y); i++){
                for (int j = ecart_x; j < nW - (taille_x-ecart_x); j++){
                    float valpx = 0;
                    for (int k = 0; k<taille_y; k++){
                        for (int l = 0; l<taille_x; l++){
                            valpx += filtre[k*taille_x+l] * ImgIn[(i+k-ecart_y)*nW + (j+l-ecart_x)];
                        }   
                    }
                    ImgRes[i*nW+j] = valpx;
                }
            }

            if(AFFICHAGE_DEBUG) {
                std::cout << "FIN TRAITEMENT CONVOLUTION" << std::endl;
            }
        }
        void convolution_miroir(int nH, int nW,int taille_x, int taille_y, OCTET *ImgIn,OCTET *ImgRes,float* filtre){
            if(AFFICHAGE_DEBUG) {
                std::cout << "TRAITEMENT CONVOLUTION MIRROIR" << std::endl;
                std::cout << "HEIGHT : " << nH << std::endl;
                std::cout << "WIDTH : " << nW << std::endl;
                std::cout << "Taille X : " << taille_x << std::endl;
                std::cout << "Taille Y : " << taille_y << std::endl;
            }

            int ecart_x = (taille_x-1)/2;
            int ecart_y = (taille_y-1)/2;

            for (int i = 0; i < nH; i++){
                for (int j = 0; j < nW; j++){
                    float valpx = 0;
                    for (int k = -ecart_y; k <= ecart_y; k++){
                        for (int l = -ecart_x; l <= ecart_x; l++){
                            int lin = abs(i+k);
                            int col = abs(j+l);
                            lin = (lin<nH)? lin : nH-1;
                            col = (col<nW)? col : nW-1;
                            valpx += filtre[(ecart_y+k)*taille_x+ecart_x+l] * ImgIn[lin*nW+col];
                        }
                    }
                    ImgRes[i*nW+j] = valpx;
                }
            }
        }

    public :
        Traitement_convo(masque masqueconv);
        virtual ~Traitement_convo();
        void appliquer();

        Traitement_convo* set_flg(int flg);
        Traitement_convo* set_flg_mixing(int flg);

        virtual void afficher(std::ostream& os);
        //friend std::ostream& operator<<(std::ostream& os, const Traitement_convo& elem);  
};

Traitement_convo::Traitement_convo(masque masqueconv) : Traitement_image(){
    masque_convo = masqueconv;
}
Traitement_convo::~Traitement_convo(){
    this->masque_convo.onDestroy();
}
void Traitement_convo::appliquer(){
    bool write_img = getTraitementEcriture();
    switch (flag){
    case FLG_NORMAL:
        if(AFFICHAGE_DEBUG) std::cout << "TRAITEMENT CONVO APPLIQUER : NORMAL" << std::endl;
        for (size_t i = 0; i < img_in.size(); i++){
            OCTET* res;
            allocation_tableau(res,OCTET,img_in[i]->nH*img_in[i]->nW);
            convolution(img_in[i]->nH,img_in[i]->nW,masque_convo.nb_x,masque_convo.nb_y,img_in[i]->img_data,res,masque_convo.masque);
            char* nom = new char[TAILLE_MAX_NAME_IMG];
            
            if(this->flg_mix == FLG_MIXING){
                if(AFFICHAGE_DEBUG) std::cout << "MIXING" << std::endl;
                sprintf(nom,"%s_convolutionMIX(%s)", img_in[i]->nom,masque_convo.name);
                OCTET* res_mix;
                allocation_tableau(res_mix,OCTET,img_in[i]->nH*img_in[i]->nW*2);
                
                for (int j = 0; j < img_in[i]->nH; j++){
                    for (int k = 0; k < img_in[i]->nW; k++){
                        res_mix[j*img_in[i]->nW*2 + k] = img_in[i]->img_data[j*img_in[i]->nW + k];
                    }
                    for (int k = 0; k < img_in[i]->nW; k++){
                        res_mix[j*img_in[i]->nW*2 + k+img_in[i]->nW] = res[j*img_in[i]->nW + k];
                    }
                }
                img_out.push_back(new IMG(this->format_res,img_in[i]->nH, img_in[i]->nW*2, res_mix, write_img,nom));
            }else{
                if(AFFICHAGE_DEBUG) std::cout << "NO MIXING" << std::endl;
                sprintf(nom,"%s_convolution(%s)", img_in[i]->nom,masque_convo.name);
                img_out.push_back(new IMG(this->format_res,img_in[i]->nH, img_in[i]->nW, res, write_img,nom));
            }
        }
        break;
    case FLG_MIROR :
        if(AFFICHAGE_DEBUG) std::cout << "TRAITEMENT CONVO APPLIQUER : MIRROR" << std::endl;
        for (size_t i = 0; i < img_in.size(); i++){
            OCTET* res;
            allocation_tableau(res,OCTET,img_in[i]->nH*img_in[i]->nW);
            convolution_miroir(img_in[i]->nH,img_in[i]->nW,masque_convo.nb_x,masque_convo.nb_y,img_in[i]->img_data,res,masque_convo.masque);
            char* nom = new char[TAILLE_MAX_NAME_IMG];
            
            if(this->flg_mix == FLG_MIXING){
                if(AFFICHAGE_DEBUG) std::cout << "MIXING" << std::endl;
                sprintf(nom,"%s_convolution_mirorMIX(%s)", img_in[i]->nom,masque_convo.name);
                OCTET* res_mix;
                allocation_tableau(res_mix,OCTET,img_in[i]->nH*img_in[i]->nW*2);
                for (int j = 0; j < img_in[i]->nH; j++){
                    for (int k = 0; k < img_in[i]->nW; k++){
                        res_mix[j*img_in[i]->nW*2 + k] = img_in[i]->img_data[j*img_in[i]->nW + k];
                    }
                    for (int k = img_in[i]->nW; k < img_in[i]->nW*2; k++){
                        res_mix[j*img_in[i]->nW*2 + k] = res[j*img_in[i]->nW + k];
                    }
                }
                img_out.push_back(new IMG(this->format_res,img_in[i]->nH, img_in[i]->nW*2, res_mix, write_img,nom));
            }else{
                if(AFFICHAGE_DEBUG) std::cout << "NO MIXING" << std::endl;
                sprintf(nom,"%s_convolution_miror(%s)", img_in[i]->nom,masque_convo.name);
                img_out.push_back(new IMG(this->format_res,img_in[i]->nH, img_in[i]->nW, res, write_img, nom));
            }
        }
        break;
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_convo::afficher(std::ostream& os){
    os << "*** TRAITEMENT CONVOLUTIF ***" << std::endl;
    Traitement_image::afficher(os);
    os << "=> MASQUE : " << std::endl;
    masque_convo.afficher(os);
    os << "=> FLAG : " << flag << std::endl;
}
Traitement_convo* Traitement_convo::set_flg(int flg){
    this->flag = flg;
    return this;
}
Traitement_convo* Traitement_convo::set_flg_mixing(int flg){
    this->flg_mix = flg;
    return this;
}


/**   ====================         
 *     TRAITEMENT MEDIANE        
 *    ====================**/
class Traitement_mediane : public Traitement_image{
    private :
        int flag = FLG_NORMAL;
        int taille_x;
        int taille_y;
        
        void mediane(int nH, int nW, OCTET *ImgIn,OCTET *ImgRes){
            if(AFFICHAGE_DEBUG) {
                std::cout << "TRAITEMENT MEDIANE" << std::endl;
                std::cout << "HEIGHT : " << nH << std::endl;
                std::cout << "WIDTH : " << nW << std::endl;
            }
            std::vector<int> valeurs_px;
            int ecart_x = (taille_x-1)/2;
            int ecart_y = (taille_y-1)/2;

            for (int i = ecart_y; i < nH- (taille_y-ecart_y); i++){
                for (int j = ecart_x; j < nW - (taille_x-ecart_x); j++){
                    for (int k = 0; k<taille_y; k++){
                        for (int l = 0; l<taille_x; l++){
                            valeurs_px.push_back(ImgIn[(i+k-ecart_y)*nW + (j+l-ecart_x)]);
                        }   
                    }
                    std::sort(valeurs_px.begin(), valeurs_px.end());
                    ImgRes[i*nW+j] = valeurs_px[(taille_x*taille_y)/2];
                    valeurs_px.clear();
                }
            }
        }
        void mediane_miroir(int nH, int nW, OCTET *ImgIn,OCTET *ImgRes){
            if(AFFICHAGE_DEBUG) {
                std::cout << "TRAITEMENT MEDIANE MIRROIR" << std::endl;
                std::cout << "HEIGHT : " << nH << std::endl;
                std::cout << "WIDTH : " << nW << std::endl;
            }
            std::vector<int> valeurs_px;
            int ecart_x = (taille_x-1)/2;
            int ecart_y = (taille_y-1)/2;

            for (int i = 0; i < nH; i++){
                for (int j = 0; j < nW; j++){
                    float valpx = 0;
                    for (int k = -ecart_y; k <= ecart_y; k++){
                        for (int l = -ecart_x; l <= ecart_x; l++){
                    
                            int lin = abs(i+k);
                            int col = abs(j+l);
                            lin = (lin<nH)? lin : nH-1;
                            col = (col<nW)? col : nW-1;
                            valeurs_px.push_back(ImgIn[lin*nW+col]);
                        }
                    }
                    std::sort(valeurs_px.begin(), valeurs_px.end());
                    ImgRes[i*nW+j] = valeurs_px[(taille_x*taille_y)/2];
                    valeurs_px.clear();
                }
            }
        }

    public :
        Traitement_mediane(int t_x, int t_y);
        virtual ~Traitement_mediane();
        void appliquer();
        virtual void afficher(std::ostream& os);
        Traitement_mediane* set_flg(int flg);
        
        //friend std::ostream& operator<<(std::ostream& os, const Traitement_mediane& elem);  
};

Traitement_mediane::Traitement_mediane(int t_x, int t_y): Traitement_image(){
    taille_x = t_x;
    taille_y = t_y;
}
Traitement_mediane::~Traitement_mediane(){}
void Traitement_mediane::appliquer(){
    bool write_img = getTraitementEcriture();
    switch (flag){
    case FLG_NORMAL:
        if(AFFICHAGE_DEBUG) std::cout << "TRAITEMENT MEDIANE APPLIQUER : NORMAL" << std::endl;
        for (size_t i = 0; i < img_in.size(); i++){
            OCTET* res;
            allocation_tableau(res,OCTET,img_in[i]->nH*img_in[i]->nW);
            mediane(img_in[i]->nH,img_in[i]->nW,img_in[i]->img_data,res);
            char* nom = new char[TAILLE_MAX_NAME_IMG];
            sprintf(nom,"%s_mediane%dx%d", img_in[i]->nom, this->taille_x,this->taille_y);
            img_out.push_back(new IMG(this->format_res,img_in[i]->nH, img_in[i]->nW, res, write_img,nom));
        }
        break;
    case FLG_MIROR :
        if(AFFICHAGE_DEBUG) std::cout << "TRAITEMENT MEDIANE APPLIQUER : MIRROR" << std::endl;
        for (size_t i = 0; i < img_in.size(); i++){
            OCTET* res;
            allocation_tableau(res,OCTET,img_in[i]->nH*img_in[i]->nW);
            mediane_miroir(img_in[i]->nH,img_in[i]->nW,img_in[i]->img_data,res);
            char* nom = new char[TAILLE_MAX_NAME_IMG];
            sprintf(nom,"%s_mediane_miror%dx%d", img_in[i]->nom, this->taille_x,this->taille_y);
            img_out.push_back(new IMG(this->format_res,img_in[i]->nH, img_in[i]->nW, res, write_img, nom));
        }
        break;
    }

    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_mediane::afficher(std::ostream& os){
    os << "*** TRAITEMENT MEDIANE ***" << std::endl;
    Traitement_image::afficher(os);
    os << "=> TAILLE X : " << taille_x<< std::endl;
    os << "=> TAILLE Y : " << taille_y<< std::endl;
    os << "=> FLAG : " << flag << std::endl;
}
Traitement_mediane* Traitement_mediane::set_flg(int flg){
    this->flag = flg;
    return this;
}

#endif