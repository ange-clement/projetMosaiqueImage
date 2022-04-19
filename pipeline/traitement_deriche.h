#ifndef DERICHE_H
#define DERICHE_H

#include "pipeline_first.h"

#define FLG_ZEROPASSING_4 0
#define FLG_ZEROPASSING_8 1

/**   ==================== ==================== ==================== ====================         
 *              ==================== ==================== ====================
 *    ==================== ==================== ==================== ====================**/
/**   ====================         
 *    TRAITEMENT DERICHE        
 *    ====================**/
class Traitement_deriche : public Traitement_image{
    private :

    float alpha;

    void horizontal_detection(const OCTET* img, float alpha, int nH, int nW,OCTET* leftToRight,OCTET* rightToLeft, OCTET* res);
    void vertical_detection(const OCTET* img, float alpha, int nH, int nW,OCTET* leftToRight,OCTET* rightToLeft, OCTET* res);


    public :
        Traitement_deriche(float alpha);
        virtual ~Traitement_deriche();
        void appliquer();
        virtual void afficher(std::ostream& os);
};

void Traitement_deriche::horizontal_detection(const OCTET* img, float a, int nH, int nW,OCTET* leftToRight,OCTET* rightToLeft, OCTET* res){
    float k = ((1.0f - std::exp(-a))*(1.0f - std::exp(-a))) / (1.0f+2.0f*a*std::exp(-a) - std::exp(-2.0f*a));
    float a1 = k;
    float a2 = k * std::exp(-a)*(a-1.0f);
    float a3 = k * std::exp(-a)*(a+1.0f);
    float a4 = -k*std::exp(- 2.0f * a);

    float b1 = 2.0f * std::exp(-a);
    float b2 = - std::exp(-a * 2.0f);

    float s1, s2, s;

    for (int j = 0; j < nH; j++){
        s2 = img[(nW - 1) + j * nW];
        s1 = s2;

        for (int i = 0; i < 2; i++){
            s = s1;
            s1 = a3 * img[nW-1 + j*nW] + a4 * img[nW-1 + j*nW]  + b1 * s1 + b2 * s2;
            s2 = s;        
        }
        for (int i = nW-3; i >= 0; i--){
            s = s1;
            s1 = a3 * img[i+(j+1)*nW] + a4 * img[i+(nH+2)*nW]  + b1 * s1 + b2 * s2;
            s2 = s; 
        }

        rightToLeft[0+ j * nW] = s1;
        rightToLeft[1+ j * nW] = s2;

        for (int i = 2; i < nW; i++){
            rightToLeft[i+j*nW] = a1*img[i+j*nW] + a2 * img[i-1+j*nW] + b1 * rightToLeft[i-1 + j * nW] + b2 * rightToLeft[i-2 + j * nW];
            s = s1;
            s1 =  a1*img[i+j*nW] + a2 * img[i-1+j*nW] + b1 * s1 + b2 * s2;
            s2 = s;
        }

        leftToRight[nW -1 +j*nW] = s1;
        leftToRight[nW -2 +j*nW] = s2;

        for (int i = nW - 3; i >= 0; i--){
            leftToRight[i+j*nW] = a3 * img[i+1+j*nW] + a4 * img[i+2+j*nW] + b1 * leftToRight[i+1+j*nW] + b2 * leftToRight[i+2+j*nW];
        }
    }
    for (int i = 0; i < nW*nH; i++){
        res[i] = rightToLeft[i] + leftToRight[i];    
    }
}
void Traitement_deriche::vertical_detection(const OCTET* img, float a, int nH, int nW,OCTET* upToBottom,OCTET* bottomToUp, OCTET* res){
    float k = ((1.0f - std::exp(-a))*(1.0f - std::exp(-a))) / (1.0f+2.0f*a*std::exp(-a) - std::exp(-2.0f*a));
    float a1 = k;
    float a2 = k * std::exp(-a)*(a-1.0f);
    float a3 = k * std::exp(-a)*(a+1.0f);
    float a4 = -k*std::exp(- 2.0f * a);

    float b1 = 2.0f * std::exp(-a);
    float b2 = - std::exp(-a * 2.0f);

    float s1, s2, s;

    for (int i = 0; i < nW; i++){
        s2 = img[i + (nH - 1)* nW];
        s1 = s2;

        for (int j = 0; j< 2; j++){
            s = s1;
            s1 = a3 * img[i + (nH -1)*nW] + a4 * img[i + (nH -1)*nW]  + b1 * s1 + b2 * s2;
            s2 = s; 
        }

        for (int j = nH-3; j >= 0; j--){
            s = s1;
            s1 = a3 * img[i+(j+1)*nW] + a4 * img[i+(j+2)*nW]  + b1 * s1 + b2 * s2;
            s2 = s; 
        }

        bottomToUp[i] = s1;
        bottomToUp[i+nW] = s2;
        
        for (int j = 2; j < nH; j++){
            bottomToUp[i+j*nW] = a1*img[i+j*nW] + a2 * img[i+(j-1)*nW] + b1 * bottomToUp[i + (j-1) * nW] + b2 * bottomToUp[i + (j-2) * nW];
            s = s1;
            s1 =  a1*img[i+j*nW] + a2 * img[i+(j-1)*nW] + b1 * s1 + b2 * s2;
            s2 = s;
        }

        upToBottom[i+(nH-1)*nW] = s1;
        upToBottom[i+(nH-2)*nW] = s2;

        for (int j = nH-1; j >= 0; j--){
            upToBottom[i+j*nW] = a3 * img[i+(j+1)*nW] + a4 * img[i+(2+j)*nW] + b1 * upToBottom[i+(1+j)*nW] + b2 * upToBottom[i+(2+j)*nW];
        }
    }
    
    for (int i = 0; i < nW*nH; i++){
        res[i] = bottomToUp[i] + upToBottom[i];    
    } 
}
Traitement_deriche::Traitement_deriche(float a): Traitement_image(){
    this->alpha = a;
}
Traitement_deriche::~Traitement_deriche(){}
void Traitement_deriche::appliquer(){
    //TODO


    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_deriche::afficher(std::ostream& os){
    os << "TRAITEMENT DERICHE" << std::endl;
    os << "ALPHA : " << this->alpha <<  std::endl;
    Traitement_image::afficher(os);
}


/**   ====================         
 TRAITEMENT ZERO PASSING DETECTION        
 *    ====================**/
//On part du principe que les deux images sont données à la suite
class Traitement_zeroPassingDetection : public Traitement_image{
    private :
    int flg;
    bool reverse = false;
    float couleur_fond;

    void zeroPassingDetection_4(int nH, int nW, const OCTET* ImgIn , const OCTET* ImgIn_laplace, OCTET* ImgOut);
    void zeroPassingDetection_8(int nH, int nW, const OCTET* ImgIn , const OCTET* ImgIn_laplace, OCTET* ImgOut);

    public :
        Traitement_zeroPassingDetection(int flg, float couleur_fond = 0);
        virtual ~Traitement_zeroPassingDetection();
        void appliquer();
        virtual void afficher(std::ostream& os);
};

//Img in : pour chercher la direction
//Img laplace : pour chercher la valeur issue de la direction
void Traitement_zeroPassingDetection::zeroPassingDetection_8(int nH, int nW, const OCTET* ImgIn , const OCTET* ImgIn_laplace, OCTET* ImgOut){
       float pi8 = PI/8.0f;
        for (int i = 0; i < nH-1; i++){
            for (int j = 0; j < nW-1; j++){
                float h = ((reverse)? ImgIn[i*nW+j] - ImgIn[(i+1)*nW+j]: ImgIn[(i+1)*nW+j] - ImgIn[i*nW+j]);
                float w = ((reverse)? ImgIn[i*nW+j] - ImgIn[i*nW+j+1] : ImgIn[i*nW+j+1] - ImgIn[i*nW+j]);
                
                float val = std::max(abs(w), abs(h));

                val = (val<0)? 0 : (val>255)? 255 : val;

                float ang = atan2(h,w);

                //std::cout << "ANG : " << ang << std::endl;
                
                int indice_i, indice_j;

                //Si on est à gauche
                if (ang < -7.0*pi8) {
                    indice_i = 0;
                    indice_j = -1;
                }
                //On est en haut à gauche
                else if (ang < -5.0*pi8) {
                    indice_i = -1;
                    indice_j = -1;
                }
                //On est en haut
                else if (ang < -3.0*pi8) {
                    indice_i = -1;
                    indice_j = 0;
                }
                //On est en haut à droite
                else if (ang < -pi8) {
                    indice_i = -1;
                    indice_j = 1;
                }
                //On est à droite
                else if (ang < pi8) {
                    indice_i = 0;
                    indice_j = 1;
                }
                //On est en bas à droite
                else if (ang < 3.0*pi8) {
                    indice_i = 1;
                    indice_j = 1;
                }
                //On est en bas
                else if (ang < 5.0*pi8) {
                    indice_i = 1;
                    indice_j = 0;
                }
                //On est en bas à gauche
                else if (ang < 7.0*pi8) {
                    indice_i = 1;
                    indice_j = -1;
                }
                //On est à gauche de nouveau
                else {
                    indice_i = 0;
                    indice_j = -1;
                }
                
                //Si on sort du cadre de l'image avec les indices des voisins
                if (indice_i+i < 0 || indice_i+i > nH){
                    indice_i = 0;
                }
                if (indice_j+j < 0 || indice_j+j > nW){
                    indice_j = 0;
                }

                float val_lp = ImgIn_laplace[i*nW+j] - 128;
                float val_lp_vois = ImgIn_laplace[(i+indice_i)*nW+j+indice_j] - 128;

            
                if (val_lp > 0 && val_lp_vois < 0) {
                    ImgOut[i*nW+j] = val;
                }
                else if (val_lp < 0 && val_lp_vois > 0){
                    ImgOut[i*nW+j] = val;
                }
                else {
                    ImgOut[i*nW+j] = couleur_fond;
                }
            }
        }
}
void Traitement_zeroPassingDetection::zeroPassingDetection_4(int nH, int nW, const OCTET* ImgIn , const OCTET* ImgIn_laplace, OCTET* ImgOut){
        float pi4 = PI/4.0f;
        for (int i = 0; i < nH-1; i++){
            for (int j = 0; j < nW-1; j++){
                float h = ((reverse)? ImgIn[i*nW+j] - ImgIn[(i+1)*nW+j]: ImgIn[(i+1)*nW+j] - ImgIn[i*nW+j]);
                float w = ((reverse)? ImgIn[i*nW+j] - ImgIn[i*nW+j+1] : ImgIn[i*nW+j+1] - ImgIn[i*nW+j]);
                float val = std::max(abs(w), abs(h));

                val = (val<0)? 0 : (val>255)? 255 : val;

                float ang = atan2(h,w);

                //std::cout << "ANG : " << ang << std::endl;
                
                int indice_i, indice_j;

                //Si on est à gauche
                if (ang < -3.0*pi4) {
                    indice_i = 0;
                    indice_j = -1;
                }
                //On est en haut
                else if (ang < -pi4) {
                    indice_i = -1;
                    indice_j = 0;
                }
                //On est à droite
                else if (ang < pi4) {
                    indice_i = 0;
                    indice_j = 1;
                }
                //On est en bas
                else if (ang < 3.0*pi4) {
                    indice_i = 1;
                    indice_j = 0;
                }
                //On est à gauche de nouveau
                else {
                    indice_i = 0;
                    indice_j = -1;
                }
                
                //Si on sort du cadre de l'image avec les indices des voisins
                if (indice_i+i < 0 || indice_i+i > nH){
                    indice_i = 0;
                }
                if (indice_j+j < 0 || indice_j+j > nW){
                    indice_j = 0;
                }

                float val_lp = ImgIn_laplace[i*nW+j] - 128;
                float val_lp_vois = ImgIn_laplace[(i+indice_i)*nW+j+indice_j] - 128;

            
                if (val_lp > 0 && val_lp_vois < 0) {
                    ImgOut[i*nW+j] = val;
                }
                else if (val_lp < 0 && val_lp_vois > 0){
                    ImgOut[i*nW+j] = val;
                }
                else {
                    ImgOut[i*nW+j] = couleur_fond;
                }
            }
        }
}
Traitement_zeroPassingDetection::Traitement_zeroPassingDetection(int flg, float couleur_fond):  Traitement_image(){
    this->flg = flg;
    this->couleur_fond = couleur_fond;
}
Traitement_zeroPassingDetection::~Traitement_zeroPassingDetection(){}
void Traitement_zeroPassingDetection::appliquer(){
    bool write_img = getTraitementEcriture();
    if(AFFICHAGE_DEBUG){
        std::cout << "TRAITEMENT ZERO APPLIQUER" << std::endl;
        std::cout << "WRITE IMG : " << write_img << std::endl;
    }

    for (size_t k = 0; k < this->img_in.size(); k++){
        if(AFFICHAGE_DEBUG){
            std::cout << "IMG IN " << std::endl;
            std::cout << *img_in[k];
        }
        IMG * img_cur = this->img_in[k];
        
        OCTET *img_gradient;
        OCTET *img_laplacien;
        OCTET *img_out;

        int nTaille = img_cur->nH*img_cur->nW/2;
        allocation_tableau(img_gradient, OCTET, nTaille);
        allocation_tableau(img_laplacien, OCTET, nTaille);
        allocation_tableau(img_out, OCTET, nTaille/2);

        for (int i = 0; i < img_cur->nH; i++){
            for (int j = 0; j < img_cur->nW/2; j++){
                img_gradient[i*img_cur->nW/2+j] = img_cur->img_data[i*img_cur->nW+j];
            }
            for (int j = 0; j < img_cur->nW/2; j++){
                img_laplacien[i*img_cur->nW/2+j] = img_cur->img_data[i*img_cur->nW+j+img_cur->nW/2];
            }
        }

        if(flg == FLG_ZEROPASSING_4){
            zeroPassingDetection_4(img_cur->nH,img_cur->nW/2,img_gradient,img_laplacien,img_out);

            char* nom_zerocross = new char[TAILLE_MAX_NAME_IMG];
            sprintf(nom_zerocross,"%s_4zerocrossing", img_cur->nom); 
            this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW/2, img_out, write_img,nom_zerocross));
        }else if(flg == FLG_ZEROPASSING_8){
            zeroPassingDetection_8(img_cur->nH,img_cur->nW/2,img_gradient,img_laplacien,img_out);

            char* nom_zerocross = new char[TAILLE_MAX_NAME_IMG];
            sprintf(nom_zerocross,"%s_8zerocrossing", img_cur->nom); 
            this->img_out.push_back(new IMG(this->format_res,img_cur->nH, img_cur->nW/2, img_out, write_img,nom_zerocross));
        }
    }
    if(AFFICHAGE_DEBUG){
        afficher(std::cout);
    }
}
void Traitement_zeroPassingDetection::afficher(std::ostream& os){
    os << "TRAITEMENT ZERO DETECTION" << std::endl;
    os << "REVERSE : " << ((reverse)? "OUI" : "NON") <<  std::endl;
    os << "COULEUR FOND : " << couleur_fond <<  std::endl;
    os << "NB DIR : " << ((flg == FLG_ZEROPASSING_4)? "4" : ((flg == FLG_ZEROPASSING_8)? "8" : "NSP")) <<  std::endl;
    Traitement_image::afficher(os);
}

#endif