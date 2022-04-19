#ifndef PIPELINEFIRST_H
#define PIPELINEFIRST_H

#include <vector>
#include <stdio.h>
#include <iostream>
#include "../image_ppm.h"
#include <algorithm>  
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

#define PI 3.14159265358979323846


#define TAILLE_MAX_NAME_IMG 250
#define AFFICHAGE_DEBUG false

/*//Pour indiquer si on doit lire l'image d'ailleurs ou ecrire
#define FLG_IMG_READ 0
#define FLG_IMG_WRITE 1*/

//Format d'image (TraitementImage.set_format_res et IMG.format)
#define FLG_IMG_PGM 2
#define FLG_IMG_PPM 3

//Indiquer si l'image doit être prise en compte pour la suite du pipeline (IMG.continu_traitement)
#define FLG_IMG_NO_SUITE 4
#define FLG_IMG_SUITE 5

//Indiquer comment doit se sauvegarder une image résultante d'un traitement (Traitement.set_ecriture_res)
#define FLG_SAVEIMG_NO_WRITE 0
#define FLG_SAVEIMG_WRITE_SEPARATED 1
#define FLG_SAVEIMG_WRITE_REGROUPED 2
#define FLG_SAVEIMG_WRITE_BOTH 3

const char* FILE_PATH_RESULT = "res";
const char* FILE_PATH_INIT = "./images";
const char* DEFAULT_FILE_NAME = ""; 
const char* DEFAULT_FILE_NAME_FUSION = "fusionPipeline"; 


/**   ==================== ==================== ==================== ====================         
 *              ==================== ==================== ====================
 *    ==================== ==================== ==================== ====================**/
/**   ====================         
 *        CLASSE IMG        
 *    ====================**/
class IMG{
    public :

    int nH;
    int nW;
    OCTET* img_data;
    bool ecrire;
    const char* nom;
    int format;
    int continu_traitement = FLG_IMG_SUITE;

    IMG(const char* name, int flg = FLG_IMG_PGM);
    IMG(int flg_img,int h, int w, OCTET* data, bool write, const char* name, int ct = FLG_IMG_SUITE);
    IMG(IMG * img);
    
    virtual ~IMG();
    
    void init(int flg);
    bool save(const char* dossier);
    bool read();

    friend std::ostream& operator<<(std::ostream& os, const IMG& elem);
};

IMG::IMG(const char* name, int flg){
    this->nH = 0;
    this->nW = 0;
    this->nom = name;
    this->format = flg;
    this->ecrire = false;
    this->img_data = 0;
}
IMG::IMG(int flg_img, int h, int w, OCTET* data, bool write, const char* name, int ct){
    this->format = flg_img;
    nH = h;
    nW = w;
    img_data = data;
    ecrire = write;
    nom = name;
    this->continu_traitement = ct;
}
IMG::IMG(IMG * img){
    this->nH = img->nH;
    this->nW = img->nW;
    allocation_tableau(this->img_data, OCTET, this->nH *this->nW);
    for (int i = 0; i < this->nH *this->nW; i++){
        this->img_data[i] = img->img_data[i];
    }

    this->ecrire = img->ecrire;

    char* n = new char[TAILLE_MAX_NAME_IMG];
    sprintf(n,"%s", img->nom);
    this->nom = n;
    
    this->format = img->format;
    this->continu_traitement = img->continu_traitement;
}
IMG::~IMG(){
    free(img_data);
    delete[] nom; 
}
bool IMG::save(const char* dossier){
    if(AFFICHAGE_DEBUG) {
        std::cout << " => IMG ECRITURE DE L'IMAGE : " << nom << std::endl;
        std::cout << " => DANS LE DOSSIER : " << dossier << std::endl;
    }
    char* nom_save = new char[TAILLE_MAX_NAME_IMG];
    bool res = false;
    switch (format){
    case FLG_IMG_PGM:
        sprintf(nom_save,"%s/%s.pgm",dossier, nom);
        ecrire_image_pgm(nom_save, img_data,  nH, nW);    
        res = true;
        break;
    case FLG_IMG_PPM:
        sprintf(nom_save,"%s/%s.ppm",dossier, nom);
        ecrire_image_ppm(nom_save, img_data,  nH, nW);  
        res = true;  
        break;
    default:
        std::cout << "IMPOSSIBLE DE SAUVEGARDER L'IMAGE" << std::endl;
        break;
    }
    return res;
}
bool IMG::read(){
    char* path_img = new char[TAILLE_MAX_NAME_IMG];
    bool res = false;
    switch (format){
        case FLG_IMG_PGM:
            sprintf(path_img, "%s/%s.pgm", FILE_PATH_INIT,this->nom);
            std::cout << "INITIALISATION DE L'IMAGE : " << path_img << std::endl;
            lire_nb_lignes_colonnes_image_pgm(path_img, &this->nH, &this->nW);
            allocation_tableau(this->img_data, OCTET, this->nH*this->nW);
            lire_image_pgm(path_img, this->img_data, this->nH * this->nW);
            res = true;
            break;
        case FLG_IMG_PPM:
            sprintf(path_img, "%s/%s.ppm", FILE_PATH_INIT,this->nom);
            std::cout << "INITIALISATION DE L'IMAGE : " << path_img << std::endl;
            lire_nb_lignes_colonnes_image_ppm(path_img, &this->nH, &this->nW);
            allocation_tableau(this->img_data, OCTET, this->nH*this->nW*3);
            lire_image_ppm(path_img, this->img_data, this->nH * this->nW);
            res = true;
            break;
        default:
            std::cout << "READ : FORMAT D'IMAGE NON SUPPORTE" << std::endl;
            break;
    }
    return res;
}
std::ostream& operator<<(std::ostream& os, const IMG& elem){
    os << "*** IMAGE : " << elem.nom << std::endl;
    os << "=> ECRIRE : " << ((elem.ecrire)? "OUI":"NON") << std::endl;
    os << "=> LARGEUR : " << elem.nW << std::endl;
    os << "=> HAUTEUR : " << elem.nH << std::endl;
    return os;
}


/**   ====================         
 *      TRAITEMENT IMG        
 *    ====================**/
class Traitement_image{
protected:
    std::vector<IMG*> img_in;
    std::vector<IMG*> img_out;
    std::vector<Traitement_image*> traitement_suite;

    int format_res = FLG_IMG_PGM;
    int ecriture_res = FLG_SAVEIMG_WRITE_REGROUPED;

public:
    Traitement_image();
    virtual ~Traitement_image();

    Traitement_image* ajouter_img_in(IMG* img);
    Traitement_image* ajouter_traitement_img(Traitement_image* img);

    Traitement_image* set_ecriture_res(int flg);
    Traitement_image* set_format_res(int flg);

    virtual void appliquer();
    virtual void afficher(std::ostream& os);
    virtual bool getTraitementEcriture();

    std::vector<IMG*> getImgOut();

    friend class Pipeline;
    //friend std::ostream& operator<<(std::ostream& os, const Traitement_image& elem);
};

bool Traitement_image::getTraitementEcriture(){
    return (ecriture_res == FLG_SAVEIMG_WRITE_SEPARATED) || (ecriture_res == FLG_SAVEIMG_WRITE_BOTH);
}
Traitement_image::Traitement_image(){}
Traitement_image::~Traitement_image(){
    img_in.clear();
    img_out.clear();
    traitement_suite.clear();
}
Traitement_image* Traitement_image::ajouter_img_in(IMG* img){
    img_in.push_back(img);
    return this;
}
Traitement_image* Traitement_image::ajouter_traitement_img(Traitement_image* img){
    traitement_suite.push_back(img);
    return this;
}
void Traitement_image::appliquer(){
    for (size_t i = 0; i < this->img_in.size(); i++){
        this->img_out.push_back(this->img_in[i]);
    }
}
void Traitement_image::afficher(std::ostream& os){
    os << "*** TRAITEMENT IMAGE ***" << std::endl;
    os << " => IMAGE(S) IN : " << std::endl;
    for (size_t i = 0; i < img_in.size(); i++){
        os << *img_in[i];
    }
    os << " => IMAGE(S) OUT : " << std::endl;
    for (size_t i = 0; i < img_out.size(); i++){
        os << *img_out[i];
    }
    os << "NB TRAITEMENTS ENSUITE : " << traitement_suite.size() << std::endl;
    os << "TYPE RES : " << format_res << std::endl;
    os << "ECRITURE RES : " << ecriture_res << std::endl;
}
Traitement_image* Traitement_image::set_ecriture_res(int flg){
    this->ecriture_res = flg;
    return this;
}
Traitement_image* Traitement_image::set_format_res(int flg){
    this->format_res = flg;
    return this;
}
std::vector<IMG*> Traitement_image::getImgOut(){
    return this->img_out;
}


/**   ====================         
 *       CLASSE FUSION        
 *    ====================**/
class Fusion{
    private :

    OCTET* data_fu;
    int nH;
    int nW;
    int id;
    const char* nom;

    public :

    Fusion(const char* nom,int h, int w, OCTET* data, int ind);
    virtual ~Fusion();
    
    friend class Pipeline;
    friend bool compareFunction(Fusion* i1,Fusion* i2);
    friend std::ostream& operator<<(std::ostream& os, const Fusion& elem);
};

Fusion::Fusion(const char* nom,int h, int w, OCTET* data, int ind){
    id = ind;
    nH = h;
    nW = w;
    data_fu = data;
    this->nom = nom;
    if(AFFICHAGE_DEBUG) {
        std::cout << *this;
    }
}
Fusion::~Fusion(){
    free(data_fu);
}
std::ostream& operator<<(std::ostream& os, const Fusion& elem){
    os << "*** FUSION ***" << std::endl;
    os << "=> HAUTEUR : " << elem.nH<< std::endl;
    os << "=> LARGEUR : " << elem.nW<< std::endl;
    os << "=> ID : " << elem.id<< std::endl;
    return os;
}


/**   ====================         
 *          PIPELINE        
 *    ====================**/
//TODO mettre en place communication entre des pipeline pour n'appliquer des traitements qu'à une partie de la sortie d'un pipeline
class Pipeline{
private:
    Traitement_image* racine;
    std::vector<Fusion*> img_fusion_data;

    std::vector<const char*> nom_fusions;

    IMG* img_fus;
    char* dossier_save;
    int cpt=0;

    void sauvegarder_donnes(IMG* img);
    void execute_rec(Traitement_image* img);
    void fusion();
    int getIndiceImg();
    bool getFusionElem(Traitement_image* img);
    bool getProchainTraitementElem(IMG* img);
    void affImgFusion();

public:
    Pipeline(const char* file_name, Traitement_image* img);
    virtual ~Pipeline();

    Pipeline* init_pipeline();
    Pipeline* set_img_fus(IMG* img_fus);
    void execute();

    friend class Traitement_image;

    friend std::ostream& operator<<(std::ostream& os, const Pipeline& elem);
};

Pipeline::Pipeline(const char* file_name, Traitement_image* img){
    this->dossier_save = new char[TAILLE_MAX_NAME_IMG];
    sprintf( this->dossier_save, "%s/%s/%s", FILE_PATH_INIT,FILE_PATH_RESULT, file_name);
    img_fus = 0;
    racine = img;

    init_pipeline();
}
Pipeline::~Pipeline(){
    delete racine;
}
int Pipeline::getIndiceImg(){
    cpt++;
    return cpt;
}
bool Pipeline::getFusionElem(Traitement_image* img){
    return img->ecriture_res==FLG_SAVEIMG_WRITE_REGROUPED || img->ecriture_res==FLG_SAVEIMG_WRITE_BOTH;
}
bool Pipeline::getProchainTraitementElem(IMG* img){
    return img->continu_traitement==FLG_IMG_SUITE;
}
void Pipeline::execute_rec(Traitement_image* img){
    if(AFFICHAGE_DEBUG){
        std::cout << "===== PIPELINE REC" << std::endl;
        //img->afficher(std::cout); 
    }

    img->appliquer();
    
    for (size_t i = 0; i < img->img_out.size(); i++){
        if(img->img_out[i]->ecrire){
            if(AFFICHAGE_DEBUG) std::cout << "PIPELINE REC ECRIRE" << std::endl;
            sauvegarder_donnes(img->img_out[i]);  
        }

        if(getFusionElem(img)){
            if(AFFICHAGE_DEBUG) std::cout << "PIPELINE REC FUSION" << std::endl;
            img_fusion_data.push_back(new Fusion(img->img_out[i]->nom,img->img_out[i]->nH,img->img_out[i]->nW, img->img_out[i]->img_data, getIndiceImg()));
        }

        for (size_t j = 0; j < img->traitement_suite.size(); j++){
            if(getProchainTraitementElem(img->img_out[i])){
                img->traitement_suite[j]->img_in.push_back(img->img_out[i]);
            }
        }
    }
    for (size_t j = 0; j < img->traitement_suite.size(); j++){
        execute_rec(img->traitement_suite[j]);
    }

    if(AFFICHAGE_DEBUG) std::cout << "===== FIN PIPELINE REC" << std::endl;
}
void Pipeline::execute(){
    if(AFFICHAGE_DEBUG){
        std::cout << *this<<std::endl;
    }
    std::cout << "===== PIPELINE EXECUTION : " << dossier_save << std::endl;
    if(racine){
        if(AFFICHAGE_DEBUG){
            std::cout << "OK !" << std::endl;
        }
        execute_rec(racine);
    }else{
        if(AFFICHAGE_DEBUG) std::cout << "PIPELINE VIDE ! " << std::endl;
    }

    fusion();
}
void Pipeline::sauvegarder_donnes(IMG* img){
    /*//TODO gérer les différents cas de type d'images
    //type dans IMG
    if(AFFICHAGE_DEBUG) std::cout << " => PIPELINE ECRITURE DE L'IMAGE : " << img->nom << std::endl;
    char* nom = new char[100];
    sprintf(nom,"%s.pgm", img->nom);
    ecrire_image_pgm(nom, img->img_data,  img->nH, img->nW);*/

    if(!img->save(this->dossier_save)){
        std::cout << "PIPELINE : ERREUR DE SAUGARDE" <<std::endl;
    }
}
Pipeline* Pipeline::set_img_fus(IMG* img_f){
    img_fus = img_f;
    return this;
}
bool compareFunction(Fusion* i1,Fusion* i2){
    return (i1->id < i2->id);
}
void Pipeline::fusion(){
    if(!img_fus){
        img_fus = new IMG(DEFAULT_FILE_NAME_FUSION,FLG_IMG_PGM);
    }
        if(AFFICHAGE_DEBUG) {
            std::cout << "===== PIPELINE : DEBUT FUSION" << std::endl;
            std::cout << *this << std::endl;
        }
    
        sort(img_fusion_data.begin(),img_fusion_data.end(),compareFunction);
        
        if(AFFICHAGE_DEBUG) std::cout << "FIN SORT" << std::endl; 

        int nw_max = 0;
        int nh_max = 0;
    
        for (size_t i = 0; i < img_fusion_data.size(); i++){
            if(AFFICHAGE_DEBUG) std::cout << *img_fusion_data[i] << std::endl;
            if(img_fusion_data[i]->nH > nh_max){
                nh_max = img_fusion_data[i]->nH;
            }
            nw_max += img_fusion_data[i]->nW;
        }

        allocation_tableau(img_fus->img_data, OCTET, nw_max*nh_max);
        img_fus->nW = nw_max;
        img_fus->nH = nh_max;

    
        if(AFFICHAGE_DEBUG){
            std::cout << "ALLOCATION " << std::endl;
            std::cout << "->MAX H " << nh_max <<  std::endl;
            std::cout << "->MAX W " << nw_max <<  std::endl;
        }

        int nw_dep = 0;
        for (int k = 0; k < img_fusion_data.size(); k++){
            for (int i = 0; i< img_fusion_data[k]->nH; i++){
                for (int j = 0; j < img_fusion_data[k]->nW; j++){
                    img_fus->img_data[i*nw_max+j+nw_dep] = img_fusion_data[k]->data_fu[i * img_fusion_data[k]->nW + j];
                }
            }
            nw_dep += img_fusion_data[k]->nW;
            nom_fusions.push_back(img_fusion_data[k]->nom);
        }
        sauvegarder_donnes(img_fus);
        affImgFusion();
        if(AFFICHAGE_DEBUG) std::cout << "FIN FUSION : " << std::endl;
}
void Pipeline::affImgFusion(){
    std::cout << "IMAGE FUSION (de gauche à droite) :" << std::endl;
    for (size_t i = 0; i < this->nom_fusions.size(); i++){
        std::cout << " - " << i+1 << " : " << this->nom_fusions[i] << std::endl;
    }
}
Pipeline* Pipeline::init_pipeline(){
    if(AFFICHAGE_DEBUG) std::cout << " => SETUP DU DOSSIER : " << this->dossier_save << std::endl;
    if(fs::create_directories(this->dossier_save)){
        if(AFFICHAGE_DEBUG) std::cout << " => DOSSIER CREE" << std::endl;
    }else{
        if(AFFICHAGE_DEBUG) std::cout << " => DOSSIER NON CREE" << std::endl;
        char* nom_sup = new char[TAILLE_MAX_NAME_IMG];
        sprintf(nom_sup, "exec rm -r %s/*", this->dossier_save);
        system(nom_sup);
        if(AFFICHAGE_DEBUG) std::cout << " => ON EXEC : " << nom_sup << std::endl;
    }
    return this;
}
std::ostream& operator<<(std::ostream& os, const Pipeline& elem){
    os << "*** PIPELINE ***" << std::endl;
    os << "=> RACINE : " << std::endl;
    elem.racine->afficher(os);
    os << "=> FIN RACINE" << std::endl;

    os << "=> FUSIONS : " << std::endl;
    for (size_t i = 0; i < elem.img_fusion_data.size(); i++){
        os << *elem.img_fusion_data[i];
    }

    if(elem.img_fus){
        os << "=> IMG FUS : " << std::endl << *elem.img_fus;
    }else{
        os << "=> IMG FUS : VIDE" << std::endl;
    }
    
    os << "DOSSIER SAVE : " << *elem.dossier_save<< std::endl;
    return os;
}

#endif