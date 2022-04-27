#include "metrics_nothread.h"

#include <opencv2/opencv.hpp>
#include <opencv2/opencv_modules.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/core/utils/filesystem.hpp>
#include <opencv2/core/core.hpp> 
#include <opencv2/highgui.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <map>


#define W_SCREEN 1000
#define H_SCREEN 600

#define NBBLOCKX 100
#define NBBLOCKY 128


#define TRAITEMENTSIFTGRILLE    false
#define AFFICHAGE_FINAL         false
#define AFFICHAGE_MATCHING      false
#define AFFICHAGE_PRESENTATION  false

#define COMPUTE_CLASSIC         true

#define COMPUTE_EGALISATION     false
#define COMPUTE_DISTANCE        false
#define COMPUTE_QUAD_TREE       false
int quad_minSize = 1;
int quad_seuil = 10;

#define COMPUTE_MOY_VAR         false
float coef_moy                  = 1.f;
float coef_var                  = 0.1f;
float seuil_var_moy             = 1.0f;


const unsigned int nb_max_imagette    = 1000;
const unsigned int nb_max_imagette_ut = -1;
const float seuil_imagette_moy        = .5;
const float seuil_imagette_distance   = .0;
const float seuil_arret_distance             = -100000;
const float seuil_arret_distance_egalisation = -100000;

//Seuil pour le choix de passer sur un algo de moyenne ou non
const int Seuil_Passage_moy = -1;


using namespace cv;
//=== Gestion sauvegarde de métriques
std::string db_file_format          = ".pgm";
bool compute_moy_db                 = false;
bool compute_var_db                 = false;
bool compute_egalisation_db         = false;

bool compute_siftdata_db            = false;
bool compute_siftdataIMGKP_db       = false;



//=== Gestion sauvegarde des images (imwrite) 
std::string dir_name_result                 = "./res";
std::string dir_name_pre_result;
std::string dir_name_pre_result_block;
std::string dir_name_final_result;
const char* nomImg_finale;

bool save_img_grid                  = true;
bool save_img_grid_blocks           = false;
bool save_img_kp_gray               = true;
bool save_img_kp_color              = true;
bool save_step_kp                   = false;
bool save_grid_kp                   = false;

const int default_width_img_save    = 800;
const int default_height_img_save   = 500;

//=== Gestion affichage des images (imshow)
bool affichage_img                  = false;
const int default_width_img_aff     = 800;
const int default_height_img_aff    = 500;


//=== Données de l'algo FIRST
const int nb_imagettes_max    = 250;
const int taille_nom_img      = 250;

const int block_width         = 128;
const int block_height        = 128;

//SIFT PARAMS
//IMAGE ORIGINALE
//Le nombre de points d'intérets que l'on veut retenir (ces points sont rangés par leur score, ou contraste local dans SIFT).
const int original_nfeatures             =   0;
//Nombre de couches que l'on va prendre en compte dans la pyramide de gaussienne (niveaux de résolution). 3 est la valeur utilisée dans le papier de D. Lowe.
const int original_nOctaveLayers         =   3;
//Seuil de contraste pour filtrer les point de plus faibles intérêts (contraste) dans les régions semi-uniformes. Plus le seuil est grand, et moins le détecteur va produire de points d'intérêt.
const double original_contrastThreshold  =   0.04;  //Entre 0.03 et 0.09
//Seuil utilisé pour filtrer les points d'intérêts sur les arètes. Plus ce seuil est grand, et moins il y aura de points d'intérêts de filtrés.
const double original_edgeThreshold      =   10;
//Le sigma du filtre Gaussien appliqué à l'imge d'entrée. Si l'image de départ est de faible qualité, il vaut mieux avoir un sigma plus faible. 
const double original_sigma              =   1.6;

//Seuil pour la mise en comparaison ensuite (choix du PPV)
const double Seuil_Comparaison_Distance = 100.0f;



typedef struct block_img {
    cv::Mat img_data;
    int imagette;
    bool sift_imagette;
    int x0;
    int y0;

    bool egal_histo = false;
    cv::Mat histo_data;
}BLOCK;

//=== Fonctions UTILS
//Pour la sauvegarde de l'image
void sauvegarde(std::string name, Mat& img, int w = default_width_img_save, int h = default_height_img_save) {
    Mat resized;
    resize(img, resized, Size(w, h), INTER_LINEAR);
    imwrite(name, img);
}

//Pour l'affichage de l'image
void affichage(std::string name, Mat& img) {
      Mat resized;
      
      int dif_w = img.cols - W_SCREEN;
      int dif_h = img.rows - H_SCREEN;

      if (dif_h > 0 || dif_w > 0){
          int max = (dif_w > dif_h) ? dif_w : dif_h;
          resize(img, resized, Size(img.cols-max, img.rows-max), INTER_LINEAR);
      }

      namedWindow(name, WINDOW_AUTOSIZE);
      imshow(name,resized);
      affichage_img = true;
}

//Pour la sauvagarde de l'image de grille (décomposition)
void sauvegarde_grid_decomp(std::vector<BLOCK>& blocks) {
    std::cout << "SAUVEGARDE DE LA GRILLE (decomp)" << std::endl;
    for (int j = 0; j < blocks.size(); j++) {
        std::string blockId = std::to_string(j);
        std::string blockImgName = dir_name_pre_result_block +"/b_" + blockId + ".jpg";
        sauvegarde(blockImgName,blocks[j].img_data, blocks[j].img_data.cols, blocks[j].img_data.rows);
    }
}


void sauvegarde_grid_gray(std::vector<BLOCK>& blocks, Mat img_init, int nbbX, int nbbY) {
    uchar color_fond = 0;
    int taille_border = 1;

    Size sz_grd = Size(img_init.cols + (nbbX* taille_border )+ taille_border, img_init.rows + (nbbY* taille_border )+ taille_border);

    std::cout << "SAUVEGARDE IMAGE GRID (unique) : (" << sz_grd.width << "," << sz_grd.height << ")" << std::endl;
    Mat one_grid = Mat::zeros(sz_grd, CV_8UC1);
    std::cout << "nbROWS :  " << one_grid.rows << std::endl;
    std::cout << "nbCOLS :  " << one_grid.cols << std::endl;


    int cpt = 0;
    for (BLOCK img : blocks){

        int offset_borderX = (cpt%nbbX)*taille_border + taille_border;
        int offset_borderY = (cpt/nbbX)*taille_border + taille_border;

        cpt++;

        for (int i = 0; i < img.img_data.rows; i++){
            
            for (int j = 0; j < img.img_data.cols; j++) {

                one_grid.at<uchar>(
                    img.y0 + i + offset_borderY,
                    img.x0 + j + offset_borderX
                    ) = img.img_data.at<uchar>(i, j);
            }
        }
    }

    std::cout << "SAUVEGARDE" << std::endl;
    sauvegarde(dir_name_pre_result + "/img_grid.jpg", one_grid,one_grid.cols, one_grid.rows);
}


void sauvegarde_grid_gray_kp(std::vector<BLOCK>& blocks, Mat img_init, int nbbX, int nbbY) {
    uchar color_fond = 0;
    int taille_border = 1;

    Size sz_grd = Size(img_init.cols + (nbbX * taille_border) + taille_border, img_init.rows + (nbbY * taille_border) + taille_border);

    std::cout << "SAUVEGARDE IMAGE GRID KP (unique) : (" << sz_grd.width << "," << sz_grd.height << ")" << std::endl;
    std::cout << "WinitIMG : " << img_init.cols << std::endl;
    std::cout << "HinitIMG : " << img_init.rows << std::endl;
    std::cout << "PinitIMG : " << img_init.depth() << std::endl;

    Mat img_init_gray;
    cvtColor(img_init, img_init_gray, COLOR_BGR2GRAY);
    Mat one_grid = Mat::zeros(sz_grd, CV_8UC1);


    std::cout << "nbROWS :  " << one_grid.rows << std::endl;
    std::cout << "nbCOLS :  " << one_grid.cols << std::endl;


    int cpt = 0;
    for (BLOCK img : blocks) {

        int offset_borderX = (cpt % nbbX) * taille_border + taille_border;
        int offset_borderY = (cpt / nbbX) * taille_border + taille_border;

        cpt++;

        for (int i = 0; i < img.img_data.rows; i++) {

            for (int j = 0; j < img.img_data.cols; j++) {

                one_grid.at<uchar>(
                    img.y0 + i + offset_borderY,
                    img.x0 + j + offset_borderX
                    ) = img_init_gray.at<uchar>(i+ img.y0, j+ img.x0);
            }
        }
    }

    std::cout << "SAUVEGARDE" << std::endl;
    sauvegarde(dir_name_pre_result + "/img_grid_kp.jpg", one_grid, one_grid.cols, one_grid.rows);
}


void setNames(std::string img_name) {
    std::cout << "SET NAMES : " << img_name.c_str() << std::endl;

    dir_name_pre_result         = dir_name_result+"/"+ img_name + "/pre_result";
    dir_name_pre_result_block   = dir_name_pre_result + "/block_res";
    dir_name_final_result       = dir_name_result + "/" + img_name + "/final_result";
}


//=== Fonctions FIRST
//Diviser une image en blocs réguliers
int blockImage(const cv::Mat& img, const int blockWidth, const int blockHeight, std::vector<BLOCK>& blocks, int& nbBlockX, int& nbBlockY){
    if (!img.data || img.empty()){
        std::wcout << "Erreur lors du chargement de l'image au moment de la diviser" << std::endl;
        return EXIT_FAILURE;
    }

    //Dimensions de l'image
    int imgWidth = img.cols;
    int imgHeight = img.rows;
    std::cout << "DIM image : " << "(" << imgWidth << "," << imgHeight << ")" << std::endl;

    //Taille des blocs
    int bloc_wSize;
    int bloc_hSize;

    int y0 = 0;
    nbBlockY = 0;

    while (y0 < imgHeight){
        //Hauteur du bloc
        bloc_hSize = ((y0 + blockHeight) > imgHeight) * (blockHeight - (y0 + blockHeight - imgHeight)) + ((y0 + blockHeight) <= imgHeight) * blockHeight;

        int x0 = 0;
        nbBlockY++;

        nbBlockX = 0;
        while (x0 < imgWidth){
            //Lageur du bloc
            bloc_wSize = ((x0 + blockWidth) > imgWidth) * (blockWidth - (x0 + blockWidth - imgWidth)) + ((x0 + blockWidth) <= imgWidth) * blockWidth;

            BLOCK b = BLOCK();
            b.img_data = img(cv::Rect(x0, y0, bloc_wSize, bloc_hSize)).clone();
            b.x0 = x0;
            b.y0 = y0;
            //On récupère la partie qui nous intéresse de l'image de base (bloc régulier)
            blocks.push_back(b);

            x0 = x0 + blockWidth;
            nbBlockX++;
        }

        y0 = y0 + blockHeight;
    }
    return EXIT_SUCCESS;
}

//Pour le calcul des voisins les plus proches par rapport à notre image d'origine
//Distance euclidienne
float euclidDistance(const Mat& vec1,const Mat& vec2) {
    //std::cout << "EUCLIANDISTANCE" << std::endl;
    float sum = 0.0;
    
    int dim = vec1.cols;
    float v1_val;
    float v2_val;
    float val;

    for (int i = 0; i < dim; i++) {
        v1_val = vec1.at<float>(0, i);
        v2_val = vec2.at<float>(0, i);
        val = v1_val - v2_val;

        sum += val*val;
    }
    return sqrt(sum);
}

//Trouver l'indexe du plus proche voisin du point fixe, à un seuil près
int nearestNeighbor(const Mat& vec,const std::vector<KeyPoint> keypoints,const Mat& descriptors) {
    int neighbor = -1;
    double minDist = 1e6;

    //std::cout << "NEAREST NEIGHBOURG" << std::endl;
    KeyPoint pt;
    Mat v;
    double d;

    for (int i = 0; i < descriptors.rows; i++) {
        pt = keypoints[i];
        v = descriptors.row(i);
        d = euclidDistance(vec, v);
        
        //printf("%d %f\n", v.cols, d);
        if (d < minDist) {
            minDist = d;
            neighbor = i;
        }

        //v.release();
        //assert(!v.data);
    }

    if (minDist < Seuil_Comparaison_Distance) {
        return neighbor;
    }

    return -1;
}

//Définition des pairs de points proches entre eux (clusters)
void findPairs(const std::vector<KeyPoint> keypoints_or,const Mat& descriptors_or,const std::vector<KeyPoint> keypoints_imagette,const Mat& descriptors_imagette, std::vector<Point2f> srcPoints, std::vector<Point2f> dstPoints) {
    //std::cout << "FIND PAIR" << std::endl;
    KeyPoint pt1;
    Mat desc1;
    KeyPoint pt2;

    for (int i = 0; i < descriptors_or.rows; i++) {
        pt1 = keypoints_or[i];
        desc1 = descriptors_or.row(i);

        //std::cout << "     -> " << i << std::endl;

        int nn = nearestNeighbor(desc1, keypoints_imagette, descriptors_imagette);
        
        if (nn >= 0) {
            pt2 = keypoints_imagette[nn];
            srcPoints.push_back(pt1.pt);
            dstPoints.push_back(pt2.pt);
        }

        //desc1.release();
        //assert(!desc1.data);
    }
}


std::string getNameFromPath(std::string path) {
    int start = path.find_last_of("/\\") + 1;
    int end = path.find_last_of('.') - start;

    return path.substr(start, end);
}


void sift_block(BLOCK* curB,std::vector<KeyPoint> kp_block, std::vector<cv::Mat> imagettes, std::vector<std::vector<KeyPoint>> imagette_kps, std::vector<Mat> imagette_descriptors, std::vector<int> utilisee) {
    Mat descriptor_block;
    std::vector<Point2f> srcPoints;
    std::vector<Point2f> dstPoints;

    
    curB->sift_imagette = true;

    
    Ptr<DescriptorExtractor> extractor = SIFT::create();
    extractor->compute(curB->img_data, kp_block, descriptor_block);
    
    float best_score_block = -1;
    int best_ind_imagette = -1;

    
    for (int ind_imagette = 0; ind_imagette < imagettes.size(); ind_imagette++) {
        
        if (/*!utilisee[ind_imagette] && */ imagette_kps[ind_imagette].size() > 0) {
            
            findPairs(kp_block, descriptor_block, imagette_kps[ind_imagette], imagette_descriptors[ind_imagette], srcPoints, dstPoints);

            float cur_score = (float(srcPoints.size()) / float(kp_block.size())) * 10.0f;

            if (best_score_block < cur_score) {
                best_ind_imagette = ind_imagette;
                best_score_block = cur_score;
            }

            srcPoints.clear();
            dstPoints.clear();
        }
    }

    if (best_ind_imagette >= 0) {
        utilisee[best_ind_imagette] = true;
        curB->imagette = best_ind_imagette;
    }
    else {
        std::cout << "PAS DE SIFT TROUVE" << std::endl;
    }
}


void sift_block_bpb(BLOCK* curB, std::vector<KeyPoint> kp_block, std::vector<cv::Mat> imagettes, std::vector<std::vector<KeyPoint>> imagette_kps, std::vector<Mat> imagette_descriptors, std::vector<int> utilisee, Mat& descriptor_block) {
    std::vector<Point2f> srcPoints;
    std::vector<Point2f> dstPoints;

    curB->sift_imagette = true;
    float best_score_block = -1;
    int best_ind_imagette = -1;


    for (int ind_imagette = 0; ind_imagette < imagettes.size(); ind_imagette++) {

        if (/*!utilisee[ind_imagette] && */ imagette_kps[ind_imagette].size() > 0) {

            findPairs(kp_block, descriptor_block, imagette_kps[ind_imagette], imagette_descriptors[ind_imagette], srcPoints, dstPoints);

            float cur_score = (float(srcPoints.size()) / float(kp_block.size())) * 10.0f;

            if (best_score_block < cur_score) {
                best_ind_imagette = ind_imagette;
                best_score_block = cur_score;
            }

            srcPoints.clear();
            dstPoints.clear();
        }
    }

    if (best_ind_imagette >= 0) {
        utilisee[best_ind_imagette] = true;
        curB->imagette = best_ind_imagette;
    }
    else {
        std::cout << "PAS DE SIFT TROUVE" << std::endl;
    }
}


void moy_block(BLOCK* curB, std::vector<float> moyennes, std::vector<int> utilisee) {
    curB->sift_imagette = false;
    
    float couleur_moy_b = 0;
    for (int y = 0; y < curB->img_data.rows; y++) {
        for (int x = 0; x < curB->img_data.cols; x++) {
            couleur_moy_b += (int)curB->img_data.at<uchar>(y, x);
        }
    }
    couleur_moy_b /= float(curB->img_data.rows * curB->img_data.cols);

    float cur_best_moy_dif = 256;
    int cur_best_moy_ind = -1;

    for (int moy_ind = 0; moy_ind < moyennes.size(); moy_ind++) {
        int dif = couleur_moy_b - moyennes[moy_ind];
        dif = (dif < 0) ? - dif : dif;
        if (dif < cur_best_moy_dif) {
            cur_best_moy_ind = moy_ind;
            cur_best_moy_dif = dif;
        }
        else {
            if (dif < cur_best_moy_dif+seuil_imagette_moy && utilisee[moy_ind]< utilisee[cur_best_moy_ind]) {
                cur_best_moy_ind = moy_ind;
                cur_best_moy_dif = dif;
            }
        }
    }

    if (cur_best_moy_ind >= 0) {
        utilisee[cur_best_moy_ind]++;
        curB->imagette = cur_best_moy_ind;
    }
    else {
        std::cout << "PAS DE MOYENNE TROUVEE" << std::endl;
    }
}


void moy_variance_block(BLOCK* curB, std::vector<float> moyennes, std::vector<float> var, std::vector<int> utilisee) {
    curB->sift_imagette = false;

    double* nbElementsLus = new double[256]{};
    int nTaille = curB->img_data.rows * curB->img_data.cols;

    for (int i = 0; i < curB->img_data.rows; i++) {
        for (int j = 0; j < curB->img_data.cols; j++) {
            int val = (int)curB->img_data.at<uchar>(i, j);
            nbElementsLus[val]++;
        }
    }
    float moyenne = 0.;
    for (int i = 0; i < 256; i++) {
        moyenne += (i * nbElementsLus[i]);
    }
    moyenne /= (float)nTaille;

    float variance = 0.;
    for (int i = 0; i < 256; i++) {
        variance += (i * i * nbElementsLus[i]);
    }
    variance = (variance / (float)nTaille) - (moyenne * moyenne);


    float cur_best_score = moyennes[0] + var[0];
    int cur_best_moy_ind = 0;

    float initial_score = coef_moy * moyenne + coef_var * variance;

    for (int moy_ind = 0; moy_ind < moyennes.size(); moy_ind++) {
        float currentScore = abs((moyennes[moy_ind] * coef_moy + coef_var * var[moy_ind]) - initial_score);
        float actual_bestscore = abs(cur_best_score - initial_score);
        if (currentScore < actual_bestscore) {
            cur_best_moy_ind = moy_ind;
            cur_best_score = currentScore;
        }else {
            if (currentScore < (actual_bestscore + seuil_var_moy) && utilisee[moy_ind] < utilisee[cur_best_moy_ind]) {
                cur_best_moy_ind = moy_ind;
                cur_best_score = currentScore;
            }
        }
    }

    if (cur_best_moy_ind >= 0) {
        utilisee[cur_best_moy_ind]++;
        curB->imagette = cur_best_moy_ind;
    }
    else {
        std::cout << "PAS DE MOYENNE VAR TROUVEE" << std::endl;
    }
}


void dist_block(BLOCK* curB, std::vector<Mat> reshaped, std::vector<int> utilisee) {
    curB->egal_histo = true;
    //std::cout << "DEBUT ALGO" << std::endl;
    //std::cout << "W bdata : " << curB->img_data.rows <<  std::endl;
    //std::cout << "H bdata : " << curB->img_data.cols << std::endl;
    //std::cout << "depth bdata : " << curB->img_data.depth() << std::endl;
    //std::cout << "W bhisto : " << curB->histo_data.rows << std::endl;
    //std::cout << "H bhisto : " << curB->histo_data.cols << std::endl;
    //std::cout << "depth bhisto : " << curB->histo_data.depth() << std::endl;

    Mat block = Mat::zeros(curB->img_data.size(), curB->img_data.type());
    //std::cout << "W block : " << block.rows << std::endl;
    //std::cout << "H block : " << block.cols << std::endl;
    //std::cout << "depth block : " << block.depth() << std::endl;

    blurlin(curB->img_data, block);
    int histo[256] = {};
    float ddp[256];
    float F[256];
    float Finv[256];
    int min, max;

    //std::cout << "cpt histo : " << std::endl;
    computeHisto(block, histo);
    //std::cout << "cpt ddp : " << std::endl;
    computeDdp(block, histo, ddp);
    //std::cout << "cpt F : " << std::endl;
    computeF(ddp, F);
    //std::cout << "cpt inv : " << std::endl;
    computeInverse(F, Finv);
    //std::cout << "cpt minmax : " << std::endl;
    computeMinMaxFromHisto(histo, &min, &max);

    //std::cout << "applyF : " << std::endl;
    applyF(block, F, block);

    //Parcourir les images et trouver celle avec une distance la plus petite
    float valeurPlusProche = distance(block, reshaped[0]);
    int indicePlusProche = 0;
    //std::cout << "parcour img" << std::endl;
    for (int k = 1; k < reshaped.size(); k++) {
        float valeur = distance(block, reshaped[k]);
        if (valeur < valeurPlusProche) {
            valeurPlusProche = valeur;
            indicePlusProche = k;
        }
        //else {
        //    if (valeur < valeurPlusProche + seuil_imagette_distance && utilisee[k] < utilisee[indicePlusProche]) {
        //        valeurPlusProche = valeur;
        //        indicePlusProche = k;
        //    }
        //}
        //if (valeurPlusProche < seuil_arret_distance_egalisation) {
        //    break;
        //}
    }

    //std::cout << "ECRITURE DES VALEURS" << std::endl;
    int valeur;
    if (indicePlusProche >= 0) {
        utilisee[indicePlusProche]++;

        for (int y = 0; y < curB->histo_data.rows; y++) {
            for (int x = 0; x < curB->histo_data.cols; x++) {
                valeur = reshaped[indicePlusProche].at<uchar>(y,x);
                curB->histo_data.at<uchar>(y,x) = (int)clamp(255.0 * Finv[(int)valeur], min, max);
            }
        }
    }
    else {
        std::cout << "PAS DE Plus Proche TROUVE" << std::endl;
    }
}


void showMatching(std::string name,std::string nom1, std::string nom2, DescriptorMatcher::MatcherType descriptorMatcherType, float ratio_thresh = 0.7, int nb_match = 2) {
    Mat img1_col = imread(nom1, IMREAD_ANYCOLOR | IMREAD_ANYDEPTH);
    Mat img2_col = imread(nom2, IMREAD_ANYCOLOR | IMREAD_ANYDEPTH);

    Mat img1 = imread(nom1, IMREAD_GRAYSCALE);
    Mat img2 = imread(nom2, IMREAD_GRAYSCALE);
    
    Ptr<SIFT> detector = SIFT::create();
    std::vector<KeyPoint> keypoints1, keypoints2;
    Mat descriptors1, descriptors2;
    
    detector->detectAndCompute(img1, noArray(), keypoints1, descriptors1);
    detector->detectAndCompute(img2, noArray(), keypoints2, descriptors2);
    
    //Matching
    Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create(descriptorMatcherType);
    std::vector< std::vector<DMatch> > knn_matches;
    matcher->knnMatch(descriptors1, descriptors2, knn_matches, nb_match);
    
    //Filtration
    std::vector<DMatch> good_matches;
    for (size_t i = 0; i < knn_matches.size(); i++){
        if (knn_matches[i][0].distance < ratio_thresh * knn_matches[i][1].distance){
            good_matches.push_back(knn_matches[i][0]);
        }
    }

    Mat img_matches;
    drawMatches(
        img1_col, 
        keypoints1, 
        img2_col, 
        keypoints2, 
        good_matches, 
        img_matches, 
        Scalar::all(-1),
        Scalar::all(-1), 
        std::vector<char>(), 
        DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS
    );

    affichage(name, img_matches);
}


void CptMatchingGrid(BLOCK* curB, std::vector<KeyPoint> kp_block, std::vector<cv::Mat> imagettes, std::vector<std::vector<KeyPoint>> imagette_kps, std::vector<Mat> imagette_descriptors, std::vector<bool> utilisee, Mat& descriptor_block, DescriptorMatcher::MatcherType descriptorMatcherType, float ratio_thresh = 0.7, float distances = 10.0f) {
    std::cout << " === MATCHING GRID === " << std::endl;
    std::cout << " ELEM : " << descriptorMatcherType << std::endl;
    std::cout << " SEUIL : " << ratio_thresh << std::endl;

    //Matching
    Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create(descriptorMatcherType);
    std::vector< std::vector<DMatch> > matches;
    std::vector<DMatch> good_matches;

    curB->sift_imagette = true;
    float best_score_block = -1;
    int best_ind_imagette = -1;


    for (int ind_imagette = 0; ind_imagette < imagettes.size(); ind_imagette++) {

        if (imagette_kps[ind_imagette].size() > 0) {
            matcher->knnMatch(descriptor_block, imagette_descriptors[ind_imagette], matches, 2);
            
            //Filtration
            for (size_t i = 0; i < matches.size(); i++) {
                if (matches[i][0].distance < ratio_thresh * matches[i][1].distance) {
                    good_matches.push_back(matches[i][0]);
                }
            }

            float cur_score = (float(good_matches.size()) / float(matches.size())) * 10.0f;

            if (best_score_block > cur_score) {
                best_ind_imagette = ind_imagette;
                best_score_block = cur_score;
            }

            good_matches.clear();
            matches.clear();
        }
    }

    if (best_ind_imagette >= 0) {
        utilisee[best_ind_imagette] = true;
        curB->imagette = best_ind_imagette;
    }
    else {
        std::cout << "PAS DE SIFT TROUVE" << std::endl;
    }
}



/**===========================================================================================================**/
struct Block {
    int startI, startJ;
    int endI, endJ;

    int width, height;

    int indexImgChoosen;

    Block** childs;
    Block* parent;

    Block(int startI, int startJ, int width, int height) {
        this->startI = startI;
        this->startJ = startJ;
        this->width = width;
        this->height = height;
        this->endI = startI + width;
        this->endJ = startJ + width;

        this->indexImgChoosen = -1;

        childs = new Block * [4];
        childs[0] = NULL;
        childs[1] = NULL;
        childs[2] = NULL;
        childs[3] = NULL;
        parent = NULL;
    }

    void setChild(int i, Block* child) {
        child->parent = this;
        childs[i] = child;
    }
};

void quadTree(Block* current, Mat& ImgIn, int seuil, int minSize, std::vector<Mat> imagettes, std::vector<BLOCK>& res) {
    std::cout << "CALL QUAD TREE" << std::endl;
    
    int startI = current->startI;
    int startJ = current->startJ;
    int endI = current->endI;
    int endJ = current->endJ;

    int bW = current->width;
    int bH = current->height;

    int bW2 = bW / 2;
    int bH2 = bH / 2;
    int bN = bW * bH;

    std::vector<Mat> imagettes_resized_2;
    imagettes_resized_2.resize(imagettes.size());
    for (size_t i = 0; i < imagettes.size(); i++) {
        resize(imagettes[i], imagettes_resized_2[i], Size(bW, bH), INTER_LINEAR);

        //egaliser(imagettes_resized_2[i], imagettes_resized_2[i]);
        //blurlin(imagettes_resized_2[i], imagettes_resized_2[i]);
    }

    Mat bloc_img = ImgIn(cv::Rect(startI, startJ, bW, bH)).clone();
    //blurlin(bloc_img, bloc_img);


    /*int histo[256] = {}; //Initialisé à 0
    float ddp[256];
    float F[256];
    float Finv[256];
    int min, max;

    computeHisto(bloc_img, histo);
    computeDdp(bloc_img, histo, ddp);
    computeF(ddp, F);
    computeInverse(F, Finv);
    computeMinMaxFromHisto(histo, &min, &max);

    applyF(bloc_img, F, bloc_img);*/

    // parcourir les images et trouver celle avec une distance la plus petite
    float value = distance(bloc_img, imagettes_resized_2[0]);
    int indicePlusProche = 0;
    float valeur;
    for (int k = 1; k < imagettes_resized_2.size(); k++) {
        valeur = distance(bloc_img, imagettes_resized_2[k]);
        if (valeur < value) {
            value = valeur;
            indicePlusProche = k;
        }
    }

    std::cout << "VALUE " << value << std::endl;

    if (value < seuil || bH <= minSize || bW <= minSize) {
        std::cout << "COPIE " << startI << " " << startJ << std::endl;
        //Création d'un block pour le rendu
        BLOCK b = BLOCK();
        b.img_data = bloc_img;
        b.x0 = startI;
        b.y0 = startJ;
        res.push_back(b);
    }
    else {
        std::cout << "RECURSION " << startI << " " << startJ << std::endl;
        Block* HG = new Block(startI, startJ, bW2, bH2);
        Block* HD = new Block(startI, startJ + bW2, bW2, bH2);
        Block* BG = new Block(startI + bH2, startJ, bW2, bH2);
        Block* BD = new Block(startI + bH2, startJ + bW2, bW2, bH2);

        current->setChild(0, HG);
        current->setChild(1, HD);
        current->setChild(2, BG);
        current->setChild(3, BD);

        quadTree(HG, ImgIn, seuil, minSize, imagettes, res);
        quadTree(HD, ImgIn, seuil, minSize, imagettes, res);
        quadTree(BG, ImgIn, seuil, minSize, imagettes, res);
        quadTree(BD, ImgIn, seuil, minSize, imagettes, res);
    }
}

void callQuadTree(Mat& ImgIn, int seuil, int minSize, std::vector<Mat> imagettes, std::vector<BLOCK>& res) {
    Block* current = new Block(0, 0, ImgIn.cols, ImgIn.rows);
    int bW2 = ImgIn.cols / 2;
    int bH2 = ImgIn.rows / 2;

    Block* HG = new Block(0, 0, bW2, bH2);
    Block* HD = new Block(0, 0 + bW2, bW2, bH2);
    Block* BG = new Block(0 + bH2, 0, bW2, bH2);
    Block* BD = new Block(0 + bH2, 0 + bW2, bW2, bH2);

    current->setChild(0, HG);
    current->setChild(1, HD);
    current->setChild(2, BG);
    current->setChild(3, BD);

    quadTree(HG, ImgIn, seuil, minSize, imagettes, res);
    quadTree(HD, ImgIn, seuil, minSize, imagettes, res);
    quadTree(BG, ImgIn, seuil, minSize, imagettes, res);
    quadTree(BD, ImgIn, seuil, minSize, imagettes, res);
}
/**===========================================================================================================**/



int main(int argc, char** argv) {
    std::cout << "=== MOSAIQUE PAR POINTS D'INTERETS ===" << std::endl;

    if (argc < 6) {
        std::cerr << "PAS ASSEZ D'ARGUMENTS POUR LA MOSAIQUE \n => usage : ImageOrigine.pgm ImageOut.pgm nbBlockX nbBlockY db_dir\n " << std::endl;
        return -1;
    }

    /* =========== =========== VARIABLES =========== =========== */
    //Nom de l'image produite finale
    char* nomImgFinale = new char[taille_nom_img];
    //Imagettes
    std::vector<cv::Mat> imagettes;
    //Blocs réguliers de l'image de base
    std::vector<BLOCK> blocks;
    int nbBlockX, nbBlockY;

    Mat originalImage_color;
    Mat originalImage_gray;


    //En premier on met à jour les noms
    setNames(getNameFromPath(argv[1]));

    //CREATION DES REPERTOIRES DE SAUVEGARDE DES FICHIERS
    {
        if (cv::utils::fs::createDirectories(dir_name_pre_result.c_str())) {
            std::cout << "CREATION DU REPERTOIRE : " << dir_name_pre_result.c_str() << std::endl;
        }
        if (cv::utils::fs::createDirectories(dir_name_final_result.c_str())) {
            std::cout << "CREATION DU REPERTOIRE : " << dir_name_final_result.c_str() << std::endl;
        }
        if (cv::utils::fs::createDirectories(dir_metric.c_str())) {
            std::cout << "CREATION DU REPERTOIRE : " << dir_metric.c_str() << std::endl;
        }
        if (cv::utils::fs::createDirectories(dir_name_pre_result_block.c_str())) {
            std::cout << "CREATION DU REPERTOIRE : " << dir_name_pre_result_block.c_str() << std::endl;
        }
    }

    /* =========== =========== LECTURES DES IMAGES =========== =========== */
    {
        std::cout << "=== LECTURE DES IMAGES ===" << std::endl;

        std::cout << "LECTURE DE L'IMAGE MOSAIQUE : ";
        originalImage_gray = imread(argv[1], IMREAD_GRAYSCALE);
        if (!originalImage_gray.data) {
            std::cout << "ERREUR" << std::endl;
            return -1;
        }
        else {
            std::cout << "OK" << std::endl;
        }

        std::cout << "LECTURE DE L'IMAGE MOSAIQUE (couleur) : ";
        originalImage_color = imread(argv[1], IMREAD_ANYCOLOR | IMREAD_ANYDEPTH);
        if (!originalImage_color.data) {
            std::cout << "ERREUR" << std::endl;
            return -1;
        }
        else {
            std::cout << "OK" << std::endl;
        }
        std::cout << "=> Width : " << originalImage_gray.cols << std::endl;
        std::cout << "=> Height: " << originalImage_gray.rows << std::endl;

        sprintf_s(nomImgFinale, taille_nom_img, "%s/%s", dir_name_final_result.c_str(), argv[2]);
        std::cout << "IMAGE FINALE STOCKEE : " << nomImgFinale << std::endl;
    }

    //EXTRACTION DES IMAGETTES + METRIQUES
    std::vector<std::string> nom_imagette_sift;
    {
        std::cout << "=== CALCUL DES METRIQUES ===" << std::endl;
        std::vector<cv::String> nom_imagette;
        char* fp = new char[taille_nom_img];

        sprintf_s(fp, taille_nom_img, "%s/*%s", ((Seuil_Passage_moy<0)? argv[5] : "./db_500"), db_file_format.c_str());
        glob(fp, nom_imagette, false);

        std::cout << "LECTURE DES IMAGETTES : " << ((nom_imagette.size() < nb_max_imagette)? nom_imagette.size() : nb_max_imagette) << std::endl;

        for (int i = 0; i < nom_imagette.size() && i < nb_max_imagette; i++) {
            std::cout << " + AVANCEMENT : " << int(((i + 1) / float(nom_imagette.size())) * 100) << " %\r" << std::flush;
            nom_imagette_sift.push_back(getNameFromPath(nom_imagette[i]));
            imagettes.push_back(imread(nom_imagette[i], IMREAD_GRAYSCALE));
        }
        std::cout << "\nFIN LECTURE DES IMAGETTES"<< std::endl;


        //On calcule les métriques
        if (compute_moy_db) {
            std::cout << "CALCUL DE LA MOYENNE" << std::endl;
            calc_moy(imagettes);
        }
        if (compute_siftdata_db) {
            std::cout << "COMPUTE KP" << std::endl;
            calc_sift(imagettes, nom_imagette_sift);
        }
        if (compute_siftdataIMGKP_db) {
            std::cout << "COMPUTE KP img" << std::endl;
            compute_kp_img(imagettes, nom_imagette_sift);
        }
        if (compute_var_db) {
            std::cout << "CALCUL DE LA VARIANCE" << std::endl;
            calc_var(imagettes);
        }
        if (compute_egalisation_db) {
            std::cout << "CALCUL DE L'EGALISATION" << std::endl;
            calc_egalisation(imagettes, nom_imagette_sift);
        }
    }


    Mat res_kp_gray;
    Mat res_kp_color;
    Ptr<DescriptorExtractor> extractor = SIFT::create();

    std::vector<KeyPoint>               originale_SIFT_keypoints;
    cv::Ptr<cv::SIFT>                   originale_detector;
    Mat                                 originale_descriptor;
    if (Seuil_Passage_moy < 0) {
        //CALCUL DES POINTS D'INTERETS DE L'IMAGE MOSAIQUE
        {
            std::cout << "=== TRAITEMENT SIFT ===" << std::endl;
            std::cout << "CREATION DU DETECTEUR SIFT" << std::endl;
            //SIFT
            originale_detector = cv::SIFT::create(
                original_nfeatures,
                original_nOctaveLayers,
                original_contrastThreshold,
                original_edgeThreshold,
                original_sigma
            );

            std::cout << "DEBUT DE LA DETECTION" << std::endl;
            originale_detector->detect(originalImage_gray, originale_SIFT_keypoints);
            printf("===> FIN DETECTION : %d points d'interet ont ete trouves.\n", (int)originale_SIFT_keypoints.size());

            if (save_step_kp) {
                Mat res_kp_gray;
                cv::drawKeypoints(originalImage_gray, originale_SIFT_keypoints, res_kp_gray, cv::Scalar(255, 255, 0));
                //TODO mettre desronds et tout et meme une image de matching ensuite
                sauvegarde(dir_name_pre_result + "/kp_detection.jpg", res_kp_gray, res_kp_gray.cols, res_kp_gray.rows);
            }

            std::cout << "DEBUT DE L'EXTRACTION" << std::endl;
            extractor->compute(originalImage_gray, originale_SIFT_keypoints, originale_descriptor);


            if (save_step_kp) {
                sauvegarde(dir_name_pre_result + "/descriptor.jpg", originale_descriptor, originale_descriptor.cols, originale_descriptor.rows);
            }
            printf("===> FIN EXTRACTION : %d points d'interet ont ete trouves.\n", (int)originale_SIFT_keypoints.size());


            if (save_grid_kp) {
                cv::drawKeypoints(originalImage_gray, originale_SIFT_keypoints, res_kp_gray, cv::Scalar::all(-1), cv::DrawMatchesFlags::DRAW_RICH_KEYPOINTS);
                if (save_img_kp_gray) {
                    sauvegarde(dir_name_pre_result + "/gray_kp.jpg", res_kp_gray, res_kp_gray.cols, res_kp_gray.rows);
                }
            }
            if (save_img_kp_color) {
                cv::drawKeypoints(originalImage_color, originale_SIFT_keypoints, res_kp_color, cv::Scalar::all(-1), cv::DrawMatchesFlags::DRAW_RICH_KEYPOINTS);
                sauvegarde(dir_name_pre_result + "/color_kp.jpg", res_kp_color, res_kp_color.cols, res_kp_color.rows);
            }
        }
    }


    int wblock;
    int hblock;
    //TRAITEMENT DE LA GRILLE
    {
        std::cout << "CALCUL DE LA GRILLE" << std::endl;
        sscanf_s(argv[3], "%d", &nbBlockX);
        sscanf_s(argv[4], "%d", &nbBlockY);

        nbBlockX = NBBLOCKX;
        nbBlockY = NBBLOCKY;

        std::cout << "=> Nombre de blocs en X : " << nbBlockX << std::endl;
        std::cout << "=> Nombre de blocs en Y : " << nbBlockY << std::endl;

        wblock = originalImage_gray.cols / nbBlockX;
        hblock = originalImage_gray.rows / nbBlockY;

        std::cout << "=> largeur bloc : " << wblock << std::endl;
        std::cout << "=> hauteur bloc : " << hblock << std::endl;

        int nb_blockX_reel = 0;
        int nb_blockY_reel = 0;
        int res_block = blockImage(originalImage_gray, wblock, hblock, blocks, nb_blockX_reel, nb_blockY_reel);
        std::cout << "==> Nombre reel de blocks : " << blocks.size() << std::endl;
        std::cout << "==> Nombre reel de blocks en X : " << nb_blockX_reel << std::endl;
        std::cout << "==> Nombre reel de blocks en Y : " << nb_blockY_reel << std::endl;

        if (save_img_grid_blocks) {
            sauvegarde_grid_decomp(blocks);
        }
        if (save_img_grid) {
            sauvegarde_grid_gray(blocks, originalImage_gray, nb_blockX_reel, nb_blockY_reel);
        }
        if (save_grid_kp) {
            sauvegarde_grid_gray_kp(blocks, res_kp_gray, nb_blockX_reel, nb_blockY_reel);
        }
    }


    std::vector<float> moyennes;

    if (COMPUTE_CLASSIC) {
        //TRAITEMENT DE SIFT SUR LA GRILLE
        if (TRAITEMENTSIFTGRILLE && (Seuil_Passage_moy >= 0)) {
            std::vector<KeyPoint>* kp_pour_b = new std::vector<KeyPoint>[blocks.size()];
            std::vector<Mat> descriptor_pour_b;
            descriptor_pour_b.resize(blocks.size());

            cv::Ptr<cv::SIFT> b_detector = cv::SIFT::create();

            //CALCUL DES SIFTS POUR CHAQUES BLOCS
            {
                std::cout << "=== TRAITEMENT SIFT (grille) ===" << std::endl;
                for (int i = 0; i < blocks.size(); i++) {
                    std::cout << " + AVANCEMENT : " << int(((i + 1) / float(blocks.size())) * 100) << " %\r" << std::flush;
                    b_detector->detect(blocks[i].img_data, kp_pour_b[i]);
                    extractor->compute(blocks[i].img_data, kp_pour_b[i], descriptor_pour_b[i]);
                }
            }

            //PHASE DE CALCUL DE L'IMAGE MOSAIQUE
            {
                std::cout << "\nCOMPARAISON ET CHOIX DES IMAGETTES BLOC PAR BLOC" << std::endl;
                std::cout << "RECUPERATION DES MOYENNES" << std::endl;

                std::vector<int> utilisee;
                get_moy(imagettes.size(), moyennes);
                for (size_t i = 0; i < moyennes.size(); i++) {
                    utilisee.push_back(0);
                }
                std::cout << " => MOYENNES : " << moyennes.size() << std::endl;
                std::cout << " => UTILISEE : " << utilisee.size() << std::endl;


                std::vector<Mat> imagette_descriptors;
                std::vector<std::vector<KeyPoint>> imagette_kps;
                imagette_descriptors.resize(imagettes.size());
                imagette_kps.resize(imagettes.size());

                for (int i = 0; i < imagettes.size(); i++) {
                    std::cout << "RECUPERATION DES DECRIPTEURS ET KP : " << int(((i + 1) / float(imagettes.size())) * 100) << " %\r" << std::flush;
                    std::string name_img = nom_imagette_sift[i];

                    //recuperer_keypoint(name_img, imagette_kps[i]);
                    recuperer_donnees_sift_kp(name_img, imagette_kps[i]);
                    recuperer_descriptor(name_img, imagette_descriptors[i]);
                }

                int nb_sift_test = 0;
                int nb_moy_test = 0;

                std::vector<KeyPoint> kp_block;
                BLOCK* curB;

                std::cout << "\nDEBUT DE L'ALGO" << std::endl;

                for (int i = 0; i < blocks.size(); i++) {
                    std::cout << " + AVANCEMENT : " << int(((i + 1) / float(blocks.size())) * 100) << " %\r" << std::flush;

                    kp_block = kp_pour_b[i];
                    curB = &blocks[i];

                    if (kp_block.size() < Seuil_Passage_moy) {
                        //Traitement moyenne
                        nb_moy_test++;
                        moy_block(curB, moyennes, utilisee);
                    }
                    else {
                        //Traitement SIFT
                        nb_sift_test++;
                        sift_block_bpb(curB, kp_block, imagettes, imagette_kps, imagette_descriptors, utilisee, descriptor_pour_b[i]);
                    }
                }

                std::cout << "\nFIN SIFT : " << nb_sift_test << std::endl;
                std::cout << "FIN MOY : " << nb_moy_test << std::endl;
            }
        }
        else {
            if (Seuil_Passage_moy >= 0) {
                //GESTION DES POINTS D'INTERETS PAR BLOCKS
                std::vector<KeyPoint>* kp_pour_b = new std::vector<KeyPoint>[blocks.size()];
                {
                    std::cout << "IMPLEMENTATION DES POINTS D'INTERETS DANS LES BLOCS" << std::endl;
                    std::cout << "->Nombre de KP : " << originale_SIFT_keypoints.size() << std::endl;
                    std::cout << "->Nombre de blocks : " << blocks.size() << std::endl;

                    std::vector<KeyPoint>* kp_pour_b_pre = new std::vector<KeyPoint>[blocks.size()];

                    //Calcul des keypoints dans les blocks
                    for (KeyPoint kp : originale_SIFT_keypoints) {
                        int x_kp = kp.pt.x;
                        int y_kp = kp.pt.y;


                        for (int i = 0; i < blocks.size(); i++) {
                            BLOCK b = blocks[i];
                            //Si le kp est dans le block
                            if (b.x0 < x_kp && x_kp <= (b.x0 + b.img_data.cols) && b.y0 < y_kp && y_kp <= (b.y0 + b.img_data.rows)) {
                                kp_pour_b[i].push_back(kp);
                                //kp_pour_b_pre[i].push_back(kp);
                                break;
                            }
                        }
                    }

                    /*std::cout << "RECADRAGE DES POINTS D'INTERETS AU SEIN DU BLOC" << std::endl;
                    {
                        for (int i = 0; i < blocks.size(); i++) {
                            std::vector<KeyPoint> kpb = kp_pour_b_pre[i];
                            for (size_t j = 0; j < kpb.size(); j++){
                                KeyPoint toclone = kpb[j];

                                KeyPoint curkp = KeyPoint(
                                    toclone.pt.x - float(blocks[i].x0),
                                    toclone.pt.y - float(blocks[i].y0),
                                    toclone.size,
                                    toclone.angle,
                                    toclone.response,
                                    toclone.octave,
                                    toclone.class_id
                                );
                                kp_pour_b[i].push_back(curkp);
                            }
                        }
                    }*/

                    int nb_max = 0;
                    for (int i = 0; i < blocks.size(); i++) {
                        int v = kp_pour_b[i].size();
                        if (v > nb_max) {
                            nb_max = v;
                        }
                    }
                    std::cout << "NOMBRE MAXIMUM DE POINTS D'INTERETS DANS UN SEUL BLOCK : " << nb_max << std::endl;
                }

                //PHASE DE CALCUL DE L'IMAGE MOSAIQUE
                {
                    std::cout << "COMPARAISON ET CHOIX DES IMAGETTES BLOC PAR BLOC" << std::endl;
                    std::cout << "RECUPERATION DES MOYENNES" << std::endl;

                    std::vector<int> utilisee;

                    get_moy(imagettes.size(), moyennes);
                    for (size_t i = 0; i < moyennes.size(); i++) {
                        utilisee.push_back(0);
                    }
                    std::cout << " => MOYENNES : " << moyennes.size() << std::endl;
                    std::cout << " => UTILISEE : " << utilisee.size() << std::endl;


                    std::vector<Mat> imagette_descriptors;
                    std::vector<std::vector<KeyPoint>> imagette_kps;
                    imagette_descriptors.resize(imagettes.size());
                    imagette_kps.resize(imagettes.size());

                    for (int i = 0; i < imagettes.size(); i++) {
                        std::cout << "RECUPERATION DES DECRIPTEURS ET KP : " << int(((i + 1) / float(imagettes.size())) * 100) << " %\r" << std::flush;
                        std::string name_img = nom_imagette_sift[i];
                        //recuperer_keypoint(name_img, imagette_kps[i]);
                        recuperer_donnees_sift_kp(name_img, imagette_kps[i]);
                        recuperer_descriptor(name_img, imagette_descriptors[i]);
                    }

                    /*for (int i = 0; i < imagette_kps.size(); i++) {
                        std::cout << "TAILLE IMAGETTE KP : " << imagette_kps[i].size() << std::endl;
                    }*/

                    int nb_sift_test = 0;
                    int nb_moy_test = 0;

                    std::cout << "\nDEBUT DE L'ALGO" << std::endl;
                    std::vector<KeyPoint> kp_block;
                    BLOCK* curB;



                    for (int i = 0; i < blocks.size(); i++) {
                        std::cout << " + AVANCEMENT : " << int(((i + 1) / float(blocks.size())) * 100) << " %\r" << std::flush;

                        kp_block = kp_pour_b[i];
                        curB = &blocks[i];

                        if (kp_block.size() < Seuil_Passage_moy) {
                            //Traitement moyenne
                            nb_moy_test++;
                            moy_block(curB, moyennes, utilisee);
                        }
                        else {
                            //Traitement SIFT
                            nb_sift_test++;
                            sift_block(curB, kp_block, imagettes, imagette_kps, imagette_descriptors, utilisee);
                        }
                    }

                    std::cout << "\nFIN SIFT : " << nb_sift_test << std::endl;
                    std::cout << "FIN MOY : " << nb_moy_test << std::endl;
                }
            }
            else {

                std::cout << "COMPARAISON ET CHOIX DES IMAGETTES BLOC PAR BLOC" << std::endl;
                std::cout << "RECUPERATION DES MOYENNES" << std::endl;

                std::vector<int> utilisee;
                get_moy(imagettes.size(), moyennes);
                for (size_t i = 0; i < moyennes.size(); i++) {
                    utilisee.push_back(0);
                }
                std::cout << " => MOYENNES : " << moyennes.size() << std::endl;
                std::cout << " => UTILISEE : " << utilisee.size() << std::endl;

                for (int i = 0; i < blocks.size(); i++) {
                    std::cout << " + AVANCEMENT : " << int(((i + 1) / float(blocks.size())) * 100) << " %\r" << std::flush;
                    BLOCK* curB = &blocks[i];
                    moy_block(curB, moyennes, utilisee);
                }
                std::cout << "\nFIN de l'algo" << std::endl;
            }
        }


        std::cout << "=================================================" << std::endl;

        //RECONSTRUCTION DE L'IMAGE MOSAIQUE A L'AIDE DE TOUTES LES IMAGETTES SELECTIONNEES
        {
            std::cout << "RECONSTRUCTION DE L'IMAGE MOSAIQUE" << std::endl;
            Mat img_mosaique = Mat::zeros(originalImage_gray.size(), CV_8U);

            float nb_sift = 0;
            float nb_moy = 0;

            for (int k = 0; k < blocks.size(); k++) {
                std::cout << " + AVANCEMENT : " << int(((k + 1) / float(blocks.size())) * 100) << " %\r" << std::flush;
                BLOCK curB = blocks[k];

                if (curB.sift_imagette) {
                    nb_sift++;
                }
                else {
                    nb_moy++;
                }

                Mat resized_blk;
                //INTER_LINEAR 2x2 voisinage
                //INTER_CUBIC 4x4 voisinnage
                //INTER_NEAREST nearest neighbour 1 voisinnage
                resize(imagettes[curB.imagette], resized_blk, Size(curB.img_data.cols, curB.img_data.rows), INTER_CUBIC);

                //On ecrit le resultat
                for (int i = 0; i < resized_blk.rows; i++) {
                    for (int j = 0; j < resized_blk.cols; j++) {
                        img_mosaique.at<uchar>(
                            curB.y0 + i,
                            curB.x0 + j
                            ) = resized_blk.at<uchar>(i, j);
                    }
                }
                resized_blk.deallocate();
            }

            imwrite((dir_name_final_result + "/mosaique" + std::to_string(nbBlockX) + "X" + std::to_string(nbBlockY) + "INTERCUBIC.jpg").c_str(), img_mosaique);

            std::cout << "\n == PROPORTIONS : " << std::endl;
            std::cout << "    -> SIFT : " << nb_sift << " (" << (nb_sift / float(blocks.size())) * 100.0f << "%)" << std::endl;
            std::cout << "    -> AUTRE TRAITEMENT : " << nb_moy << " (" << (nb_moy / float(blocks.size())) * 100.0f << "%)" << std::endl;

            if (AFFICHAGE_FINAL) {
                affichage("IMAGE points d'intérets", res_kp_color);
                affichage("IMAGE initiale", originalImage_color);
                affichage("IMAGE initiale (gray)", originalImage_gray);
                affichage("IMAGE finale", img_mosaique);
            }
        }

        if (affichage_img) {
            waitKey(0);
        }

        //Affichage d'images additionnel
        if (AFFICHAGE_MATCHING) {
            std::cout << "SHOW MATCHING : FLANNBASED" << std::endl;
            showMatching("Comparaison : FLANNBASED", "./images/Facile.ppm", "./images/Difficile.ppm", DescriptorMatcher::FLANNBASED);
            std::cout << "SHOW MATCHING : BRUTEFORCE" << std::endl;
            showMatching("Comparaison : BRUTEFORCE", "./images/Facile.ppm", "./images/Difficile.ppm", DescriptorMatcher::BRUTEFORCE);
            std::cout << "SHOW MATCHING : BRUTEFORCE_L1" << std::endl;
            showMatching("Comparaison : BRUTEFORCE_L1", "./images/baboon.jpg", "./images/ImgB.jpg", DescriptorMatcher::BRUTEFORCE_L1);

            //Only for binaries value descriptors (ORB, FREAK ...)
            //std::cout << "SHOW MATCHING : BRUTEFORCE_HAMMING" << std::endl;
            //showMatching("Comparaison : BRUTEFORCE_HAMMING","./images/baboon.jpg", "./images/ImgB.jpg",DescriptorMatcher::BRUTEFORCE_HAMMING);
            //std::cout << "SHOW MATCHING : BRUTEFORCE_HAMMINGLUT" << std::endl;
            //showMatching("Comparaison : BRUTEFORCE_HAMMINGLUT","./images/baboon.jpg", "./images/ImgB.jpg",DescriptorMatcher::BRUTEFORCE_HAMMINGLUT);
            std::cout << "SHOW MATCHING : BRUTEFORCE_SL2" << std::endl;
            showMatching("Comparaison : BRUTEFORCE_SL2", "./images/baboon.jpg", "./images/ImgB.jpg", DescriptorMatcher::BRUTEFORCE_SL2);
        }
    }

   if (affichage_img) {
       waitKey();
   }



   int taille_theorique_block = wblock * hblock;
   if (COMPUTE_DISTANCE) {
       std::vector<Mat> reshaped_img_dist;
       reshaped_img_dist.resize(imagettes.size());
       {
           std::cout << "=== DEBUT ALGORITHME DISTANCE ===" << std::endl;
           int cpt_moy = 0;
           int cpt_dist = 0;

           if(moyennes.size() == 0) get_moy(imagettes.size(), moyennes);

           std::vector<int> utilisee;
           for (size_t i = 0; i < imagettes.size(); i++) {
               utilisee.push_back(0);
           }

           std::cout << "=== RESHAPE DE TOUTES LES IMAGETTES ===" << std::endl;
           for (int i = 0; i < imagettes.size(); i++) {
               std::cout << " + AVANCEMENT : " << int(((i + 1) / float(imagettes.size())) * 100) << " %\r" << std::flush;
               resize(imagettes[i], reshaped_img_dist[i], Size(wblock, hblock), INTER_LINEAR);
           }

           std::cout << "\n=== DEBUT DE L'ALGO ===" << std::endl;
           for (size_t i = 0; i < blocks.size(); i++) {
               std::cout << " + AVANCEMENT : " << int(((i + 1) / float(blocks.size())) * 100) << " %\r" << std::flush;
               BLOCK* curB = &blocks[i];
               int taille_b = curB->img_data.cols * curB->img_data.rows;

               if (taille_b < taille_theorique_block) {
                   //Traitement moyenne
                   moy_block(curB, moyennes, utilisee);
                   cpt_moy++;
               }
               else {
                   //Traitement distance
                   //Parcourir les images et trouver celle avec une distance la plus petite
                   float valeurPlusProche = distance(curB->img_data, reshaped_img_dist[0]);
                   int indicePlusProche = 0;

                   for (int k = 1; k < reshaped_img_dist.size(); k++) {
                       float valeur = distance(curB->img_data, reshaped_img_dist[k]);
                       if (valeur < valeurPlusProche) {
                           valeurPlusProche = valeur;
                           indicePlusProche = k;
                       }
                       else {
                           if (valeur < valeurPlusProche + seuil_imagette_distance && utilisee[k] < utilisee[indicePlusProche]) {
                               valeurPlusProche = valeur;
                               indicePlusProche = k;
                           }
                       }
                       if (valeurPlusProche < seuil_arret_distance) {
                           break;
                       }
                   }
                   cpt_dist++;
               }
           }

           std::cout << "\n=== FIN DE L'ALGO ===" << std::endl;
           std::cout << "RECONSTRUCTION DE L'IMAGE MOSAIQUE PAR DISTANCE ET EGALISATION" << std::endl;
           Mat img_mosaique = Mat::zeros(originalImage_gray.size(), CV_8SC1);

           for (int k = 0; k < blocks.size(); k++) {
               std::cout << " + AVANCEMENT : " << int(((k + 1) / float(blocks.size())) * 100) << " %\r" << std::flush;
               BLOCK curB = blocks[k];


               Mat resized_blk;
               resize(imagettes[curB.imagette], resized_blk, Size(curB.img_data.cols, curB.img_data.rows), INTER_LINEAR);

               //On ecrit le resultat
               for (int i = 0; i < resized_blk.rows; i++) {
                   for (int j = 0; j < resized_blk.cols; j++) {
                       img_mosaique.at<uchar>(
                           curB.y0 + i,
                           curB.x0 + j
                           ) = resized_blk.at<uchar>(i, j);
                   }
               }

           }

           imwrite((dir_name_final_result + "/mosaique" + std::to_string(nbBlockX) + "X" + std::to_string(nbBlockY) + "dist.jpg").c_str(), img_mosaique);
       }
   }
   

   if (COMPUTE_EGALISATION) {
       //Calcul de l'image avec specification d'histogramme
       const std::vector<std::string> names;
       std::vector<Mat> egalisation;
       std::vector<Mat> egalisation_blured;

       std::vector<Mat> reshaped_img;
       reshaped_img.resize(imagettes.size());
       {
           std::cout << "=== DEBUT ALGORITHME SPECIFICATION HISTOGRAMME ===" << std::endl;
           std::cout << "=== RECUPERATION EGALISATION ===" << std::endl;
           get_egalisation(nom_imagette_sift, egalisation, egalisation_blured);

           std::cout << "\nEGALISATION :" << egalisation.size() << std::endl;
           std::cout << "EGALISATION BLURED : " << egalisation_blured.size() << std::endl;

           int cpt_moy = 0;
           int cpt_egal = 0;

           if (moyennes.size() == 0) get_moy(imagettes.size(), moyennes);
           std::vector<int> utilisee;
           for (size_t i = 0; i < imagettes.size(); i++) {
               utilisee.push_back(0);
           }

           std::cout << "=== RESHAPE DE TOUTES LES IMAGETTES ===" << std::endl;
           for (int i = 0; i < imagettes.size(); i++) {
               std::cout << " + AVANCEMENT : " << int(((i + 1) / float(imagettes.size())) * 100) << " %\r" << std::flush;
               //resize(egalisation[i], reshapped_img[i], Size(wblock, hblock), INTER_LINEAR);
               resize(egalisation_blured[i], reshaped_img[i], Size(wblock, hblock), INTER_LINEAR);
           }

           std::cout << "\n=== DEBUT DE L'ALGO ===" << std::endl;
           for (size_t i = 0; i < blocks.size(); i++) {
               std::cout << " + AVANCEMENT : " << int(((i + 1) / float(blocks.size())) * 100) << " %\r" << std::flush;
               BLOCK* curB = &blocks[i];
               int taille_b = curB->img_data.cols * curB->img_data.rows;

               if (taille_b < taille_theorique_block) {
                   //Traitement moyenne
                   moy_block(curB, moyennes, utilisee);
                   cpt_moy++;
               }
               else {
                   //Traitement égalisation
                   curB->histo_data = Mat::zeros(curB->img_data.size(), CV_8U);
                   dist_block(curB, reshaped_img, utilisee);
                   cpt_egal++;
               }
           }
           std::cout << "\n=== FIN DE L'ALGO ===" << std::endl;
           std::cout << "NB egal : " << cpt_egal << std::endl;
           std::cout << "NB Moy : " << cpt_moy << std::endl;


           float nb_egal_test = 0;
           float nb_moy_test = 0;

           std::cout << "RECONSTRUCTION DE L'IMAGE MOSAIQUE PAR DISTANCE ET EGALISATION" << std::endl;
           Mat img_mosaique = Mat::zeros(originalImage_gray.size(), CV_8SC1);

           for (int k = 0; k < blocks.size(); k++) {
               std::cout << " + AVANCEMENT : " << int(((k + 1) / float(blocks.size())) * 100) << " %\r" << std::flush;
               BLOCK curB = blocks[k];

               if (curB.egal_histo) {
                   nb_egal_test++;
                   for (int i = 0; i < curB.histo_data.rows; i++) {
                       for (int j = 0; j < curB.histo_data.cols; j++) {
                           img_mosaique.at<uchar>(
                               curB.y0 + i,
                               curB.x0 + j
                               ) = curB.histo_data.at<uchar>(i, j);
                       }
                   }
               }
               else {
                   nb_moy_test++;
                   Mat resized_blk;
                   resize(imagettes[curB.imagette], resized_blk, Size(curB.img_data.cols, curB.img_data.rows), INTER_LINEAR);

                   //On ecrit le resultat
                   for (int i = 0; i < resized_blk.rows; i++) {
                       for (int j = 0; j < resized_blk.cols; j++) {
                           img_mosaique.at<uchar>(
                               curB.y0 + i,
                               curB.x0 + j
                               ) = resized_blk.at<uchar>(i, j);
                       }
                   }
               }
           }

           imwrite((dir_name_final_result + "/mosaique" + std::to_string(nbBlockX) + "X" + std::to_string(nbBlockY) + "EGAL_dist.jpg").c_str(), img_mosaique);

           std::cout << "\n == PROPORTIONS : " << std::endl;
           std::cout << "    -> EGALISATION : " << nb_egal_test << " (" << (nb_egal_test / float(blocks.size())) * 100.0f << "%)" << std::endl;
           std::cout << "    -> MOY : " << nb_moy_test << " (" << (nb_moy_test / float(blocks.size())) * 100.0f << "%)" << std::endl;
       }
   }

   if (COMPUTE_QUAD_TREE) {
       std::vector<BLOCK> res;
       std::cout << "RECONSTRUCTION DE L'IMAGE MOSAIQUE PAR Quad tree" << std::endl;

       callQuadTree(originalImage_gray, quad_seuil, quad_minSize,imagettes, res);

       Mat img_mosaique = Mat::zeros(originalImage_gray.size(), CV_8SC1);

       for (int k = 0; k < res.size(); k++) {
           std::cout << " + AVANCEMENT : " << int(((k + 1) / float(res.size())) * 100) << " %\r" << std::flush;
           BLOCK curB = res[k];

               //On ecrit le resultat
               for (int i = 0; i < curB.img_data.rows; i++) {
                   for (int j = 0; j < curB.img_data.cols; j++) {
                       img_mosaique.at<uchar>(
                           curB.y0 + i,
                           curB.x0 + j
                           ) = curB.img_data.at<uchar>(i, j);
                   }
               }
       }

       imwrite((dir_name_final_result + "/mosaique" + std::to_string(nbBlockX) + "X" + std::to_string(nbBlockY) + "Quad" + std::to_string(quad_minSize) + "_" + std::to_string(quad_seuil) + ".jpg").c_str(), img_mosaique);

   }


   if (COMPUTE_MOY_VAR) {

       std::vector<float> variances;
       std::cout << "RECONSTRUCTION DE L'IMAGE MOSAIQUE PAR Moyenne et variance" << std::endl;

       get_var(imagettes.size(), variances);
       if (moyennes.size() == 0) get_moy(imagettes.size(), moyennes);
       std::vector<int> utilisee;
       for (size_t i = 0; i < imagettes.size(); i++) {
           utilisee.push_back(0);
       }
       std::cout << " => MOYENNE : " << moyennes.size() << std::endl;
       std::cout << " => VARIANCE : " << variances.size() << std::endl;
       std::cout << " => UTILISEE : " << utilisee.size() << std::endl;


       for (int i = 0; i < blocks.size(); i++) {
           std::cout << " + AVANCEMENT : " << int(((i + 1) / float(blocks.size())) * 100) << " %\r" << std::flush;
           BLOCK* curB = &blocks[i];
           moy_variance_block(curB, moyennes, variances, utilisee);
       }
       std::cout << "\nFIN de l'algo" << std::endl;


       Mat img_mosaique = Mat::zeros(originalImage_gray.size(), CV_8SC1);

       std::cout << "CONSTRUCTION DE L'IMAGE MOSAIQUE" << std::endl;
       for (int k = 0; k < blocks.size(); k++) {
           std::cout << " + AVANCEMENT : " << int(((k + 1) / float(blocks.size())) * 100) << " %\r" << std::flush;
           BLOCK curB = blocks[k];

           Mat resized_blk;
           //INTER_LINEAR 2x2 voisinage
           //INTER_CUBIC 4x4 voisinnage
           //INTER_NEAREST nearest neighbour 1 voisinnage
           resize(imagettes[curB.imagette], resized_blk, Size(curB.img_data.cols, curB.img_data.rows), INTER_AREA);

           //On ecrit le resultat
           for (int i = 0; i < resized_blk.rows; i++) {
               for (int j = 0; j < resized_blk.cols; j++) {
                   img_mosaique.at<uchar>(
                       curB.y0 + i,
                       curB.x0 + j
                       ) = resized_blk.at<uchar>(i, j);
               }
           }
       }

       imwrite((dir_name_final_result + "/mosaique" + std::to_string(nbBlockX) + "X" + std::to_string(nbBlockY) + "Moy" + std::to_string(coef_moy) + "Var" + std::to_string(coef_var) + ".jpg").c_str(), img_mosaique);

   }

    std::cout << "=== FIN DE L'ALGORITHME ==="<< std::endl;
    
    return 0;
}
