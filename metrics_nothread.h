#ifndef H_METRICS
#define	H_METRICS

#include <opencv2/opencv.hpp>
#include <opencv2/core/utils/filesystem.hpp>
#include <opencv2/core/persistence.hpp>
#include <opencv2/features2d.hpp>

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <sys/stat.h>

#include <fstream>

using namespace cv;

std::string dir_metric          = "./metric";
std::string file_metric_moy     = dir_metric + "/moy.dat";
std::string file_metric_var     = dir_metric + "/var.dat";

//IMAGETTES

const int imagette_nfeatures = 10;
const int imagette_nOctaveLayers = 3;
const double imagette_contrastThreshold = 0.04;  //Entre 0.03 et 0.09
const double imagette_edgeThreshold = 10;
const double imagette_sigma = 1.6;



/**   ====================
 *     TRAITEMENT MOYENNE
 *    ====================**/
int calc_moy(std::vector<Mat> imgs) {
    int nbel = imgs.size();
    float* moyenne = new float[nbel] {};

    std::cout << "START" << std::endl;


    for (size_t l = 0; l < nbel; l++) {
        for (int i = 0; i < imgs[l].cols; i++) {
            for (int j = 0; j < imgs[l].rows; j++) {
                moyenne[l] += (int) imgs[l].at<uchar>(j, i);
            }   
        }

        float nTaille = imgs[l].cols * imgs[l].rows;
        moyenne[l] /= nTaille;
    }

    std::ofstream monFlux;
    monFlux.open(file_metric_moy);
    std::cout << "OUVERTURE DE : " << file_metric_moy.c_str() << std::endl;

    for (int i = 0; i < nbel; i++) {
        monFlux << i << " " << moyenne[i] << std::endl;
    }
    monFlux.close();
    return 0;
}

void get_moy(int taille, std::vector<float>& res) {
    res.resize(taille);
    std::ifstream monFlux;
    monFlux.open(file_metric_moy);
    std::cout << "OUVERTURE DE : " << file_metric_moy.c_str() << std::endl;

    for (int i = 0; i < taille; i++){
        float m;
        int ind;
        monFlux >>ind >> m;
        res[ind] = m;
    }
    monFlux.close();
}

int calc_var(std::vector<Mat> imgs) {
    int nbel = imgs.size();
    
    double* variance = new double[nbel] {};

    
    std::cout << "START" << std::endl;

    Mat el;
    for (size_t l = 0; l < nbel; l++) {
        double* nbElementsLus = new double[256]{};
        el = imgs[l];
        int nTaille = el.rows * el.cols;

        for (int i = 0; i < el.rows; i++) {
            for (int j = 0; j < el.cols; j++) {
                int val = (int)imgs[l].at<uchar>(i, j);
                nbElementsLus[val]++;
            }
        }

        double moyenne = 0.;
        for (int i = 0; i < 256; i++) {
            moyenne += (i * nbElementsLus[i]);
        }
        moyenne /= (double)nTaille;

        for (int i = 0; i < 256; i++) {
            variance[l] += (i * i * nbElementsLus[i]);
        }
        variance[l] = (variance[l] / (double)nTaille) - (moyenne * moyenne);
    }


    std::ofstream monFlux;
    monFlux.open(file_metric_var);
    std::cout << "OUVERTURE DE : " << file_metric_var.c_str() << std::endl;

    for (int i = 0; i < nbel; i++) {
        monFlux << i << " " << variance[i] << std::endl;
    }
    monFlux.close();
    return 0;
}

void get_var(int taille, std::vector<float>& res) {
    res.resize(taille);
    std::ifstream monFlux;
    monFlux.open(file_metric_var);
    std::cout << "OUVERTURE DE : " << file_metric_var.c_str() << std::endl;

    for (int i = 0; i < taille; i++) {
        double m;
        int ind;
        monFlux >> ind >> m;
        res[ind] = (float) m;
    }
    monFlux.close();
}

void get_egalisation(const std::vector<std::string> names, std::vector<Mat>& egalisation, std::vector<Mat>& egalisation_blured) {
    egalisation.clear();
    egalisation_blured.clear();

    for (int i = 0; i < names.size(); i++) {
        std::cout << " + AVANCEMENT : " << int(((i + 1) / float(names.size())) * 100) << " %\r" << std::flush;

        std::string read = dir_metric + "/" + names[i] + "/egalisation.pgm";
        egalisation.push_back(imread(read, IMREAD_ANYCOLOR | IMREAD_ANYDEPTH));
        
        read = dir_metric + "/" + names[i] + "/egalisation_blured.pgm";
        egalisation_blured.push_back(imread(read, IMREAD_ANYCOLOR | IMREAD_ANYDEPTH));
    }
}



/**   ====================
 *      TRAITEMENT SIFT
 *    ====================**/
void recuperer_keypoint(std::string id, std::vector<KeyPoint>& kps) {
    //std::cout << "RECUPERATION KP : " << (dir_metric + "/" + id + "/kps.dat").c_str() << std::endl;

    kps.clear();
    std::ifstream monFlux(dir_metric + "/" + id + "/kps.dat");

    if (!monFlux.is_open()) {
        std::cout << "ERREUR RECUPERATION DE : " << (dir_metric + "/" + id + "/kps.dat").c_str() << std::endl;
    }

    int nbel;
    monFlux >> nbel;

    for (int i = 0; i < nbel; i++) {
        float 	x;
        float 	y;
        float 	size;
        float 	angle;
        float 	response;
        int 	octave;
        int 	class_id;

        monFlux >> x >> y >> size >> angle >> response >> octave >> class_id;
        kps.push_back(KeyPoint(x, y, size, angle, response, octave, class_id));
    }
    monFlux.close();
}

void sauvegarde_keypoint(std::string id,std::vector<KeyPoint>& kps) {
    //std::cout << "SAUVEGARDE KP : " << id.c_str() << std::endl;
    std::ofstream monFlux;
    monFlux.open(dir_metric+"/"+id + "/kps.dat");
    monFlux << kps.size() << std::endl;
    
    for (KeyPoint kp : kps) {
        float 	x = kp.pt.x;
        float 	y = kp.pt.y;
        float 	size = kp.size;
        float 	angle = kp.angle;
        float 	response = kp.response;
        int 	octave = kp.octave;
        int 	class_id = kp.class_id;
        
        monFlux << x <<" "<< y << " " << size << " " << angle << " " << response << " " << octave << " " << class_id << std::endl;
    }
    monFlux.close();
}

void recuperer_descriptor(std::string id, Mat& descriptor) {
    std::string read = dir_metric + "/" + id + "/descriptor.hdr";
    //std::cout << "READ DESCRIPTORS : " << read.c_str() << std::endl;
    descriptor = imread(read, IMREAD_ANYCOLOR | IMREAD_ANYDEPTH);
}

void sauvegarde_descriptor(std::string id, Mat& descriptor) {
    std::string save = dir_metric + "/" + id + "/descriptor.hdr";
    //std::cout << "SAVE DESCRIPTORS : " << save.c_str() << std::endl;
    if (!descriptor.empty()) {
        imwrite(save, descriptor);
    }
}

void sauvegarder_donnees_sift_kp(std::string id,const std::vector<KeyPoint> kps) {
    //std::cout  << "SAVE SKP : " << (dir_metric + "/" + id + "/siftKP.yml").c_str() << std::endl;
    FileStorage fs(dir_metric + "/" + id + "/siftKP.yml", FileStorage::WRITE);
    std::ostringstream oss;

    //std::cout << "KP : " << kps.size() << std::endl;

    for (size_t i = 0; i < kps.size(); ++i) {
        oss << "kp" << i;
        fs << oss.str() << kps[i];
        oss.str("");
    }
    fs.release();
}

void recuperer_donnees_sift_kp(std::string id, std::vector<KeyPoint>& kps) {
    kps.clear();
    FileStorage fs(dir_metric + "/" + id + "/siftKP.yml", FileStorage::READ);
    std::ostringstream oss;
    KeyPoint aKeypoint;
    int cpt = 0;
    while (1) {
        oss <<"kp" << cpt;
        fs[oss.str()] >> aKeypoint;
        if (fs[oss.str()].isNone() == 1) {
            break;
        }

        oss.str("");
        kps.push_back(aKeypoint);
        cpt++;
    }
    fs.release();
}

int calc_sift(std::vector<Mat> imgs, std::vector<std::string> names) {
    int nbel = imgs.size();

    std::cout << "START SIFT" << std::endl;

    float cpt = 0.0f;
    cv::Ptr<cv::SIFT>    detector;
    std::vector<KeyPoint> SIFT_keypoints;
    Mat descriptor;
    
    detector = cv::SIFT::create(
        imagette_nfeatures,
        imagette_nOctaveLayers,
        imagette_contrastThreshold,
        imagette_edgeThreshold,
        imagette_sigma
    );

    cv::Ptr<cv::DescriptorExtractor> extractor;
    extractor = cv::SiftDescriptorExtractor::create();

    for (size_t l = 0; l < nbel; l++) {
        std::cout << " + AVANCEMENT : " << int(((cpt+1) / float(nbel))*100) << " %\r" << std::flush;

        std::string nm = names[l];
        std::string dir = dir_metric + "/" + nm;
        struct stat buffer;
        if (stat(dir.c_str(), &buffer) != 0) {
            cv::utils::fs::createDirectories(dir.c_str());
        }

        Mat img = imgs[l];

        detector->detect(img, SIFT_keypoints);
        extractor->compute(img, SIFT_keypoints, descriptor);

        //Sauvegarde
        //sauvegarde_keypoint(nm, SIFT_keypoints);
        sauvegarder_donnees_sift_kp(nm, SIFT_keypoints);
        sauvegarde_descriptor(nm, descriptor);

        cpt++;
        SIFT_keypoints.clear();
        descriptor.release();
    }
    return 0;
}

void compute_kp_img(std::vector<Mat> imgs, std::vector<std::string> names) {
    int cpt = 0;
    for (std::string str : names) {
        std::vector<KeyPoint> kp;
        recuperer_keypoint(str, kp);
        Mat res;
        cv::drawKeypoints(imgs[cpt], kp, res, cv::Scalar::all(-1), cv::DrawMatchesFlags::DRAW_RICH_KEYPOINTS);
        
        std::string save = dir_metric + "/" + str + "/img_kp.jpg";
        //std::cout << "SAVE DESCRIPTORS : " << save.c_str() << std::endl;
        if (!res.empty()) {
            imwrite(save, res);
        }

        kp.clear();
        res.release();
        cpt++;
    }
}


/**   ====================
 *TRAITEMENT EGALISATION GAUSSIENNE
 *    ====================**/

float clamp(float val, float min, float max) {
    return (val < min) ? min : (val > max) ? max : val;
}

float distanceEQM(const cv::Mat& Img1, const cv::Mat& Img2) {
    //Distance par erreur quadratique moyenne
    float sum = 0;
    for (int i = 0; i < Img1.cols; i++){
        for (int j = 0; j < Img1.rows; j++) {
            sum += pow(Img1.at<uchar>(j, i) - Img2.at<uchar>(j, i), 2);
        }
    }
    
    return sum /(Img1.cols+ Img1.rows);
}

float distanceGrad(const cv::Mat& Img1, const cv::Mat& Img2) {
    // On compare les gradiants globaux vers le haut et vers la droite :
    float img1Bottom, img1Top, img1Left, img1Right;
    float img2Bottom, img2Top, img2Left, img2Right;

    img1Bottom = img1Top = img1Left = img1Right = img2Bottom = img2Top = img2Left = img2Right = 0;
    for (int i = 0; i < Img1.rows; i++) {
        for (int j = 0; j < Img1.cols; j++) {
            if (i < Img1.rows / 2) {
                img1Bottom += Img1.at<uchar>(i, j);
                img2Bottom += Img2.at<uchar>(i, j);
            }
            else {
                img1Top += Img1.at<uchar>(i, j);
                img2Top += Img2.at<uchar>(i, j);
            }

            if (j < Img1.cols / 2) {
                img1Left += Img1.at<uchar>(i, j);
                img2Left += Img2.at<uchar>(i, j);
            }
            else {
                img1Right += Img1.at<uchar>(i, j);
                img2Right += Img2.at<uchar>(i, j);
            }
        }
    }
    float gradX1 = img1Top - img1Bottom;
    float gradY1 = img1Left - img1Right;
    float gradX2 = img2Top - img2Bottom;
    float gradY2 = img2Left - img2Right;

    return abs(gradX2 - gradX1) + abs(gradY2 - gradY1);
}

float distance(const cv::Mat& img1,const cv::Mat& img2) {
    return distanceEQM(img1, img2);
}

void computeHisto(const cv::Mat& m, int* histo) {
    for (int i = 0; i < m.cols; i++) {
        for (int j = 0; j < m.rows; j++) {
            int val = m.at<uchar>(j, i);
            histo[val]++;
        }
    }
}

void computeDdp(const cv::Mat& Img, const int* histo, float* ddp) {
    int nTaille = Img.cols * Img.rows;
    for (int v = 0; v < 256; v++) {
        ddp[v] = (double)histo[v] / (double)nTaille;
    }
}

void computeF(const float* ddp, float* F) {
    F[0] = 0;
    for (int v = 1; v < 256; v++) {
        F[v] = F[v - 1] + ddp[v];
    }
}

void computeInverse(const float* F, float* Finv) {
    double axisX = 1;
    double axisY = 1.0 / 255.0;
    double vectLength = sqrt(axisX * axisX + axisY * axisY);
    axisX /= vectLength;
    axisY /= vectLength;
    double x, y, xProj, yProj, x2, y2, dot;
    int xFinal, yFinal;
    float FinvPre[255];
    bool set[255];
    for (int v = 0; v < 255; v++) {
        set[v] = false;
        FinvPre[v] = 0;
    }
    for (int v = 0; v < 255; v++) {
        x = v;
        y = F[v];
        dot = axisX * x + axisY * y;
        xProj = axisX * dot;
        yProj = axisY * dot;
        x2 = x + (xProj - x) * 2.0;
        y2 = y + (yProj - y) * 2.0;
        xFinal = x2;
        FinvPre[xFinal] = y2;
        set[xFinal] = true;
    }

    FinvPre[0] = 0;
    float inv = 1.0 / 2.0;
    Finv[0] = FinvPre[0];
    for (int v = 1; v < 255; v++) {
        if (!set[v]) {
            FinvPre[v] = FinvPre[v - 1];
        }
        Finv[v] = inv * (FinvPre[v - 1] + FinvPre[v]);
    }
}

int computeMaxFromHisto(const int* histo) {
    for (int i = 255; i >= 0; i--) {
        if (histo[i] > 0) {
            return i;
        }
    }
    return 0;
}

int computeMinFromHisto(const int* histo) {
    for (int i = 0; i < 255; i++) {
        if (histo[i] > 0) {
            return i;
        }
    }
    return 255;
}

void computeMinMaxFromHisto(const int* histo, int* min, int* max) {
    *min = computeMinFromHisto(histo);
    *max = computeMaxFromHisto(histo);
}

void applyF(const cv::Mat& ImgIn, const float* F,cv::Mat& ImgOut) {
    for (int i = 0; i < ImgIn.cols; i++) {
        for (int j = 0; j < ImgIn.rows; j++) {
            int val_img_in = ImgIn.at<uchar>(j, i);
            ImgOut.at<uchar>(j, i) = clamp(255.0 * F[val_img_in], 0, 255);
        }
    }
}

void egaliser(const cv::Mat& ImgIn, cv::Mat& ImgOut){
    int histo[256] = {};
    float ddp[256];
    float F[256];

    computeHisto(ImgIn, histo);
    computeDdp(ImgIn, histo, ddp);
    computeF(ddp, F);

    applyF(ImgIn, F, ImgOut);
}

void blurlin(const cv::Mat& ImgIn,cv::Mat& ImgOut) {
    float coeffs[9] = {
        1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
        1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
        1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0
    };
    float h;
    int indI, indJ;

    for (int i = 0; i < ImgIn.rows; i++) {
        for (int j = 0; j < ImgIn.cols; j++) {
            h = 0;
            for (int m = -1; m < 2; m++) {
                for (int n = -1; n < 2; n++) {
                    indI = abs(i + m);
                    indJ = abs(j + n);
                    if (indI >= ImgIn.rows) indI = i - 1;
                    if (indJ >= ImgIn.cols) indJ = j - 1;
                    h += coeffs[(m + 1) * 3 + (n + 1)] * ImgIn.at<uchar>(indI,indJ);
                }
            }
            ImgOut.at<uchar>(i, j) = clamp(h, 0, 255);
        }
    }
}

void calc_egalisation(std::vector<Mat> imgs, const std::vector<std::string> names) {
    int nbel = imgs.size();
    std::cout << "START" << std::endl;

    Mat curM;
    std::string curN;
    for (int l = 0; l < nbel; l++) {
        std::cout << " + AVANCEMENT : " << int(((l + 1) / float(nbel)) * 100) << " %\r" << std::flush;

        std::string nm = names[l];
        std::string dir = dir_metric + "/" + nm;
        struct stat buffer;
        if (stat(dir.c_str(), &buffer) != 0) {
            cv::utils::fs::createDirectories(dir.c_str());
        }


        curM = imgs[l];
        curN = names[l];

        Mat resToSave = Mat::zeros(curM.size(), curM.type());
        egaliser(curM, resToSave);

        std::string save = dir_metric + "/" + curN + "/egalisation.pgm";
        if (!resToSave.empty()) {
            imwrite(save, resToSave);
        }

        blurlin(resToSave, resToSave);
        save = dir_metric + "/" + curN + "/egalisation_blured.pgm";
        if (!resToSave.empty()) {
            imwrite(save, resToSave);
        }
    }
}

#endif // !H_METRICS
