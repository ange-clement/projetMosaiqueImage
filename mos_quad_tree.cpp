#include <stdio.h>
#include <cmath>
#include <iostream>
#include "../lib/image_ppm.h"

float clamp(float val, float min, float max) {
    return (val < min) ? min : (val > max) ? max : val;
}

float distanceEQM(const float* Img1, const float* Img2, size_t nW, size_t nH) {
    //distance par erreur quadratique moyenne
    size_t nTaille = nW*nH;
    float sum = 0;
    for (int i=0; i < nTaille; i+=3) {
        sum += pow(Img1[i] - Img2[i], 2);
    }
    return sum / nTaille;;
}

float distanceGrad(const float* Img1, const float* Img2, size_t nW, size_t nH) {
    // On compare les gradiants globaux vers le haut et vers la droite :
    float img1Bottom, img1Top, img1Left, img1Right;
    float img2Bottom, img2Top, img2Left, img2Right;

    img1Bottom = img1Top = img1Left = img1Right = img2Bottom = img2Top = img2Left = img2Right = 0;
    for (int i=0; i < nH; i++) {
        for (int j=0; j < nW; j++) {
            if (i < nH/2) {
                img1Bottom += Img1[i*nW+j];
                img2Bottom += Img2[i*nW+j];
            }
            else {
                img1Top += Img1[i*nW+j];
                img2Top += Img2[i*nW+j];
            }

            if (j < nW/2) {
                img1Left += Img1[i*nW+j];
                img2Left += Img2[i*nW+j];
            }
            else {
                img1Right += Img1[i*nW+j];
                img2Right += Img2[i*nW+j];
            }
        }
    }
    float gradX1 = img1Top - img1Bottom;
    float gradY1 = img1Left - img1Right;
    float gradX2 = img2Top - img2Bottom;
    float gradY2 = img2Left - img2Right;

    return abs(gradX2 - gradX1) + abs(gradY2 - gradY1);
}

float distance(const float* Img1, const float* Img2, size_t nW, size_t nH) {
    return distanceEQM(Img1, Img2, nW, nH);
}

void computeHisto(const float* Img, int* histo, size_t nW, size_t nH) {
    for (int i=0; i < nH; i++) {
        for (int j=0; j < nW; j++) {
            histo[(int)Img[i*nW+j]]++;
        }
    }
}

void computeDdp(const float* Img, const int* histo, float* ddp, size_t nW, size_t nH) {
    int nTaille = nW*nH;
    for (int v = 0; v < 256; v++) {
        ddp[v] = (double)histo[v] / (double)nTaille;
    }
}

void computeF(const float* Img, const float* ddp, float* F, size_t nW, size_t nH) {
    F[0] = 0;

    for (int v = 1; v < 256; v++) {
        F[v] = F[v-1] + ddp[v];
    }
}

void computeInverse(const float* F, float* Finv) {
    double axisX = 1;
    double axisY = 1.0/255.0;
    double vectLength = sqrt(axisX*axisX + axisY*axisY);
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
        x2 = x + (xProj - x)*2.0;
        y2 = y + (yProj - y)*2.0;
        xFinal = x2;
        FinvPre[xFinal] = y2;
        set[xFinal] = true;
    }

    FinvPre[0] = 0;
    Finv[0] = FinvPre[0];
    for (int v = 1; v < 255; v++) {
        if (!set[v]) {
            FinvPre[v] = FinvPre[v-1];
        }
        Finv[v] = FinvPre[v];
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

void computeMinMaxFromHisto(const int* histo, int * min, int * max) {
    *min = computeMinFromHisto(histo);
    *max = computeMaxFromHisto(histo);
}

void applyF(const float* ImgIn, const float* F, float* ImgOut, size_t nW, size_t nH) {
    int nTaille = nW*nH;
    for (int i=0; i<nTaille; i++) {
        ImgOut[i] = clamp(255.0*F[(int)ImgIn[i]], 0, 255);
    }
}

void egaliser(const float* ImgIn, float* ImgOut, size_t nW, size_t nH) {
    int histo[256] = {}; //Initialisé à 0
    float ddp[256];
    float F[256];

    computeHisto(ImgIn, histo, nW, nH);
    computeDdp(ImgIn, histo, ddp, nW, nH);
    computeF(ImgIn, ddp, F, nW, nH);

    applyF(ImgIn, F, ImgOut, nW, nH);
}

void blur(const float* ImgIn, float* ImgOut, size_t nW, size_t nH) {
    float coeffs[9] = {
        1.0/9.0, 1.0/9.0, 1.0/9.0,
        1.0/9.0, 1.0/9.0, 1.0/9.0,
        1.0/9.0, 1.0/9.0, 1.0/9.0
    };
    float h;
    int indI, indJ;
    for (int i=0; i < nH; i++) {
        for (int j=0; j < nW; j++) {
            h = 0;
            for (int m = -1; m < 2; m++) {
                for (int n = -1; n < 2; n++) {
                    indI = abs(i+m);
                    indJ = abs(j+n);
                    if (indI >= nH) indI = i-1;
                    if (indJ >= nW) indJ = j-1;
                    h += coeffs[(m+1)*3+(n+1)] * ImgIn[indI*nW + indJ];
                }
            }
            ImgOut[i*nW+j] = clamp(h, 0, 255);
        }
    }
}

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
        this->endI = startI+width;
        this->endJ = startJ+width;

        this->indexImgChoosen = -1;

        childs = new Block*[4];
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

void quadTree(Block* current, const OCTET* ImgIn, OCTET* ImgOut, int nbImgLue, char** cNomImgLue, int seuil, int minSize, int nW, int nH) {
    int startI = current->startI;
    int startJ = current->startJ;
    int endI = current->endI;
    int endJ = current->endJ;
    int bW = current->width;
    int bH = current->height;
    int bW2 = bW/2;
    int bH2 = bH/2;

    int bN = bW*bH;


    float* moyenne = new float[nbImgLue*bN];
    float* moyenneImg = new float[bN];
    float bbH, bbW, bbN;
    int bstartI, bstartJ;
    int bendI, bendJ;
    int indBlock;
    int nHLu, nWLu, nTaille;
    OCTET* ImgTmp;

    // lecture images
    // 1. découpage des images d'entrés en blocks représentant les pixels de l'image de sortie (moyenne)
    //TODO FONCTIONS
    for (int i = 0; i < nbImgLue; i++) {
        lire_nb_lignes_colonnes_image_pgm(cNomImgLue[i], &nHLu, &nWLu);
        nTaille = nHLu * nWLu;
        allocation_tableau(ImgTmp, OCTET, nTaille);

        lire_image_pgm(cNomImgLue[i], ImgTmp, nHLu * nWLu);

        bbH = nHLu / bH;
        bbW = nWLu / bW;
        bbN = bbH * bbW;

        for (int bi = 0; bi < bH; bi++) {
            for (int bj = 0; bj < bW; bj++) {
                indBlock = bi*bW+bj;
                bstartI = bi*bbH;
                bstartJ = bj*bbH;
                bendI = bstartI + bbH;
                bendJ = bstartJ + bbW;

                if (bbH >= 1 && bbW >= 1) {
                    moyenneImg[indBlock] = 0;
                    for (int y = bstartI; y < bendI; y++) {
                        for (int x = bstartJ; x < bendJ; x++) {
                            moyenneImg[indBlock] += ImgTmp[y*nWLu+x];
                        }
                    }
                    moyenneImg[indBlock] /= bbN;
                }
                else {
                    moyenneImg[indBlock] = ImgTmp[bstartI*nWLu+bstartJ];
                }
            }
        }
        egaliser(moyenneImg, moyenneImg, bW, bH);
        blur(moyenneImg, &moyenne[i*bN], bW, bH);

        free(ImgTmp);
    }
    free(moyenneImg);

    float* blockLu = new float[bN];
    float* block = new float[bN];

    float value;

    int ligneImage, colonneImage;
    
    // mettre les valeurs de l'image d'origine dans le block
    for (int y = startI; y < endI; y++) {
        for (int x = startJ; x < endJ; x++) {
            ligneImage   = y-startI;
            colonneImage = x-startJ;
            blockLu[ligneImage*bW+colonneImage] = ImgIn[y*nW+x];
        }
    }
    blur(blockLu, block, bW, bH);
    free(blockLu);


    int histo[256] = {}; //Initialisé à 0
    float ddp[256];
    float F[256];
    float Finv[256];
    int min, max;

    computeHisto(block, histo, bW, bH);
    computeDdp(block, histo, ddp, bW, bH);
    computeF(block, ddp, F, bW, bH);
    computeInverse(F, Finv);
    computeMinMaxFromHisto(histo, &min, &max);

    applyF(block, F, block, bW, bH);

    // parcourir les images et trouver celle avec une distance la plus petite
    value = distance(block, &moyenne[0], bW, bH);
    int indicePlusProche = 0;
    float valeur;
    for (int k = 1; k < nbImgLue; k++) {
        valeur = distance(block, &moyenne[k*bN], bW, bH);
        if (valeur < value) {
            value = valeur;
            indicePlusProche = k;
        }
    }
    free(block);

    std::cout << "VALUE " << value << std::endl;

    if (value < seuil || bH <= minSize || bW <= minSize) {
        std::cout << "COPIE " << startI << " " << startJ << std::endl;
        // copier l'image choisie dans l'image de sortie
        for (int y = startI; y < endI; y++) {
            for (int x = startJ; x < endJ; x++) {
                ligneImage   = y-startI;
                colonneImage = x-startJ;
                valeur = moyenne[indicePlusProche * bN + ligneImage * bW + colonneImage];
                ImgOut[y*nW+x] = (int)clamp(255.0*Finv[(int)valeur], min, max);
            }
        }
    }
    else {
        std::cout << "RECURSION " << startI << " " << startJ << std::endl;
        Block* HG = new Block(startI    , startJ    , bW2, bH2);
        Block* HD = new Block(startI    , startJ+bW2, bW2, bH2);
        Block* BG = new Block(startI+bH2, startJ    , bW2, bH2);
        Block* BD = new Block(startI+bH2, startJ+bW2, bW2, bH2);
        
        current->setChild(0, HG);
        current->setChild(1, HD);
        current->setChild(2, BG);
        current->setChild(3, BD);

        quadTree(HG, ImgIn, ImgOut, nbImgLue, cNomImgLue, seuil, minSize, nW, nH);
        quadTree(HD, ImgIn, ImgOut, nbImgLue, cNomImgLue, seuil, minSize, nW, nH);
        quadTree(BG, ImgIn, ImgOut, nbImgLue, cNomImgLue, seuil, minSize, nW, nH);
        quadTree(BD, ImgIn, ImgOut, nbImgLue, cNomImgLue, seuil, minSize, nW, nH);
    }
    free(moyenne);
}

void callQuadTree(const OCTET* ImgIn, OCTET* ImgOut, int nbImgLue, char** cNomImgLue, int seuil, int minSize, int nW, int nH) {
    Block* current = new Block(0, 0, nW, nH);
    int bW2 = nW/2;
    int bH2 = nH/2;
    
    Block* HG = new Block(0    , 0    , bW2, bH2);
    Block* HD = new Block(0    , 0+bW2, bW2, bH2);
    Block* BG = new Block(0+bH2, 0    , bW2, bH2);
    Block* BD = new Block(0+bH2, 0+bW2, bW2, bH2);

    current->setChild(0, HG);
    current->setChild(1, HD);
    current->setChild(2, BG);
    current->setChild(3, BD);

    quadTree(HG, ImgIn, ImgOut, nbImgLue, cNomImgLue, seuil, minSize, nW, nH);
    quadTree(HD, ImgIn, ImgOut, nbImgLue, cNomImgLue, seuil, minSize, nW, nH);
    quadTree(BG, ImgIn, ImgOut, nbImgLue, cNomImgLue, seuil, minSize, nW, nH);
    quadTree(BD, ImgIn, ImgOut, nbImgLue, cNomImgLue, seuil, minSize, nW, nH);
}

int main(int argc, char* argv[])
{
    char cNomImgOrigine[250], cNomImgEcrite[250];
    char** cNomImgLue;
    cNomImgLue = new char*[250];
    for (int i = 0; i < 250; i++) {
        cNomImgLue[i] = new char[250];
    }
    int nbImgLue, seuil, minSize;

    if (argc < 6) {
        printf("Usage: ImageOrigine.pgm ImageOut.pgm seuil minSize I1.pgm I2.pgm ...\n"); 
        exit (1) ;
    }

    sscanf (argv[1],"%s",cNomImgOrigine);
    sscanf (argv[2],"%s",cNomImgEcrite);
    sscanf (argv[3],"%d",&seuil);
    sscanf (argv[4],"%d",&minSize);
    int startIndexImage = 5;
    for (int i = startIndexImage; i < argc; i++) {
        sscanf(argv[i],"%s", cNomImgLue[i-startIndexImage]);
    }


    nbImgLue = argc - startIndexImage;

    OCTET *ImgOr, *ImgOut;

    int nHO, nWO, nTailleO;

    lire_nb_lignes_colonnes_image_pgm(cNomImgOrigine, &nHO, &nWO);
    nTailleO = nWO * nHO;
    
    allocation_tableau(ImgOr, OCTET, nTailleO);
    lire_image_pgm(cNomImgOrigine, ImgOr, nHO * nWO);

    allocation_tableau(ImgOut, OCTET, nTailleO);

    callQuadTree(ImgOr, ImgOut, nbImgLue, cNomImgLue, seuil, minSize, nWO, nHO);

    ecrire_image_pgm(cNomImgEcrite, ImgOut, nHO, nWO);
    free(ImgOr);
    return 1;
}
