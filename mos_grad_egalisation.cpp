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
    float inv = 1.0/2.0;
    Finv[0] = FinvPre[0];
    for (int v = 1; v < 255; v++) {
        if (!set[v]) {
            FinvPre[v] = FinvPre[v-1];
        }
        Finv[v] = inv * (FinvPre[v-1] + FinvPre[v]);
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

int main(int argc, char* argv[])
{
    char cNomImgOrigine[250], cNomImgLue[250][250], cNomImgEcrite[250];
    int nbImgLue, nbI;

    if (argc < 5) {
        printf("Usage: ImageOrigine.pgm ImageOut.pgm nbI I1.pgm I2.pgm ...\n"); 
        exit (1) ;
    }

    sscanf (argv[1],"%s",cNomImgOrigine);
    sscanf (argv[2],"%s",cNomImgEcrite);
    sscanf (argv[3],"%d",&nbI);
    for (int i = 4; i < argc; i++) {
        sscanf(argv[i],"%s", cNomImgLue[i-4]);
    }


    nbImgLue = argc - 4;

    OCTET *ImgOr, *ImgIn, *ImgOut;
    int nH, nW, nTaille;

    int nHO, nWO, nTailleO;

    lire_nb_lignes_colonnes_image_pgm(cNomImgOrigine, &nHO, &nWO);
    nTailleO = nWO * nHO;
    
    allocation_tableau(ImgOr, OCTET, nTailleO);
    lire_image_pgm(cNomImgOrigine, ImgOr, nHO * nWO);

    allocation_tableau(ImgOut, OCTET, nTailleO);

    int tailleBlockW = nWO / nbI;
    int tailleBlockH = nHO / nbI;
    int tailleBlock = tailleBlockW * tailleBlockH;

    int startI, startJ, endI, endJ;

    float* moyenne = new float[nbImgLue*tailleBlock];
    float* moyenneImg = new float[tailleBlock];
    int tailleTailleBlockH;
    int tailleTailleBlockW;
    int tailleTailleBlock;
    int indBlock;

    // lecture images
    // 1. découpage des images d'entrés en blocks représentant les pixels de l'image de sortie (moyenne)
    for (int i = 0; i < nbImgLue; i++) {
        lire_nb_lignes_colonnes_image_pgm(cNomImgLue[i], &nH, &nW);
        nTaille = nH * nW;
        allocation_tableau(ImgIn, OCTET, nTaille);

        lire_image_pgm(cNomImgLue[i], ImgIn, nH * nW);

        tailleTailleBlockH = nH / tailleBlockH;
        tailleTailleBlockW = nW / tailleBlockW;
        tailleTailleBlock = tailleTailleBlockH * tailleTailleBlockW;

        for (int bi = 0; bi < tailleBlockH; bi++) {
            for (int bj = 0; bj < tailleBlockW; bj++) {
                startI = bi*tailleTailleBlockH;
                startJ = bj*tailleTailleBlockH;
                endI = startI + tailleTailleBlockH;
                endJ = startJ + tailleTailleBlockW;

                indBlock = bi*tailleBlockW+bj;

                moyenneImg[indBlock] = 0;
                for (int y = startI; y < endI; y++) {
                    for (int x = startJ; x < endJ; x++) {
                        moyenneImg[indBlock] += ImgIn[y*nW+x];
                    }
                }
                moyenneImg[indBlock] /= tailleTailleBlock;
            }
        }
        
        egaliser(moyenneImg, moyenneImg, tailleBlockW, tailleBlockH);
        blur(moyenneImg, &moyenne[i*tailleBlock], tailleBlockW, tailleBlockH);

        free(ImgIn);
    }

    float* blockLu = new float[tailleBlock];
    float* block = new float[tailleBlock];

    float valeur;
    float valeurPlusProche;
    int indicePlusProche;

    int ligneImage, colonneImage;

    // pour chaque block :
    for (int i = 0; i < nbI; i++) {
        for (int j = 0; j < nbI; j++) {
            startI = i*tailleBlockH;
            startJ = j*tailleBlockW;
            endI = startI + tailleBlockH;
            endJ = startJ + tailleBlockW;

            // mettre les valeurs de l'image d'origine dans le block
            for (int y = startI; y < endI; y++) {
                for (int x = startJ; x < endJ; x++) {
                    ligneImage   = y-startI;
                    colonneImage = x-startJ;
                    blockLu[ligneImage*tailleBlockW+colonneImage] = ImgOr[y*nWO+x];
                }
            }
            blur(blockLu, block, tailleBlockW, tailleBlockH);

            int histo[256] = {}; //Initialisé à 0
            float ddp[256];
            float F[256];
            float Finv[256];
            int min, max;

            computeHisto(block, histo, tailleBlockW, tailleBlockH);
            computeDdp(block, histo, ddp, tailleBlockW, tailleBlockH);
            computeF(block, ddp, F, tailleBlockW, tailleBlockH);
            computeInverse(F, Finv);
            computeMinMaxFromHisto(histo, &min, &max);

            applyF(block, F, block, tailleBlockW, tailleBlockH);

            // parcourir les images et trouver celle avec une distance la plus petite
            valeurPlusProche = distance(block, &moyenne[0], tailleBlockW, tailleBlockH);
            indicePlusProche = 0;
            for (int k = 1; k < nbImgLue; k++) {
                valeur = distance(block, &moyenne[k*tailleBlock], tailleBlockW, tailleBlockH);
                if (valeur < valeurPlusProche) {
                    valeurPlusProche = valeur;
                    indicePlusProche = k;
                }
            }

            float valeur;
            // copier l'image choisie dans l'image de sortie
            for (int y = startI; y < endI; y++) {
                for (int x = startJ; x < endJ; x++) {
                    ligneImage   = y-startI;
                    colonneImage = x-startJ;
                    valeur = moyenne[indicePlusProche * tailleBlock + ligneImage*tailleBlockW + colonneImage];
                    //ImgOut[y*nWO+x] = (int)(ImgOr[y*nWO+x]*0.3 + 0.7*clamp(255.0*Finv[(int)valeur], min, max));
                    ImgOut[y*nWO+x] = (int)clamp(255.0*Finv[(int)valeur], min, max);
                }
            }
        }
    }


    ecrire_image_pgm(cNomImgEcrite, ImgOut, nHO, nWO);
    free(ImgOr);
    return 1;
}
