#include <stdio.h>
#include <iostream>
#include "../lib/image_ppm.h"

float clamp(float val, float min, float max) {
    return (val < min) ? min : (val > max) ? max : val;
}

float distanceEQM(const float* Img1, const float* Img2, size_t nW, size_t nH) {
    size_t nTaille = nW*nH;
    float sum = 0;
    for (int i=0; i < nTaille; i+=3) {
        sum += pow(Img1[i] - Img2[i], 2);
    }
    return sum / nTaille;
}


float distance(const float* Img1, const float* Img2, size_t nW, size_t nH) {
    return distanceEQM(Img1, Img2, nW, nH);
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
        
        blur(moyenneImg, &moyenne[i*tailleBlock], tailleBlockW, tailleBlockH);
            if (i == 0) {
                int bi = tailleBlockH-1;
                    for (int bj = 0; bj < tailleBlockW; bj++) {
                        std::cout << bi*tailleBlockW + bj << " ";
                        std::cout << moyenne[i*tailleBlock + bi*tailleBlockW + bj] << " ";
                        std::cout << moyenneImg[bi*tailleBlockW + bj] << std::endl;
                    }
            }
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

            // copier l'image choisie dans l'image de sortie
            for (int y = startI; y < endI; y++) {
                for (int x = startJ; x < endJ; x++) {
                    ligneImage   = y-startI;
                    colonneImage = x-startJ;
                    //ImgOut[y*nWO+x] = ImgOr[y*nWO+x]*0.3 + moyenne[indicePlusProche][ligneImage*tailleBlockW + colonneImage]*0.7;
                    ImgOut[y*nWO+x] = moyenne[indicePlusProche * tailleBlock + ligneImage*tailleBlockW + colonneImage];
                }
            }
        }
    }


    ecrire_image_pgm(cNomImgEcrite, ImgOut, nHO, nWO);
    free(ImgOr);
    return 1;
}
