#include <stdio.h>
#include <iostream>
#include "../lib/image_ppm.h"

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
    return distanceGrad(Img1, Img2, nW, nH) * distanceEQM(Img1, Img2, nW, nH);
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

    OCTET *ImgOr, *ImgIn[nbImgLue], *ImgOut;
    int nH[nbImgLue], nW[nbImgLue], nTaille[nbImgLue];

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

    float moyenne[nbImgLue][tailleBlock];
    int tailleTailleBlockH;
    int tailleTailleBlockW;
    int tailleTailleBlock;
    int indBlock;

    // lecture images
    // 1. découpage des images d'entrés en blocks représentant les pixels de l'image de sortie (moyenne)
    for (int i = 0; i < nbImgLue; i++) {
        lire_nb_lignes_colonnes_image_pgm(cNomImgLue[i], &nH[i], &nW[i]);
        nTaille[i] = nH[i] * nW[i];
        allocation_tableau(ImgIn[i], OCTET, nTaille[i]);

        lire_image_pgm(cNomImgLue[i], ImgIn[i], nH[i] * nW[i]);

        tailleTailleBlockH = nH[i] / tailleBlockH;
        tailleTailleBlockW = nW[i] / tailleBlockW;
        tailleTailleBlock = tailleTailleBlockH * tailleTailleBlockW;

        for (int bi = 0; bi < tailleBlockH; bi++) {
            for (int bj = 0; bj < tailleBlockW; bj++) {
                startI = bi*tailleTailleBlockH;
                startJ = bj*tailleTailleBlockH;
                endI = startI + tailleTailleBlockH;
                endJ = startJ + tailleTailleBlockW;

                indBlock = bi*tailleBlockW+bj;

                moyenne[i][indBlock] = 0;
                for (int y = startI; y < endI; y++) {
                    for (int x = startJ; x < endJ; x++) {
                        moyenne[i][indBlock] += ImgIn[i][y*nW[i]+x];
                    }
                }
                moyenne[i][indBlock] /= tailleTailleBlock;
            }
        }

        free(ImgIn[i]);
    }

    float block[tailleBlock];

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
                    block[ligneImage*tailleBlockW+colonneImage] = ImgOr[y*nWO+x];
                }
            }

            // parcourir les images et trouver celle avec une distance la plus petite
            valeurPlusProche = distance(block, moyenne[0], tailleBlockW, tailleBlockH);
            indicePlusProche = 0;
            for (int k = 1; k < nbImgLue; k++) {
                valeur = distance(block, moyenne[k], tailleBlockW, tailleBlockH);
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
                    ImgOut[y*nWO+x] = moyenne[indicePlusProche][ligneImage*tailleBlockW + colonneImage];
                }
            }
        }
    }


    ecrire_image_pgm(cNomImgEcrite, ImgOut, nHO, nWO);
    free(ImgOr);
    return 1;
}
