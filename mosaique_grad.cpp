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

void computeGradAngle(const float* Img1, const float* Img2, float* GradImg1, float* GradImg2, size_t nW, size_t nH) {
    float ndg1, ndg1Gauche, ndg1Haut;
    float ndg2, ndg2Gauche, ndg2Haut;
    float vertical1, horizontal1;
    float vertical2, horizontal2;

    for (int i=0; i < nH; i++) {
        for (int j=0; j < nW; j++) {
            ndg1 = Img1[i*nW+j];
            ndg1Gauche = (j > 0) ? Img1[i*nW+j-1]   : Img1[i*nW+j];
            ndg1Haut   = (i > 0) ? Img1[(i-1)*nW+j] : Img1[i*nW+j];
            ndg2 = Img2[i*nW+j];
            ndg2Gauche = (j > 0) ? Img2[i*nW+j-1]   : Img2[i*nW+j];
            ndg2Haut   = (i > 0) ? Img2[(i-1)*nW+j] : Img2[i*nW+j];

            vertical1   = ndg1Haut   - ndg1;
            horizontal1 = ndg1Gauche - ndg1;
            vertical2   = ndg2Haut   - ndg2;
            horizontal2 = ndg2Gauche - ndg2;

            GradImg1[i*nW+j] = atan2(horizontal1, vertical1);
            GradImg2[i*nW+j] = atan2(horizontal2, vertical2);
        }
    }
}

float distanceGrad(const float* Img1, const float* Img2, size_t nW, size_t nH) {
    //distanceEQM entre les gradients de l'image
    size_t nTaille = nW*nH;
    float* GradImg1 = new float[nTaille];
    float* GradImg2 = new float[nTaille];
    computeGradAngle(Img1, Img2, GradImg1, GradImg2, nW, nH);
    return distanceEQM(GradImg1, GradImg2, nW, nH);
}

float distance(const float* Img1, const float* Img2, size_t nW, size_t nH) {
    // distanceEQM des gradients + distanceEQM des pixels
    return distanceGrad(Img1, Img2, nW, nH)*5000 + distanceEQM(Img1, Img2, nW, nH);
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
                    ImgOut[y*nWO+x] = moyenne[indicePlusProche][ligneImage*tailleBlockW + colonneImage];
                }
            }
        }
    }


    ecrire_image_pgm(cNomImgEcrite, ImgOut, nHO, nWO);
    free(ImgOr);
    return 1;
}
