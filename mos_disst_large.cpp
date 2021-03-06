#include <stdio.h>
#include <iostream>
#include "../lib/image_ppm.h"

#define maxNbLu = 1024

float distance(float* Img1, float* Img2, size_t nTaille) {
    //distance par erreur quadratique moyenne
    float sum = 0;
    for (int i=0; i < nTaille; i+=3) {
        sum += pow(Img1[i] - Img2[i], 2);
    }

    float EQM = sum / nTaille;
    //float PSNR = 10 * log10((255*255)/EQM);
    return EQM;
}

int main(int argc, char* argv[])
{
    char cNomImgOrigine[250], cNomImgLue[250][250], cNomImgEcrite[250];
    int nbImgLueTotal, nbImgLue, nbI;

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


    nbImgLueTotal = argc - 4;
    if (nbImgLueTotal > maxNbLu) {
        nbImgLue = maxNbLu;
    }
    else {
        nbImgLue = nbImgLueTotal;
    }

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

    for (int i = 0; i < nbI; i++) {
        for (int j = 0; j < nbI; j++) {
            startI = i*tailleBlockH;
            startJ = j*tailleBlockW;
            endI = startI + tailleBlockH;
            endJ = startJ + tailleBlockH;

            for (int y = startI; y < endI; y++) {
                for (int x = startJ; x < endJ; x++) {
                    ligneImage   = y-startI;
                    colonneImage = x-startJ;
                    block[ligneImage*tailleBlockW+colonneImage] = ImgOr[y*nWO+x];
                }
            }

            valeurPlusProche = distance(block, moyenne[0], tailleBlock);
            indicePlusProche = 0;
            for (int k = 1; k < nbImgLue; k++) {
                valeur = distance(block, moyenne[k], tailleBlock);
                if (valeur < valeurPlusProche) {
                    valeurPlusProche = valeur;
                    indicePlusProche = k;
                }
            }


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
