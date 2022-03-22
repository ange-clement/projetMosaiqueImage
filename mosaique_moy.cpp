#include <stdio.h>
#include <iostream>
#include "../lib/image_ppm.h"

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
    float moyenne[nbImgLue];

    for (int i = 0; i < nbImgLue; i++) {

        lire_nb_lignes_colonnes_image_pgm(cNomImgLue[i], &nH[i], &nW[i]);
        nTaille[i] = nH[i] * nW[i];
        allocation_tableau(ImgIn[i], OCTET, nTaille[i]);

        lire_image_pgm(cNomImgLue[i], ImgIn[i], nH[i] * nW[i]);

        moyenne[i] = 0;
        for (int k = 0; k < nTaille[i]; k++) {
            moyenne[i] += ImgIn[i][k];
        }
        moyenne[i] /= nTaille[i];

    }

    int nHO, nWO, nTailleO;
    float moyenneO;

    lire_nb_lignes_colonnes_image_pgm(cNomImgOrigine, &nHO, &nWO);
    nTailleO = nWO * nHO;
    
    allocation_tableau(ImgOr, OCTET, nTailleO);
    lire_image_pgm(cNomImgOrigine, ImgOr, nHO * nWO);

    allocation_tableau(ImgOut, OCTET, nTailleO);

    int tailleBlockW = nWO / nbI;
    int tailleBlockH = nHO / nbI;
    int tailleBlock = tailleBlockW * tailleBlockH;

    int startI, startJ, endI, endJ;

    float moyennePlusProche;
    int indicePlusProche;

    int ligneImage, colonneImage;

    for (int i = 0; i < nbI; i++) {
        for (int j = 0; j < nbI; j++) {
            startI = i*tailleBlockH;
            startJ = j*tailleBlockW;
            endI = startI + tailleBlockH;
            endJ = startJ + tailleBlockH;

            moyenneO = 0;

            for (int y = startI; y < endI; y++) {
                for (int x = startJ; x < endJ; x++) {
                    moyenneO += ImgOr[y*nWO+x];
                }
            }

            moyenneO /= tailleBlock;

            moyennePlusProche = moyenne[0];
            indicePlusProche = 0;
            for (int k = 1; k < nbImgLue; k++) {
                if (abs(moyenne[k] - moyenneO) < abs(moyennePlusProche - moyenneO)) {
                    moyennePlusProche = moyenne[k];
                    indicePlusProche = k;
                }
            }


            // std::cout << indicePlusProche << " : " << cNomImgLue[indicePlusProche] << " " << nW[indicePlusProche] << " " << nH[indicePlusProche] << std::endl;
            // std::cout << startI << " " << endI << std::endl;
            // std::cout << startJ << " " << endJ << std::endl;
            // std::cout << tailleBlockH << " " << tailleBlockW << std::endl;
            // std::cout << (startI-startI) /(float) tailleBlockH * nH[indicePlusProche] << " " << (endI-startI) /(float) tailleBlockH * nH[indicePlusProche] << std::endl;
            // std::cout << (startJ-startJ) /(float) tailleBlockW * nW[indicePlusProche] << " " << (endJ-startJ) /(float) tailleBlockW * nW[indicePlusProche] << std::endl;


            for (int y = startI; y < endI; y++) {
                for (int x = startJ; x < endJ; x++) {
                    ligneImage   = (y-startI) /(float) tailleBlockH * nH[indicePlusProche];
                    colonneImage = (x-startJ) /(float) tailleBlockW * nW[indicePlusProche];
                    ImgOut[y*nWO+x] = ImgIn[indicePlusProche][ligneImage*nW[indicePlusProche] + colonneImage];
                }
            }
        }
    }


    ecrire_image_pgm(cNomImgEcrite, ImgOut, nHO, nWO);

    for (int i = 0; i < nbImgLue; i++) {
        free(ImgIn[i]);
    }
    free(ImgOr);
    return 1;
}
