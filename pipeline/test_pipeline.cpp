#include "pipeline.h"

const char* img_depart = "baboon";
const char* img_bis = "fusionPipeline";

const char* img_res_fusion = "fusion";


//TEST DES ¨PIPELINES DE TRAITEMENT
int main(int argc, char* argv[]){
    std::cout << "==================================================" << std::endl;
    std::cout << "            ===== TEST PIPELINE =====            " << std::endl;
    std::cout << "==================================================" << std::endl;

    
    IMG* image_traitement = new IMG(img_depart);
    image_traitement->read();
    
    IMG* image_traitement_bis = new IMG(img_bis);
    image_traitement_bis->read();

    /*(new Pipeline("test_convolution",
        (new Traitement_image())
        ->ajouter_img_in(image_traitement)
        ->ajouter_traitement_img((new Traitement_convo(masque().moyenneur(10))))
        ->ajouter_traitement_img((new Traitement_convo(masque().moyenneur(10)))->set_flg(FLG_MIROR))
        ->ajouter_traitement_img((new Traitement_mediane(10,1)))
        ->ajouter_traitement_img((new Traitement_mediane(10,10))->set_flg(FLG_MIROR))

        ->ajouter_traitement_img((new Traitement_convo(masque().gaussian(5, 0.01))))
        ->ajouter_traitement_img((new Traitement_convo(masque().gaussian(5, 0.01)))->set_flg(FLG_MIROR))

        ->ajouter_traitement_img((new Traitement_convo(masque().gaussian(5))))
        ->ajouter_traitement_img((new Traitement_convo(masque().gaussian(5)))->set_flg(FLG_MIROR))

        ->ajouter_traitement_img((new Traitement_convo(masque().binomial_newton_ordre())))
        ->ajouter_traitement_img((new Traitement_convo(masque().sobel_prewit_y(false, 2.0))))
        ->ajouter_traitement_img((new Traitement_convo(masque().sobel_prewit_x(false, 2.0)))->set_flg(FLG_MIROR))

        ->ajouter_traitement_img((new Traitement_convo(masque().sobel_prewit_y(true, 2.0)))->set_flg(FLG_MIROR))
        ->ajouter_traitement_img((new Traitement_convo(masque().sobel_prewit_x(true, 2.0))))

        ->ajouter_traitement_img((new Traitement_convo(masque().laplacien_4connexes())))
        ->ajouter_traitement_img((new Traitement_convo(masque().laplacien_8connexes())))
        ->ajouter_traitement_img((new Traitement_convo(masque().kirsh())))
        ->ajouter_traitement_img((new Traitement_convo(masque().robinson())))
        ->ajouter_traitement_img((new Traitement_convo(masque().robinson_laplacien())))

    ))->execute();*/

    /*(new Pipeline("test_gradient",
        (new Traitement_image())
        ->ajouter_img_in(image_traitement)
        ->ajouter_traitement_img((new Traitement_gradient(FLG_GRADIENT_RES_BOTH))->set_ecriture_res(FLG_SAVEIMG_NO_WRITE)
        ->ajouter_traitement_img((new Traitement_hessien(FLG_HESSIEN_DROITE_HAUT))->set_ecriture_res(FLG_SAVEIMG_WRITE_REGROUPED))
    )))->execute();*/

    /*(new Pipeline("test_autre",
        (new Traitement_image())
        ->ajouter_img_in(image_traitement)
        ->ajouter_traitement_img((new Traitement_convo(masque().moyenneur(10)))->set_flg_mixing(FLG_MIXING)->set_ecriture_res(FLG_SAVEIMG_WRITE_SEPARATED)
            ->ajouter_traitement_img((new Traitement_soustraction()))
        )
        ->ajouter_traitement_img((new Traitement_uniforme(10))->set_ecriture_res(FLG_SAVEIMG_WRITE_SEPARATED))       
        ->ajouter_traitement_img((new Traitement_seuillage(128,FLG_SEUILLAGE_SUP))
            ->ajouter_traitement_img((new Traitement_erosiondilatation(0,5,FLG_DILATATION,FLG_EROSIONDILATATION_ORDRE_DEUX)))
            ->ajouter_traitement_img((new Traitement_erosiondilatation(255,5,FLG_EROSION,FLG_EROSIONDILATATION_ORDRE_DEUX)))
            ->ajouter_traitement_img((new Traitement_erosiondilatation(0,5,FLG_DILATATION,FLG_EROSIONDILATATION_CROIX)))
            ->ajouter_traitement_img((new Traitement_erosiondilatation(255,5,FLG_EROSION,FLG_EROSIONDILATATION_CROIX)))
        )
    ))->execute();*/

    /*(new Pipeline("test_seuil",
        (new Traitement_image())
        ->ajouter_img_in(image_traitement)
        ->ajouter_traitement_img((new Traitement_seuillage(128,FLG_SEUILLAGE_INF)))
        ->ajouter_traitement_img((new Traitement_seuillage(128,FLG_SEUILLAGE_SUP)))

        ->ajouter_traitement_img((new Traitement_hysteresis(50,150)))
        ->ajouter_traitement_img((new Traitement_hysteresis(50,150,255,FLG_HISTERESIS_4V)))
        ->ajouter_traitement_img((new Traitement_hysteresis(50,150,0)))
        ->ajouter_traitement_img((new Traitement_hysteresis(50,150,0,FLG_HISTERESIS_4V)))

        ->ajouter_traitement_img((new Traitement_seuillageExtremums(50,150)))

        ->ajouter_traitement_img((new Traitement_seuillageMoy(FLG_SEUILLAGE_SUP)))
        ->ajouter_traitement_img((new Traitement_seuillageMoy(FLG_SEUILLAGE_INF)))

        ->ajouter_traitement_img((new Traitement_seuillageMedian(5)))
        ->ajouter_traitement_img((new Traitement_seuillageMedian(10)))
        ->ajouter_traitement_img((new Traitement_seuillageMedian(10))->ajouterProportion(0,.5))

        ->ajouter_traitement_img((new Traitement_seuillageMinLocaux(100)))
    ))->execute();*/

    /*(new Pipeline("test_histogramme",
        (new Traitement_image())
        ->ajouter_img_in(image_traitement)
        ->ajouter_traitement_img((new Traitement_egalisationHistogramme()))

        ->ajouter_traitement_img((new Traitement_expensionDynamique()))
        ->ajouter_traitement_img((new Traitement_specificationHistogramme(image_traitement_bis)))

    ))->execute();*/

    /*(new Pipeline("test_compression",
        (new Traitement_image())
        ->ajouter_img_in(image_traitement)
        ->ajouter_traitement_img((new Traitement_decorrelation()))
        ->ajouter_traitement_img((new Traitement_ondelette(1,FLG_ONDELETTE_RES)))
        ->ajouter_traitement_img((new Traitement_ondelette(1,FLG_ONDELETTE_4RES)))    
    ))->execute();*/

    /*(new Pipeline("test_contours",
        (new Traitement_image())
        ->ajouter_img_in(image_traitement)
        ->ajouter_traitement_img((new Traitement_gradient(FLG_GRADIENT_RES_CARTE))->set_ecriture_res(FLG_SAVEIMG_NO_WRITE)
            ->ajouter_traitement_img((new Traitement_hysteresis(0,5)))
            ->ajouter_traitement_img((new Traitement_hysteresis(5,10)))
            ->ajouter_traitement_img((new Traitement_hysteresis(10,15)))
            ->ajouter_traitement_img((new Traitement_hysteresis(15,20)))
            ->ajouter_traitement_img((new Traitement_hysteresis(20,25)))
        )
        
        ->ajouter_traitement_img((new Traitement_gradient(FLG_GRADIENT_RES_GRADIENT))
            ->set_ecriture_res(FLG_SAVEIMG_NO_WRITE)
            ->ajouter_traitement_img((new Traitement_convo(masque().laplacien_4connexes()))
                ->set_flg_mixing(FLG_MIXING)
                ->set_ecriture_res(FLG_SAVEIMG_WRITE_SEPARATED)
                ->ajouter_traitement_img((new Traitement_zeroPassingDetection(FLG_ZEROPASSING_4))
                    ->set_ecriture_res(FLG_SAVEIMG_WRITE_SEPARATED)
                    ->ajouter_traitement_img((new Traitement_hysteresis(5,10))
                        ->set_ecriture_res(FLG_SAVEIMG_WRITE_SEPARATED)
                    )
                )
                ->ajouter_traitement_img((new Traitement_zeroPassingDetection(FLG_ZEROPASSING_8))
                    ->set_ecriture_res(FLG_SAVEIMG_WRITE_SEPARATED)
                    ->ajouter_traitement_img((new Traitement_hysteresis(5,10))
                        ->set_ecriture_res(FLG_SAVEIMG_WRITE_SEPARATED)
                    )
                )
            )

            ->ajouter_traitement_img((new Traitement_convo(masque().laplacien_8connexes()))
                ->set_flg_mixing(FLG_MIXING)
                ->set_ecriture_res(FLG_SAVEIMG_WRITE_SEPARATED)
                ->ajouter_traitement_img((new Traitement_zeroPassingDetection(FLG_ZEROPASSING_4))
                    ->set_ecriture_res(FLG_SAVEIMG_WRITE_SEPARATED)
                    ->ajouter_traitement_img((new Traitement_hysteresis(5,10))
                        ->set_ecriture_res(FLG_SAVEIMG_WRITE_SEPARATED)
                    )
                )
                ->ajouter_traitement_img((new Traitement_zeroPassingDetection(FLG_ZEROPASSING_8))
                    ->set_ecriture_res(FLG_SAVEIMG_WRITE_SEPARATED)
                    ->ajouter_traitement_img((new Traitement_hysteresis(5,10))
                        ->set_ecriture_res(FLG_SAVEIMG_WRITE_SEPARATED)
                    )
                )
            )
        )
        ->ajouter_traitement_img((new Traitement_convo(masque().laplacien_8connexes())))

    ))->execute();*/

    Traitement_rossenfeld * rossenfled_1 = new Traitement_rossenfeld();
    Traitement_rossenfeld * rossenfled_2 = new Traitement_rossenfeld();
    Traitement_rossenfeld * rossenfled_3 = new Traitement_rossenfeld();
    Traitement_rossenfeld * rossenfled_4 = new Traitement_rossenfeld();

    (new Pipeline("test_pointsFixes",
        (new Traitement_image())
        ->ajouter_img_in(image_traitement)
        ->ajouter_traitement_img((new Traitement_gradient(FLG_GRADIENT_RES_GRADIENT))
            ->set_ecriture_res(FLG_SAVEIMG_NO_WRITE)
            ->ajouter_traitement_img((new Traitement_hessien(FLG_HESSIEN_GAUCHE_BAS))
                ->set_ecriture_res(FLG_SAVEIMG_WRITE_REGROUPED)
                ->ajouter_traitement_img(rossenfled_1)
            )
            ->ajouter_traitement_img(rossenfled_1)

            ->ajouter_traitement_img((new Traitement_hessien(FLG_HESSIEN_GAUCHE_HAUT))
                ->set_ecriture_res(FLG_SAVEIMG_WRITE_REGROUPED)
                ->ajouter_traitement_img(rossenfled_2)
            )
            ->ajouter_traitement_img(rossenfled_2)

            ->ajouter_traitement_img((new Traitement_hessien(FLG_HESSIEN_DROITE_BAS))
                ->set_ecriture_res(FLG_SAVEIMG_WRITE_REGROUPED)
                ->ajouter_traitement_img(rossenfled_3)
            )
            ->ajouter_traitement_img(rossenfled_3)

            ->ajouter_traitement_img((new Traitement_hessien(FLG_HESSIEN_DROITE_HAUT))
                ->set_ecriture_res(FLG_SAVEIMG_WRITE_REGROUPED)
                ->ajouter_traitement_img(rossenfled_4)
            )
            ->ajouter_traitement_img(rossenfled_4)
        )

        ->ajouter_traitement_img((new Traitement_gradient(FLG_GRADIENT_RES_GRADIENT))
            ->set_ecriture_res(FLG_SAVEIMG_NO_WRITE)
            ->ajouter_traitement_img((new Traitement_beaudet())
                ->ajouter_traitement_img((new Traitement_seuillage(250, FLG_SEUILLAGE_SUP)))
                ->ajouter_traitement_img((new Traitement_hysteresis(240,250)))
            )
        )
        ->ajouter_traitement_img((new Traitement_convo(masque().gaussian(5, 0.01)))
            ->set_ecriture_res(FLG_SAVEIMG_NO_WRITE)
            ->ajouter_traitement_img((new Traitement_HarrisStephan(.09))
                ->ajouter_traitement_img((new Traitement_seuillage(250, FLG_SEUILLAGE_SUP)))
                ->ajouter_traitement_img((new Traitement_hysteresis(240,250)))
            )
        )
    ))->execute();
}