# Fr-PPIChem
Scripts related to the Fr-PPIChem project for the PPI-HitProfiler models.

The dataset used to train the models is not available at the moment.

For more details about the methodology, see: https://pubs.acs.org/doi/10.1021/acschembio.0c00179

## QSAR
Knime workflow to train J48, Jrip and RF models. SVM R model is simply loaded here for simplicity.
Worklow also performs model testing and Y-scrambling

## QSAR_CV
Knime worfklow to perform 5-fold CV for J48, Jrip, RF and SVM models

## QSAR_zinc_prediction
Knime workflow to predict every chunck of the ZINC db

## QSAR_SVM_CV_radial
R script that performs grid search to tune the SVM model and train the tuned model

## SVM_radial_optimized.rds
Optimised SVM model (R object)

## variable_selection
Knime script to perform variable selection as described in: https://doi.org/10.1021/acs.jcim.7b00435

## LICENSE
All the files in this repository are covered by the licence in the file LICENSE.

Under the -BY clause, we request attribution for subsequent use of these files.

For publications using these files, the primary current citation is:
(article submitted)

