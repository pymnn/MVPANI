function [RANKED, WEIGHT] = reliefF( X, Y, K )
%   [RANKED,WEIGHT] = relieff(X,Y,K) computes ranks and weights of
%     attributes (predictors) for input data matrix X and response vector Y
%     using ReliefF algorithm for classification or RReliefF for regression
%     with K nearest neighbors. For classification, relieff uses K nearest
%     neighbors per class. RANKED are indices of columns in X ordered by
%     attribute importance, meaning RANKED(1) is the index of the most
%     important predictor. WEIGHT are attribute weights ranging from -1 to 1
%     with large positive weights assigned to important attributes.
%  
%     If Y is numeric, relieff by default performs RReliefF analysis for
%     regression. If Y is categorical, logical, a character array, or a cell
%     array of strings, relieff by default performs ReliefF analysis for
%     classification.
%  
%     Attribute ranks and weights computed by relieff usually depend on K. If
%     you set K to 1, the estimates computed by relieff can be unreliable for
%     noisy data. If you set K to a value comparable with the number of
%     observations (rows) in X, relieff can fail to find important
%     attributes. You can start with K=10 and investigate the stability and
%     reliability of relieff ranks and weights for various values of K.
%
% Matlab Code-Library for Feature Selection
%  Version 5.0 August 2017
%  Support: Giorgio Roffo
%  E-mail: giorgio.roffo@glasgow.ac.uk

%  If you use our toolbox please cite our supporting papers:
% 
%  BibTex
%  ------------------------------------------------------------------------
% @InProceedings{RoffoICCV17, 
% author={Giorgio Roffo and Simone Melzi and Umberto Castellani and Alessandro Vinciarelli}, 
% booktitle={2017 IEEE International Conference on Computer Vision (ICCV)}, 
% title={Infinite Latent Feature Selection: A Probabilistic Latent Graph-Based Ranking Approach}, 
% year={2017}, 
% month={Oct}}
%  ------------------------------------------------------------------------
% @Inbook{Roffo2017,
% author="Roffo, Giorgio
% and Melzi, Simone",
% editor="Appice, Annalisa
% and Ceci, Michelangelo
% and Loglisci, Corrado
% and Masciari, Elio
% and Ra{\'{s}}, Zbigniew W.",
% title="Ranking to Learn:",
% bookTitle="New Frontiers in Mining Complex Patterns: 5th International Workshop, NFMCP 2016, Held in Conjunction with ECML-PKDD 2016, Riva del Garda, Italy, September 19, 2016, Revised Selected Papers",
% year="2017",
% publisher="Springer International Publishing",
% address="Cham",
% pages="19--35",
% isbn="978-3-319-61461-8",
% doi="10.1007/978-3-319-61461-8_2",
% url="https://doi.org/10.1007/978-3-319-61461-8_2"
% }
%  ------------------------------------------------------------------------
% @InProceedings{RoffoICCV15, 
% author={G. Roffo and S. Melzi and M. Cristani}, 
% booktitle={2015 IEEE International Conference on Computer Vision (ICCV)}, 
% title={Infinite Feature Selection}, 
% year={2015}, 
% pages={4202-4210}, 
% doi={10.1109/ICCV.2015.478}, 
% month={Dec}}
%  ------------------------------------------------------------------------

fprintf('\n+ Feature selection method: Relief-F \n');
%% Wrapper: use Matlab implementation
[RANKED,WEIGHT] = relieff(X,Y,K);
% Matlab Code-Library for Feature Selection
% Contact: Giorgio Roffo email: giorgio.roffo@univr.it
