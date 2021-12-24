function ranking = UDFS(X,nClass)
% run UDFS feature selection algorithm
% REF:
% @inproceedings{Yang:2011:LNR:2283516.2283660,
%  author = {Yang, Yi and Shen, Heng Tao and Ma, Zhigang and Huang, Zi and Zhou, Xiaofang},
%  title = {L2,1-norm Regularized Discriminative Feature Selection for Unsupervised Learning},
%  booktitle = {Proceedings of the Twenty-Second International Joint Conference on Artificial Intelligence - Volume Volume Two},
%  series = {IJCAI'11},
%  year = {2011},
%  isbn = {978-1-57735-514-4},
%  location = {Barcelona, Catalonia, Spain},
%  pages = {1589--1594},
%  numpages = {6},
%  url = {http://dx.doi.org/10.5591/978-1-57735-516-8/IJCAI11-267},
%  doi = {10.5591/978-1-57735-516-8/IJCAI11-267},
%  acmid = {2283660},
%  publisher = {AAAI Press},
% } 
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

%======================setup===========================
gammaCandi = 10.^(-5);
lamdaCandi = 10.^(-5);
knnCandi = 1;
paramCell = fs_unsup_udfs_build_param(knnCandi, gammaCandi, lamdaCandi);
%======================================================
X = double(X);

disp('UDFS: Regularized Discriminative Feature Selection for Unsupervised Learning');
param = paramCell{1};
L = LocalDisAna(X', param);
if sum(sum(isnan(L)))>0
    ranking = [1:size(X,2)];
else
    A = X'*L*X;
    W = fs_unsup_udfs(A, nClass, param.gamma);
    [~, ranking] = sort(sum(W.*W,2),'descend');
end

end