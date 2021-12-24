function [fea, score] = mRMR(X_train, Y_train, K)
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
% MID scheme according to MRMR
%
% By Hanchuan Peng
% April 16, 2003
%
% fprintf('\n+ Feature selection method: mRMR \n');

bdisp=0;

nd = size(X_train,2);
nc = size(X_train,1);

t1=cputime;
for i=1:nd, 
   t(i) = mutualinfo(X_train(:,i), Y_train);
end; 
% fprintf('calculate the marginal dmi costs %5.1fs.\n', cputime-t1);

[tmp, idxs] = sort(-t);
fea_base = idxs(1:K);

fea(1) = idxs(1);

KMAX = min(1000,nd); %500

idxleft = idxs(2:KMAX);

k=1;
% if bdisp==1,
% % fprintf('k=1 cost_time=(N/A) cur_fea=%X_train #left_cand=%X_train\n', ...
% %       fea(k), length(idxleft));
% end;

for k=2:K,
   t1=cputime;
   ncand = length(idxleft);
   curlastfea = length(fea);
   for i=1:ncand,
      t_mi(i) = mutualinfo(X_train(:,idxleft(i)), Y_train); 
      mi_array(idxleft(i),curlastfea) = getmultimi(X_train(:,fea(curlastfea)), X_train(:,idxleft(i)));
      c_mi(i) = mean(mi_array(idxleft(i), :)); 
   end;

   [score(k), fea(k)] = max(t_mi(1:ncand) - c_mi(1:ncand));

   tmpidx = fea(k); fea(k) = idxleft(tmpidx); idxleft(tmpidx) = [];
   
%    if bdisp==1,
% %    fprintf('k=%X_train cost_time=%5.4f cur_fea=%X_train #left_cand=%X_train\n', ...
%       k, cputime-t1, fea(k), length(idxleft));
%    end;
end;

return;

%===================================== 
function c = getmultimi(da, dt) 
for i=1:size(da,2), 
   c(i) = mutualinfo(da(:,i), dt);
end; 
    
