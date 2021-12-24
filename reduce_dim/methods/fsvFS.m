function [ rank , w] = fsvFS( X,Y,numF )
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
fprintf('\n+ Feature selection method: FSV \n');
 
loop=0;
finished=0;
alpha = 5; % Default
[m,n] = size(X);
v = zeros(n,1);

%% Main LOOP
while (~finished),
    loop=loop+1;    
    scale = alpha*exp(-alpha*v);
    
    A=[diag(Y)*X, -diag(Y)*X, Y, -Y, -eye(m)];
    Obj = [scale',scale',0,0,zeros(1,m)];
    b = ones(m,1);
    x = slinearsolve(Obj',A,b,Inf);
    w = x(1:n)-x(n+1:2*n);
    b0 = x(2*n+1)-x(2*n+2);
    vnew=abs(w);
    
    if (norm(vnew-v,1)<10^(-5)*norm(v,1)),
        finished=1;
    else
        v=vnew;
    end;
    if (loop>2),
        finished=1;
    end;
    nfeat=length(find(vnew>100*eps));
    
    disp(['Iter ' num2str(loop) ' - feat ' num2str(nfeat)]);
    
    if nfeat<numF,
        finished=1;
    end;
end;

[~, ind] = sort(-abs(w));

w=w;
rank=ind;


end

