function [p dv] = nsvm_test(dd,AA,model)

    % dd: labels
    % AA: test data vector
    % model: result of training
    
    dv = AA*model.w-model.gamma; % decision values
    p=sign(dv); % predicted labels
    % not necessary for our toolbox, but in case you want the accuracy too
    % activate this and pass as third output
%     corr=length(find(p==dd))/size(AA,1)*100;

end