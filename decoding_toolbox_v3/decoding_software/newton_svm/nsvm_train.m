function model = nsvm_train(d,A,nu,wp1,arm)

    % d: labels
    % A: training data vector
    % nu: weighted parameter
    %     -1 - easy estimation
    %     0  - hard estimation
    %     any other value - used as nu by the algorithm
    %     default - 0
    % wp1: percentage of weight for class 1, wp must be 0<=wp1<=1  
    %      the precentage for class -1 is 1-wp1 (1 is more important,
    %      0 is less important		    
    % arm: 1 - use armijo, 0 - otherwise, default is 0
    
    if nargin<5,arm=0;end
    if nargin<4,wp1=.5;end
    if     nargin<3 || nu==0,   nu = EstNuLong(A,d);
    elseif nu==-1,              nu = EstNuShort(A,d);end

    [m,n]=size(A);

    if m>=n
        [model.w,model.gamma,model.iter]=nsvm_with_smw(A,d,nu,arm,wp1);
    else
        [model.w,model.gamma,model.iter]=nsvm_without_smw(A,d,nu,arm,wp1);
    end
    

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        NSVM for m>=n                                           %    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,gamma,iter]=nsvm_with_smw(A,d,nu,arm,wp1)

    % with armijo and with SMW
    maxIter=100;
    [m,n]=size(A);
    iter=0;
    u=zeros(m,1);e=ones(m,1);

    H=[spdiags(d,0,m,m)*A -d]; 

    % When balance is on this only approximates the norm
    alpha=1.1*((1/nu)+(norm(H',2)^2));

    if wp1==0.5
       v=u/nu+H*(H'*u)-e;
    else
       vt=ones(m,1);
       vt(d==1)=(1-wp1)*ones(length(find(d==1)),1);
       vt(d==-1)=wp1*ones(length(find(d==-1)),1);
       v=vt.*u/nu+H*(H'*u)-e;
    end  
    hu=-max((v-alpha*u),0)+v;




    while norm(hu)>10^(-3) && (iter < maxIter)  
        iter=iter+1;
        E=sign(max(v-alpha*u,0));
        if wp1==.5
           temp=(1./((alpha-(1/nu))*E+(1/nu)));    
        else
           temp=(1./((alpha*E)-E.*vt*(1/nu)+vt.*(1/nu)));
       end
        G=spdiags(temp.*(1-E),0,m,m);
        FHU=temp.*hu;
        LD=H'*FHU;
        GH=G*H;
        SM=speye(n+1)+H'*GH;
        R=chol(SM);
        sol=R'\LD;
        sol=R\sol;
        delta=FHU-GH*sol;

        if arm==1 %armijo step without the balance yet
            lambda=1;
            v1=v+e;
            v2=v;
            su=0.5*(u'*v1)-sum(u)+(1/(2*alpha))*(norm(max(-alpha*u+v2,0))^2-norm(v2)^2);        
            unew=u-lambda*delta;
            v1=u/nu+H*(H'*unew);
            v2=v1-e;
            sunew=0.5*(unew'*v1)-sum(unew)+(1/(2*alpha))*(norm(max(-alpha*unew+v2,0))^2-norm(v2)^2);
            while su-sunew < -(0.25)*lambda*hu

                lambda=0.1*lambda;           
                unew=u-lambda*delta;
                v1=u/nu+H*(H'*unew);
                v2=v1-e;
                sunew=0.5*(unew'*v1)-sum(unew)+(1/(2*alpha))*(norm(max(-alpha*unew+v2,0))^2-norm(v2)^2);
            end

        else
            unew=u-delta;
        end    

        u=unew;

        if wp1==.5
           v=u/nu+H*(H'*u)-e;
       else
           v=vt.*u/nu+H*(H'*u)-e;
       end
        hu=-max((v-alpha*u),0)+v;



    end

    w=A'*(d.*u);gamma=-sum(d.*u);

    return


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        NSVM for m<n                                           %    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [w,gamma,iter]=nsvm_without_smw(A,d,nu,arm,wp1)
    maxIter=100;
    [m,n]=size(A);
    iter=0;
    u=zeros(m,1);e=ones(m,1);
    HH=diag(d)*(A*A'+1)*diag(d);
    if wp1==.5
       Q=speye(m)/nu+HH;
    else
       vt=ones(m,1);
       vt(d==1)=(1-w1)*ones(length(find(d==1)));
       vt(d==-1)=w1*ones(length(find(d==-1)));
       Q=spdiags(vt,0,m,m)/nu+HH;
    end
    alpha=1.1*((1/nu)+norm(HH));
    v=Q*u-e;
    hu=-max((v-alpha*u),0)+v;

    while norm(hu)>10^(-3) && (iter < maxIter)
        iter=iter+1;
        dhu=sign(max(((Q-alpha*eye(m))*u-e),0));% the 1/2 thing
        dhu=sparse((diag(1-dhu))*Q+diag(alpha*dhu));
        delta=dhu\hu;

        if arm==1 %armijo step
            lambda=1;
            v=Q*u;
            ve=v-e;
            su=0.5*(u'*v)-sum(u)+(1/(2*alpha))*(norm(max(-alpha*u+ve,0))^2-norm(ve)^2);           
            unew=u-lambda*delta;
            v=Q*unew;
            ve=v-e;
            sunew=0.5*(unew'*v)-sum(unew)+(1/(2*alpha))*(norm(max(-alpha*unew+ve,0))^2-norm(ve)^2);
            at=0;
            while (su-sunew < -(0.25)*lambda*hu) && (at<5)
                at=at+1;
                disp('armijo');
                lambda=0.5*lambda;           
                unew=u-lambda*delta;
                v=Q*unew;
                ve=v-e;
                sunew=0.5*(unew'*v)-sum(unew)+(1/(2*alpha))*(norm(max(-alpha*unew+ve,0))^2-norm(ve)^2);
            end
            if at==5
                unew=u-delta;
            end    
        else
            unew=u-delta;
        end    


        u=unew;
        v=Q*unew-e;
        hu=-max((v-alpha*u),0)+v;


    end

    w=A'*(d.*u);gamma=-sum(d.*u);

    return

end




%%%%%%%%%%%%%%EstNuLong%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hard way to estimate nu if not specified by the user
function value = EstNuLong(C,d)

    [m,n]=size(C);e=ones(m,1);
    H=[C -e];
    if m<201
    H2=H;d2=d;
    else
    r=rand(m,1);
    [s1,s2]=sort(r);
    H2=H(s2(1:200),:);
    d2=d(s2(1:200));
    end

    lamda=1;
    [vu,u]=eig(H2*H2');u=diag(u);p=length(u);

    yt=d2'*vu;
    lamdaO=lamda+1;

    cnt=0;
    while (abs(lamdaO-lamda)>10e-4) && (cnt<100)
         cnt=cnt+1;
         nu1=0;pr=0;ee=0;waw=0;
         lamdaO=lamda;
         for i=1:p
             nu1= nu1 + lamda/(u(i)+lamda);
             pr= pr + u(i)/(u(i)+lamda)^2;
             ee= ee + u(i)*yt(i)^2/(u(i)+lamda)^3;
             waw= waw + lamda^2*yt(i)^2/(u(i)+lamda)^2;
         end

         lamda=nu1*ee/(pr*waw);
    end

    value =lamda;
    if cnt==100
        value=1;
    end    

    return

end


%%%%%%%%%%%%%%%%%EstNuShort%%%%%%%%%%%%%%%%%%%%%%%

% easy way to estimate nu if not specified by the user
function value = EstNuShort(C,d)

    value = 1/(sum(sum(C.^2))/size(C,2));
    return
    
end