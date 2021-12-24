function [neg1,pos1]=PlotROC(Predict_test_Label,strClassAlgorithm)
value=Predict_test_Label(:,3);
Label=Predict_test_Label(:,1);
ReLabel=unique(Label);
Label(find(Label==max(ReLabel)))=max(ReLabel);
Label(find(Label==min(ReLabel)))=min(ReLabel);
switch strClassAlgorithm
    case 'LDA(Linear Discriminant Analysis)'
        value=-Predict_test_Label(:,3);
    otherwise 
        value=Predict_test_Label(:,3);
end
MinLabel=min(Label);
MaxLabel=max(Label);
xaxis(1)=0;
yaxis(1)=0;
for i=1:length(value)
    for j=1:length(Label)
        if value(j)<value(i)
            prelabel(j,1)=MinLabel;
        else
            prelabel(j,1)=MaxLabel;
        end
    end
    PredictionPos=find(prelabel==MaxLabel);%predict positive value
    Pos_num=find(Label==MaxLabel); %true positive value
    Neg_num=find(Label==MinLabel);%true negative value
    TP=length(intersect(PredictionPos,Pos_num));
    FP=length(intersect(PredictionPos,Neg_num));
    xaxis(i+1)=FP/length(find(Label==MinLabel));
    yaxis(i+1)=TP/length(find(Label==MaxLabel));

end
[neg,Index]=sort(xaxis);
pos=yaxis(Index);
[pos1,Index]=sort(pos);
neg1=neg(Index);



