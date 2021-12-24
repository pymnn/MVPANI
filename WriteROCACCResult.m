function AUC=WriteROCACCResult(meanaccuracy,Predict_test_Label,Cfg,strClassAlgorithm)
sel_checkboxROC=Cfg.sel_checkboxROC;
OutputDir=Cfg.OutputDir;
SelectFeatureNum=Cfg.SelectFeatureNum;
AUC=0;
if length(meanaccuracy)>1 %FeatureSelect
    h=figure;
    plot(SelectFeatureNum,meanaccuracy,'bo-', 'LineWidth',3)
    grid on
    xlabel('Number of Features','FontSize',15)
    ylabel('Accuracy (%)','FontSize',15)
    set(gca,'FontSize',15)
    set(gca, 'LineWidth',2)
    saveas(h,[OutputDir,filesep,'MeanAccuracy.bmp'])
    saveas(h,[OutputDir,filesep,'MeanAccuracy.fig'])
end

if sel_checkboxROC
    if length(size(Predict_test_Label))<3
        [neg1,pos1]=PlotROC(Predict_test_Label,strClassAlgorithm);
        
        h=figure;
        switch strClassAlgorithm
            case {'SVM (Support Vector Machine)','LR (Logistic Regression)'}
                plot(neg1,pos1,'r','LineWidth',3)
                AUC=trapz(neg1,pos1);
            otherwise
                plot(pos1,neg1,'r','LineWidth',3)
                AUC=1-trapz(neg1,pos1);
        end
        
        axis([0 1 0 1])
        xlabel('False Positive Rate','FontSize',15)
        ylabel('True Positive Rate','FontSize',15)
        set(gca, 'LineWidth',2)
        set(gca,'FontSize',15)
        title(['AUC=',num2str(AUC)]);
        saveas(h,[OutputDir,filesep,'ROC.bmp'])
        saveas(h,[OutputDir,filesep,'ROC.fig'])
        
    else
        h=figure;
        c=jet(size(Predict_test_Label,3));
        for jPreTes=1:size(Predict_test_Label,3)
            [neg1,pos1]=PlotROC(Predict_test_Label(:,:,jPreTes),strClassAlgorithm);
            
            hold on
            switch strClassAlgorithm
                case {'SVM (Support Vector Machine)','LDA (Linear Discriminant Analysis)','LR (Logistic Regression)'}
                    lh=plot(neg1,pos1,'color',c(jPreTes,:),'LineWidth',3);
                    AUC(jPreTes)=trapz(neg1,pos1);
                otherwise
                    lh=plot(pos1,neg1,'color',c(jPreTes,:),'LineWidth',3);
                    AUC(jPreTes)=1-trapz(neg1,pos1);
            end
        end
        axis([0 1 0 1])
        xlabel('False Positive Rate','FontSize',15)
        ylabel('True Positive Rate','FontSize',15)
        set(gca,'FontSize',15)
        set(gca, 'LineWidth',2)
        legendtext=[];
        for jPreTes=1:size(Predict_test_Label,3)
            legendtext=[legendtext; ['SelFeature',sprintf('%03d',jPreTes),'AUC',sprintf('%.2f',AUC(jPreTes))]];
        end
        legend(legendtext)
        saveas(h,[OutputDir,filesep,'ROC.bmp'])
        saveas(h,[OutputDir,filesep,'ROC.fig'])
    end
end