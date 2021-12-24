function WriteCoeffResult(r,Cfg)
SelectFeatureNum=Cfg.SelectFeatureNum;
OutputDir=Cfg.OutputDir;
if length(r)>1 %FeatureSelect
    h=figure;
    plot(SelectFeatureNum,r,'bo-', 'LineWidth',3)
    xlabel('Number of Features','FontSize',15)
    ylabel('Correlation Coefficent','FontSize',15)
    set(gca,'FontSize',15)
    set(gca, 'LineWidth',2)
    saveas(h,[OutputDir,filesep,'CorrCoefficent.bmp'])
    saveas(h,[OutputDir,filesep,'CorrCoefficent.fig'])
end