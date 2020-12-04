    function IPM_plot_data_orBound(IPMObject,IPMDesign,X,Y,PlotType,IndexCondition,Ci)
 
    if PlotType==1 % scatter all data
        scatter(X,Y,'o','MarkerFaceColor','b','MarkerFaceAlpha',0.3);
        box on; grid on; hold on;

    elseif PlotType==2 % plot IPM bounds
        X_predict=linspace(min(X),max(X),IPMObject.IntegrationKnots);
        I_predict=IPMObject.Predict(X_predict,IPMDesign.OptThetalow,IPMDesign.OptThetaup);% predict method
        I_predict= IPMObject.De_Normalize(I_predict,min(Y),max(Y));
        [~,Xindex]=sort(X_predict);
        plot(X_predict(Xindex), I_predict(:,Xindex),'k','LineWidth',1.5);
        box on; grid on; hold on;

    elseif PlotType==3 % scatter data with condition
        scatter(X(IndexCondition),Y(IndexCondition),'xr','MarkerEdgeAlpha',0.1')
        box on; grid on; hold on;
  
    elseif PlotType==4 % plot bounds minmax
        X_predict=linspace(min(X),max(X),IPMObject.IntegrationKnots);
        I_predict=IPMObject.Predict(X_predict,IPMDesign.OptThetalow,IPMDesign.OptThetaup);% predict method
        [~,Xindex]=sort(X_predict);
        plot(X_predict(Xindex), IPMObject.De_Normalize(I_predict,min(Y),max(Y)),'k','LineWidth',1.5);
        hold on
        plot(X_predict(Xindex), IPMObject.De_Normalize(I_predict+IPMDesign.EmpiricalCosts(Ci),min(Y),max(Y)),':k','LineWidth',1.5);
        plot(X_predict(Xindex), IPMObject.De_Normalize(I_predict-IPMDesign.EmpiricalCosts(Ci),min(Y),max(Y)),':k','LineWidth',1.5);
    end