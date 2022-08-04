
function [Rt,Rt_S,Rt_P,Tt,Tt_S,Tt_P] = f_plot_RT_scatt(N, D, lcoher, wl, theta,Rexp,Texp,alpha)
    
   

    Tt_S=zeros(length(wl),length(theta));
    Tt_P=zeros(length(wl),length(theta));
    Rt_S=zeros(length(wl),length(theta));
    Rt_P=zeros(length(wl),length(theta));
    Rt=zeros(length(wl),length(theta));
    Tt=zeros(length(wl),length(theta));
    Ae=zeros(length(wl),length(theta));
    
    for k1=1:length(theta)
        for k2=1:length(wl)
            
            [Rt_S(k2,k1), Rt_P(k2,k1), Tt_S(k2,k1), Tt_P(k2,k1), ~, ~, ~, ~] = RTF_Abeles_F(N(k2,:), D', wl(k2),theta(k1)*pi/180,[],lcoher,30);
            
        end
        
        Rt(:,k1) = 0.5*(Rt_S(:,k1) + Rt_P(:,k1));
        Tt(:,k1) = 0.5*(Tt_S(:,k1) + Tt_P(:,k1));
        Ae(:,k1) = 1-Rexp(:,k1)-Texp(:,k1);

    end

    aux = 0;
    figure(2)
    clf
    for jj=1:length(theta)
        subplot(length(theta),2,jj+aux)
        plot(wl,Texp(:,jj),'LineWidth',1.5,'color','k','LineStyle','-')
        hold on
        plot(wl,Tt(:,jj)-alpha(jj)*Ae(:,jj),'LineWidth',1.5,'color','r','LineStyle','--')
        xlabel("\lambda (nm)")
        ylabel("T")
        xlim([min(wl),max(wl)])

        aux = aux+1;
        subplot(length(theta),2,jj+aux)
        plot(wl,Rexp(:,jj),'LineWidth',1.5,'color','k','LineStyle','-')
        hold on
        plot(wl,Rt(:,jj)-(1-alpha(jj))*Ae(:,jj),'LineWidth',1.5,'color','r','LineStyle','--')
        xlabel("\lambda (nm)")
        ylabel("R")
        xlim([min(wl),max(wl)])
    end
    
    
end