%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   This work is licensed under the Creative Commons Reconocimiento-NoComercial-Compapsi_tirIgual 4.0 Internacional License.
%   To view a copy of this license, visit hdelta_tp://creativecommons.org/licenses/by-nc-sa/4.0/.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [psi_t,delta_t] = f_plot_pd(N, D,s00, lcoher, wl, theta,psi_exp,delta_exp,onlyplot)
    
   

    psi_t=zeros(length(wl),length(theta));
    delta_t=zeros(length(wl),length(theta));
    
    for k1=1:length(theta)
        for k2=1:length(wl)
            [psi_t(k2,k1),delta_t(k2,k1)] = Abeles_pd(N(k2,:), D',s00, wl(k2),theta(k1)*pi/180,[],lcoher,30);
            
        end

    end


    if length(theta)==1
        psi_exp = psi_exp';
        delta_exp = delta_exp';
    end

    aux = 0;
    figure(2)
    clf
    for jj=1:length(theta)
        subplot(2,length(theta),jj)
        if onlyplot == false
            plot(wl,delta_exp(:,jj),'LineWidth',1.5,'color','k','LineStyle','-')
            hold on
        end
        plot(wl,delta_t(:,jj),'LineWidth',1.5,'color','r','LineStyle','--')
        xlabel("\lambda (nm)")
        ylabel("\Delta (º)")
        xlim([min(wl),max(wl)])

        aux = aux+1;
%         subplot(2,length(theta),jj+aux)
%         if onlyplot == false
%             plot(wl,psi_exp(:,jj),'LineWidth',1.5,'color','k','LineStyle','-')
%             hold on
%         end
%         plot(wl,psi_t(:,jj),'LineWidth',1.5,'color','r','LineStyle','--')
%         xlabel("\lambda (nm)")
%         ylabel("\psi (º)")
%         xlim([min(wl),max(wl)])
    end

    for jj=1:length(theta)
%         subplot(2,length(theta),jj+aux)
%         if onlyplot == false
%             plot(wl,delta_exp(:,jj),'LineWidth',1.5,'color','k','LineStyle','-')
%             hold on
%         end
%         plot(wl,delta_t(:,jj),'LineWidth',1.5,'color','r','LineStyle','--')
%         xlabel("\lambda (nm)")
%         ylabel("\Delta (º)")
%         xlim([min(wl),max(wl)])

        
        subplot(2,length(theta),jj+aux)
        if onlyplot == false
            plot(wl,psi_exp(:,jj),'LineWidth',1.5,'color','k','LineStyle','-')
            hold on
        end
        plot(wl,psi_t(:,jj),'LineWidth',1.5,'color','r','LineStyle','--')
        xlabel("\lambda (nm)")
        ylabel("\psi (º)")
        xlim([min(wl),max(wl)])

    end
    
end