%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   This work is licensed under the Creative Commons Reconocimiento-NoComercial-CompartirIgual 4.0 Internacional License.
%   To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Rt,Rt_S,Rt_P,Tt,Tt_S,Tt_P] = f_plot_RT(N, D, lcoher, wl, theta,Rexp,Texp,onlyplot)
    
   

    Tt_S=zeros(length(wl),length(theta));
    Tt_P=zeros(length(wl),length(theta));
    Rt_S=zeros(length(wl),length(theta));
    Rt_P=zeros(length(wl),length(theta));
    Rt=zeros(length(wl),length(theta));
    Tt=zeros(length(wl),length(theta));
    
    for k1=1:length(theta)
        for k2=1:length(wl)
            
            [Rt_S(k2,k1), Rt_P(k2,k1), Tt_S(k2,k1), Tt_P(k2,k1), ~, ~, ~, ~] = RTF_Abeles_F(N(k2,:), D', wl(k2),theta(k1)*pi/180,[],lcoher,30);
            
        end
        
        Rt(:,k1) = 0.5*(Rt_S(:,k1) + Rt_P(:,k1));
        Tt(:,k1) = 0.5*(Tt_S(:,k1) + Tt_P(:,k1));


    end


    if length(theta)==1
        Rexp = Rexp';
        Texp = Texp';
    end

    aux = 0;
    figure(2)
    clf
    for jj=1:length(theta)
        subplot(length(theta),2,jj+aux)
        if onlyplot == false
            plot(wl,Texp(:,jj),'LineWidth',1.5,'color','k','LineStyle','-')
            hold on
        end
        plot(wl,Tt(:,jj),'LineWidth',1.5,'color','r','LineStyle','--')
        xlabel("\lambda (nm)")
        ylabel("T")
        xlim([min(wl),max(wl)])

        aux = aux+1;
        subplot(length(theta),2,jj+aux)
        if onlyplot == false
            plot(wl,Rexp(:,jj),'LineWidth',1.5,'color','k','LineStyle','-')
            hold on
        end
        plot(wl,Rt(:,jj),'LineWidth',1.5,'color','r','LineStyle','--')
        xlabel("\lambda (nm)")
        ylabel("R")
        xlim([min(wl),max(wl)])
    end

    
    
end