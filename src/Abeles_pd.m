%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   This work is licensed under the Creative Commons Reconocimiento-NoComercial-CompartirIgual 4.0 Internacional License.
%   To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ajimenez 08.05.21 
% cbujalance/ecabello 03.12.21 (coherence average in RT)
% mromero 13.01.22 (negative z field & absorbed power)

% Ohta & Ishida 10.1364/AO.29.001952
% input: 
 % n: refractive index (array)
 % d: thickness (array) // size(d)=size(n)-2
 % wl: wavelength (scalar)
 % ang_inc: incident angle in rad (scalar)
 % z: depth (array)
% output
 % Rs (scalar)
 % Rp (scalar)
 % Ts (scalar)
 % Tp (scalar)
 % Fs (scalar)
 % Fs (scalar)
    
% output:

function [psi, delta] = Abeles_pd(n00,d00,s00,wl,ang_inc,z,lcoher,pp)

    %s00 = zeros(length(n00)-1,1);

%     s00 = [0,15,0,0];

pp=0;

    %Averages R, T and E with phase increments wwhen d>lcoher
        phyt = 0;
        auxi = find((d00-lcoher)>0);

     if isempty(auxi) == 0
         
         if pp == 0
             phyt = 0;
         else
             phyt = ([0:pp:360-pp])*pi/180;
         end

     else
         clear auxi
        auxi = 0;
     end
     
     len = length(phyt);
     
    %introduces a 0th air layer because the Abeles algorithm starts
    %calculating the field in the 1st layer
    d = [0 d00  0];
    n = [n00];
    nn = length(n);

    if isempty(s00) == true
        s00 = zeros(length(d)-1);
    end
    
    if(length(d00)~=(length(n00)-2))
    
        error('error in D or N dimensions');
    end

    
    % initialization
    Rs = 0;
    Rp = 0;
    Ts = 0;
    Tp = 0;
        
    t = zeros(size(n));
    c = zeros(size(n));
    f = zeros([1 nn]);
    rs = zeros([1 nn-1]);
    rp = zeros([1 nn-1]);
    ts = zeros([1 nn-1]);
    tp = zeros([1 nn-1]);

    t(1) = ang_inc;
    c(1) = cos(t(1));

    Cs0 = zeros([2 2 nn-1 len]);
    Cp0 = zeros([2 2 nn-1 len]);
    Cs = repmat(eye([2 2]),[1 1 len]);
    Cp = repmat(eye([2 2]),[1 1 len]);

    Fs = zeros(size(z));
    Fp = zeros(size(z));
    
    for k2 = 1:len
        
        for k1=2:nn
            % angle
            t(k1) = asin((n(k1-1)*sin(t(k1-1)))/n(k1));
            c(k1) = cos(t(k1));        

            % fresnel coeff capa (k1-1) - (k1) + scattering correction
            
            rs(k1-1) = (n(k1-1).*c(k1-1) - n(k1)*c(k1))   /  (n(k1-1).*c(k1-1) + n(k1)*c(k1)) *(exp(-2*(2*pi*s00(k1-1)/wl)^2*(n(k1-1))^2)) ;
            rp(k1-1) = (n(k1-1).*c(k1)   - n(k1)*c(k1-1)) /  (n(k1-1).*c(k1)   + n(k1)*c(k1-1)) *(exp(-2*(2*pi*s00(k1-1)/wl)^2*(n(k1-1))^2));

            ts(k1-1) = 2*n(k1-1)*c(k1-1) / (n(k1-1)*c(k1-1)+n(k1)*c(k1)) *(exp(-(2*pi*s00(k1-1)/wl)^2*(n(k1)-n(k1-1))^2/2));
            tp(k1-1) = 2*n(k1-1)*c(k1-1) / (n(k1-1)*c(k1)+n(k1)*c(k1-1)) *(exp(-(2*pi*s00(k1-1)/wl)^2*(n(k1)-n(k1-1))^2/2));    

                 
            % change of phase       

            phy=0;
            
            if ismember(k1-1,auxi+1)
                phy = phyt(k2);
            else 
                phy = 0;
            end

            f(k1-1) = 2*pi/wl*n(k1-1)*c(k1-1)*d(k1-1) + phy;
            
          % propagation matrix  
          
            Cs0(:,:,k1-1,k2) =  [         exp(-1i*(f(k1-1)))  rs(k1-1)*exp(-1i*(f(k1-1))); ...
                              rs(k1-1)*exp(+1i*(f(k1-1)))           exp(+1i*(f(k1-1)))];        
            Cp0(:,:,k1-1,k2) =  [         exp(-1i*(f(k1-1)))  rp(k1-1)*exp(-1i*(f(k1-1))); ...
                              rp(k1-1)*exp(+1i*(f(k1-1)))           exp(+1i*(f(k1-1)))];              
                    

            Cs(:,:,k2) = Cs(:,:,k2)*Cs0(:,:,k1-1,k2);
            Cp(:,:,k2) = Cp(:,:,k2)*Cp0(:,:,k1-1,k2);

            % r and t

            Rs0 = Cs(2,1,k2)./Cs(1,1,k2);
            Rp0 = Cp(2,1,k2)./Cp(1,1,k2);

        end
        
        % R and T
        Rs = Rs+Rs0;
        Rp = Rp+Rp0;
    
    end
    
   
   Rs = Rs./len;
   Rp = Rp./len;

   c_ratio = Rp0/Rs0;

   psi = atand(abs(c_ratio));
    
   % Elipso 1
   delta = abs(180/pi*phase(-c_ratio));

%    % Elipso 2
%    delta = (180/pi*phase(-c_ratio)) ;

   % Elipso 3
%    delta = (angle(-c_ratio)+pi)*180/pi/2 ;

%    if delta<0
%         delta = delta+360;
%    end
   

   
end
    



