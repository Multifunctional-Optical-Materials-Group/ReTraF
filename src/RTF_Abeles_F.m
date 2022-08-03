 
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

function [Rs, Rp, Ts, Tp, Fs, Fp, Pabs_s, Pabs_p] = RTF_Abeles_F(n00,d00,wl,ang_inc,z,lcoher,pp)

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

            % fresnel coeff capa (k1-1) - (k1)
            
            rs(k1-1) = (n(k1-1).*c(k1-1) - n(k1)*c(k1))   /  (n(k1-1).*c(k1-1) + n(k1)*c(k1)) ;
            rp(k1-1) = (n(k1-1).*c(k1)   - n(k1)*c(k1-1)) /  (n(k1-1).*c(k1)   + n(k1)*c(k1-1));

            ts(k1-1) = 2*n(k1-1)*c(k1-1) / (n(k1-1)*c(k1-1)+n(k1)*c(k1));
            tp(k1-1) = 2*n(k1-1)*c(k1-1) / (n(k1-1)*c(k1)+n(k1)*c(k1-1));    

                 
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

            Ts0 = prod(ts)./Cs(1,1,k2);
            Tp0 = prod(tp)./Cp(1,1,k2);

        end
        
        % R and T
        Rs = Rs+abs(Rs0)^2;
        Rp = Rp+abs(Rp0)^2;

        Ts = Ts+real(n(end)/n(1)*c(end)/c(1))*abs(Ts0).^2;
        Tp = Tp+real(conj(n(end))/conj(n(1))*c(end)/c(1))*abs(Tp0).^2;
    
    
    %Field

       if isempty(z) == 1
            Fs = [];
            Fp = [];
       else

        Ds = repmat(eye([2 2]),[1 1 nn len]);
        Dp = repmat(eye([2 2]),[1 1 nn len]);

        Esp = ones([1 nn-2]);
        Esn = ones([1 nn-2]);

        Epp = ones([1 nn-2]);
        Epn = ones([1 nn-2]);    

        Eszp = zeros(size(z));
        Eszn = zeros(size(z));
        Epzp = zeros(size(z));
        Epzn = zeros(size(z));

        Ex = zeros(size(z));
        Ey = zeros(size(z));
        Ez = zeros(size(z));

        E0s = 1./sqrt(2);
        E0p = 1./sqrt(2);


    
        % field coef
        Ds(:,:,1,k2) = Cs(:,:,k2);
        Dp(:,:,1,k2) = Cp(:,:,k2);

        for k1=2:nn
            for k3=k1:nn-1
                Ds(:,:,k1,k2) = squeeze(Ds(:,:,k1,k2))*squeeze(Cs0(:,:,k3,k2));
                Dp(:,:,k1,k2) = squeeze(Dp(:,:,k1,k2))*squeeze(Cp0(:,:,k3,k2));
            end
            Esp(k2,k1-1) = prod(ts(1:k1-1))*squeeze(Ds(1,1,k1,k2))/Cs(1,1,k2)*E0s;
            Epp(k2,k1-1) = prod(tp(1:k1-1))*squeeze(Dp(1,1,k1,k2))/Cp(1,1,k2)*E0p;

            Esn(k2,k1-1) = prod(ts(1:k1-1))*squeeze(Ds(2,1,k1,k2))/Cs(1,1,k2)*E0s;
            Epn(k2,k1-1) = prod(tp(1:k1-1))*squeeze(Dp(2,1,k1,k2))/Cp(1,1,k2)*E0p;
        end
    


       %Field inside the film
     
       if z(1) < d(1)
           
           h = [cumsum(d)] ;
           h(end) = inf;
           ini = 1;
        
        z0 = z(z<h(1));
        
        fin = ini+length(z0)-1;
        ff = 2*pi/wl*n(1)*c(1);
        
        Eszp(ini:fin) = E0s.*exp(+1i*ff*z0);

        Eszn(ini:fin) = Rs0*E0s.*exp(-1i*ff*z0);

        Epzp(ini:fin) = E0p.*exp(+1i*ff*z0);

         Epzn(ini:fin) = Rp0*E0p.*exp(-1i*ff*z0);
        
         Ex(ini:fin) = (Epzp(ini:fin) + Epzn(ini:fin))*c(1,1);
        Ey(ini:fin) = (Eszp(ini:fin) + Eszn(ini:fin));
        Ez(ini:fin) = (Epzp(ini:fin) - Epzn(ini:fin))*sin(t(1,1));
        
        ini = fin+1;
        
         
       else
           
           h = [cumsum(d)];
           h(end) = inf;
           ini = 1;
        
           
       end
       

            for i=1:length(h)-1
         
                if ismember(i,auxi+1)
                    phy = phyt(k2);
                else 
                    phy = 0;
                end
       
      
                z0 = z(z>h(i) & z<=h(i+1)) - h(i);
                fin = ini+length(z0)-1;
                
                ff = 2*pi/wl*n(i+1)*c(i+1);
                            
                Eszp(ini:fin) = Esp(k2,i).*exp(+1i*ff*z0 + 1i*phy);
                Eszn(ini:fin) = Esn(k2,i).*exp(-1i*ff*z0 - 1i*phy);
                Epzp(ini:fin) = Epp(k2,i).*exp(+1i*ff*z0 + 1i*phy);
                Epzn(ini:fin) = Epn(k2,i).*exp(-1i*ff*z0 - 1i*phy);

                 Ex(ini:fin) = (Epzp(ini:fin) + Epzn(ini:fin))*c(1,i+1);
                Ey(ini:fin) = (Eszp(ini:fin) + Eszn(ini:fin));
                 Ez(ini:fin) = (Epzp(ini:fin) - Epzn(ini:fin))*sin(t(1,i+1));

                ini = fin+1;
            end    

           

        Fx = abs(Ex).^2./abs(E0p).^2;
        Fy = abs(Ey).^2./abs(E0s).^2;
        Fz = abs(Ez).^2./abs(E0p).^2;

        Fs = Fs + Fy;
        Fp = Fp + Fx + Fz;

   
       end
    end
    
   
   Fs = Fs./len;
   Fp = Fp./len;
   
   Ts = Ts./len;
   Tp = Tp./len;
   Rs = Rs./len;
   Rp = Rp./len;

   h = [cumsum(d)];
   h(end) = inf;

   Pabs_s=zeros(length(z),1);
   Pabs_p=zeros(length(z),1);

   if isempty(z)==0
       aux=find(z<h(1));
       Pabs_s(aux) = 4*pi.*real(n00(1)).*imag(n00(1)).* Fs(aux)./(cos(ang_inc)*wl*n00(1));
       Pabs_p(aux) = 4*pi.*real(n00(1)).*imag(n00(1)).* Fp(aux)./(cos(ang_inc)*wl*n00(1));
       for k=1:length(h)-1
           aux=find(z>h(k) & z<=h(k+1));
           Pabs_s(aux) = 4*pi.*real(n00(k+1)).*imag(n00(k+1)).* Fs(aux)./(cos(ang_inc)*wl.*n00(1));
           Pabs_p(aux) = 4*pi.*real(n00(k+1)).*imag(n00(k+1)).* Fp(aux)./(cos(ang_inc)*wl.*n00(1));
       end
       % Total power absorbed (s pol) = trapz(z,Pabs_s);
       % Total power absorbed (p pol) = trapz(z,Pabs_p);
       %Pabs_s=trapz(z,Pabs_s);
       %Pabs_p=trapz(z,Pabs_p);
   else
       Pabs_s=[];
       Pabs_p=[];
   end

   
end
    



