function [Rs] = trnsfr_mm_dspless_spol(num, d, n, wv, phi)
%num = number of total layers including air and substrate (layers+2)
%d = matrix of thicknesses 
%n = matrix of refractive indices corresponding to matrix of thicknesses
%wv = vector of wavelengths
%phi = vector of angles in degrees
theta  = zeros(1,num);
kn = zeros(1,num-2);
Ip = zeros(2,2,num-2);
P = zeros(2,2,num-2);
Rp = zeros(length(wv),length(phi));

for n_wv = 1:length(wv)%wavelength loop for dispersionless materials
    
    k0 = 2*pi/wv(n_wv); % free spave wave-vector
    
    
    for n_phi = 1:length(phi)%AOI loop
        
        
        theta(1) = phi(n_phi)*pi/180; %from degree to radian
        
        for k=2: num
            theta(k) = asin(n(1)*sin(theta(1))/n(k)); %Snell's law
        end
        
        
        % This calculates the Tansfer matrix, I, for reflection and
        % transmission at an interface between materials with complex dielectric
        % constant n(i) and n(i+1).
        
        for i=1:num-1
            
            rs=(n(i)*cos(theta(i))-n(i+1)*cos(theta(i+1)))/(n(i)*cos(theta(i))+n(i+1)*cos(theta(i+1)));
            ts=2*n(i)*cos(theta(i))/(n(i)*cos(theta(i))+n(i+1)*cos(theta(i+1)));
            Is(:,:,i)=[1 rs; rs 1]/ts;
            
            
            
%             rp=(n(i)*cos(theta(i+1))-n(i+1)*cos(theta(i)))/(n(i)*cos(theta(i+1))+n(i+1)*cos(theta(i)));
%             tp=2*n(i)*cos(theta(i))/(n(i)*cos(theta(i+1))+n(i+1)*cos(theta(i)));
%             Ip(:,:,i)=[1 rp; rp 1]/tp;
            
        end
        
        %
        % This calculates the propagation matrix, P, through a material of
        % complex dielectric constant n and thickness d for the wavelength wv.
        
        for i=1:num-2
            kn(i) = k0*n(i+1)*cos(theta(i+1));
            P(:,:,i) = [exp(-1i*kn(i)*d(i)) 0; 0 exp(1i*kn(i)*d(i))];
        end
        
%         Sp=Ip(:,:,1);
         Ss=Is(:,:,1);
        
        for i=1:num-2
                 Ss=Ss*P(:,:,i)*Is(:,:,i+1);
%             Sp=Sp*P(:,:,i)*Ip(:,:,i+1);
        end
        
        %Reflection
%         Rp(n_wv, n_phi) = (abs(Sp(2,1)/Sp(1,1)))^2;
         Rs(n_wv, n_phi) = (abs(Ss(2,1)/Ss(1,1)))^2;
        
        %Transmission
        %Tp(n_wv, n_phi) = abs(1/Sp(1,1))^2*real(n(num)*cos(conj(theta(num))))/real(n(1)*cos(conj(theta(1))));
        %   Ts(n_wv, n_phi) = abs(1/Ss(1,1))^2*real(n(num)*cos(conj(theta(num))))/real(n(1)*cos(conj(theta(1))));
        
        % %   Absorption
        %   Ap(n_phi) = 1-Tp(n_phi)-Rp(n_phi);
        %   As(n_phi) = 1-Ts(n_phi)-Rs(n_phi);
    end
end
end
