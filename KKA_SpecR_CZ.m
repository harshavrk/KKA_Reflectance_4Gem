%author: Harshavardhan Reddy Kalluru Date: 22nd March 2026
rfs=dlmread('CZ_specR.dat');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%read data from tranmission spectrum saves as ASCII .dat file
d=size(rfs,1);
l0=rfs(1,1)-1;
l=zeros(d,1);R=zeros(d,1);theta=zeros(d,1);nfin=zeros(d,1);ls=zeros(d,1);Rs=zeros(d,1);%%%%Declaration of Matrices
igm=zeros(d,d);s=zeros(d,1);sg=zeros(d,1);sh=zeros(d,1);sf=zeros(d,1);thetalo=zeros(d,1);thetam=zeros(d,1);thetaup=zeros(d,1); IGg=zeros(d,1);IGh=zeros(d,1);%%%%Declaration of Matrices
thetafin=zeros(d,1);u=zeros(d,1);v=zeros(d,1);p=zeros(d,1);q=zeros(d,1);IGgI=zeros(d,1);IGhI=zeros(d,1);%%%%Declaration of Matrices
en=zeros(d,1);PPg=zeros(d,1);PPh=zeros(d,1);n=zeros(d,1);dnf=zeros(d,d);kl=zeros(d,1);%%%%Declaration of Matrices

for i=1:d
l(i,1) = rfs(i,1);%%%%Declaration of wavelength row vector
en(i,1)=1239.9219/l(i,1);%%%%Declaration of enrgy row vector
R(i,1) = rfs(i,2)/100;%%%%Declaration of reflectance row vector
end

ds=1;%%%%Declaration of wavelength resolution
a=199;%%%%Declaration of lower limit of KK analysis
b=1001;%%%%Declaration of upper limit of KK analysis
g=950;%%%%Declaration of first wavelength point (g) in saturation region of Reflection spectrum
h=980;%%%%Declaration of second wavelength point (h) in saturation region of Reflection spectrum

for i=1:d
    for j=1:d
        if (l0+j*ds)==g
            IGgI(i,j) = 0;%%%%Declaration of Integrand Matrix of KK analysis for Phase retrieval at wavelength g
        else
            IGgI(i,j) = 1e9*(1/(pi))*(log(R(j,1))*l(i,1))/(g^2-(l0+j*ds).^2);%%%%Declaration of Integrand Matrix of KK analysis for Phase retrieval at wavelength g
        end
    end
end
sg=sum(IGgI,2);

for i=1:d
   PPg(i,1)=(1e-9)*(sg(i,1)-IGgI(i,1)/2-IGgI(i,d)/2);%%%%Trapezoidal numerical integration for Phase retrieval
end
Pg=(PPg(751,1));%%%%Selecting the Phase value at first wavelength point (g) 

for i=1:d
    for j=1:d
        if (l0+j*ds)==h
            IGhI(i,j) = 0;%%%%Declaration of Integrand Matrix of KK analysis for Phase retrieval at wavelength h
        else
            IGhI(i,j) = 1e9*(1/(pi))*(log(R(j,1))*l(i,1))/(h^2-(l0+j*ds).^2);%%%%Declaration of Integrand Matrix of KK analysis for Phase retrieval at wavelength h
        end
    end
end
sh=sum(IGhI,2);

for i=1:d
   PPh(i,1)=(1e-9)*(sh(i,1)-IGhI(i,1)/2-IGhI(i,d)/2);%%%%Trapezoidal numerical integration for Phase retrieval
end
Ph=(PPh(781,1));%%%%Selecting the Phase value at second wavelength point (h) 

ag = log(abs((a+g)/(a-g)));bg=log(abs((b+g)/(b-g)));
ah = log(abs((a+h)/(a-h)));bh=log(abs((b+h)/(b-h)));
A=(Ph*bg-Pg*bh)/(ag*bh-ah*bg);%%%%Evaluating the coefficient A 
B=(Pg*ah-ag*Ph)/(ag*bh-ah*bg);%%%%Evaluating the coefficient B

for i=1:d
    for j=1:d
        if l(i,1)==l0+j*ds
            igm(i,j) = 0;%%%%Declaration of Integrand Matrix of KK analysis for Phase retrieval over the range (a,b)
        else
            igm(i,j) = 1e9*(1/(pi))*(log(R(j,1))*l(i,1))/(l(i,1).^2-(l0+j*ds).^2);%%%%Declaration of Integrand Matrix of KK analysis for Phase retrieval over the range (a,b)
        end
    end
end
s=sum(igm,2);

for i=1:d
   thetam(i,1)=(1e-9)*(s(i,1)-igm(i,1)/2-igm(i,d)/2);%%%%Trapezoidal numerical integration for Phase retrieval
end

for i=1:d
   thetafin(i,1)=(A)*log(abs((a+l(i,1))/(a-l(i,1)))) + thetam(i,1) + (B)*log(abs((b+l(i,1))/(b-l(i,1))));%%%%Stitching Phase function from zero to infinity
end

u(:,1)=(1-R(:,1))./(1+R(:,1)-2*cos(thetafin(:,1)).*sqrt(R(:,1)));%%%%Evaluating bulk refractive index from Phase and Reflectance
v(:,1)=(-2*sqrt(R(:,1)).*sin(thetafin(:,1)))./(1+R(:,1)-2*cos(thetafin(:,1)).*sqrt(R(:,1)));%%%%Evaluating bulk extinction constant from Phase and Reflectance

%Doubly subtractive KK anlysis intiation

n2=2.1180;%known refractive index at wavelength of 800.99 nm 
n1=2.4591;%known refractive index at wavelength of 208.00 nm 

wl2=800.99;%known pivot point wavelength
wl1=208.00;%known pivot point wavelength

for i=1:d
kl(i,1)= (v(i,1)*l(i,1)*(1e-9))/(4*pi*2*6*1e-3);%Evaluating the extinction coefficient for a 1 mm thick slide
end

for i=1:d
    for j=1:d
        if l(i,1)==l0+j*ds||l0+j*ds==wl1||l0+j*ds==wl2
            dnf(i,j) = 0;
        else
            dnf(i,j) = (kl(j,1)/(l(i,1)^2))*((l0+j*ds)^3/((l(i,1)^2-(l0+j*ds)^2)*(wl2^2-(l0+j*ds)^2)*(wl1^2-(l0+j*ds)^2)));%%%%Declaration of Integrand Matrix of DSKK analysis for refractive index retrieval
        end
    end
end

sf=sum(dnf,2);

for i=1:d
n(i,1)= 1 + (n1-1)*(wl1^2/l(i,1)^2)*((wl2^2-l(i,1)^2)/(wl2^2 - wl1^2)) + (n2-1)*((wl2^2/l(i,1)^2)*((wl1^2-l(i,1)^2)/(wl1^2 - wl2^2))) + (2/pi)*(l(i,1)^2-wl2^2)*(l(i,1)^2-wl1^2)*(l(d,1)-l(d-1,1))*((sf(i,1)-dnf(i,2)/2-dnf(i,d-1)/2));%%%%Trapezoidal numerical integration for DSKK refractive index evaluation
end

figure; yyaxis left; plot(l,thetafin,'b');ylabel('Phase (radian)');yyaxis right;plot(l,R,'r'); ylabel('Specular Reflectance'); axis square; grid on; xlabel('wavelength');legend('Cubic Zirconia (CZ)');
figure; plot(l,n,'b');axis square; legend('DSKK estimate for Cubic Zirconia (CZ)');grid on; xlabel('wavelength'); ylabel('Refractive Index');

for i=1:d
    if i>201 && i<501
    nfin(i,1)= n(i,1);
    else
    nfin(i,1)=0;
    end
end

navg=sum(nfin)/300;%Mean value of DSKK estimated refractive index over visible range (400 nm to 700 nm)
