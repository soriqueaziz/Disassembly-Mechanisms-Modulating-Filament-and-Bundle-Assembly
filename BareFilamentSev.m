%Author: Sorique Aziz, Email: sxasps@rit.edu

close all;
clear all;

MaxS=100000; 
t1=5000; 
MaxTraj=1000; 
num = 1; 
r = 0.3; 
s = 0.0075; 
Ntot=1000;

p1= zeros(1, Ntot*1000); 
ptot = zeros(1,Ntot*1000);
plong = zeros(1,Ntot*1000); 
hbin_vals = zeros(1,MaxTraj); 
NumLags = 40; 

for j=1:MaxTraj
    
    m = nan(num,MaxS);
    mtot = nan(1,MaxS);
    m(:,1) = 1;
    mtot(1) = sum(m(:,1));

    mlong = nan(1,MaxS); 
    mlong(1) = max(m(:,1));

    monomers=Ntot; 
    T=zeros(1,MaxS); 
    T(1)=0; 
    
    for i=1:MaxS

        fil = randi([1 num]);
        s1 = s*m(fil,i); 
        k2 = r*Ntot; 
        k1 = s1;
                
        k0=k1+k2;
        CoinFlip1=rand;
        tau(i)=(1/k0)*log(1/CoinFlip1); 
        T(i+1)= T(i)+tau(i);
            
        CoinFlip2=rand;

        if CoinFlip2<=(k1/k0)
            if m(fil,i)==1
                m(fil,i+1)= m(fil,i);
            else
                s1 = floor((CoinFlip2/(k1/k0))*m(fil,i));    
                m(fil,i+1)=m(fil,i) - s1;
                monomers = monomers + s1;
            end
        else
            m(fil,i+1)=m(fil,i)+1;
            monomers = monomers - 1;        
        end

        for newind = 1:num
            if newind ~= fil
                m(newind,i+1)=m(newind,i);
            end
        end

        mtot(i+1)=sum(m(:,i+1),1);
        mlong(i+1)=max(m(:,i+1));

        if T(i+1)>=t1
            break;
        end
    end

    p1(m(1,i+1)) = p1(m(1,i+1))+1; 
    ptot(mtot(i+1))= ptot(mtot(i+1))+1; 
    plong(mlong(i+1)) = plong(mlong(i+1))+1; 

    m1 = m(1,:); 

    ta = T(10000:i);
    ma = m1(10000:i);
    mtota = mtot(10000:i);

    hbin = 0.3;
    hbin_vals(j) = hbin;
    tb = ta(1):hbin:ta(end);
    mb = interp1(ta,ma,tb);
    mtotb = interp1(ta,mtota,tb);

    [acf, lag] = autocorr(mb,NumLags);
    acf_vals(j,:) = acf;
end

acf_avg = sum(acf_vals,1)/MaxTraj;
acf_stddev = std(acf_vals,1,1);
acf_sqr = acf_avg.^2;

p1 = p1/sum(p1);
x = 1:1:Ntot*1000;
Avg = sum(x.*p1);
variance= sum((x.^2).*p1)-(sum(x.*p1))^2;

T1= 0:0.1:t1;
L = (2^(1/2)*Ntot^(1/2)*r^(1/2)*tanh((2^(1/2)*Ntot^(1/2)*T1*r^(1/2)*s^(1/2))/2))/s^(1/2);

x_theory = linspace(0, Ntot, 1000);
P_theory = (s .* x_theory) ./ (r * Ntot) .* exp(- (s .* x_theory.^2) ./ (2 * r * Ntot));

acf_error = acf_stddev/sqrt(MaxTraj);
lagtime = lag*hbin; 
Lss = sqrt(2*r*Ntot/s);
alpha = sqrt((8 * s * r * Ntot) / pi);
AC = exp(-alpha.*lagtime);

% === Figures ===

figure(1);
plot(T,m1*0.004,'-','LineWidth', 2, 'Color', 'r', 'DisplayName', 'Sev Sim');
hold on
plot(T1,L*0.004,'LineWidth', 2, 'Color', 'k', 'DisplayName', 'Sev The');
xlim([0 100])
ylim([0 5])
xlabel('Time (s)');
ylabel('Length (micro meter)');
legend('Location', 'best');

figure(2)
plot(x*0.004,p1,'-','LineWidth', 2, 'Color', 'r', 'DisplayName', 'Sev Sim');
hold on
plot(x_theory*0.004, P_theory, 'k-', 'LineWidth', 2, 'DisplayName', 'Sev The');
xlim([0 3])
xlabel('Length (micro meter)');
ylabel('Probability distribution');
legend('Location', 'best');

figure(3)
errorbar(lagtime, acf_avg, acf_error,".",'MarkerSize',30,'LineWidth',2,'Color', 'r', 'DisplayName', 'Sev (Simulated Avg \pm SE)');
hold on
plot(lagtime,AC,'LineWidth',2, 'Color', 'k', 'DisplayName', 'Sev The');
xlim([0 5])
xlabel('Lagtime (s)');
ylabel('Autocorrelation');
legend('Location', 'best');

