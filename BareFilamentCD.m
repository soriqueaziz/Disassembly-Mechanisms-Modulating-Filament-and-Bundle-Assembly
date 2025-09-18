%Author: Dr.Sorique Aziz, Email: sxasps@rit.edu
close all;
clear all;

num = 1; 
r = 0.3; 
d = 225; 
Ntot=1000;
MaxS=1000050; 
t1=2000; 
MaxTraj=1000; 

p1= zeros(1, Ntot); 
ptot = zeros(1,Ntot);
hbin_vals = zeros(1,MaxTraj); 
NumLags = 10; 
mb_all_trajectories = cell(MaxTraj, 1);

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
          
        k2 = r*(monomers-mtot(1)); 
        k1 = d;                         

        k0=k1+k2;
        CoinFlip1=rand;
        tau(i)=(1/k0)*log(1/CoinFlip1);
        T(i+1)= T(i)+tau(i);
            
        CoinFlip2=rand;

        if CoinFlip2<=(k1/k0)
            if m(fil,i)==1
                m(fil,i+1)= m(fil,i);
            else
                m(fil,i+1) = m(fil,i)-1;
                monomers = monomers + 1;
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
    plong = zeros(1,Ntot); 

    m1 = m(1,:);

    ta = T(5000:i);
    ma = m1(5000:i);
    ma1 = m1(i:1000);
    mtota = mtot(5000:i);

    hbin=1;
    hbin_vals(j) = hbin;
    tb = ta(1):hbin:ta(end);
    mb = interp1(ta,ma,tb);
    mtotb = interp1(ta,mtota,tb);

    [acf, lag] = autocorr(mb, NumLags);
    acf_vals(j,:) = acf;
end

acf_avg = sum(acf_vals,1)/MaxTraj;
acf_stddev = std(acf_vals,1,1);
acf_sqr = acf_avg.^2;

p1 = p1/sum(p1);
x = 1:1:Ntot;
Avg= sum(x.*p1);
variance= sum((x.^2).*p1)-(sum(x.*p1))^2;

T1= 0:0.1:t1;
Lss = Ntot-d/r;
L = Lss*(1- exp(-r.*T1));

x_theory = 0:Ntot;

log_P_0 = -Ntot * log(r / d) ...
          - d / r ...
          - log(gammainc(d / r, Ntot + 1, 'upper')) ...
          - gammaln(Ntot + 1);

log_P = x_theory * log(r / d) ...
        + gammaln(Ntot + 1) ...
        - gammaln(Ntot - x_theory + 1) ...
        + log_P_0;

P_theory = exp(log_P);

fig=figure(1);
plot(T,m1*0.004,'-','LineWidth', 2, 'Color', 'b', 'DisplayName', 'CD Sim');
hold on
plot(T1,L*0.004,'LineWidth', 2, 'Color', 'k', 'DisplayName', 'CD The');
xlim([0 100])
ylim([0 3])
xlabel('Time (s)');
ylabel('Length (micro meter)');
legend('Location', 'best');

figure(2)
plot(x*0.004,p1,'-','LineWidth', 2, 'Color', 'b', 'DisplayName', 'CD Sim');
hold on
plot(x_theory*0.004, P_theory, 'k-', 'LineWidth', 2, 'DisplayName', 'CD Theory');
xlabel('Length (micro meter)');
ylabel('Probability distribution');
legend('Location', 'best');

figure(3)
acf_error = acf_stddev/sqrt(MaxTraj);
lagtime = lag*hbin; 
errorbar(lagtime, acf_avg, acf_error,".",'MarkerSize',30,'LineWidth',2, 'Color', 'b', 'DisplayName', 'CD (Simulated Avg \pm SE)');
hold on
alpha = num*r;
AC = exp(-alpha.*lagtime);
lastColor = get(gca, 'ColorOrder');
plot(lagtime,AC,'LineWidth',2, 'Color', 'k', 'DisplayName', 'CD The');
xlim([0 10])
xlabel('Lagtime (s)');
ylabel('Autocorrelation');
legend('Location', 'best');



