% Author: Sorique Aziz, Email: sxasps@rit.edu

clear all;
close all;

% Parameters
num = 6; 
r = 0.01; 
s = 0.074; 
Ntot = 100000; 
MaxS = 100000; 
t1 = 4000; 
MaxTraj = 10000; 
pMax = zeros(1, Ntot); 
p1 = zeros(1, Ntot); 
NumLags = 30; 

figure(1)

for j = 1:MaxTraj
    m = nan(num, MaxS);
    mtot = nan(1, MaxS);
    mmax = nan(1, MaxS);
    m(:, 1) = 1;
    mtot(1) = sum(m(:, 1));
    mmax(1) = max(m(:, 1));

    monomers = Ntot;
    T = zeros(1, MaxS);
    T(1) = 0;
    
    for i = 1:MaxS
        fil = randi([1 num]);
       
        s1 = s * m(fil, i);
        k2 = r*Ntot;
        k1 = s1;
    
        k0 = k1 + k2;
        CoinFlip1 = rand;
        tau(i) = (1 / k0) * log(1 / CoinFlip1);   
        T(i + 1) = T(i) + tau(i);

        CoinFlip2 = rand;
        if CoinFlip2 <= (k1 / k0)
            if m(fil, i) == 1
                m(fil, i + 1) = m(fil, i);
            else
                s1 = min(floor((CoinFlip2 / (k1 / k0)) * m(fil, i)), m(fil, i) - 1);
                m(fil, i + 1) = m(fil, i) - s1;
                monomers = min(monomers + s1, Ntot);
            end
        else
            growth_units = min(1, monomers);
            m(:, i + 1) = m(:, i);
            m(fil, i + 1) = m(fil, i + 1) + growth_units;
            monomers = max(monomers - growth_units, 0);
        end
    
        for newind = 1:num
            if newind ~= fil
                m(newind, i + 1) = m(newind, i);
            end
        end

        mmax(i + 1) = max(m(:, i + 1));

        if T(i + 1) >= t1
            break;
        end
    end

    p1(m(1,i+1)) = p1(m(1,i+1))+1;
    pMax(mmax(1,i+1)) = pMax(mmax(1,i+1))+1; 
    
    ta = T(1000:i);
    hbin = 0.2;
    tb = ta(1):hbin:ta(end);
    m1a = m(1,1000:i);
    mmaxa = mmax(1000:i);
    m1b = interp1(ta,m1a,tb);
    mmaxb = interp1(ta,mmaxa,tb);

    [acf1, lag1] = autocorr(m1b, NumLags);
    [acf, lag] = autocorr(mmaxb, NumLags);

    acf_vals1(j,:) = acf1;
    acf_vals(j,:) = acf;
end

plot(T, mmax*0.004, '-r', 'LineWidth', 2);
xlim([0 50])
xlabel('Time (s)');
ylabel('Bundle Length (micro meter)');
legend('Sev Sim','Location','Best');
box on;
hold off;

pMax = pMax/sum(pMax);
p1 = p1/sum(p1);
x = 1:1:Ntot;

L = 1:1000;
p = ((L .* s .* num) / (r .* Ntot)) .* (1 - exp(-(s/(2*r*Ntot)) * L.^2)).^(num - 1) .* exp(-(s/(2*r*Ntot)) * L.^2);

figure(2)
plot(x*0.004,pMax,'-','LineWidth', 2, 'Color', 'r','DisplayName', 'Sev Sim');
hold on
plot(L*0.004, p, 'LineWidth', 2, 'Color', 'k','DisplayName', 'Sev The');
xlim([0 3])
xlabel('Length (micro meter)');
ylabel('Probability distribution');
legend('Location', 'best');

acf_avg1 = sum(acf_vals1,1)/MaxTraj;
acf_stddev1 = std(acf_vals1,1,1);
acf_error1 = acf_stddev1/sqrt(MaxTraj);
lagtime = lag*hbin;

figure(3)
errorbar(lagtime, acf_avg1, acf_error1,".",'MarkerSize',30,'LineWidth',2, 'Color', 'r', 'DisplayName', 'Sev Sim');
hold on;
alpha_s = (1/num)*sqrt((8 * s * r * Ntot) / pi);
t = linspace(0, 10, 100); 
AC = exp(-alpha_s * t);
plot(t,AC,'LineWidth',2, 'Color', 'k', 'DisplayName', 'Sev The');
xlabel('Lagtime (s)');
ylabel('Autocorrelation');
legend('show','Location', 'best');
xlim([0 3])
