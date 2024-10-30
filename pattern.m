clear all,
filename = "patternfile.txt";
data = readmatrix(filename);
pattern_flatten = squeeze(data(:,3));
pattern_grid = 10*log10(reshape(pattern_flatten,181,360));
pattern_grid(182:360,:) = flip(pattern_grid(2:end-1,:),1);

theta = linspace(-pi,pi,360);
phi = deg2rad(linspace(0,359,360));

cut_phi = pattern_grid(91,:);
cut_theta = pattern_grid(:,180);

figure(1)
polarplot(phi,cut_phi)
rmin = min(cut_phi); rmax = max(cut_phi);
rlim([rmin rmax])
figure(2)
polarplot(theta,cut_theta)
rmin = min(cut_theta); rmax = max(cut_theta);
rlim([rmin rmax])
