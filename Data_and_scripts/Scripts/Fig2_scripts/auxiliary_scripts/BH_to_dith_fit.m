function err=BH_to_dith_fit(l)

global nm_y dit_y BH_y ox_y

x00 = [350 700];
y00 = [l(1) l(2)];
c00 = [[1; 1]  x00(:)]\y00(:);                        % Calculate Parameter Vector
slope_m = c00(2);
intercept_b = c00(1);
c=l(3)*1e-3;
kox=l(4);
c=c/(1-kox);
BH1=(BH_y-(nm_y*slope_m + intercept_b)-kox*ox_y)/(1-kox);
BH1=BH1/c;

err=sqrt(sum((BH1-dit_y).^2));