function [H, W]=SCIhcv(betaTran,alphaTran,para)
for j=1:2
for l=1:2
H(j,l)=betaTran(j,l).*para(l,1).*para(j,2)./(para(j,3).*(para(l,4)+para(l,5)).*para(l,6));
W(j,l)=alphaTran(j,l).*(para(l,7)+para(l,3)+para(l,8))./((para(l,7)+para(l,3)).*(para(l,7)+para(l,3)+para(l,8)+para(l,9)));
end
end