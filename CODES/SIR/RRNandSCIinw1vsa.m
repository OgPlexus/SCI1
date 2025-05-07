clc;clear all; close all;
N1=200;
N2=200; 
w2=0.8;
w0=0.9;
%w0=linspace(-1,1,100);%If you need to test different parameter spaces.
a=linspace(-1,1,1000);
w1=linspace(0,1,1000);
beta=1;gam=1; % You can choose different values, but they will cancel out.
N=N1+N2;
  for i=1:length(w1)
     for j=1:length(a)
         A1= (1+w0).*(1+a(j))./N1;
         A2= (1-w0).*(1-a(j))./N2;
         sqPartC=sqrt((A1.*w1(i)+A2.*w2).^2+4*A1.*A2.*(1-w1(i)-w2));
         A(i,j)=(N/4)*((A1*w1(i)+A2*w2)+sqPartC);
         if A(i,j)>1
             flag1(i,j)=2;
         elseif A(i,j)==1
             flag1(i,j)=1;
         else
             flag1(i,j)=0;
         end
          C1(i,j)=(N*beta/(4*gam))*((A1*w1(i)-A2*w2)+sqPartC)./A(i,j);
          C2(i,j)=(N*beta/(4*gam))*((-A1*w1(i)+A2*w2)+sqPartC)./A(i,j);
          if C1(i,j)>C2(i,j)
             flag2(i,j)=2;
         elseif C1(i,j)==C2(i,j)
             flag2(i,j)=1;
         else
             flag2(i,j)=0;
         end
     end
 end
 
 % for i=1:length(w0)
 %     for j=1:length(a)
 %         A1= (1+w0(i)).*(1+a(j))./N1;
 %         A2= (1-w0(i)).*(1-a(j))./N2;
 %         sqPartC=sqrt((A1.*w1+A2.*w2).^2+4*A1.*A2.*(1-w1-w2));
 %         A(i,j)=(N/4)*((A1*w1+A2*w2)+sqPartC);
 %         if A(i,j)>1
 %             flag1(i,j)=1;
 %         else
 %             flag(i,j)=0;
 %         end
 %         C1(i,j)=(N*beta/(4*gam))*((A1*w1-A2*w2)+sqPartC);
 %         C2(i,j)=(N*beta/(4*gam))*((-A1*w1+A2*w2)+sqPartC);
 %     end
 % end


colormap_custom = [
    0.8 0.8 0.8; % Color for values <= 1 (e.g., light gray)
    1.0 0.0 0.0; % Color for values > 1 (e.g., red)
];


% Define the color points
cmapPoints = [-2, -0.2,-0.0001, 0.0001,0.2, 2];
colorValues = [1, 0.4, 0;   % Orange for values <= -2
                1, 0.9, 0;   % Orange for values <= -2
                1, 1, 1;   % Orange for values <= -2
               1, 1, 1;     % Gray for -2 < values <= -1
               0.9, 0.7, 1;  % Light purple for 0 < values <= 1
               0.4, 0, 0.4]; % Purple for values > 1

% Set the number of colors for the colormap
numColors = 256;

% Create the custom colormap
cmap = zeros(numColors, 3);
for i = 1:3
    cmap(:, i) = interp1(cmapPoints, colorValues(:, i), linspace(min(cmapPoints), max(cmapPoints), numColors), 'pchip');
end

% Ensure colormap values are within [0, 1] range
cmap = max(0, min(1, cmap));

%% 
close all; % RRN values in the \{(w_1, a)\} parameter space
figure;
imagesc(w1, a, real(A)');hold on;
contour(w1,a,real(A)','k','linew',1.5)
  set(gca, 'YDir', 'normal');
xlabel('$w_1$');
ylabel('$a$');
  colormap(colormap_custom);
 colorbar
  clim([0,2])
hfig = gcf;
%clc;
fname = 'FlagAw0_0p1W20p8N1_350';
  %%
  close all; % C1-C2 (SCI difference) values in the \{(w_1, a)\} parameter space
figure;
  imagesc(w1,a,real((C1-C2))');hold on;
  contour(w1,a,real((C1-C2))','k','linew',1.5)
  set(gca, 'YDir', 'normal');
xlabel('$w_1$');
ylabel('$a$');
  colormap(cmap);
  %clim([-0.01,0.01])
   colorbar
hfig = gcf;
clc;
fname = 'C1mC2w0_0p9W20p8';
%%

%

% Mark regions with rectangles

text(0.3, 0.5, '$\mathbf{\Omega_1}$', 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'k', 'Interpreter', 'latex');
text(0.5, -0.7, '$\mathbf{\Omega_2}$', 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'k', 'Interpreter', 'latex');
text(0.9, 0.6, '$\mathbf{\Omega_3}$', 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'k', 'Interpreter', 'latex');


