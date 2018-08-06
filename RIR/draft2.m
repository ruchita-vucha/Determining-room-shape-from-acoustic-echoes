c = 340;                    % Sound velocity (m/s)
fs = 96000;                 % Sample frequency (samples/s)
r = [3 3 3 ; 1 1 1 ;  5 4 6 ; 2 8 4 ; 4 6 9];    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
s = [3 5 6];              % Source position [x y z] (m)
L = [6 10 12];  %max distance 32.3              % Room dimensions [x y z] (m)
beta = [0.4 0.4 0.4 0.4 0.4 0.4];                 % Reverberation time (s)
n = 4096;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = 1;                 % -1 equals maximum reflection order!
dim = 3;                    % Room dimension
orientation = 0;            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter

h = rir_generator(c, fs, r, s, L, beta, n, mtype, order, dim, orientation, hp_filter);

h1 = h(1, :);
h2 = h(2, :);
h3 = h(3, :);
h4 = h(4, :);
h5 = h(5, :);

%[pks1,locs1] = findpeaks(h(1,4096));

%disp(locs1);

plot(linspace(0,4096,4096),h1);


%p = 

%prev = 0;
%data = [25 8 15 5 6 10 10 3 1 20 7];
%pks = findpeaks(data);

%pks1 = findpeaks(h1,fs);
%pks2 = findpeaks(h2,fs);
%pks3 = findpeaks(h3,fs);
%pks4 = findpeaks(h4,fs);
%pks5 = findpeaks(h5,fs);
[pks1,locs1] = findpeaks(h1);
[pks2,locs2] = findpeaks(h2);
[pks3,locs3] = findpeaks(h3);
[pks4,locs4] = findpeaks(h4);
[pks5,locs5] = findpeaks(h5);


%plot(linspace(0,100,length(locs1)),locs1);


  [B1,I1] = sort(pks1,'descend');
  [B2,I2] = sort(pks2,'descend');
  [B3,I3] = sort(pks3,'descend');
  [B4,I4] = sort(pks4,'descend');
  [B5,I5] = sort(pks5,'descend');
  
  
  
  K1 = locs1(:,I1(:,1:7));
  K2 = locs2(:,I2(:,1:7));
  K3 = locs3(:,I3(:,1:7));
  K4 = locs4(:,I4(:,1:7));
  K5 = locs5(:,I5(:,1:7));
  
  K1 = c*K1*1/fs;
  K2 = c*K2*1/fs;
  K3 = c*K3*1/fs;
  K4 = c*K4*1/fs;
  K5 = c*K5*1/fs;
 
 
K1 = K1.*K1;
K2 = K2.*K2;
K3 = K3.*K3;
K4 = K4.*K4;
K5 = K5.*K5;


%K1 = K1(:,2:end);
%K2 = K2(:,2:end);
%K3 = K3(:,2:end);
%K4 = K4(:,2:end);
%K5 = K5(:,2:end);

%  x = v(7:end);

D = pdist(r,'euclidean');  % euclidean distance
D_Matrix = squareform(D);

 E = zeros(5,1);
 Y1 = zeros(6,1);
 Y = zeros(5,1);
%this is b{P,sstressasic
% 
 W = zeros(5,1);
 Z = zeros(1,1);
% 
% %loop
 for f1=1:7
 E(1,1) = K1(1,f1);
     %loop
 for f2=1:7
 E(2,1) = K2(1,f2);
     %loop
 for f3=1:7
 E(3,1) = K3(1,f3);
     %loop
 for f4=1:7
 E(4,1) = K4(1,f4);
     %loop
 for f5=1:7
 E(5,1) = K5(1,f5);
% 
% 
 W = horzcat(W,E);
% 
 C = horzcat(D_Matrix,E);
 Zr = zeros(1,1);
% 
 T = horzcat((E)',Zr);
 K = vertcat(C,T);
% 
% %[P,sstress] = mdscale(K,1);
% %Z  =  horzcat(Z,sstress);
% 
[P,sstress] = mdscale(K,3,'criterion','metricsstress');

Z = horzcat(Z,sstress);

Y1 = horzcat(Y1,P);

% 
% 
% R = rank(K);
% % 
% disp(R);
% % 
% if R<=5
%  Y = vertcat(Y,E);
% end
% 
% 
% 
[B0,I0] = sort(Z,'ascend');

% %loop
 end
% %loop
 end
% %loop
 end
% %loop
 end
% %loop
 end
% 
% %except first column of Y rest are candidates
% 
W = W(:,2:end);
Z = Z(:,2:end);
Y1 = Y1(:,2:end);

IO = I0(:,1:80);

WAR = abs(sqrt(W(:,IO)));