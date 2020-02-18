clc;
clear;
close all;
% Read Input Image

InputImage=imread('download.jpg');

%Resize the Image

InputImage=imresize(InputImage,[256 256]);

% Display the Image

imshow(InputImage);

% Get Inputs from Mouse,Select 4 Seed Points in Image

[Col Row]=ginput(4);

c =Col;
r =Row;

% Select polygonal region of interest
BinaryMask = roipoly(InputImage,c,r);
figure, imshow(BinaryMask);title('Selected Region of Interest');

%Create Buffer for ROI
ROI=zeros(256,256);

%Create Buffer for NONROI
NONROI=zeros(256,256);
for i=1:256

for j=1:256

if BinaryMask(i,j)==1
ROI(i,j)=InputImage(i,j);

else
NONROI(i,j)=InputImage(i,j);
end

end

end

%Display ROI and Non ROI
figure;
subplot(1,2,1);imshow(ROI,[]);title('ROI');
subplot(1,2,2);imshow(NONROI,[]);title('NON ROI');


%---------------------------------------------------------
% Subtracting each image pixel value by clc;

I =uint8(ROI);
I1=I;
[row coln]= size(I);
I= double(I);
%---------------------------------------------------------
% Subtracting each image pixel value by 128
%--------------------------------------------------------
I = I - (128*ones(row,coln)); % uint8 [0,255] --> int8 [-128,127]

quality = input('Enter Quality of Compression[0->100] for ROI: ');

%----------------------------------------------------------
% Quality Matrix Formulation
%----------------------------------------------------------
Q50 = [ 16 11 10 16 24 40 51 61; % Quantization matrix [8X8]
12 12 14 19 26 58 60 55;
14 13 16 24 40 57 69 56;
14 17 22 29 51 87 80 62;
18 22 37 56 68 109 103 77;
24 35 55 64 81 104 113 92;
49 64 78 87 103 121 120 101;
72 92 95 98 112 100 103 99];
if quality > 50
QX = round(Q50.*(ones(8)*((100-quality)/50)));
QX = uint8(QX);
elseif quality < 50
QX = round(Q50.*(ones(8)*(50/quality)));
QX = uint8(QX);
elseif quality == 50
QX = Q50;
end
%----------------------------------------------------------
% Formulation of forward DCT Matrix and inverse DCT matrix
%----------------------------------------------
DCT_matrix8 = dct(eye(8));
iDCT_matrix8 = DCT_matrix8'; %inv(DCT_matrix8);

%----------------------------------------------------------
% Jpeg Compression
%----------------------------------------------------------
dct_restored = zeros(row,coln);
QX = double(QX);

%----------------------------------------------------------
% Jpeg Encoding
%----------------------------------------------------------
% Forward Discret Cosine Transform
%----------------------------------------------------------
for i1=[1:8:row]
for i2=[1:8:coln]
zBLOCK=I(i1:i1+7,i2:i2+7);
win1=DCT_matrix8*zBLOCK*iDCT_matrix8;
dct_domain(i1:i1+7,i2:i2+7)=win1;
end
end
%-----------------------------------------------------------
% Quantization of the DCT coefficients
%-----------------------------------------------------------
for i1=[1:8:row]
for i2=[1:8:coln]
win1 = dct_domain(i1:i1+7,i2:i2+7);
win2=round(win1./QX);
dct_quantized(i1:i1+7,i2:i2+7)=win2;
end
end


%-----------------------------------------------------------
% Jpeg Decoding
%-----------------------------------------------------------
% Dequantization of DCT Coefficients
%-----------------------------------------------------------
for i1=[1:8:row]
for i2=[1:8:coln]
win2 = dct_quantized(i1:i1+7,i2:i2+7);
win3 = win2.*QX;
dct_dequantized(i1:i1+7,i2:i2+7) = win3;
    end
end
%-----------------------------------------------------------
% Inverse DISCRETE COSINE TRANSFORM
%-----------------------------------------------------------
for i1=[1:8:row]
    for i2=[1:8:coln]
        win3 = dct_dequantized(i1:i1+7,i2:i2+7);
        win4=iDCT_matrix8*win3*DCT_matrix8;
        dct_restored(i1:i1+7,i2:i2+7)=win4;
    end
end
I2=dct_restored;


% ---------------------------------------------------------
% Conversion of Image Matrix to Intensity image
%----------------------------------------------------------
KROI=mat2gray(I2);


%----------------------------------------------------------
%Display of Results
%----------------------------------------------------------
figure(4);imshow(I1);title('original image ROI');
figure(5);imshow(KROI);title('Compressed Image ROI');

%BG work start
I =uint8(NONROI);
I1=I;
[row coln]= size(I);
I= double(I);
%---------------------------------------------------------
% Subtracting each image pixel value by 128
%--------------------------------------------------------
I = I - (128*ones(row,coln)); % uint8 [0,255] --> int8 [-128,127]

quality = input('Enter Quality of Compression for background[0->100]: ');

%----------------------------------------------------------
% Quality Matrix Formulation
%----------------------------------------------------------
Q50 = [ 16 11 10 16 24 40 51 61; % Quantization matrix [8X8]
12 12 14 19 26 58 60 55;
14 13 16 24 40 57 69 56;
14 17 22 29 51 87 80 62;
18 22 37 56 68 109 103 77;
24 35 55 64 81 104 113 92;
49 64 78 87 103 121 120 101;
72 92 95 98 112 100 103 99];
if quality > 50
QX = round(Q50.*(ones(8)*((100-quality)/50)));
QX = uint8(QX);
elseif quality < 50
QX = round(Q50.*(ones(8)*(50/quality)));
QX = uint8(QX);
elseif quality == 50
QX = Q50;
end
%----------------------------------------------------------
% Formulation of forward DCT Matrix and inverse DCT matrix
%----------------------------------------------
DCT_matrix8 = dct(eye(8));
iDCT_matrix8 = DCT_matrix8'; %inv(DCT_matrix8);

%----------------------------------------------------------
% Jpeg Compression
%----------------------------------------------------------
dct_restored = zeros(row,coln);
QX = double(QX);

%----------------------------------------------------------
% Jpeg Encoding
%----------------------------------------------------------
% Forward Discret Cosine Transform
%----------------------------------------------------------
for i1=[1:8:row]
for i2=[1:8:coln]
zBLOCK=I(i1:i1+7,i2:i2+7);
win1=DCT_matrix8*zBLOCK*iDCT_matrix8;
dct_domain(i1:i1+7,i2:i2+7)=win1;
end
end
%-----------------------------------------------------------
% Quantization of the DCT coefficients
%-----------------------------------------------------------
for i1=[1:8:row]
for i2=[1:8:coln]
win1 = dct_domain(i1:i1+7,i2:i2+7);
win2=round(win1./QX);
dct_quantized(i1:i1+7,i2:i2+7)=win2;
end
end


%-----------------------------------------------------------
% Jpeg Decoding
%-----------------------------------------------------------
% Dequantization of DCT Coefficients
%-----------------------------------------------------------
for i1=[1:8:row]
for i2=[1:8:coln]
win2 = dct_quantized(i1:i1+7,i2:i2+7);
win3 = win2.*QX;
dct_dequantized(i1:i1+7,i2:i2+7) = win3;
    end
end
%-----------------------------------------------------------
% Inverse DISCRETE COSINE TRANSFORM
%-----------------------------------------------------------
for i1=[1:8:row]
    for i2=[1:8:coln]
        win3 = dct_dequantized(i1:i1+7,i2:i2+7);
        win4=iDCT_matrix8*win3*DCT_matrix8;
        dct_restored(i1:i1+7,i2:i2+7)=win4;
    end
end
I2=dct_restored;


% ---------------------------------------------------------
% Conversion of Image Matrix to Intensity image
%----------------------------------------------------------
KNONROI=mat2gray(I2);


%----------------------------------------------------------
%Display of Results
%----------------------------------------------------------
figure(6);imshow(I1);title('original image NON ROI');
figure(7);imshow(KNONROI);title('Compressed Image NON ROI');

%COMBINEING OF ROI AND NONROI
C = imfuse(KNONROI,KROI,'blend','Scaling','joint');
%C=imfuse(KROI,KNONROI);
final=mat2gray(C);
figure(8);imshow(final);title('combined images');

%Calculate MSE and PSNR
n=size(I1);
M=n(1);
N=n(2);
I1=double(I1);
MSE = sum(sum((I1-final).^2))/(M*N);
PSNR = 10*log10(N*M/MSE);
fprintf('\nMSE: %7.2f ', MSE);
fprintf('\nPSNR: %9.7f dB', PSNR);