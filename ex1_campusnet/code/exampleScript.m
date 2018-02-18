% Dette script er et eksemple p� brug af funktioner til f�rste opgave i
% matematisk modllering - 02526
% 
% Anders Bjorholm Dahl
% Februar 2015
% 

% Folder hvor mine datafiler ligger
dirIn = '/Users/abd/Documents/Teaching/02526_matematisk_modellering/Exercises/data/toCampusnet/';


%% Eksempel p� at loade multispektralt billede

[multiIm, annotationIm] = loadMulti('multispectral_day20.mat' , 'annotation_day20.png', dirIn);

% multiIm er et multispektralt billede - her ses dimensioner
size(multiIm)

%% her vises det billed ved spektralb�nd 7

figure
imagesc(multiIm(:,:,7))

%% annotationIm er et bin�rt billede med tre dimensioner

size(annotationIm)

%% I hver dimesion er et bin�rt billede 1 - baggrund med p�lse, 2 - fedtannotering, 
% 3 - k�dannotering. Her vises k�dannoteringen

figure
imagesc(annotationIm(:,:,2))

%% Funktionen getPix udtr�kker multispektrale pixels for annoteringerne 

% Her er et eksemple p� k�d- og fedtannoteringerne
[fatPix, fatR, fatC] = getPix(multiIm, annotationIm(:,:,2));
[meatPix, meatR, meatC] = getPix(multiIm, annotationIm(:,:,3));

% Her plottes middelv�rdierne for pixels med k�d hhv. fedt
figure
plot(mean(meatPix), 'b');
hold on
plot(mean(fatPix), 'r');


%% Funktionen showHistogram laver et histogram som den returnerer som en
% vektor

% Her er et eksempel - sidste argument g�r at funktionen plotter
% histogrmmet for k�d og fedt
h = showHistograms(multiIm, annotationIm(:,:,2:3), 3, 1);
axis square

%% histogrammet er ogs� givet i h - h er altid 256 dimensioner, og da der
% kun er lave pixelv�rdier viser jeg kun de f�rste 50

figure
plot(h(1:50,:))

%% Funktionen setImagePix giver et farvebillede, hvor pixel-koordinater
% gives som argument.

% RGB billede loades
imRGB = imread([dirIn 'color_day20.png']);

% Pixel-koordinater for fedtannoteringerne f�s med getPix
[fatPix, fatR, fatC] = getPix(multiIm, annotationIm(:,:,2));

% Pixel-koordinater s�ttes sammen til en matrix
pixId = [fatR, fatC];

% Det nye billede laves
rgbOut = setImagePix(imRGB, pixId);
figure
imagesc(rgbOut)


