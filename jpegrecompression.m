clc
clear all

testImagesDir = 'KodakTest_Images\BMPs';
jpegExeDir = 'libjpegturbo';
outputDir = 'temp_dir';

test_image = 'kodim03.bmp';
origBmpFile = fullfile(testImagesDir, test_image);
imageName  = test_image(1:(length(test_image)-4));

XY_Qvalues = 60:10:100;
m = length(XY_Qvalues);

YZ_Qvalues = 30:5:100;
n = length(YZ_Qvalues); % number  of quality_factors values


bppXY = zeros(1,m);
ssimXY = zeros(1,m);
psnrXY = zeros(1,m);
msssimXY = zeros(1,m);

bppL2 = zeros(m, 1,n);
ssimXZ = zeros(m, 1,n);
psnrXZ = zeros(m, 1,n);
msssimXZ = zeros(m,1,n);

ssimYZ = zeros(m, 1,n);
psnrYZ = zeros(m, 1,n);
msssimYZ = zeros(m, 1,n);
X_raw = imread(fullfile(testImagesDir, test_image));
[rows, cols] = size(X_raw);


jpg_encoder = fullfile(jpegExeDir, "cjpeg.exe");
jpg_decoder = fullfile(jpegExeDir, "djpeg.exe");


legendsXZ = strings(1,n);
legendsYZ = strings(1,n);

markersXZ = {'-o', '-+', '-*', '-s', '-x', '-d', ':^', '-v', '-->', ':<','-p','-h', '--o', ':+', '-.*'};
markersYZ = {':o', ':+', ':*', ':s', ':x', ':d', ':^', '-v', '-->', ':<','-p','-h', '--o', ':+', '-.*'};
figure(1)
hold on
for xy = 1:m
    XYq = XY_Qvalues(xy);
    L1_jpeg_file = fullfile(outputDir, [imageName '_XYq' num2str(XYq) '.jpg']);
    L1_decoded_bmp = fullfile(outputDir, [imageName '_XYq' num2str(XYq) '.bmp']);
    jpg_cmd = sprintf("%s -quality %d -outfile %s %s", jpg_encoder, XYq, L1_jpeg_file, origBmpFile);
    system(jpg_cmd);
    de_jpg_cmd = sprintf("%s -outfile %s %s", jpg_decoder, L1_decoded_bmp, L1_jpeg_file);
    system(de_jpg_cmd);
    Y_raw = imread(L1_decoded_bmp);
    sL1 = dir(L1_jpeg_file);
    bppXY(1,xy) = (sL1.bytes*8)/(rows*cols);
    psnrXY(1,xy) = psnr(Y_raw, X_raw);
    ssimXY(1,xy) = ssim(Y_raw, X_raw);
    msssimXY(1,xy) = mean(squeeze(multissim(Y_raw, X_raw)));
    disp([...
        'Image : ' imageName ...
        '    X-Y Q=> ' num2str(XYq)  ...
        '    bpp-X-Y =>' num2str(bppXY(1,xy))  ...
        '    psnr-X-Y => ' num2str( psnrXY(1,xy)) ...
        '    ssim-X-Y =>  ' num2str(msssimXY(1,xy))])
    
    for yz = 1:n
        YZq = YZ_Qvalues(yz);
        L2_jpeg_file = fullfile(outputDir, [imageName '_XYq' num2str(XYq) '_YZq' num2str(YZq) '.jpg']);
        L2_decoded_bmp = fullfile(outputDir, [imageName '_XYq' num2str(XYq) '_YZq' num2str(YZq) '.bmp']);
        jpg_cmd_L2 = sprintf("%s -quality %d -outfile %s %s", jpg_encoder, YZq, L2_jpeg_file, L1_decoded_bmp);
        system(jpg_cmd_L2);
        de_jpg_cmd_L2 = sprintf("%s -outfile %s %s", jpg_decoder, L2_decoded_bmp, L2_jpeg_file);
        system(de_jpg_cmd_L2);
        
        Z_raw = imread(L2_decoded_bmp);
        sL1 = dir(L2_jpeg_file);
        bppL2(xy, 1,yz) = (sL1.bytes*8)/(rows*cols);
        
        psnrXZ(xy, 1,yz) = psnr(Z_raw, X_raw);
        ssimXZ(xy, 1,yz) = ssim(Z_raw, X_raw);
        msssimXZ(xy, 1,yz) = mean(squeeze(multissim(Z_raw, X_raw)));
        
        psnrYZ(xy, 1,yz) = psnr(Z_raw, Y_raw);
        ssimYZ(xy, 1,yz) = ssim(Z_raw, Y_raw);
        msssimYZ(xy, 1,yz) = mean(squeeze(multissim(Z_raw, Y_raw)));
        disp([...
            'Image : ' imageName ...
            '    Y-Z Q  => ' num2str(YZq)  ...
            '    bpp-Y-Z  =>' num2str(bppL2(1,yz))  ...
            '    psnr-X-Z => ' num2str( psnrXZ(xy, 1,yz)) ...
            '    msssim-X-Z =>  ' num2str(msssimXZ(xy, 1,yz))...
            '    psnr-X-Y =>  ' num2str(psnrYZ(xy, 1, yz)) ...
            '    ms-ssim-X-Y =>  ' num2str(msssimYZ(xy, 1, yz)) ...
            ])
    end
    legendsXZ = ['X-Z--XYq' num2str(XYq)];
    legendsYZ = ['Y-Z -XYq' num2str(XYq)];
    
    bpp = squeeze(bppL2(xy,1,:));
    psnrXYp = squeeze(psnrXZ(xy,1,:));
    psnrYZp = squeeze(psnrYZ(xy,1,:));
    plot(bpp, psnrXYp, markersYZ{xy},'MarkerSize',8, 'LineWidth',2, 'DisplayName',legendsXZ)
    plot(bpp, psnrYZp, markersXZ{xy}, 'MarkerSize',8, 'LineWidth',2, 'DisplayName',legendsYZ)
end
bppXYp = squeeze(bppXY);
psnrXYp = squeeze(psnrXY);
plot(bppXYp, psnrXYp, '-p' ,'MarkerSize',8, 'LineWidth',2, 'DisplayName','X-Y')

legend
grid on
grid minor
title(imageName)
xlabel('bpp','FontSize',14,'FontWeight','bold')
ylabel('psnr (dB)', 'FontSize',14,'FontWeight','bold')
hold off

