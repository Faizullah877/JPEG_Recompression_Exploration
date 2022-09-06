clc
clear all

testImagesDir = 'KodakTest_Images\BMPs';
jpegExeDir = 'libjpegturbo';
outputDir = 'temp_dir';

test_image = 'kodim02.bmp';
origBmpFile = fullfile(testImagesDir, test_image);
imageName  = test_image(1:(length(test_image)-4));

L1Qvalues = 70:10:100;
m = length(L1Qvalues);

qualityFactor = 30:10:100;
n = length(qualityFactor); % number  of quality_factors values

ssimL1 = zeros(1,m);
psnrL1 = zeros(1,m);
bppL1 = zeros(1,m);
msssimL1 = zeros(1,m);

bppL2 = zeros(m, 1,n);
ssimGT = zeros(m, 1,n);
psnrGT = zeros(m, 1,n);
msssimGT = zeros(m,1,n);

ssimRC = zeros(m, 1,n);
psnrRC = zeros(m, 1,n);
msssimRC = zeros(m, 1,n);
GT_img = imread(fullfile(testImagesDir, test_image));
[rows, cols] = size(GT_img);


jpg_encoder = fullfile(jpegExeDir, "cjpeg.exe");
jpg_decoder = fullfile(jpegExeDir, "djpeg.exe");


legendsGT = strings(1,n);
legendsRC = strings(1,n);

markers = {'-o', '-+', '-*', '-..', '-x',':s', '--d', ':^', '-v', '-->', ':<','-p','-h', '--o', ':+', '-.*'};
figure(1)
hold on
for QL1 = 1:m
    L1q = L1Qvalues(QL1);
    L1_jpeg_file = fullfile(outputDir, [imageName '_L1q' num2str(L1q) '.jpg']);
    L1_decoded_bmp = fullfile(outputDir, [imageName '_L1q' num2str(L1q) '.bmp']);
    jpg_cmd = sprintf("%s -quality %d -outfile %s %s", jpg_encoder, L1q, L1_jpeg_file, origBmpFile);
    system(jpg_cmd);
    de_jpg_cmd = sprintf("%s -outfile %s %s", jpg_decoder, L1_decoded_bmp, L1_jpeg_file);
    system(de_jpg_cmd);
    L1_recImg = imread(L1_decoded_bmp);
    sL1 = dir(L1_jpeg_file);
    bppL1(1,QL1) = (sL1.bytes*8)/(rows*cols);
    psnrL1(1,QL1) = psnr(L1_recImg, GT_img);
    ssimL1(1,QL1) = ssim(L1_recImg, GT_img);
    msssimL1(1,QL1) = mean(squeeze(multissim(L1_recImg, GT_img)));
    disp([...
        'Image : ' imageName ...
        '    L1Q => ' num2str(L1q)  ...
        '    bpp  =>' num2str(bppL1(1,QL1))  ...
        '    psnr => ' num2str( psnrL1(1,QL1)) ...
        '    ssim =>  ' num2str(msssimL1(1,QL1))])
    
    for QL2 = 1:n
        L2q = qualityFactor(QL2);
        L2_jpeg_file = fullfile(outputDir, [imageName '_L1q' num2str(L1q) '_L2q' num2str(L2q) '.jpg']);
        L2_decoded_bmp = fullfile(outputDir, [imageName '_L1q' num2str(L1q) '_L2q' num2str(L2q) '.bmp']);
        jpg_cmd_L2 = sprintf("%s -quality %d -outfile %s %s", jpg_encoder, L2q, L2_jpeg_file, L1_decoded_bmp);
        system(jpg_cmd_L2);
        de_jpg_cmd_L2 = sprintf("%s -outfile %s %s", jpg_decoder, L2_decoded_bmp, L2_jpeg_file);
        system(de_jpg_cmd_L2);
        
        L2_recImg = imread(L2_decoded_bmp);
        sL1 = dir(L2_jpeg_file);
        bppL2(QL1, 1,QL2) = (sL1.bytes*8)/(rows*cols);
        
        psnrGT(QL1, 1,QL2) = psnr(L2_recImg, GT_img);
        ssimGT(QL1, 1,QL2) = ssim(L2_recImg, GT_img);
        msssimGT(QL1, 1,QL2) = mean(squeeze(multissim(L2_recImg, GT_img)));
        
        psnrRC(QL1, 1,QL2) = psnr(L2_recImg, L1_recImg);
        ssimRC(QL1, 1,QL2) = ssim(L2_recImg, L1_recImg);
        msssimRC(QL1, 1,QL2) = mean(squeeze(multissim(L2_recImg, L1_recImg)));
        disp([...
            'Image : ' imageName ...
            '    L2Q => ' num2str(L2q)  ...
            '    bpp  =>' num2str(bppL2(1,QL2))  ...
            '    psnrGT => ' num2str( psnrGT(QL1, 1,QL2)) ...
            '    msssimGT =>  ' num2str(msssimGT(QL1, 1,QL2))...
            '    psnrRC =>  ' num2str(psnrRC(QL1, 1, QL2)) ...
            '    ms-ssimRC =>  ' num2str(msssimRC(QL1, 1, QL2)) ...
            ])
    end
    legendsGT = [imageName '-recQ' num2str(L1q) '-GT'];
    legendsRC = [imageName '-recQ' num2str(L1q) '-rec'];
    
    bpp = squeeze(bppL2(QL1,1,:));
    psnrGTa = squeeze(psnrGT(QL1,1,:));
    psnrRCa = squeeze(psnrRC(QL1,1,:));
    plot(bpp, psnrGTa, markers{QL1}, 'DisplayName',legendsGT)
    plot(bpp, psnrRCa, markers{QL1}, 'DisplayName',legendsRC)
%     set(leg, 'Interpreter', 'none')
end
bppl1aa = squeeze(bppL1);
psnrl1a = squeeze(psnrL1);
plot(bppl1aa, psnrl1a, markers{1}, 'DisplayName',imageName)
hold off

