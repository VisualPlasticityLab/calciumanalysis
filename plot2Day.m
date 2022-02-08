function plot2Day(Day1,Day2,xx,yy)


p1=load('D:\Box\Enhancement-Andrea\1847\329&404\1847_329_000_1&1847_329_003_1&&1847_404_000_1&1847_404_001_1_plane1_selected142.mat');


[SI1.baselineR,SI1.baselineSDR,SI1.peakR,SI1.errorR,...
            SI1.baselineS,SI1.baselineSDS,SI1.peakS,SI1.errorS,cal_ER_option]= sigFcmp(Day1.sigF,Day1.win_sig,Day1.matrix);  % sigR:seg,1,Var,ncell
        [SI1.OSI_R,SI1.gOSI_R,SI1.DSI_R, SI1.gDSI_R,SI1.pref_ori_R,SI1.pref_dir_R]=calOSI(SI1.peakR);
        [SI1.OSI_S,SI1.gOSI_S,SI1.DSI_S, SI1.gDSI_S,SI1.pref_ori_S,SI1.pref_dir_S]=calOSI(SI1.peakS);
        
        [SI2.baselineR,SI2.baselineSDR,SI2.peakR,SI2.errorR,...
            SI2.baselineS,SI2.baselineSDS,SI2.peakS,SI2.errorS,cal_ER_option]= sigFcmp(Day2.sigF,Day2.win_sig,Day2.matrix);  % sigR:seg,1,Var,ncell
        [SI2.OSI_R,SI2.gOSI_R,SI2.DSI_R, SI2.gDSI_R,SI2.pref_ori_R,SI2.pref_dir_R]=calOSI(SI2.peakR);
        [SI2.OSI_S,SI2.gOSI_S,SI2.DSI_S, SI2.gDSI_S,SI2.pref_ori_S,SI2.pref_dir_S]=calOSI(SI2.peakS);
        
        
               
Cor1=(1:171);
Cor2=(1:455);
day1_picked = [53 39 80 76 57 146 11 162 44 149];
day2_picked = [4 27 24 57 68 67 78 119 96 201];
        if ~figskip
            hsigF1=sigFplt(Day1.sigF,Day1.run.matrix,Day1.win_sig,Cor1(day1_picked),folder1);
            hsigF2=sigFplt(Day2.sigF,Day2.run.matrix,Day2.win_sig,Cor2(day2_picked),folder2);  % sigF:seg,rep,Var,ncell
        end

        OSI_change = SI2.OSI_R(pair2_final)- SI1.OSI_R(pair1_final); %OSI difference day5-day1
        pref_ori_change = SI2.pref_ori_R(pair2_final)- SI1.pref_ori_R(pair1_final); %prefered orinetation change from day5 to day1
        pref_dir_change = SI2.pref_dir_R(pair2_final)- SI1.pref_dir_R(pair1_final); %prefered direction change from day5 to day1
        SI2_mag=SI2.peakR(pair2_final); %magnitude common cells day5
        SI1_mag=SI1.peakR(pair1_final); %magnitude common cells day1
        mag_change = SI2_mag-SI1_mag; % magnitude change from day5 to day1

        %select cells with positive or negative dOSI
        OSI_change_positive = find(OSI_change>0); 
        OSI_change_negative = find(OSI_change<0);

        %select cells with positive or negative dmagnitude
        mag_change_positive = find(mag_change>0);
        mag_change_negative = find(mag_change<0);

        %select cells with positive or negative dprefered
        %orientaion/direction
        pref_ori_change_positive =  find(pref_ori_change>0);
        d0_pref_dir =  find(pref_dir_change==0); %prefered direction stays the same
        dorth_pref_dir = find(pref_dir_change==2);%prefered direction changes to orthogonal direction
        d1_pref_dir = find(pref_dir_change==1); 
        d3_pref_dir = find(pref_dir_change==3);
        dsmall_pref_dir = [d1_pref_dir,d3_pref_dir];

        pref_ori_change_negative =  find(pref_ori_change<0);
        pref_dir_change_negative =  find(pref_dir_change<0);

        [OSI_mag_positive]=intersect(OSI_change_positive,mag_change_positive); %cells with increased magnitude and OSI

        [pref_dir0_mag_positive]=intersect(d0_pref_dir,mag_change_positive); %cells with increased magnitude and pref dir
        
        [OSI_pref_ori_positive]=intersect(OSI_change_positive,pref_ori_change_positive);
        [OSI_pref_dir_positive]=intersect(OSI_change_positive,d0_pref_dir);

        [OSI_mag_negative]=intersect(OSI_change_negative,mag_change_negative);%cells with decreased magnitude and OSI
        
        [OSI_pref_ori_negative]=intersect(OSI_change_negative,pref_ori_change_negative);
        [OSI_pref_dir_negative]=intersect(OSI_change_negative,pref_dir_change_negative);

        [OSI_positive_mag_negative]=intersect(OSI_change_positive,mag_change_negative);%cells with increased OSI and decreased magnitude 

        [OSI_positive_pref_ori_negative]=intersect(OSI_change_positive,pref_ori_change_negative);
        [OSI_positive_pref_dir_negative]=intersect(OSI_change_positive,pref_dir_change_negative);

        [OSI_negative_mag_positive]=intersect(OSI_change_negative,mag_change_positive);%cells with decreased OSI and increased magnitude

        [OSI_negative_pref_ori_positive]=intersect(OSI_change_negative,pref_ori_change_positive);
        [OSI_negative_pref_dir_positive]=intersect(OSI_change_negative,d0_pref_dir);

        SI1_cat= cat(1,SI1.peakR(:,:,pair1_final));
        SI2_cat= cat(1,SI2.peakR(:,:,pair2_final));
        all_cat= [SI1_cat; SI2_cat];
        Cor= [1:142];
