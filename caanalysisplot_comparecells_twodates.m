
%%
% compare the reponses in awake mice and awake mice between the 2
% dates

%compare 329_000 with 404_000
pair1_data_folder = 'D:\Box\Enhancement-Andrea\1847\032917\000\data\1847_329_000_1\1847_329_000_1\picked171_sp+ca'; %Peak for 03/29 000                                        03
pair2_data_folder = 'D:\Box\Enhancement-Andrea\1847\040417\000\1847_404_000_1\picked455_sp+ca'; %Peak for 04/04 000

%for the lab PC
p1=load('D:\Box\Enhancement-Andrea\1847\329&404\1847_329_000_1&1847_329_003_1&&1847_404_000_1&1847_404_001_1_plane1_selected142.mat');
%for the home PC
%p1=load('C:\Users\UCL053324\Box\Enhancement-Andrea\1847\329&404\1847_404_000&1847_404_001&&1847_329_000&1847_329_003_plane0_selected142.mat');

pair1_final = p1.pair1; % 329 000&003 pair 1
pair2_final=p1.pair2; %404 000&001 pair2



% pickedi = [9 19 20 35 40 43 50,60 65 70 84 96 109 130]
% pickedm = [2 22 26 39 63 64 80 86 116 133]
% pickedd = [14 17 69 75 76 92 121]
% pickedb = [52 62 139]

%prefered orientaion 
caanalysisplot_selectedcells(pair1_data_folder,pair1_final(OSI_pref_ori_positive),1); %file 329 000/3 selected cells
caanalysisplot_selectedcells(pair2_data_folder,pair2_final(OSI_pref_ori_positive),1);% file 404 000/1 selected cells

%prefered direction positive both
caanalysisplot_selectedcells(pair1_data_folder,pair1_final(OSI_pref_dir_positive),1); %file 329 000/3 selected cells
caanalysisplot_selectedcells(pair2_data_folder,pair2_final(OSI_pref_dir_positive),1);% file 404 000/1 selected cells

%magnitude change positive both
caanalysisplot_selectedcells(pair1_data_folder,pair1_final(OSI_mag_positive),1); %file 329 000/3 selected cells
caanalysisplot_selectedcells(pair2_data_folder,pair2_final(OSI_mag_positive),1);% file 404 000/1 selected cells

%magnitude change negative both
caanalysisplot_selectedcells(pair1_data_folder,pair1_final(OSI_mag_negative),1); %file 329 000/3 selected cells
caanalysisplot_selectedcells(pair2_data_folder,pair2_final(OSI_mag_negative),1);% file 404 000/1 selected cells

%positive OSI negative magnitude
caanalysisplot_selectedcells(pair1_data_folder,pair1_final(OSI_positive_mag_negative),1); %file 329 000/3 selected cells
caanalysisplot_selectedcells(pair2_data_folder,pair2_final(OSI_positive_mag_negative),1);% file 404 000/1 selected cells

%negative OSI positive magnitude
caanalysisplot_selectedcells(pair1_data_folder,pair1_final(OSI_negative_mag_positive),1); %file 329 000/3 selected cells
caanalysisplot_selectedcells(pair2_data_folder,pair2_final(OSI_negative_mag_positive),1);% file 404 000/1 selected cells

