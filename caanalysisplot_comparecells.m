pair1_data_folder = 'D:\Box\Enhancement-Andrea\1847\040417\000\1847_404_000_plane0\picked424_sp\';
pair12_data1_folder = 'D:\Box\Enhancement-Andrea\1847\040417\000\1847_404_000_plane0\picked433_sp\sp_wo_baseadj\sp_wo_baseadj';
pair12_data2_folder ='D:\Box\Enhancement-Andrea\1847\040417\001\1847_404_001_1_plane0\picked433_sp\sp_wo_baseadj';
pair2_data_folder = 'D:\Box\Enhancement-Andrea\1847\040417\001\1847_404_001_plane0\picked143_sp\sp_wo_baseadj';

p1=load('D:\Box\Enhancement-Andrea\1847\040417\000&001\plane 0\1847_404_000&1847_404_001_plane0_selected110.mat');
p2=load('D:\Box\Enhancement-Andrea\1847\040417\000&001+000\plane 0\1847_404_000_plane0&1847_404_000_1&1847_404_001_1_selected382.mat');


[pair1_final,idx1,idx2] = intersect(p1.pair1,p2.pair1);

pair2_final=p1.pair2(idx1);
pair12_final=p2.pair2(idx2);

picked =[2 3 5 8:10 11 13 14 15 17 20];
picked = [5 17];
caanalysisplot(pair1_data_folder,pair1_final(picked),1);% file 000 selected cells

caanalysisplot(pair12_data1_folder,pair12_final(picked),0);% first part of file 000&001 selected cells
caanalysisplot(pair2_data_folder,pair2_final(picked),0); % file 001 selected cells
caanalysisplot(pair12_data2_folder,pair12_final(picked),0); % 2nd part of file 000&001 selected cells

suite2p_000_001 = load('D:\Box\Enhancement-Andrea&Jennifer\1847\040417\suite2p\plane0\Fall.mat');
suite2p_001 =load('D:\Box\Enhancement-Andrea&Jennifer\1847\040417\001\suite2p\plane0\Fall.mat');

caanalysisplot(pair1_data_folder,pair1_final(1:10:100),1);
caanalysisplot(pair12_data1_folder,pair12_final(1:10:100),1);

caanalysisplot(pair2_data_folder,pair2_final(1:10:100),1);
caanalysisplot(pair12_data2_folder,pair12_final(1:10:100),1);



