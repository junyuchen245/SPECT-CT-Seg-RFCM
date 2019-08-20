import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy import ndimage, misc
import ctypes
import numpy.ctypeslib as ctl
from sys import exit, argv

################### Python functions #######################
def get_duration(mask, direction):
    Begin = 0; End = 0; flag = 0
    if(direction == 'z'):
        tmp_shape = np.shape(mask)
        for i in range(0,tmp_shape[2]):
            tmp = np.reshape(mask[:,:,i],[tmp_shape[0], tmp_shape[1]])
            if flag == 0 and np.sum(tmp) != 0:
                flag = 1
                Begin = i
            if flag == 1 and np.sum(tmp) == 0:
                End = i
                break
    if(direction == 'y'):
        tmp_shape = np.shape(mask)
        for i in range(0,tmp_shape[1]):
            tmp = np.reshape(mask[:,i,:],[tmp_shape[0], tmp_shape[2]])
            if flag == 0 and np.sum(tmp) != 0:
                flag = 1
                Begin = i
            if flag == 1 and np.sum(tmp) == 0:
                End = i
                break
    if(direction == 'x'):
        tmp_shape = np.shape(mask)
        for i in range(0,tmp_shape[0]):
            tmp = np.reshape(mask[i,:,:],[tmp_shape[1], tmp_shape[2]])
            if flag == 0 and np.sum(tmp) != 0:
                flag = 1
                Begin = i
            if flag == 1 and np.sum(tmp) == 0:
                End = i
                break
    if End == 0:
        tmp_shape = np.shape(mask)
        End = tmp_shape[0]
    return Begin, End

def extract_reg_w_same_label(label_map_in, seedx, seedy, seedz):
    label_map = np.zeros(label_map_in.shape)
    #print(label_map_in.shape)
    #print([seedx, seedy, seedz])
    label_map[label_map_in == label_map_in[seedx, seedy, seedz]] = 1
    return label_map

def calc_vol(truth):
    vol = np.sum(truth)*(np.power(2.2,3)) #mm^3
    return vol

def dice_eval(seg_result,gt_bw):
    # get ground truth bw image within the ROI
    dice_coef = 2 * np.sum(np.multiply(gt_bw,seg_result)) / (np.sum(gt_bw) + np.sum(seg_result))
    return dice_coef

def uptake_ratio(img, bone_mask, tumor_mask):
    theta_ratio = uptake_mean(img, tumor_mask)/uptake_mean(img, bone_mask)
    return theta_ratio

def uptake_mean(img, mask):
    img = np.multiply(img,mask)
    img = img.flatten()
    theta = np.mean(img[np.nonzero(img)])
    return theta

def uptake_eval(img, seg_result, gt_bw, string):
    theta_gt  = uptake_mean(img, gt_bw)
    theta_seg = uptake_mean(img, seg_result)
    if   string == 'BIAS':
        out = abs(theta_seg - theta_gt)
    elif string == 'VAR':
        out = theta_seg
    elif string == 'MSE':
        out = np.power(abs(theta_seg - theta_gt),2)
    return out

# eval_seg_3D_cluster_v1.py beta no_noise_realizations_in
beta_in = float(argv[1])
no_noise_realization_in = float(argv[3])
gamma   = float(argv[2])

################ Import C functions ################
priorlib=ctl.ctypes_load_library('libsegment', './')
# MRF_EM
c_func_FCM = priorlib.FCM
c_func_FCM.restype=ctypes.c_int
c_func_FCM.argtypes= \
	[ctl.ndpointer(np.float32,flags='c_contiguous'), #voxels
	 ctl.ndpointer(np.int16, flags='c_contiguous'), #output labels
	 ctypes.c_int, #xdim
	 ctypes.c_int, #ydim
	 ctypes.c_int, #zdim
	 ctypes.c_int, #nclasses
	 ctypes.c_double, # beta
     ctypes.c_double, # q
    ]
#print(c_func_MRF_EM.argtypes)

# MRF_EM_ctinfo
priorlib=ctl.ctypes_load_library('libsegment', './')
c_func_FCM_ctinfo = priorlib.FCM_ctinfo
c_func_FCM_ctinfo.restype=ctypes.c_int
c_func_FCM_ctinfo.argtypes= \
	[ctl.ndpointer(np.float32, flags='c_contiguous'), #voxels
	 ctl.ndpointer(np.int16, flags='c_contiguous'), #output labels
     ctl.ndpointer(np.int16, flags='c_contiguous'), #CT prior labels
	 ctypes.c_int, #xdim
	 ctypes.c_int, #ydim
	 ctypes.c_int, #zdim
	 ctypes.c_int, #nclasses
	 ctypes.c_double, # beta
     ctypes.c_double, # q
     ctypes.c_double, # gamma
    ]
#print(c_func_MRF_EM_ctinfo.argtypes)

# region growing
c_func_region_grow_cropped3Dimg = priorlib.region_grow_cropped3Dimg
c_func_region_grow_cropped3Dimg.restype=ctypes.c_int
c_func_region_grow_cropped3Dimg.argtypes= \
	[ctl.ndpointer(np.int16, flags='c_contiguous'), #input labels
	 ctypes.c_int, #xdim
	 ctypes.c_int, #ydim
	 ctypes.c_int, #zdim
	 ctypes.c_double, # seedx
     ctypes.c_double, # seedy
     ctypes.c_double, # seedz
     ctypes.c_int, #nclasses
    ]
#print(c_func_region_grow_cropped3Dimg.argtypes)

################ Evaluation Starts Here #######################
#cluster directory:

directory          = '/scratch/jchen/'
ex_drive_directory = '/scratch/jchen/eval_results/'

#mac directory:
#directory          = '/Volumes/My Passport/Research/Seg_eval/DMIP/'
#ex_drive_directory = '/Volumes/My Passport/Research/Seg_eval/'

################ load truth images ###################
# lesion
les_truth = sio.loadmat(directory+'truth/voimask_d.mat')
les_truth = les_truth['out']
#lt.figure(1); plt.imshow(les_truth[:,:,50]); plt.show()

# cort bone
cort_truth = sio.loadmat(directory+'truth/cortbone_d.mat')
cort_truth = cort_truth['out']

# trab bone
trab_truth = sio.loadmat(directory+'truth/trab_truth_d.mat')
trab_truth = trab_truth['out']

# bone
bone_truth = trab_truth + cort_truth - les_truth
bone_truth[bone_truth <= 0] = 0
bone_truth[bone_truth != 0] = 1

################ load CT images ###################
CTimg = sio.loadmat(directory+'CTImages/simulatedCT.5.mat')
CTimg = CTimg['out']

################ begin evaluation ###################
# initialzation
no_lesion = 170
fuzziness = 2
#beta_in   = 1
#no_noise_realization_in = 10

for beta_i in range(0,1):
    beta = beta_in
    print('beta = '+str(beta_in))
    for nr_i in range(0, 1):
       nr_i = int(no_noise_realization_in)
       print('noise realization num = '+str(no_noise_realization_in))
       nr_directory = 'imgs/n'+str(nr_i)
       for recon_i in range(2,3):
            print('recon num = '+str(recon_i))
           ################ load SPECT images ###################
            try:
                #print(directory+nr_directory+'/osads.wl.n'+str(nr_i)+'.'+str(recon_i)+'.mat')
                SPECTimg = sio.loadmat(directory+nr_directory+'/osads.wl.n'+str(nr_i)+'.'+str(recon_i)+'.mat')
                SPECTimg = SPECTimg['out']
                tmp_shape= np.shape(SPECTimg)
                if(tmp_shape[0] > 250):
                    continue
            except:
            #    print('error loading SPECTimg')
               continue
                
            
            # variable initialization
            dice_out_les          = []
            dice_out_bone         = []
            uptake_bias_out_les   = [] 
            uptake_bias_out_bone  = []
            uptake_bias_out_ratio = []
            uptake_out_les        = []
            uptake_out_bone       = []
            uptake_MSE_out_les    = [] 
            uptake_MSE_out_bone   = []
            uptake_MSE_out_ratio  = []

            for les_i in range(1,no_lesion+1):
                print('leision num = '+str(les_i))
                ################ load seeds & ROI ###################
                try:
                    ROIimg = sio.loadmat(directory+'masks_3D/masks_les'+str(les_i)+'.mat')
                    ROIimg = ROIimg['ROIimgNew']
                    seeds  = sio.loadmat(directory+'seeds_3D/seeds_les'+str(les_i)+'.mat')
                    seeds  = seeds['seeds']
                    seed_les  = seeds[0][0]
                    seed_les  = [seed_les[0,0]-1, seed_les[0,1]-1, seed_les[0,2]-1]
                    seed_bone = seeds[0][1]
                    seed_bone = [seed_bone[0,0]-1, seed_bone[1,0]-1, seed_bone[2,0]-1]
                    seed_cort = seeds[0][2]
                    seed_cort = [seed_cort[0,0]-1, seed_cort[1,0]-1, seed_cort[2,0]-1]
                    seed_trab = seeds[0][3]
                    seed_trab = [seed_trab[0,0]-1, seed_trab[1,0]-1, seed_trab[2,0]-1]
                    #print(seed_les)
                    #print(seed_bone)
                    #print(seed_cort)
                    #print(seed_trab)
                    #print('////')
                    
                except:
                    #print('error')
                    continue
                    
                ################ obtain ROI images ###################
                xBegin, xEnd = get_duration(ROIimg, 'x')
                yBegin, yEnd = get_duration(ROIimg, 'y')
                zBegin, zEnd = get_duration(ROIimg, 'z')
                #print(xBegin)
                #print(xEnd)
                CT_ROI         = CTimg[xBegin:xEnd, yBegin:yEnd, zBegin:zEnd]
                CT_ROI         = np.divide(CT_ROI, np.max(CT_ROI))
                CT_ROI         = np.ascontiguousarray(CT_ROI, np.float32)
                SPECT_ROI      = SPECTimg[xBegin:xEnd, yBegin:yEnd, zBegin:zEnd]
                SPECT_ROI      = np.divide(SPECT_ROI,np.max(SPECT_ROI))
                SPECT_ROI      = np.ascontiguousarray(SPECT_ROI, np.float32) 
                bone_truth_ROI = bone_truth[xBegin:xEnd, yBegin:yEnd, zBegin:zEnd]
                #plt.figure(1); plt.imshow(bone_truth_ROI[:,:,seed_les[2]]); plt.show()
                les_truth_ROI  = les_truth[xBegin:xEnd, yBegin:yEnd, zBegin:zEnd]
                #plt.figure(1); plt.imshow(les_truth_ROI[:,:,seed_les[2]]); plt.show()

                #plt.figure(1); plt.imshow(CT_ROI[:,:,seed_les[2]]); plt.show()

                ################ cluster CT images ###################
                num_clus_ct = 4
                # clustering method goes here:
                zdim,ydim,xdim = CT_ROI.shape

                #plt.figure(1); plt.imshow(CT_ROI[:,:,seed_les[2]]); plt.show()
                CT_labeled = np.zeros(CT_ROI.shape)
                CT_labeled = np.ascontiguousarray(CT_labeled, np.int16)
                c_func_FCM(CT_ROI, CT_labeled, xdim, ydim, zdim, num_clus_ct, 0.0064, fuzziness)
                #plt.figure(1); plt.imshow(CT_labeled[:,:,seed_les[2]]); plt.show()
                ################ region growing on CT images ###################
                # cort bone
                CT_cort_labeled = extract_reg_w_same_label(CT_labeled, seed_cort[0], seed_cort[1], seed_cort[2])
                #plt.figure(1); plt.imshow(CT_cort_labeled[:,:,seed_les[2]]); plt.show()
                # trab bone
                CT_trab_labeled = extract_reg_w_same_label(CT_labeled, seed_trab[0], seed_trab[1], seed_trab[2])
                #plt.figure(1); plt.imshow(CT_trab_labeled[:,:,seed_les[2]]); plt.show()
                # merge regions
                CT_labeled = np.zeros(CT_labeled.shape)
                CT_labeled = CT_labeled + (CT_cort_labeled + CT_trab_labeled)
                CT_labeled[CT_labeled >= 1] = 1
                #plt.figure(1); plt.imshow(CT_labeled[:,:,seed_les[2]]); plt.show()

                ################ cluster SPECT images ###################
                num_clus_spect = 3
                # clustering method goes here:
                SPECT_labeled = np.zeros(SPECT_ROI.shape)
                SPECT_labeled = np.ascontiguousarray(SPECT_labeled, np.int16)
                CT_labeled    = np.ascontiguousarray(CT_labeled, np.int16)
                c_func_FCM_ctinfo(SPECT_ROI, SPECT_labeled, CT_labeled, xdim, ydim, zdim, num_clus_spect, beta, fuzziness, gamma)
                #plt.figure(2); plt.imshow(SPECT_labeled[:,:,seed_les[2]]); plt.show()
                
                ################ region growing on SPECT images ###################
                # lesion region growing
                SPECT_les_labeled = np.ascontiguousarray(np.asarray(list(SPECT_labeled)), np.int16)
                c_func_region_grow_cropped3Dimg(SPECT_les_labeled, xdim, ydim, zdim, seed_les[2], seed_les[1], seed_les[0], num_clus_spect)
                #plt.figure(1); plt.imshow(SPECT_les_labeled[:,:,seed_les[2]]); plt.title('seg results'); plt.show()
                #plt.figure(2); plt.imshow(SPECT_les_labeled[:,:,seed_les[2]]); plt.show()

                # bone region growing
                #SPECT_bone_labeled = np.ascontiguousarray(np.asarray(list(SPECT_labeled)), np.int16)
                #c_func_region_grow_cropped3Dimg(SPECT_bone_labeled, xdim, ydim, zdim, seed_bone[2], seed_bone[1], seed_bone[0], num_clus_spect)
                SPECT_bone_labeled = extract_reg_w_same_label(SPECT_labeled, seed_bone[0], seed_bone[1], seed_bone[2])

                ################ region growing on truth images ###################
                les_truth_ROI = np.ascontiguousarray(les_truth_ROI, np.int16)
                c_func_region_grow_cropped3Dimg(les_truth_ROI, xdim, ydim, zdim, seed_les[2], seed_les[1], seed_les[0], 2)
                #plt.figure(2); plt.imshow(les_truth_ROI[:,:,seed_les[2]]); plt.title('truth');plt.show()

                ################ calc volume ###################
                vol_les  = calc_vol(les_truth_ROI)
                vol_bone = calc_vol(bone_truth_ROI)

                dice_coef_les  = dice_eval(SPECT_les_labeled,les_truth_ROI)
                dice_coef_bone = dice_eval(SPECT_bone_labeled, bone_truth_ROI)

                dice_out_les.append([les_i, dice_coef_les, vol_les, vol_bone])
                dice_out_bone.append([les_i, dice_coef_bone, vol_les, vol_bone])
                print(dice_coef_les)
                print(dice_coef_bone)

                ################ calc uptake ratio ###################
                gt_ratio = uptake_ratio(SPECT_ROI, bone_truth_ROI, les_truth_ROI)
                seg_ratio = uptake_ratio(SPECT_ROI, SPECT_bone_labeled, SPECT_les_labeled)

                ################ calc uptake bias ###################
                uptake_bias_les  = uptake_eval(SPECT_ROI, SPECT_les_labeled, les_truth_ROI, 'BIAS')
                uptake_bias_bone = uptake_eval(SPECT_ROI, SPECT_bone_labeled, bone_truth_ROI, 'BIAS')
                uptake_bias_out_les.append([les_i, uptake_bias_les, vol_les, vol_bone])
                uptake_bias_out_bone.append([les_i, uptake_bias_bone, vol_les, vol_bone])

                uptake_bias_ratio = abs(gt_ratio - seg_ratio)
                uptake_bias_out_ratio.append([les_i, uptake_bias_ratio, vol_les, vol_bone])

                ################ calc uptake variance ###################
                uptake_var_les = uptake_eval(SPECT_ROI, SPECT_les_labeled, les_truth_ROI, 'VAR')
                uptake_out_les.append([les_i, uptake_var_les, vol_les, vol_bone])

                uptake_var_bone = uptake_eval(SPECT_ROI, SPECT_bone_labeled, bone_truth_ROI, 'VAR')
                uptake_out_bone.append([les_i, uptake_var_bone, vol_les, vol_bone])

                ################ calc uptake MSE ###################
                uptake_MSE_les = uptake_eval(SPECT_ROI, SPECT_les_labeled, les_truth_ROI, 'MSE')
                uptake_MSE_out_les.append([les_i, uptake_MSE_les, vol_les, vol_bone])

                uptake_MSE_bone = uptake_eval(SPECT_ROI, SPECT_bone_labeled, bone_truth_ROI, 'MSE')
                uptake_MSE_out_bone.append([les_i, uptake_MSE_bone, vol_les, vol_bone])

                uptake_MSE_ratio = np.power(abs(gt_ratio - seg_ratio),2)
                uptake_MSE_out_ratio.append([les_i, uptake_MSE_ratio, vol_les, vol_bone])
            

            ################ saving results ###################
            dice_out_filename = 'dice_beta'+str(beta)+'_gamma'+str(gamma)+'_n'+str(nr_i)+'_iter'+str(recon_i)+'.mat'
            with open((directory+'seg_results/dice_3D/'+dice_out_filename),'wb') as f:
                sio.savemat(f, {'dice_out_les' :dice_out_les},do_compression=True)
                sio.savemat(f, {'dice_out_bone':dice_out_bone},do_compression=True)
            
            bias_out_filename = 'bias_beta'+str(beta)+'_gamma'+str(gamma)+'_n'+str(nr_i)+'_iter'+str(recon_i)+'.mat'
            with open((directory+'seg_results/bias_3D/'+bias_out_filename),'wb') as f:
                sio.savemat(f, {'uptake_bias_out_les' :uptake_bias_out_les},do_compression=True)
                sio.savemat(f, {'uptake_bias_out_bone':uptake_bias_out_bone},do_compression=True)
                sio.savemat(f, {'uptake_bias_out_ratio':uptake_bias_out_ratio},do_compression=True)
            
            var_out_filename = 'var_beta'+str(beta)+'_gamma'+str(gamma)+'_n'+str(nr_i)+'_iter'+str(recon_i)+'.mat'
            with open((directory+'seg_results/variance_3D/'+var_out_filename),'wb') as f:
                sio.savemat(f, {'uptake_out_les' :uptake_out_les},do_compression=True)
                sio.savemat(f, {'uptake_out_bone':uptake_out_bone},do_compression=True)

            MSE_out_filename = 'MSE_beta'+str(beta)+'_gamma'+str(gamma)+'_n'+str(nr_i)+'_iter'+str(recon_i)+'.mat'
            with open((directory+'seg_results/MSE_3D/'+MSE_out_filename),'wb') as f:
                sio.savemat(f, {'uptake_MSE_out_les' :uptake_MSE_out_les},do_compression=True)
                sio.savemat(f, {'uptake_MSE_out_bone':uptake_MSE_out_bone},do_compression=True)
                sio.savemat(f, {'uptake_MSE_out_ratio':uptake_MSE_out_ratio},do_compression=True)

print('evaluation complete!')






                




