#new imports
import analysisNOGUI
import metadata_trimmer
import parameters_class
from time import sleep
from random import random
import threading
import sys
#all else
import numpy as np
import cv2
import math
import os
import time
from skimage.measure import regionprops, regionprops_table
from skimage.io import imread
from scipy import ndimage
from PIL import Image
import pandas as pd
import matplotlib.image as mpimg
from PIL import Image, ImageQt
from scipy.ndimage import label, distance_transform_edt
import multiprocessing
from joblib import Parallel, delayed
from multiprocessing import Pool, Process, Manager, Value, Lock
#from GUI_parameters import Gui_Params
import btrack
import imageio
from tifffile import imwrite
from skimage.transform import rotate
from matplotlib import pyplot as plt
from skimage.color import label2rgb, gray2rgb
from aicsimageio import AICSImage
from aicsimageio.writers import OmeTiffWriter
from deepcell.applications import CellTracking
from skimage.feature import peak_local_max
from skimage.metrics import structural_similarity as ssim
from skimage.metrics import adapted_rand_error, contingency_table, mean_squared_error, peak_signal_noise_ratio, variation_of_information
from skimage.filters import threshold_otsu, median
from sklearn import mixture
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
from skimage.morphology import disk, binary_closing, skeletonize, binary_opening, binary_erosion, white_tophat

WELL_PLATE_ROWS = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"]

#class for making a process-safe counter
class SafeCounter():
    # constructor
    def __init__(self):
        # initialize counter
        self._counter = Value('i', 0)
 
    # increment the counter
    def increment(self):
        # get the lock
        with self._counter.get_lock():
            self._counter.value += 1
 
    # get the counter value
    def value(self):
        # get the lock
        with self._counter.get_lock():
            return self._counter.value
    
class BatchAnalysis(object):
    output_folder = []
    manager = Manager()
    cell_pd_list, ch1_spot_df_list, ch2_spot_df_list = manager.list(), manager.list(), manager.list()
    ch3_spot_df_list, ch4_spot_df_list, ch5_spot_df_list = manager.list(), manager.list(), manager.list()
    def __init__(self, is_gui, input_params = None, analysisgui = None, image_analyzer = None, inout_resource_gui = None, displaygui = None, ImDisplay = None):
        self.is_gui = is_gui
        if is_gui == False:
            #intiialize nogui version
            self.ch1_spot_df, self.ch2_spot_df, self.ch3_spot_df = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
            self.ch4_spot_df, self.ch5_spot_df, self.cell_df = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
            self.df_checker = pd.DataFrame()
            self.Meta_Data_df = pd.DataFrame()
            self.spot_distances = {}
            self.experiment_name = []
            self.output_prefix = []
            self.input_params = input_params
            #initialize the gui_params class
            self.gui_params = parameters_class.Gui_Params(is_gui, input_params=input_params)
            #initialize with my analysisNOGUI class
            self.ImageAnalyzer = image_analyzer
            #initialize counters
            self.counter = SafeCounter() #masks
            self.distance_counter = SafeCounter() #spot distances
            #stuff for the progress saver
            self.trimmer = metadata_trimmer.metadata_trimmer(input_params["metadata_path"], input_params["metadata_name"])
            self.lock = Lock()
            
        else:
            #initialize gui version
            self.inout_resource_gui = inout_resource_gui
            self.AnalysisGui = analysisgui
            self.displaygui = displaygui
            self.ImageAnalyzer = image_analyzer
            self.ImDisplay = ImDisplay
            self.ch1_spot_df, self.ch2_spot_df, self.ch3_spot_df = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
            self.ch4_spot_df, self.ch5_spot_df, self.cell_df = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
            self.df_checker = pd.DataFrame()
            self.Meta_Data_df = pd.DataFrame()
            self.spot_distances = {}
            self.experiment_name = []
            self.output_prefix = []
            #initialize the gui_params class
            self.gui_params = parameters_class.Gui_Params(is_gui, analysisgui=self.AnalysisGui, inout_resource_gui=self.inout_resource_gui)

    def ON_APPLYBUTTON(self, Meta_Data_df):
        
        
        seconds1 = time.time()
        #set output directory
        #if it's the headless version, we need to pull the output dir from the gui_params class
        if self.is_gui == False:
            self.output_dir = self.gui_params.output_dir
            if os.path.isdir(self.output_dir) == False:
                os.mkdir(self.output_dir) 
            self.Meta_Data_df = Meta_Data_df
            path_list = os.path.split(self.Meta_Data_df["ImageName"][0])[0].split(r'/')
            self.experiment_name = path_list[path_list.__len__()-1]
            self.output_prefix  = path_list[path_list.__len__()-1]
            self.output_folder = os.path.join(self.output_dir,self.experiment_name)
        #if it's the gui version, just call the same thing as the original
        else:
            while self.inout_resource_gui.Output_dir ==[]:
                self.inout_resource_gui.OUTPUT_FOLDER_LOADBTN()
            self.Meta_Data_df = Meta_Data_df
            path_list = os.path.split(self.Meta_Data_df["ImageName"][0])[0].split(r'/')
            self.experiment_name = path_list[path_list.__len__()-2]
            self.output_prefix  = path_list[path_list.__len__()-1]
            self.output_folder = os.path.join(self.inout_resource_gui.Output_dir,self.experiment_name)
        
        #after we get the path, gui or not, we make sure it exists
        if os.path.isdir(self.output_folder) == False:
            os.mkdir(self.output_folder) 
            
            #saves the configuration to a csv in the folder that it creates
        csv_config_folder = os.path.join(self.output_folder, 'configuration_files')
        if os.path.isdir(csv_config_folder) == False:
            os.mkdir(csv_config_folder) 
        self.config_file = os.path.join(csv_config_folder, 'analysis_configuration.csv')

        columns = np.unique(np.asarray(self.Meta_Data_df['column'], dtype=int))
        rows = np.unique(np.asarray(self.Meta_Data_df['row'], dtype=int))
        fovs = np.unique(np.asarray(self.Meta_Data_df['field_index'], dtype=int))
#         fovs = np.array([2])
        time_points = np.unique(np.asarray(self.Meta_Data_df['time_point'], dtype=int))
        #used for mask progress bars later 
        #CHANGED
        self.masks = len(columns)*len(rows)*len(fovs)*len(time_points)
        actionindices = np.unique(np.asarray(self.Meta_Data_df['action_index'], dtype=int))
        
        jobs_number=self.gui_params.NumCPUsSpinBox_value
        

#         Parallel(n_jobs=jobs_number)(delayed(self.BATCH_ANALYZER)(col,row,fov,t) for t in time_points for fov in fovs for row in rows for col in columns)
        func_args=np.zeros((0,4),dtype=int)
        for t in time_points:
            for fov in fovs:
                for row in rows:
                    for col in columns:
                        df_parallel = self.Meta_Data_df.loc[(self.Meta_Data_df['column'] == str(col)) & 
                                                              (self.Meta_Data_df['row'] == str(row)) & 
                                                              (self.Meta_Data_df['field_index'] == str(fov)) & 
                                                              (self.Meta_Data_df['time_point'] == str(t))]
                        if df_parallel.empty == False:
                            
                            func_args=np.append(func_args,np.array([col,row,fov,t]).reshape(1,4),axis=0)

        arg_len,_=func_args.shape
        start_pos=np.arange(0,arg_len,jobs_number)
        for st in start_pos:
            
            if st+jobs_number-1 >= arg_len:
                
                data_ind=np.arange(st, arg_len)
            else:
                data_ind=np.arange(st, st+jobs_number)
            
            processes=[]
            for ind1 in data_ind:
                process_args= np.array(func_args[ind1,:],dtype=int)
                if self.is_gui == False:
                    #if it's the headless version and I want to have the resume from interrupt work
                    #need to add a lock and a sleep timer to the args
                    process_args = process_args.tolist()
                    process_args.append(self.lock)
                    process_args.append(random())
                    #CHANGED
                    processes.append(Process(target=self.BATCH_ANALYZER, args=process_args))
                else:
                    processes.append(Process(target=self.BATCH_ANALYZER, args=process_args))

            # kick them off 
            for process in processes:
                process.start()
            # now wait for them to finish
            for process in processes:
                process.join()
        #if it all ran well, delete the temp file
        self.trimmer.del_temp()    
        xlsx_output_folder = os.path.join(self.output_folder, 'whole_plate_results')
        if os.path.isdir(xlsx_output_folder) == False:
            os.mkdir(xlsx_output_folder) 
        
        #get spot locations for each channel                                                                                       
        if self.gui_params.NucInfoChkBox_check_status == True: #nucinfochkbox
            
            self.cell_df = pd.concat(self.cell_pd_list)
            
            xlsx_name = ['Nuclei_Information.csv']
            xlsx_full_name = os.path.join(xlsx_output_folder, xlsx_name[0])
            self.cell_df.rename(columns={ "label":"cell_index"}, inplace = True)
            self.cell_df.to_csv(xlsx_full_name)
            
            well_nuc_folder = os.path.join(self.output_folder, 'well_nuclei_results')
            if os.path.isdir(well_nuc_folder) == False:
                os.mkdir(well_nuc_folder)
            for col in columns:
                for row in rows:

                    well_nuc_df = self.cell_df.loc[(self.cell_df['column'] == col) & (self.cell_df['row'] == row)]
                    if well_nuc_df.empty == False:
                        well_nuc_filename = self.output_prefix + '_nuclei_information_well_' + WELL_PLATE_ROWS[row-1] + str(col) + r'.csv'
                        nuc_well_csv_full_name = os.path.join(well_nuc_folder, well_nuc_filename)
                        well_nuc_df.to_csv(path_or_buf=nuc_well_csv_full_name, encoding='utf8')
        
        
        if self.ch1_spot_df_list.__len__() > 0:
            self.ch1_spot_df = pd.concat(self.ch1_spot_df_list)
            #not sure if this should be NucInfoChkBox or SpotsLocation - nt
            if self.gui_params.NucInfoChkBox_check_status == True:
                coordinates_method = self.gui_params.SpotLocationCbox_currentText
                xlsx_name = ['Ch1_Spot_Locations_' + coordinates_method + r'.csv']
                xlsx_full_name = os.path.join(xlsx_output_folder, xlsx_name[0])
                self.ch1_spot_df.to_csv(xlsx_full_name)
                
                well_spot_loc_folder = os.path.join(self.output_folder, 'well_spots_locations')
                if os.path.isdir(well_spot_loc_folder) == False:
                    os.mkdir(well_spot_loc_folder)
                for col in columns:
                    for row in rows:
                        spot_loc_df = self.ch1_spot_df.loc[(self.ch1_spot_df['column'] == col) & (self.ch1_spot_df['row'] == row)]
                        if spot_loc_df.empty == False:
                            spot_loc_filename = self.output_prefix + '_ch1_spots_locations_well_' + WELL_PLATE_ROWS[row-1] + str(col) + r'.csv'
                            spot_loc_well_csv_full_name = os.path.join(well_spot_loc_folder, spot_loc_filename)
                            spot_loc_df.to_csv(path_or_buf=spot_loc_well_csv_full_name, encoding='utf8')
                
        if self.ch2_spot_df_list.__len__() > 0:
            self.ch2_spot_df = pd.concat(self.ch2_spot_df_list)
            
            if self.gui_params.SpotsLocation_check_status == True:
                coordinates_method = self.gui_params.SpotLocationCbox_currentText
                xlsx_name = ['Ch2_Spot_Locations_' + coordinates_method + r'.csv']
                xlsx_full_name = os.path.join(xlsx_output_folder, xlsx_name[0])
                self.ch2_spot_df.to_csv(xlsx_full_name)   
                
                well_spot_loc_folder = os.path.join(self.output_folder, 'well_spots_locations')
                if os.path.isdir(well_spot_loc_folder) == False:
                    os.mkdir(well_spot_loc_folder)
                for col in columns:
                    for row in rows:
                        spot_loc_df = self.ch2_spot_df.loc[(self.ch2_spot_df['column'] == col) & (self.ch2_spot_df['row'] == row)]
                        if spot_loc_df.empty == False:
                            spot_loc_filename = self.output_prefix + '_ch2_spots_locations_well_' + WELL_PLATE_ROWS[row-1] + str(col) + r'.csv'
                            spot_loc_well_csv_full_name = os.path.join(well_spot_loc_folder, spot_loc_filename)
                            spot_loc_df.to_csv(path_or_buf=spot_loc_well_csv_full_name, encoding='utf8')
        
        if self.ch3_spot_df_list.__len__() > 0:
            self.ch3_spot_df = pd.concat(self.ch3_spot_df_list)
            
            if self.gui_params.SpotsLocation_check_status == True:
                coordinates_method = self.gui_params.SpotLocationCbox_currentText
                xlsx_name = ['Ch3_Spot_Locations_' + coordinates_method + r'.csv']
                xlsx_full_name = os.path.join(xlsx_output_folder, xlsx_name[0])
                self.ch3_spot_df.to_csv(xlsx_full_name)   
                
                well_spot_loc_folder = os.path.join(self.output_folder, 'well_spots_locations')
                if os.path.isdir(well_spot_loc_folder) == False:
                    os.mkdir(well_spot_loc_folder)
                for col in columns:
                    for row in rows:
                        spot_loc_df = self.ch3_spot_df.loc[(self.ch3_spot_df['column'] == col) & (self.ch3_spot_df['row'] == row)]
                        if spot_loc_df.empty == False:
                            spot_loc_filename = self.output_prefix + '_ch3_spots_locations_well_' + WELL_PLATE_ROWS[row-1] + str(col) + r'.csv'
                            spot_loc_well_csv_full_name = os.path.join(well_spot_loc_folder, spot_loc_filename)
                            spot_loc_df.to_csv(path_or_buf=spot_loc_well_csv_full_name, encoding='utf8')
        
        if self.ch4_spot_df_list.__len__() > 0:
            self.ch4_spot_df = pd.concat(self.ch4_spot_df_list)
            
            if self.gui_params.SpotsLocation_check_status == True:
                coordinates_method = self.gui_params.SpotLocationCbox_currentText
                xlsx_name = ['Ch4_Spot_Locations_' + coordinates_method + r'.csv']
                xlsx_full_name = os.path.join(xlsx_output_folder, xlsx_name[0])
                self.ch4_spot_df.to_csv(xlsx_full_name) 
                
                well_spot_loc_folder = os.path.join(self.output_folder, 'well_spots_locations')
                if os.path.isdir(well_spot_loc_folder) == False:
                    os.mkdir(well_spot_loc_folder)
                for col in columns:
                    for row in rows:
                        spot_loc_df = self.ch4_spot_df.loc[(self.ch4_spot_df['column'] == col) & (self.ch4_spot_df['row'] == row)]
                        if spot_loc_df.empty == False:
                            spot_loc_filename = self.output_prefix + '_ch4_spots_locations_well_' + WELL_PLATE_ROWS[row-1] + str(col) + r'.csv'
                            spot_loc_well_csv_full_name = os.path.join(well_spot_loc_folder, spot_loc_filename)
                            spot_loc_df.to_csv(path_or_buf=spot_loc_well_csv_full_name, encoding='utf8')
        
        if self.ch5_spot_df_list.__len__() > 0:
            self.ch5_spot_df = pd.concat(self.ch5_spot_df_list)
            
            if self.gui_params.SpotsLocation_check_status == True:
                coordinates_method = self.gui_params.SpotLocationCbox_currentText
                xlsx_name = ['Ch5_Spot_Locations_' + coordinates_method + r'.csv']
                xlsx_full_name = os.path.join(xlsx_output_folder, xlsx_name[0])
                self.ch5_spot_df.to_csv(xlsx_full_name)  
                
                well_spot_loc_folder = os.path.join(self.output_folder, 'well_spots_locations')
                if os.path.isdir(well_spot_loc_folder) == False:
                    os.mkdir(well_spot_loc_folder)
                for col in columns:
                    for row in rows:
                        spot_loc_df = self.ch5_spot_df.loc[(self.ch5_spot_df['column'] == col) & (self.ch5_spot_df['row'] == row)]
                        if spot_loc_df.empty == False:
                            spot_loc_filename = self.output_prefix + '_ch5_spots_locations_well_' + WELL_PLATE_ROWS[row-1] + str(col) + r'.csv'
                            spot_loc_well_csv_full_name = os.path.join(well_spot_loc_folder, spot_loc_filename)
                            spot_loc_df.to_csv(path_or_buf=spot_loc_well_csv_full_name, encoding='utf8')
        
        
        #a resume from interrupt here would have to be threadsafe, not process-safe and so would be different than the nuc mask one
        #this section just checks what kinds of results we want to have
        if self.gui_params.SpotsDistance_check_status == True:
            print("\nStarting Spot Distance calculation")
            columns = np.unique(np.asarray(self.cell_df['column'], dtype=int))
            rows = np.unique(np.asarray(self.cell_df['row'], dtype=int))
            self.distances = len(rows)*len(columns)
            Parallel(n_jobs=jobs_number, prefer="threads")(delayed(self.Calculate_Spot_Distances)(row, col) for row in rows for col in columns)
            print("\nCompleted Spot Distance Tracking")
        if self.gui_params.Cell_Tracking_check_status == True:
            
            cell_tracking_folder = os.path.join(self.output_folder, 'cell_tracking')
            if os.path.isdir(cell_tracking_folder) == False:
                os.mkdir(cell_tracking_folder)

            #### get all the unique channels containing spots and load all the images related to the spots
            columns = np.unique(np.asarray(self.cell_df['column'], dtype=int))
            rows = np.unique(np.asarray(self.cell_df['row'], dtype=int))
            cr_arr=[]
            for row in rows:
                for col in columns:
                    subcell_df=self.cell_df.loc[(self.cell_df['column']==col)&(self.cell_df['row']==row)]
                    if subcell_df.empty==False:
                        cr_arr.append([row,col])
                        
            fovs = np.unique(np.asarray(self.Meta_Data_df['field_index'], dtype=int))
#             fovs = np.array([2])
            time_points = np.unique(np.asarray(self.Meta_Data_df['time_point'], dtype=int))
            actionindices = np.unique(np.asarray(self.Meta_Data_df['action_index'], dtype=int))
            
            for cr_elem in cr_arr:
                    row=int(cr_elem[0])
                    col=int(cr_elem[1])
                    #### load and concatenate all the spot information csv files for all the wells and fields. 
                    #### For every field the csv file for spots in each channel will be loaded and concatenated 
                    spots_loc = pd.DataFrame()
                    spot_loc_substring = '_well_' + WELL_PLATE_ROWS[row-1] + str(col) + r'.csv'
                    well_spots_output_folder = os.path.join(self.output_folder, 'well_spots_locations')
                    if os.path.isdir(well_spot_loc_folder) == False:
                        os.mkdir(well_spot_loc_folder)
                    for fname in os.listdir(well_spots_output_folder):    # change directory as needed
                        if spot_loc_substring in fname:    # search for string
                            spots_loc = pd.concat([spots_loc,pd.read_csv(os.path.join(well_spots_output_folder,fname))])
                    ### get all the channels containing the spots
                    if spots_loc.empty:
                        spot_channels=np.array([])
                    else:
                        spot_channels = spots_loc['channel'].unique()
                        
                    for fov in fovs:
                        print(fovs)
                        spot_images = {}
                        spot_dict_stack={}
                        for sp_ch in spot_channels:
                            spot_images[str(int(sp_ch))] = []
                            
                        masks=[]
                        lbl_imgs=[]
                        nuc_imgs = []
                        for t in time_points:
                            self.df_checker = self.Meta_Data_df.loc[(self.Meta_Data_df['column'] == str(col)) & 
                                              (self.Meta_Data_df['row'] == str(row)) & 
                                              (self.Meta_Data_df['field_index'] == str(fov)) & 
                                              (self.Meta_Data_df['time_point'] == str(t))]
                            if self.df_checker.empty:
                                continue
                                #raw image loader takes a number as the channel value
                            ImageForNucMask = self.RAW_IMAGE_LOADER(self.gui_params.NucleiChannel_index)
                            normalized_nuc_img = cv2.normalize(ImageForNucMask, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
                            nuc_imgs.append(normalized_nuc_img)
                
                            for sp_ch in spot_channels:
                                spot_images[str(int(sp_ch))].append(self.RAW_IMAGE_LOADER(str(int(sp_ch))))
                            
                            mask_file_name = ['Nuclei_Mask_for_Col' + str(col) + r'_row' + str(row)+
                                          r'_Time' + str(t) + r'_Field' + str(fov) + r'.tif']
                            
                            nuc_mask_output_folder = os.path.join(self.output_folder, 'nuclei_masks')
                            mask_full_name = os.path.join(nuc_mask_output_folder, mask_file_name[0])
                            #print(mask_full_name)
                            mask_image = imread(mask_full_name,0)
                            #print(mask_image)
                            lbl, ncc = label(mask_image)
                            #print("lbl")
                            #print(lbl)
                            masks.append(mask_image)
                            lbl_imgs.append(lbl)
                        
                        label_stack=np.stack(lbl_imgs,axis=0)
                        
                        masks_stack=np.stack(masks,axis=0)
                        t_stack_nuc=np.stack(nuc_imgs,axis=0)
                        for sp_ch in spot_channels:
                            spot_dict_stack[str(int(sp_ch))] = np.stack(spot_images[str(int(sp_ch))],axis=2)
                        
                        #### select spots related to the current field only 
                        selected_spots = spots_loc.loc[(spots_loc['field_index'] == fov)]
                        #### get all the unique channels contating spots 
                        
                        
                        if selected_spots.empty:
                            unique_ch=np.array([])
                        else:
                            unique_ch = selected_spots['channel'].unique()
                            
                        #### run nuclei tracking algorithm based the selected method
                        if self.gui_params.NucTrackingMethod_currentText == "Bayesian":
                            print("\nStarting Bayesian cell tracking")
                            tracks_pd = self.RUN_BTRACK(label_stack, self.gui_params)
                        
                        if self.gui_params.NucTrackingMethod_currentText == "DeepCell":
                            print("\nStarting Deepcell cell tracking")
                            tracks_pd = self.deepcell_tracking(t_stack_nuc,label_stack,self.gui_params)
                        
                        ##########save whole field track images
                        font = cv2.FONT_HERSHEY_SIMPLEX
                        fontScale = int(1)
                        fontColor = (255,255,255)
                        lineType  = int(2)
                        t,ww,hh=t_stack_nuc.shape
                        rgb_t_stack = np.zeros((t,1,3,ww,hh), dtype=t_stack_nuc.dtype)
                        for i in tracks_pd['t'].unique():

                            All_Channels = gray2rgb(t_stack_nuc[i,:,:])
                            time_tracks = tracks_pd.loc[tracks_pd['t']==i]
                            for ind, pd_row in time_tracks.iterrows():

                                txt=str(pd_row["ID"])
                                bottomLeftCornerOfText = (int(round(pd_row["x"])), int(round(pd_row["y"])))
                                cv2.putText(All_Channels,txt, bottomLeftCornerOfText, font, fontScale, fontColor, lineType)

                            rgb_t_stack[i,:,:,:,:] = np.rollaxis(All_Channels, 2, 0)  
                        
                        Raw_file_name = 'track_image_col' + str(col) + r'_row' + str(row)+ r'_field' + str(fov) + r'.ome.tif'
                        whole_field_track_folder = os.path.join(cell_tracking_folder, 'whole_field_track_images')
                        if os.path.isdir(whole_field_track_folder) == False:
                            os.mkdir(whole_field_track_folder)
                        whole_field_full_name = os.path.join(whole_field_track_folder, Raw_file_name)
                        writer = OmeTiffWriter()
                        writer.save(rgb_t_stack[:,0,0,:,:],whole_field_full_name, dimension_order="TYX")

#                         writer.save(spot_dict_stack[str(int(sp_ch))],whole_field_full_name, dimension_order="TYX")
                        print("field number: " , fov)
                        #####################
                        
                        IDs = tracks_pd['ID'].unique()

                        for id_ in IDs:
                            #### create a copy of each single track df containing all the nuclei informaion
                            single_track1 = tracks_pd.loc[tracks_pd['ID']==id_]
                            single_track = single_track1[single_track1['dummy']==False]
                            single_track_copy = single_track.copy()
                            
                            if len(single_track_copy)< self.gui_params.MintrackLengthSpinbox_current_value:
                                 continue
                            
                            for chnl in unique_ch:
                                col_name1='ch'+str(int(chnl))+'_spots_number'
                                single_track_copy[col_name1]=np.zeros(len(single_track),dtype=int)
                                col_name2='ch'+str(int(chnl))+'_spots_locations'
                                single_track_copy[col_name2]=[[]]*len(single_track)
                            
                            nuc_patches = []
                            rotated_nuc = []
                            spot_coor_dict = {}
                            spots_mask = {}
                            spot_patches = {}
                            rotated_spot_patches = {}
                            rotated_spot_patches_mean = {}
                            index_to_time = {}
                            spot_channels = selected_spots['channel'].unique().astype(int)
                            for sp_ch in spot_channels:
                                spot_patches[str(int(sp_ch))] = []
                                rotated_spot_patches[str(int(sp_ch))] = []
                                rotated_spot_patches_mean[str(int(sp_ch))] = []
                                spots_mask[str(int(sp_ch))] = []
                            time_ind=0
                            no_rot_metrics = {"cosine_similarity":[], "mutual_info":[], "ssim":[], "mse":[], "inv_adapt_error":[], "psnr":[], "vi1":[], "vi2":[]}
                            rot_metrics = {"cosine_similarity":[], "mutual_info":[], "ssim":[], "mse":[],"inv_adapt_error":[], "psnr":[], "vi1":[], "vi2":[]}
                            final_metrics = {"cosine_similarity":[], "mutual_info":[], "ssim":[], "mse":[], "inv_adapt_error":[], "psnr":[], "vi1":[], "vi2":[]}
                            for ind, pd_row in single_track.iterrows():
                                time_spots = selected_spots[selected_spots['time_point']==pd_row['t']]
                                patch_spots = time_spots[(time_spots['x_location'] > int(pd_row['bbox-0'])) & 
                                                         (time_spots['x_location'] < int(pd_row['bbox-2'])) & 
                                                         (time_spots['y_location'] > int(pd_row['bbox-1'])) & 
                                                         (time_spots['y_location'] < int(pd_row['bbox-3'])) ]
                                
                                min_row = max(int(pd_row['bbox-0'])-5, 0)
                                max_row = min(int(pd_row['bbox-2'])+5, ww)
                                min_col = max(int(pd_row['bbox-1'])-5, 0)
                                max_col = min(int(pd_row['bbox-3'])+5, hh)
                                
                                img_patch1 = t_stack_nuc[int(pd_row['t']-1), min_row:max_row, min_col:max_col]
        
                                img_patch = self.zero_pad_patch(img_patch1)

                                ##### this part finds the best rotation orientation for aligned nuclei to prevent flipping 
                                angle_1 = (-180*pd_row['orientation']/np.pi)+90

                                im1_rot = rotate(img_patch, angle_1, center=np.array([128, 128]),  preserve_range=True)


                                nuc_patches.append(rotate(img_patch, 0, center=np.array([int(pd_row['y'] -min_row), int(pd_row['x'] - min_col)]), preserve_range=True))

                                if len(rotated_nuc)>0:
                                    
                                    img1 = median(rotated_nuc[-1]*255/rotated_nuc[-1].max(), disk(3))
            
                                    img2 = median(im1_rot*255/im1_rot.max(), disk(3))

                                    angle_2 = angle_1 + 180
                                    im2_rot = rotate(img_patch, angle_2, center=np.array([128, 128]), preserve_range=True)
                                    img2_rot = median(im2_rot*255/im2_rot.max() ,disk(3))

                                    im1_rot1 = self.zero_pad_patch(img2)
                                    im2_rot1 = self.zero_pad_patch(img2_rot)

                                    csi1 = np.dot(np.ravel(img1),np.ravel(im1_rot1))/(np.linalg.norm(np.ravel(img1))*np.linalg.norm(np.ravel(im1_rot1)))
                                    csi2 = np.dot(np.ravel(img1),np.ravel(im2_rot1))/(np.linalg.norm(np.ravel(img1))*np.linalg.norm(np.ravel(im2_rot1)))

                                    hist_2d_1, x_edges, y_edges = np.histogram2d(im1_rot1.ravel(), img1.ravel(), bins=20)
                                    mi1 = self.mutual_information(hist_2d_1)

                                    hist_2d_2, x_edges, y_edges = np.histogram2d(im2_rot1.ravel(), img1.ravel(), bins=20)
                                    mi2 = self.mutual_information(hist_2d_2)


                                    ssim_ind1 = ssim(img1, im1_rot1, data_range=img1.max() - img1.min())
                                    ssim_ind2 = ssim(img1, im2_rot1, data_range=img1.max() - img1.min())

                                    mse = -mean_squared_error(img1, im1_rot1)
                                    rot_mse = -mean_squared_error(img1, im2_rot1)

                                    thresh1 = threshold_otsu(img1)
                                    thresh2 = threshold_otsu(im1_rot1)
                                    thresh_rot = threshold_otsu(im2_rot1)

                                    ct = contingency_table(img1>thresh1, im1_rot1>thresh2)
                                    are = adapted_rand_error(image_true=img1>thresh1, image_test=im1_rot1>thresh2, table=ct)[0]
                                    vi1, vi2 = variation_of_information(image0=img1>thresh1, image1=im1_rot1>thresh2 , table=ct)

                                    rot_ct = contingency_table(img1>thresh1, im2_rot1>thresh_rot)
                                    rot_are = adapted_rand_error(image_true=img1>thresh1, image_test=im2_rot1>thresh_rot, table=rot_ct)[0]
                                    rot_vi1, rot_vi2 = variation_of_information(image0=img1>thresh1, image1=im2_rot1>thresh_rot , table=rot_ct)
                                    eps_are=0.000000001

                                    psnr = peak_signal_noise_ratio(img1.astype("uint8"),im1_rot1.astype("uint8"))
                                    rot_psnr = peak_signal_noise_ratio(img1.astype("uint8"),im2_rot1.astype("uint8"))


                                    first_metrics = np.array([csi1, mi1, ssim_ind1, mse, 1/(are+eps_are), psnr, vi1, vi2])
                                    sec_metrics = np.array([csi2, mi2, ssim_ind2, rot_mse, 1/(rot_are+eps_are), rot_psnr, rot_vi1, rot_vi2])

                                    for ind, key1 in enumerate(no_rot_metrics):

                                        no_rot_metrics[key1].append(first_metrics[ind])
                                        rot_metrics[key1].append(sec_metrics[ind])

                                    if sum(first_metrics > sec_metrics) > 6: 

                                        final_angle = angle_1
                                        rotated_nuc.append(im1_rot)
                                        for ind, key1 in enumerate(final_metrics):
                                            final_metrics[key1].append(first_metrics[ind])

                                    else:
                                        final_angle = angle_2
                                        rotated_nuc.append(im2_rot)
                                        for ind, key1 in enumerate(final_metrics):
                                            final_metrics[key1].append(sec_metrics[ind])

                                else:

                                    final_angle = angle_1
                                    rotated_nuc.append(im1_rot)
                                    
                                #### this part is for training for frame to frame nuclei registration
                                
                                ### convert rotation metrics to dataframes and normalize them
#                                 notrot_df=pd.DataFrame(no_rot_metrics)
#                                 norm_notrot = (notrot_df-notrot_df.mean())/notrot_df.std()
#                                 rot_df=pd.DataFrame(rot_metrics)
#                                 norm_rot = (rot_df-rot_df.mean())/rot_df.std()
#                                 final_df = pd.DataFrame(rot_metrics)
#                                 norm_final = (final_df-final_df.mean())/final_df.std()

#                                 cell_tracking_folder = "/data1/nadya/jupyter_test/new"
#                                 registration_metrics_tables = os.path.join(cell_tracking_folder, 'registration_metrics')
#                                 if os.path.isdir(registration_metrics_tables) == False:
#                                     os.mkdir(registration_metrics_tables)
#                                 csv_file_name = 'rotated_table_col' + str(col) + r'_row' + str(row)+ r'_field' + str(fov) + r'_Cell' + str(id_) + r'_Ch' + str(int(sp_ch))+ r'_spot' + str(key1) + r'.csv'
#                                 rot_df.to_csv(os.path.join(registration_metrics_tables, csv_file_name))
#                                 csv_file_name = 'notrotated_table_col' + str(col) + r'_row' + str(row)+ r'_field' + str(fov) + r'_Cell' + str(id_) + r'_Ch' + str(int(sp_ch))+ r'_spot' + str(key1) + r'.csv'
#                                 notrot_df.to_csv(os.path.join(registration_metrics_tables, csv_file_name))
#                                 csv_file_name = 'final_table_col' + str(col) + r'_row' + str(row)+ r'_field' + str(fov) + r'_Cell' + str(id_) + r'_Ch' + str(int(sp_ch))+ r'_spot' + str(key1) + r'.csv'
#                                 final_df.to_csv(os.path.join(registration_metrics_tables, csv_file_name))

                                    
                                only_spots_patch={}
                                coordinates={}
                                for chnl in unique_ch:

                                    only_spots_patch[str(int(chnl))] = np.zeros(img_patch1.shape)
                                    patch_ch_spots = patch_spots.loc[patch_spots['channel']==int(chnl)]
                                    col_name1='ch'+str(int(chnl))+'_spots_number'
                                    single_track_copy.loc[single_track_copy['t']==pd_row['t'], col_name1] = len(patch_ch_spots)
                                    col_name2='ch'+str(int(chnl))+'_spots_locations'
                                    coordinates[str(int(chnl))] = patch_ch_spots[['x_location','y_location','z_location']].to_numpy()
                                    single_track_copy.loc[single_track_copy['t']==pd_row['t'],
                                                          col_name2] = [patch_ch_spots[['x_location','y_location','z_location']].to_numpy()]
                                    nuc_rel_coordinates={}
                                    patch_rel_coordinates={}
                                    if list(coordinates[str(int(chnl))]):
                                        nuc_rel_coordinates[str(int(chnl))] = coordinates[str(int(chnl))] - pd_row[['x','y','z']].to_numpy()
                                        patch_rel_coordinates[str(int(chnl))] = coordinates[str(int(chnl))][:,:2]-np.array([min_row, min_col])

                                        coor_for_cir = np.floor(patch_rel_coordinates[str(int(chnl))].astype(np.double)).astype(int)
#                                         patch_spots_img = self.ImageAnalyzer.COORDINATES_TO_CIRCLE(coor_for_cir,img_patch1)
#                                         img_patch[patch_spots_img != 0] = [255,0,0]
                                        only_spots_patch[str(int(chnl))][coor_for_cir[:,0],coor_for_cir[:,1]] = 255
                                        spot_coor_dict[str(pd_row['t'])]=coor_for_cir        

                                    spot_patch = self.zero_pad_patch(spot_dict_stack[str(int(chnl))][ min_row:max_row, min_col:max_col, int(pd_row['t']-1)])
            
                                    spot_patch = cv2.normalize(spot_patch, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
                                    
                                    rot_spot_patch = rotate(spot_patch, final_angle, center=np.array([128 , 128]), preserve_range=True)
                                    
                                    rotated_spot_patches[str(int(chnl))].append(rot_spot_patch)
                                    rotated_spot_patches_mean[str(int(chnl))].append(rot_spot_patch.mean())

                                    rot_spot_mask = rotate(self.zero_pad_patch(only_spots_patch[str(int(chnl))]), final_angle, center=np.array([128, 128]), preserve_range=True)
                                    
                                    spots_mask[str(int(chnl))].append(rot_spot_mask) 

                                
                                
                                index_to_time[str(time_ind)] = pd_row['t']
                                time_ind = time_ind+1

                            large_patch_size=256
                            large_rotcell_stack = np.stack(rotated_nuc, axis=2)
                            large_cell_stack = np.stack(nuc_patches, axis=2)
                            large_spot_mask_stack={}
                            large_rot_spot_img_stack={}
                            for sp_ch in spot_channels:
                                large_spot_mask_stack[str(int(sp_ch))] = np.stack(spots_mask[str(int(sp_ch))], axis=2)
                                large_rot_spot_img_stack[str(int(sp_ch))] = np.stack(rotated_spot_patches[str(int(sp_ch))], axis=2)
                                
                                
                            rot_spot_coor_dict = {}
                            
                            for sp_ch in spot_channels:
                                
                                anotated_large_rot_spot_img_stack = large_rot_spot_img_stack[str(int(sp_ch))].copy()
                                unanotated_large_rot_spot_img_stack = large_rot_spot_img_stack[str(int(sp_ch))].copy()
                                w,h,c = large_spot_mask_stack[str(int(sp_ch))].shape
                                all_spots_for_gmm=np.zeros((0,2))
                                for i in range(c):
                                    plm = peak_local_max(large_spot_mask_stack[str(int(sp_ch))][:,:,i], min_distance=1)

                                    if plm.any():
                                        temp_sp_patch = np.zeros((large_patch_size ,large_patch_size ))
                                        temp_sp_patch[plm[:,0], plm[:,1]] = 255
                                        large_spot_mask_stack[str(int(sp_ch))][:,:,i]=temp_sp_patch
                                        rot_spot_coor_dict[index_to_time[str(i)]] = plm
                                        all_spots_for_gmm = np.concatenate((all_spots_for_gmm,plm), axis=0)
                                
                                # project all the spots on one plane and detect them again
                                spot_mask_max = large_spot_mask_stack[str(int(sp_ch))].max(axis=2)
                                plm = peak_local_max(spot_mask_max, min_distance=1)
                                ## define the radius to consider the spots the same 
                                spot_circles=self.ImageAnalyzer.COORDINATES_TO_CIRCLE(plm,spot_mask_max, 
                                                                        circ_radius = self.gui_params.SpotSearchRadiusSpinbox_current_value)
                                filled1 = ndimage.binary_fill_holes(spot_circles)
                                labeled_sp, number_sp = label(filled1)
                                spot_lbl_dict = {}
                                if number_sp>0:
            
                                    if (number_sp>1)&(len(all_spots_for_gmm)>number_sp+2):
#                                         silhouette_avg_list={}
#                                         for n_clusters in np.arange(2,number_sp+1):

#                                             clusterer = KMeans(n_clusters=n_clusters+1, random_state=100)
#                                             cluster_labels = clusterer.fit_predict(all_spots_for_gmm)
#                                             silhouette_avg_list[n_clusters] =silhouette_score(all_spots_for_gmm, cluster_labels)

#                                         number_of_clusters = max(silhouette_avg_list, key=silhouette_avg_list.get)
                                        number_of_clusters = number_sp
                                        dpgmm = mixture.BayesianGaussianMixture(n_components=number_of_clusters, covariance_type="full").fit(all_spots_for_gmm)

                                        ### define a dictionary for every neighborhood a.k.a assign all the spots within the neighborhood

                                        for lbl1 in np.unique(dpgmm.predict(all_spots_for_gmm)):
                                            spot_lbl_dict[lbl1] = all_spots_for_gmm[np.where(dpgmm.predict(all_spots_for_gmm)==lbl1)]
                                    else:

                                        spot_lbl_dict[1] = all_spots_for_gmm


                                    df_cols = ["column", "row", "time_point", "field_index", "channel", "x", "y", "intensity_mean",
                                               "gaussian_intensity_mean", "spot_intensity", "gaussian_spot_intensity",
                                               "normalized_spot_intensity", "gaussian_normalized_spot_intensity"]


                                    labeled_spot_intensity_arr={}
                                    for key1 in spot_lbl_dict.keys():
                                        rows = []
                                        labeled_spot_intensity_arr[key1]=[]
                                        same_spots_time = []
                                        same_spots_time_loc = {}
                                        spots_list = spot_lbl_dict[key1]
                                        for spot_coor in spots_list:

                                            for key2 in rot_spot_coor_dict.keys():

                                                 for key_arr in rot_spot_coor_dict[key2]:
                                                    if (spot_coor == key_arr).sum()==2:

                                                        same_spots_time.append(int(key2))
                                                        same_spots_time_loc[str(key2)]=spot_coor

                                        sorted_time_points= np.sort(np.array(same_spots_time))
                                        current_time_point=sorted_time_points[0]
                                        current_spot_loc=same_spots_time_loc[str(current_time_point)]
                                        for stack_ind in np.array(list(index_to_time.keys()), dtype=int):

                                            if index_to_time[str(stack_ind)] in sorted_time_points:
                                                current_time_point = index_to_time[str(stack_ind)]
                                                current_spot_loc=same_spots_time_loc[str(current_time_point)]

                                            patch_intensity_mean = rotated_spot_patches_mean[str(int(sp_ch))][stack_ind]   
                                            patch_gaussian = ndimage.gaussian_filter(large_rot_spot_img_stack[str(int(sp_ch))][:,:,stack_ind], 
                                                                                      sigma=3)
                                            patch_gaussian_intensity_mean = np.true_divide(patch_gaussian.sum(),(patch_gaussian!=0).sum())
                                            patch_spot_intensity = large_rot_spot_img_stack[str(int(sp_ch))][int(current_spot_loc[0]),
                                                                                                             int(current_spot_loc[1]),
                                                                                                             stack_ind]
                                            patch_gaussian_spot_intensity = patch_gaussian[int(current_spot_loc[0]),
                                                                                           int(current_spot_loc[1])]
                                            norm_spot_intensity = patch_spot_intensity/(patch_intensity_mean+0.000001)
                                            gaussian_norm_spot_intensity = patch_gaussian_spot_intensity/(patch_gaussian_intensity_mean+0.000001)
                                            labeled_spot_intensity_arr[key1].append(norm_spot_intensity)

                                            rows.append({

                                                         "column": col, 
                                                         "row": row, 
                                                         "time_point": index_to_time[str(stack_ind)], 
                                                         "field_index": fov, 
                                                         "channel": sp_ch,
                                                         "x": int(current_spot_loc[0]),
                                                         "y": int(current_spot_loc[1]),
                                                         "intensity_mean": patch_intensity_mean, 
                                                         "gaussian_intensity_mean": patch_gaussian_intensity_mean, 
                                                         "spot_intensity": patch_spot_intensity,
                                                         "gaussian_spot_intensity": patch_gaussian_spot_intensity,
                                                         "normalized_spot_intensity": norm_spot_intensity,
                                                         "gaussian_normalized_spot_intensity": gaussian_norm_spot_intensity
                                                    })
                                            
                                        spot_signal_df = pd.DataFrame(rows, columns = df_cols)

                                        spot_intensity_tables = os.path.join(cell_tracking_folder, 'spot_intensity_tables')
                                        if os.path.isdir(spot_intensity_tables) == False:
                                            os.mkdir(spot_intensity_tables)
                                        ######################## csv tables
                                        spot_csv_tables = os.path.join(spot_intensity_tables, 'complete_tables')
                                        spot_csv_file_name = 'spot_intensity_col' + str(col) + r'_row' + str(row)+ r'_field' + str(fov) + r'_Cell' + str(id_) + r'_Ch' + str(int(sp_ch))+ r'_spot' + str(key1) + r'.csv'
                                        if os.path.isdir(spot_csv_tables) == False:
                                            os.mkdir(spot_csv_tables)
                                        spot_intensity_full_name = os.path.join(spot_csv_tables, spot_csv_file_name)
                                        spot_signal_df.to_csv(spot_intensity_full_name)

                                        ######################## max intensity tables
                                        max_intensity_tables = os.path.join(spot_intensity_tables, 'max_intensity_tables')
                                        max_intensity_file_name = 'max_intensity_col' + str(col) + r'_row' + str(row)+ r'_field' + str(fov) + r'_Cell' + str(id_) + r'_Ch' + str(int(sp_ch))+ r'_spot' + str(key1) + r'.trk'
                                        if os.path.isdir(max_intensity_tables) == False:
                                            os.mkdir(max_intensity_tables)
                                        max_intensity_full_name = os.path.join(max_intensity_tables, max_intensity_file_name)
                                        sub_spot_df = spot_signal_df[["x", "y", "spot_intensity", "time_point"]]
                                        sub_spot_df["spot_label"] = np.zeros(len(sub_spot_df), dtype=float)
                                        sub_spot_df.to_csv(max_intensity_full_name, sep='\t', index=False, header=False)

                                        #########################normalized intensity tables

                                        norm_intensity_tables = os.path.join(spot_intensity_tables, 'norm_intensity_tables')
                                        norm_intensity_file_name = 'max_intensity_col' + str(col) + r'_row' + str(row)+ r'_field' + str(fov) + r'_Cell' + str(id_) + r'_Ch' + str(int(sp_ch))+ r'_spot' + str(key1) + r'.trk'
                                        if os.path.isdir(norm_intensity_tables) == False:
                                            os.mkdir(norm_intensity_tables)
                                        norm_intensity_full_name = os.path.join(norm_intensity_tables, norm_intensity_file_name)
                                        sub_spot_df = spot_signal_df[["x", "y", "normalized_spot_intensity", "time_point"]]
                                        sub_spot_df["spot_label"] = np.zeros(len(sub_spot_df), dtype=float)
                                        sub_spot_df.to_csv(norm_intensity_full_name, sep='\t', index=False, header=False)
                                        ####################
                                            

                                        for st_ind in range(anotated_large_rot_spot_img_stack.shape[2]):

                                            temp_img = anotated_large_rot_spot_img_stack[:, :,st_ind]
                                            circ_for_draw = self.ImageAnalyzer.COORDINATES_TO_CIRCLE(spot_signal_df[["x", "y"]].to_numpy()[st_ind].reshape(1,2),
                                                                                                     temp_img,circ_radius = 5)
                                            
                                            temp_img[circ_for_draw != 0] = 255
                                            anotated_large_rot_spot_img_stack[:, :,st_ind]=temp_img

                                    ########## save annotated and unannotated spot patches
                                    spot_img_file_name = 'spot_img_col' + str(col) + r'_row' + str(row)+ r'_field' + str(fov) + r'_Cell' + str(id_) + r'_Ch' + str(int(sp_ch)) + r'.ome.tif'
                                    spot_img_folder = os.path.join(cell_tracking_folder, 'spot_image_patches')
                                    if os.path.isdir(spot_img_folder) == False:
                                        os.mkdir(spot_img_folder)
                                    spot_img_full_name = os.path.join(spot_img_folder, spot_img_file_name)
                                    writer = OmeTiffWriter()
                                    writer.save(np.rollaxis(anotated_large_rot_spot_img_stack, 2, 0),spot_img_full_name,
                                                dimension_order="TYX")
                                    
                                    spot_img_file_name = 'unannotated_spot_img_col' + str(col) + r'_row' + str(row)+ r'_field' + str(fov) + r'_Cell' + str(id_) + r'_Ch' + str(int(sp_ch)) + r'.ome.tif'
                                    spot_img_folder = os.path.join(cell_tracking_folder, 'unannotated_spot_image_patches')
                                    if os.path.isdir(spot_img_folder) == False:
                                        os.mkdir(spot_img_folder)
                                    spot_img_full_name = os.path.join(spot_img_folder, spot_img_file_name)
                                    writer = OmeTiffWriter()
                                    writer.save(np.rollaxis(unanotated_large_rot_spot_img_stack, 2, 0),spot_img_full_name, dimension_order="TYX")
                            
                            
                            
                            single_patch_name = 'nuclei_image_for_col' + str(col) + r'_row' + str(row) + r'_field' + str(fov) + r'_cell' + str(id_) + r'.ome.tif'
                            aligned_single_patch_name = 'aligned_nuclei_image_for_col' + str(col) + r'_row' + str(row) + r'_field' + str(fov) + r'_cell' + str(id_) + r'.ome.tif'
                            track_gifs_folder = os.path.join(cell_tracking_folder, 'single_track_images')
                            if os.path.isdir(track_gifs_folder) == False:
                                os.mkdir(track_gifs_folder)
                            single_patches_full_name = os.path.join(track_gifs_folder, single_patch_name)
                            writer = OmeTiffWriter()
                            writer.save(np.rollaxis(large_cell_stack, 2, 0),single_patches_full_name,
                                        dimension_order="TYX")
                            
                            single_aligned_patches_full_name = os.path.join(track_gifs_folder, aligned_single_patch_name)
                            writer = OmeTiffWriter()
                            writer.save(np.rollaxis(large_rotcell_stack, 2, 0),single_aligned_patches_full_name,
                                        dimension_order="TYX")
                            
#                             imageio.mimsave(single_patches_full_name, 
#                                             [self.save_nuc_patches(large_cell_stack[:,:,:,i], 
#                                             large_rotcell_stack[:,:,:,i]) for i in range(len(single_track))])
#                             imageio.mimsave(Raw_full_name,nuc_patches)
#                             imageio.mimsave(Aligned_full_name, rotated_nuc)
                            trackcsv_file_name = 'track_table_for_Col' + str(col) + r'_row' + str(5)+r'_Field' + str(fov) + r'_Cell' + str(id_) + r'.csv'
                            track_table_folder = os.path.join(cell_tracking_folder, 'single_track_tables')
                            if os.path.isdir(track_table_folder) == False:
                                os.mkdir(track_table_folder)
                            trackcsv_full_name = os.path.join(track_table_folder, trackcsv_file_name)    
                            single_track_copy.to_csv(trackcsv_full_name)    
                                
        seconds2 = time.time()
        
        diff=seconds2-seconds1
        print('\nTotal Processing Time (Minutes):',diff/60)
    def save_nuc_patches(self,im1,im2):
        
        fig, ax = plt.subplots(1, 2, figsize=(12, 6))
        ax[0].imshow(im1)
        ax[0].axis('off')
        ax[0].set_title('Original')
        ax[1].imshow(im2)
        ax[1].set_title('Aligned')
        ax[1].axis('off')

        fig.canvas.draw()  # draw the canvas, cache the renderer
        image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
        image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        plt.close(fig)

        return image


    #the main nuclear mask function
    #calls the nuclear mask maker on each image at a given col/row/time/fov
    def BATCH_ANALYZER(self, col,row,fov,t, lock=None, value=None): 
        ai = 1
        self.counter.increment()    
        #process safe counter, will update before the actual images finish (ex. will say 16/50, but only the first 8 are actually done and the other 8 are still in progress)       
        progress_message = "Progress: %d%% [%d / %d] masks" % (self.counter.value() / self.masks * 100, self.counter.value(), self.masks)
        # Don't use print() as it will print in new line every time.
        sys.stdout.write("\r" + progress_message)
        sys.stdout.flush()
        self.df_checker = self.Meta_Data_df.loc[(self.Meta_Data_df['column'] == str(col)) & 
                                              (self.Meta_Data_df['row'] == str(row)) & 
                                              (self.Meta_Data_df['field_index'] == str(fov)) & 
                                              (self.Meta_Data_df['time_point'] == str(t))]
        if self.df_checker.empty == False:
        
            loadedimg_formask = self.IMG_FOR_NUC_MASK()
            if loadedimg_formask.ndim ==2:
                
                ImageForNucMask = cv2.normalize(loadedimg_formask, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
                nuc_bndry, nuc_mask = self.ImageAnalyzer.nuclei_segmenter(ImageForNucMask, self.Meta_Data_df["PixPerMic"].iloc[0])
                
                labeled_nuc, number_nuc = label(nuc_mask)
                nuc_labels = np.unique(labeled_nuc)
                if nuc_labels.max()>0:
                    dist_img = distance_transform_edt(nuc_mask)
                    dist_props = regionprops_table(labeled_nuc, dist_img, properties=('label', 'max_intensity'))
                    radial_dist_df = pd.DataFrame(dist_props)

                    
                    data = { "Experiment":[self.experiment_name]*number_nuc,
                             "column": [col]*number_nuc, "row": [row]*number_nuc, 
                             "time_point": [t]*number_nuc, "field_index": [fov]*number_nuc,
                             "z_slice": ["max_project"]*number_nuc, "action_index":[ai]*number_nuc}
                    df = pd.DataFrame(data)
                    regions = regionprops(labeled_nuc, loadedimg_formask)
                    props = regionprops_table(labeled_nuc, loadedimg_formask, properties=(
                                            'label', 'centroid', 'orientation', 'major_axis_length', 'minor_axis_length',
                                            'area', 'max_intensity', 'min_intensity', 'mean_intensity',
                                            'orientation', 'perimeter', 'solidity'))
                    pixpermicron = np.asarray(self.Meta_Data_df["PixPerMic"].iloc[0]).astype(float)
                    props_df = pd.DataFrame(props)
                    props_df['major_axis_length'] = pixpermicron*props_df['major_axis_length']
                    props_df['minor_axis_length'] = pixpermicron*props_df['minor_axis_length']
                    props_df['area'] = pixpermicron*pixpermicron*props_df['area']
                    props_df['perimeter'] = pixpermicron*props_df['perimeter']
                    image_cells_df = pd.concat([df,props_df], axis=1)
#                     self.cell_pd_list = pd.concat([self.cell_pd_list, image_cells_df], axis=0, ignore_index=True)
                    self.cell_pd_list.append(image_cells_df)

                if self.gui_params.NucMaskCheckBox_status_check == True:
                    nuc_mask_output_folder = os.path.join(self.output_folder, 'nuclei_masks')
                    if os.path.isdir(nuc_mask_output_folder) == False:
                        os.mkdir(nuc_mask_output_folder) 

                    mask_file_name = ['Nuclei_Mask_for_Col' + str(col) + r'_row' + str(row)+
                                      r'_Time' + str(t) + r'_Field' + str(fov) + r'.tif']
                    mask_full_name = os.path.join(nuc_mask_output_folder, mask_file_name[0])
#                         cv2.imwrite(mask_full_name,nuc_mask)

                    imwrite(mask_full_name,nuc_mask)
                    #CHANGED
                    #resume from interrupt trimmer that trims the related lines from the metadata file
                    #only runs for the headless version
                    if self.is_gui == False:
                        with lock:
                            self.trimmer.trim_metadata(col, row, t, fov)
                            #print(f'> sleeping for {value}')
                            sleep(value)
                    
                     

            #right now, this shouldn't run ever unless there's a bug
            else:

                nuc_bndry, nuc_mask = self.Z_STACK_NUC_SEGMENTER(ImageForNucMask)
                label_nuc_stack = self.Z_STACK_NUC_LABLER(ImageForLabel)
    
            ch1_xyz, ch1_xyz_3D, ch1_final_spots, ch2_xyz, ch2_xyz_3D, ch2_final_spots, ch3_xyz, ch3_xyz_3D, ch3_final_spots, ch4_xyz, ch4_xyz_3D, ch4_final_spots, ch5_xyz, ch5_xyz_3D, ch5_final_spots  = self.IMAGE_FOR_SPOT_DETECTION( ImageForNucMask)

            if self.gui_params.NucMaxZprojectCheckBox_status_check == True:

                if self.gui_params.SpotMaxZProject_status_check == True:

                        if self.gui_params.SpotCh1CheckBox_status_check == True:
                            if ch1_xyz!=[]:
                                ch1_xyz_round = np.round(np.asarray(ch1_xyz)).astype('int')
                                ch1_spot_nuc_labels = labeled_nuc[ch1_xyz_round[:,0], ch1_xyz_round[:,1]]

                                ch1_num_spots = ch1_xyz.__len__()
                                if labeled_nuc.max()>0:
                                    ch1_radial = self.RADIAL_DIST_CALC(ch1_xyz_round,ch1_spot_nuc_labels,
                                                                       radial_dist_df, dist_img)
                                else:
                                    ch1_radial=np.nan
                                data = { "Experiment":[self.experiment_name]*ch1_num_spots,
                                         "column": [col]*ch1_num_spots, "row": [row]*ch1_num_spots, 
                                         "time_point": [t]*ch1_num_spots, "field_index": [fov]*ch1_num_spots,
                                         "z_slice": ["max_project"]*ch1_num_spots, "channel": [1]*ch1_num_spots,
                                         "action_index": [ai]*ch1_num_spots, "cell_index": ch1_spot_nuc_labels,
                                         "x_location": ch1_xyz[:,0], "y_location": ch1_xyz[:,1],
                                         "z_location": ch1_xyz[:,2], 
                                         "radial_distance":ch1_radial
                                    
                                       }
                                df_ch1 = pd.concat([pd.DataFrame(data), ch1_final_spots], axis=1)



#                                 self.ch1_spot_df=pd.concat([self.ch1_spot_df,df_ch1],ignore_index=True)
                                self.ch1_spot_df_list.append(df_ch1)
    
#                                 cell_indices=cell_indices[cell_indices>0]
        
#                                 for ii in range(len(self.cell_pd_list)):
                                    
#                                     if ((self.cell_pd_list[ii]["column"].iloc[0] == col)&(self.cell_pd_list[ii]["row"].iloc[0] == row)&
#                                         (self.cell_pd_list[ii]["time_point"].iloc[0] == t)&(self.cell_pd_list[ii]["field_index"].iloc[0] == fov)):
                                    
#                                         df_index = ii
#                                         break
                
#                                 for ci in cell_indices:
                
#                                     cell_specific_spots=df_ch1.loc[df_ch1['cell_index']==ci]
#                                     num_spots=cell_specific_spots.__len__()
#                                     column_name = 'ch1_number_of_spots'
#                                     if column_name not in self.cell_pd_list[df_index].columns:
#                                         self.cell_pd_list[df_index][column_name]=np.zeros((len(self.cell_pd_list[df_index]),1),dtype=int)

#                                     row_index=self.cell_pd_list[df_index].loc[(self.cell_pd_list[df_index]['column']==col) & 
#                                                                   (self.cell_pd_list[df_index]['row']==row)&
#                                                                   (self.cell_pd_list[df_index]['field_index']==fov)&
#                                                                   (self.cell_pd_list[df_index]['label']==ci)&
#                                                                   (self.cell_pd_list[df_index]['time_point']==t)&
#                                                                   (self.cell_pd_list[df_index]['action_index']==ai)].index[0]
#                                     self.cell_pd_list[df_index].loc[row_index,column_name]=num_spots
                                    
                        if self.gui_params.SpotCh2CheckBox_status_check == True:
                            if ch2_xyz!=[]:
                                
                                ch2_xyz_round = np.round(np.asarray(ch2_xyz)).astype('int')
                                ch2_spot_nuc_labels = labeled_nuc[ch2_xyz_round[:,0], ch2_xyz_round[:,1]]
                                ch2_num_spots = ch2_xyz.__len__()
                                if labeled_nuc.max()>0:
                                    ch2_radial = self.RADIAL_DIST_CALC(ch2_xyz_round,ch2_spot_nuc_labels,
                                                                       radial_dist_df, dist_img)
                                else:
                                    ch2_radial=np.nan
                                data = { "Experiment":[self.experiment_name]*ch2_num_spots,
                                         "column": [col]*ch2_num_spots, "row": [row]*ch2_num_spots, 
                                         "time_point": [t]*ch2_num_spots, "field_index": [fov]*ch2_num_spots,
                                         "z_slice": ["max_project"]*ch2_num_spots, "channel": [2]*ch2_num_spots,
                                         "action_index": [ai]*ch2_num_spots, "cell_index": ch2_spot_nuc_labels,
                                         "x_location": ch2_xyz[:,0], "y_location": ch2_xyz[:,1],
                                         "z_location": ch2_xyz[:,2],
                                         "radial_distance":ch2_radial
                                       }
#                                 df_ch2 = pd.DataFrame(data)
                                df_ch2 = pd.concat([pd.DataFrame(data), ch2_final_spots], axis=1)
    

                                self.ch2_spot_df_list.append(df_ch2)
    
#                                 cell_indices=cell_indices[cell_indices>0]
        
#                                 for ii in range(len(self.cell_pd_list)):
                                    
#                                     if ((self.cell_pd_list[ii]["column"].iloc[0] == col)&(self.cell_pd_list[ii]["row"].iloc[0] == row)&
#                                         (self.cell_pd_list[ii]["time_point"].iloc[0] == t)&(self.cell_pd_list[ii]["field_index"].iloc[0] == fov)):
                                    
#                                         df_index = ii
                                        
#                                         break
                
#                                 for ci in cell_indices:
                
#                                     cell_specific_spots=df_ch2.loc[df_ch2['cell_index']==ci]
#                                     num_spots=cell_specific_spots.__len__()
#                                     temp_df_ch1 = self.cell_pd_list[df_index]
#                                     column_name = 'ch2_number_of_spots'
#                                     if column_name not in self.cell_pd_list[df_index].columns:
#                                         self.cell_pd_list[df_index][column_name]=np.zeros((len(self.cell_pd_list[df_index]),1),dtype=int)

#                                     row_index=self.cell_pd_list[df_index].loc[(self.cell_pd_list[df_index]['column']==col) & 
#                                                                   (self.cell_pd_list[df_index]['row']==row)&
#                                                                   (self.cell_pd_list[df_index]['field_index']==fov)&
#                                                                   (self.cell_pd_list[df_index]['label']==ci)&
#                                                                   (self.cell_pd_list[df_index]['time_point']==t)&
#                                                                   (self.cell_pd_list[df_index]['action_index']==ai)].index[0]
#                                     print("row index",row_index)
                                    
#                                 self.cell_pd_list[df_index].loc[row_index,column_name]=num_spots
#                                     print("DF",self.cell_pd_list[df_index].columns)
                        if self.gui_params.SpotCh3CheckBox_status_check == True:
                            if ch3_xyz!=[]:

                                ch3_xyz_round = np.round(np.asarray(ch3_xyz)).astype('int')
                                ch3_spot_nuc_labels = labeled_nuc[ch3_xyz_round[:,0], ch3_xyz_round[:,1]]
                                ch3_num_spots = ch3_xyz.__len__()
                                if labeled_nuc.max()>0:
                                    ch3_radial = self.RADIAL_DIST_CALC(ch3_xyz_round,ch3_spot_nuc_labels,
                                                                       radial_dist_df, dist_img)
                                else:
                                    ch3_radial=np.nan
                                data = { "Experiment":[self.experiment_name]*ch3_num_spots,
                                         "column": [col]*ch3_num_spots, "row": [row]*ch3_num_spots, 
                                         "time_point": [t]*ch3_num_spots, "field_index": [fov]*ch3_num_spots,
                                         "z_slice": ["max_project"]*ch3_num_spots, "channel": [3]*ch3_num_spots,
                                         "action_index": [ai]*ch3_num_spots, "cell_index": ch3_spot_nuc_labels,
                                         "x_location": ch3_xyz[:,0], "y_location": ch3_xyz[:,1],
                                         "z_location": ch3_xyz[:,2],
                                         "radial_distance": ch3_radial
                                       }

#                                 df_ch3 = pd.DataFrame(data)
                                df_ch3 = pd.concat([pd.DataFrame(data), ch3_final_spots], axis=1)
    
                                
                                self.ch3_spot_df_list.append(df_ch3)
    
#                                 cell_indices=cell_indices[cell_indices>0]
#                                 for ii in range(len(self.cell_pd_list)):
                                    
#                                     if ((self.cell_pd_list[ii]["column"].iloc[0] == col)&(self.cell_pd_list[ii]["row"].iloc[0] == row)&
#                                         (self.cell_pd_list[ii]["time_point"].iloc[0] == t)&(self.cell_pd_list[ii]["field_index"].iloc[0] == fov)):
                                    
#                                         df_index = ii
#                                         break
#                                 for ci in cell_indices:
                
#                                     cell_specific_spots=df_ch3.loc[df_ch3['cell_index']==ci]
#                                     num_spots=cell_specific_spots.__len__()
#                                     column_name = 'ch3_number_of_spots'
#                                     if column_name not in self.cell_pd_list[df_index].columns:
#                                         self.cell_pd_list[df_index][column_name]=np.zeros((len(self.cell_pd_list[df_index]),1),dtype=int)

#                                     row_index=self.cell_pd_list[df_index].loc[(self.cell_pd_list[df_index]['column']==col) & 
#                                                                   (self.cell_pd_list[df_index]['row']==row)&
#                                                                   (self.cell_pd_list[df_index]['field_index']==fov)&
#                                                                   (self.cell_pd_list[df_index]['label']==ci)&
#                                                                   (self.cell_pd_list[df_index]['time_point']==t)&
#                                                                   (self.cell_pd_list[df_index]['action_index']==ai)].index[0]
#                                     self.cell_pd_list[df_index].loc[row_index,column_name]=num_spots

                        if self.gui_params.SpotCh4CheckBox_status_check == True:
                            if ch4_xyz!=[]:

                                ch4_xyz_round = np.round(np.asarray(ch4_xyz)).astype('int')
                                ch4_spot_nuc_labels = labeled_nuc[ch4_xyz_round[:,0], ch4_xyz_round[:,1]]
                                ch4_num_spots = ch4_xyz.__len__()
                                if labeled_nuc.max()>0:
                                    ch4_radial = self.RADIAL_DIST_CALC(ch4_xyz_round,ch4_spot_nuc_labels,
                                                                       radial_dist_df, dist_img)
                                else:
                                    ch4_radial=np.nan
                                data = { "Experiment":[self.experiment_name]*ch4_num_spots,
                                         "column": [col]*ch4_num_spots, "row": [row]*ch4_num_spots, 
                                         "time_point": [t]*ch4_num_spots, "field_index": [fov]*ch4_num_spots,
                                         "z_slice": ["max_project"]*ch4_num_spots, "channel": [4]*ch4_num_spots,
                                         "action_index": [ai]*ch4_num_spots, "cell_index": ch4_spot_nuc_labels,
                                         "x_location": ch4_xyz[:,0], "y_location": ch4_xyz[:,1],
                                         "z_location": ch4_xyz[:,2],
                                         "radial_distance":ch4_radial
                                       }

#                                 df_ch4 = pd.DataFrame(data)
                                df_ch4 = pd.concat([pd.DataFrame(data), ch4_final_spots], axis=1)

                               
                                self.ch4_spot_df_list.append(df_ch4)

#                                 cell_indices=cell_indices[cell_indices>0]
#                                 for ii in range(len(self.cell_pd_list)):
                                    
#                                     if ((self.cell_pd_list[ii]["column"].iloc[0] == col)&(self.cell_pd_list[ii]["row"].iloc[0] == row)&
#                                         (self.cell_pd_list[ii]["time_point"].iloc[0] == t)&(self.cell_pd_list[ii]["field_index"].iloc[0] == fov)):
                                    
#                                         df_index = ii
#                                         break
#                                 for ci in cell_indices:
                
#                                     cell_specific_spots=df_ch4.loc[df_ch4['cell_index']==ci]
#                                     num_spots=cell_specific_spots.__len__()
#                                     column_name = 'ch4_number_of_spots'
#                                     if column_name not in self.cell_pd_list[df_index].columns:
#                                         self.cell_pd_list[df_index][column_name]=np.zeros((len(self.cell_pd_list[df_index]),1),dtype=int)

#                                     row_index=self.cell_pd_list[df_index].loc[(self.cell_pd_list[df_index]['column']==col) & 
#                                                                   (self.cell_pd_list[df_index]['row']==row)&
#                                                                   (self.cell_pd_list[df_index]['field_index']==fov)&
#                                                                   (self.cell_pd_list[df_index]['label']==ci)&
#                                                                   (self.cell_pd_list[df_index]['time_point']==t)&
#                                                                   (self.cell_pd_list[df_index]['action_index']==ai)].index[0]
#                                     self.cell_pd_list[df_index].loc[row_index,column_name]=num_spots
                                    
                        if self.gui_params.SpotCh5CheckBox_status_check == True:
                            if ch5_xyz!=[]:

                                ch5_xyz_round = np.round(np.asarray(ch5_xyz)).astype('int')
                                ch5_spot_nuc_labels = labeled_nuc[ch5_xyz_round[:,0], ch5_xyz_round[:,1]]
                                ch5_num_spots = ch5_xyz.__len__()
                                if labeled_nuc.max()>0:
                                    ch5_radial = self.RADIAL_DIST_CALC(ch5_xyz_round,ch5_spot_nuc_labels,
                                                                       radial_dist_df, dist_img)
                                else:
                                    ch5_radial=np.nan


                                data = { "Experiment":[self.experiment_name]*ch5_num_spots,
                                         "column": [col]*ch5_num_spots, "row": [row]*ch5_num_spots, 
                                         "time_point": [t]*ch5_num_spots, "field_index": [fov]*ch5_num_spots,
                                         "z_slice": ["max_project"]*ch5_num_spots, "channel": [5]*ch5_num_spots,
                                         "action_index": [ai]*ch5_num_spots, "cell_index": ch5_spot_nuc_labels,
                                         "x_location": ch5_xyz[:,0], "y_location": ch5_xyz[:,1],
                                         "z_location": ch5_xyz[:,2],
                                         "radial_distance": ch5_radial
                                       }

#                                 df_ch5 = pd.DataFrame(data)
                                df_ch5 = pd.concat([pd.DataFrame(data), ch5_final_spots], axis=1)
                                self.ch5_spot_df_list.append(df_ch5)
            
                                
#                                 cell_indices=cell_indices[cell_indices>0]
#                                 for ii in range(len(self.cell_pd_list)):
                                    
#                                     if ((self.cell_pd_list[ii]["column"].iloc[0] == col)&(self.cell_pd_list[ii]["row"].iloc[0] == row)&
#                                         (self.cell_pd_list[ii]["time_point"].iloc[0] == t)&(self.cell_pd_list[ii]["field_index"].iloc[0] == fov)):
                                    
#                                         df_index = ii
#                                         break
#                                 for ci in cell_indices:
                
#                                     cell_specific_spots=df_ch5.loc[df_ch5['cell_index']==ci]
#                                     num_spots=cell_specific_spots.__len__()
#                                     column_name = 'ch5_number_of_spots'
#                                     if column_name not in self.cell_pd_list[df_index].columns:
#                                         self.cell_pd_list[df_index][column_name]=np.zeros((len(self.cell_pd_list[df_index]),1),dtype=int)

#                                     row_index=self.cell_pd_list[df_index].loc[(self.cell_pd_list[df_index]['column']==col) & 
#                                                                   (self.cell_pd_list[df_index]['row']==row)&
#                                                                   (self.cell_pd_list[df_index]['field_index']==fov)&
#                                                                   (self.cell_pd_list[df_index]['label']==ci)&
#                                                                   (self.cell_pd_list[df_index]['time_point']==t)&
#                                                                   (self.cell_pd_list[df_index]['action_index']==ai)].index[0]
#                                     self.cell_pd_list[df_index].loc[row_index,column_name]=num_spots

    def Calculate_Spot_Distances(self, row, col):
        self.distance_counter.increment()    
        #print(self.counter.value())
        # #process safe counter for the spot distances        
        progress_message = "Progress: %d%% [%d / %d] files" % (self.distance_counter.value() / self.distances * 100, self.distance_counter.value(), self.distances)
        # Don't use print() as it will print in new line every time.
        sys.stdout.write("\r" + progress_message)
        sys.stdout.flush()
        select_cell = self.cell_df.loc[(self.cell_df['column'] == col) & (self.cell_df['row'] == row)]    

        spot_pd_dict = {}

        if self.ch1_spot_df.empty == False:

            cell_spots_ch1 = self.ch1_spot_df.loc[(self.ch1_spot_df['column'] == col) & (self.ch1_spot_df['row'] == row)] 

            spot_pd_dict['ch1'] = cell_spots_ch1

        if self.ch2_spot_df.empty == False:

            cell_spots_ch2 = self.ch2_spot_df.loc[(self.ch2_spot_df['column'] == col) & (self.ch2_spot_df['row'] == row)] 

            spot_pd_dict['ch2'] = cell_spots_ch2

        if self.ch3_spot_df.empty == False:

            cell_spots_ch3 = self.ch3_spot_df.loc[(self.ch3_spot_df['column'] == col) & (self.ch3_spot_df['row'] == row) ] 

            spot_pd_dict['ch3'] = cell_spots_ch3

        if self.ch4_spot_df.empty == False:

            cell_spots_ch4 = self.ch4_spot_df.loc[(self.ch4_spot_df['column'] == col) & (self.ch4_spot_df['row'] == row)] 

            spot_pd_dict['ch4'] = cell_spots_ch4

        if self.ch5_spot_df.empty == False:

            cell_spots_ch5 = self.ch5_spot_df.loc[(self.ch5_spot_df['column'] == col) & (self.ch5_spot_df['row'] == row)] 

            spot_pd_dict['ch5'] = cell_spots_ch5

        ch_distances = []

        for key1 in spot_pd_dict.keys():
            for key2 in spot_pd_dict.keys():

                ch_distance1 = key2 + r'_' + key1
                if ch_distance1 in ch_distances:
                    pass
                else:
                    ch_distance = key1 + r'_' + key2
                    ch_distances.append(ch_distance)
                    
                    dist_pd = self.DISTANCE_calculator(spot_pd_dict,key1,key2,select_cell, row, col)
                    well_spot_dist_folder = os.path.join(self.output_folder, 'well_spots_distances')
                    if os.path.isdir(well_spot_dist_folder) == False:
                        os.mkdir(well_spot_dist_folder)

                    spot_dist_filename = 'SpotDistances_' + key1 +'_' +key2 +'_' + WELL_PLATE_ROWS[row-1] + str(col) + r'.csv'
                    spot_dist_well_csv_full_name = os.path.join(well_spot_dist_folder, spot_dist_filename)
                    dist_pd.to_csv(path_or_buf=spot_dist_well_csv_full_name, encoding='utf8')
                    #print(spot_dist_filename)

    def DISTANCE_calculator(self, spot_pd_dict, key1, key2,select_cell, row, col):
    
        fovs = np.unique(np.asarray(select_cell['field_index'], dtype=int))
        time_points = np.unique(np.asarray(select_cell['time_point'], dtype=int))
        actionindices = np.unique(np.asarray(select_cell['action_index'], dtype=int))
        z_slices = np.unique(np.asarray(select_cell['z_slice']))
        dist_pd = pd.DataFrame()
        for f in fovs:
            for t in time_points:
                for a in actionindices:
                    for z in z_slices:
                        
                        # cells_in_field = select_cell.loc[ 
                        #                       (select_cell['time_point'] == t) & 
                        #                       (select_cell['field_index'] == f)& 
                        #                       (select_cell['action_index'] == a)&
                        #                       (select_cell['z_slice'] == z)] 
                        cells_in_field = np.unique(np.asarray(select_cell['cell_index'], dtype=int))

                        for c in cells_in_field:
                            loci1 = spot_pd_dict[key1].loc[ 
                                                  (spot_pd_dict[key1]['time_point'] == t) & 
                                                  (spot_pd_dict[key1]['field_index'] == f)& 
                                                  (spot_pd_dict[key1]['action_index'] == a)&
                                                  (spot_pd_dict[key1]['z_slice'] == z) &
                                                  (spot_pd_dict[key1]['cell_index'] == c)] 
                            loci2 = spot_pd_dict[key2].loc[ 
                                                  (spot_pd_dict[key2]['time_point'] == t) & 
                                                  (spot_pd_dict[key2]['field_index'] == f)& 
                                                  (spot_pd_dict[key2]['action_index'] == a)&
                                                  (spot_pd_dict[key2]['z_slice'] == z) &
                                                  (spot_pd_dict[key2]['cell_index'] == c)] 

                            for ind1, locus1 in loci1.iterrows():
                                for ind2, locus2 in loci2.iterrows():
                                    if ind1!=ind2:
                                        dist_2d =  math.sqrt(
                                                            math.pow((locus1['x_location'] - locus2['x_location']),2) +
                                                            math.pow((locus1['y_location'] - locus2['y_location']),2) 
                                                            )

                                        dist_3d =  math.sqrt(
                                                            math.pow((locus1['x_location'] - locus2['x_location']),2) +
                                                            math.pow((locus1['y_location'] - locus2['y_location']),2) +
                                                            math.pow((locus1['z_location'] - locus2['z_location']),2)
                                                            )

                                        s1 = key1 + '_spot_index(1)'
                                        s2 = key2 + '_spot_index(2)'
                                        data = {"Experiment":[self.experiment_name],
                                                'column':col, 'row': row, 'time_point': [t], 'field_index': [f], 
                                                'action_index': [a], 'z_slice':[z], 'cell_index': [c],
                                                str(s1): [ind1], str(s2): [ind2], 
                                                'XY-Distance(pixels)': [dist_2d], 'XYZ-Distance(pixels)': [dist_3d],
                                                'XY-Distance(micron)': [dist_2d*np.asarray(self.Meta_Data_df["PixPerMic"].iloc[0]).astype(float)], 
                                                'XYZ-Distance(micron)':[dist_3d*np.asarray(self.Meta_Data_df["PixPerMic"].iloc[0]).astype(float)]}
                                        temp_df = pd.DataFrame(data)
                                        dist_pd = pd.concat([dist_pd, temp_df], ignore_index=True)

        return dist_pd

                                                        
    def IMG_FOR_NUC_MASK(self):
        
        if self.df_checker.empty == False:

            if self.gui_params.NucMaxZprojectCheckBox_status_check == True:
                maskchannel = self.gui_params.NucleiChannel_index
                imgformask = self.df_checker.loc[(self.df_checker['channel'] == maskchannel)]
                loadedimg_formask = self.ImageAnalyzer.max_z_project(imgformask)
                

            else:
                z_imglist=[]
                maskchannel = self.gui_params.NucleiChannel_index
                imgformask = self.df_checker.loc[(self.df_checker['channel'] == maskchannel)]

                for index, row in imgformask.iterrows():
                    
                    if row['ImageName']=="dask_array":
                        im = row["Type"].compute()
                    else: 
                        im = mpimg.imread(row['ImageName'])
                    im_uint8 = cv2.normalize(im, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
                    z_imglist.append( np.asarray(im_uint8))
                loadedimg_formask = np.stack(z_imglist, axis=2)
        
        return loadedimg_formask
    
    def RAW_IMAGE_LOADER(self, maskchannel):
        
        if self.df_checker.empty == False:

            if self.gui_params.NucMaxZprojectCheckBox_status_check == True:

                imgformask = self.df_checker.loc[(self.df_checker['channel'] == maskchannel)]
                loadedimg_formask = self.ImageAnalyzer.max_z_project(imgformask)
                loaded_image = loadedimg_formask

            else:
                z_imglist=[]
                imgformask = self.df_checker.loc[(self.df_checker['channel'] == maskchannel)]

                for index, row in imgformask.iterrows():
                    
                    if row['ImageName']=="dask_array":
                        im = row["Type"].compute()
                    else: 
                        im = mpimg.imread(row['ImageName'])
                    z_imglist.append( np.asarray(im))
                loaded_image = np.stack(z_imglist, axis=2)
        
        return loaded_image

    def Z_STACK_NUC_SEGMENTER(self, ImageForNucMask):
        
        nuc_bndry_list, nuc_mask_list = [], []
        w, h, z_plane = ImageForNucMask.shape

        for z in range(z_plane):

            single_z_img = ImageForNucMask[:,:,z]
            nuc_bndry_single, nuc_mask_single = self.ImageAnalyzer.nuclei_segmenter(single_z_img,
                                                                                     self.Meta_Data_df["PixPerMic"].iloc[0])

            nuc_bndry_list.append( np.asarray(nuc_bndry_single))
            nuc_mask_list.append( np.asarray(nuc_mask_single))
        nuc_bndry = np.stack(nuc_bndry_list, axis=2)
        nuc_mask = np.stack(nuc_mask_list, axis=2)
    
        return nuc_bndry, nuc_mask
    
    def Z_STACK_NUC_LABLER(self, ImageForLabel):
        
        label_nuc_list = []
        w, h, z_plane = ImageForLabel.shape

        for z in range(z_plane):

            single_z_img = ImageForLabel[:,:,z]
            labeled_nuc, number_nuc = label(single_z_img)

            label_nuc_list.append( np.asarray(labeled_nuc))
            
        label_nuc_stack = np.stack(label_nuc_list, axis=2)
    
        return label_nuc_stack
    
    
    def IMAGE_FOR_SPOT_DETECTION(self, ImageForNucMask):
        
        ch1_xyz, ch2_xyz, ch3_xyz, ch4_xyz, ch5_xyz = [],[],[],[],[]
        ch1_xyz_3D, ch2_xyz_3D, ch3_xyz_3D, ch4_xyz_3D, ch5_xyz_3D = [],[],[],[],[]
        ch1_final_spots, ch2_final_spots, ch3_final_spots, ch4_final_spots, ch5_final_spots = pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        
        if self.gui_params.SpotCh1CheckBox_status_check == True:
            
            imgforspot = self.df_checker.loc[(self.df_checker['channel'] == '1')]
            ch1_xyz, ch1_xyz_3D, ch1_final_spots = self.XYZ_SPOT_COORDINATES( imgforspot, ImageForNucMask, 'Ch1')
                
        if self.gui_params.SpotCh2CheckBox_status_check == True:
            
            imgforspot = self.df_checker.loc[(self.df_checker['channel'] == '2')]
            ch2_xyz, ch2_xyz_3D, ch2_final_spots = self.XYZ_SPOT_COORDINATES( imgforspot, ImageForNucMask, 'Ch2')
                
        if self.gui_params.SpotCh3CheckBox_status_check == True:
            
            imgforspot = self.df_checker.loc[(self.df_checker['channel'] == '3')]
            ch3_xyz, ch3_xyz_3D, ch3_final_spots = self.XYZ_SPOT_COORDINATES( imgforspot, ImageForNucMask, 'Ch3')
                    
        if self.gui_params.SpotCh4CheckBox_status_check == True:
             
            imgforspot = self.df_checker.loc[(self.df_checker['channel'] == '4')]
            ch4_xyz, ch4_xyz_3D, ch4_final_spots = self.XYZ_SPOT_COORDINATES( imgforspot, ImageForNucMask, 'Ch4')
        
        if self.gui_params.SpotCh5CheckBox_status_check == True:
             
            imgforspot = self.df_checker.loc[(self.df_checker['channel'] == '5')]
            ch5_xyz, ch5_xyz_3D, ch5_final_spots = self.XYZ_SPOT_COORDINATES( imgforspot, ImageForNucMask, 'Ch5')
            
        return ch1_xyz, ch1_xyz_3D, ch1_final_spots, ch2_xyz, ch2_xyz_3D, ch2_final_spots, ch3_xyz, ch3_xyz_3D, ch3_final_spots, ch4_xyz, ch4_xyz_3D, ch4_final_spots, ch5_xyz, ch5_xyz_3D, ch5_final_spots 

    
    def XYZ_SPOT_COORDINATES(self, images_pd_df, ImageForNucMask, spot_channel):
        
        z_imglist = []
        coordinates_stack = np.ones((0,3), dtype='float')
        xyz_coordinates = []        
            
        for index, row in images_pd_df.iterrows():
            
            if row['ImageName']=="dask_array":
                im = row["Type"].compute()
            else: 
                im = mpimg.imread(row['ImageName'])
            z_imglist.append( np.asarray(im))
            _z_coordinates1 = np.asarray(row['z_slice']).astype('float')
            
            
            coordinates, final_spots = self.ImageAnalyzer.SpotDetector(im, ImageForNucMask, spot_channel)
            
            if coordinates.__len__()>0:
                _z_coordinates = np.ones((coordinates.__len__(),1), dtype='float')*_z_coordinates1
            else:
                coordinates=np.ones((0,2), dtype='float')
                _z_coordinates=np.ones((0,1), dtype='float')
            
            xyz_3d_coordinates = np.append(np.asarray(coordinates).astype('float'), _z_coordinates, 1)
            coordinates_stack = np.append(coordinates_stack, xyz_3d_coordinates,0)
            

        if z_imglist.__len__()>0:
            #print(row['ImageName'])
            image_stack = np.stack(z_imglist, axis=2)
            max_project = image_stack.max(axis=2)
            coordinates_max_project, final_spots = self.ImageAnalyzer.SpotDetector(max_project, ImageForNucMask, spot_channel)
        
            if coordinates_max_project.__len__()>0:

                coordinates_max_project_round = np.round(np.asarray(coordinates_max_project)).astype('int')
                coordinates_max_project = np.array(coordinates_max_project)[coordinates_max_project_round.min(axis=1)>=0,:].tolist()
                coordinates_max_project_round = coordinates_max_project_round[coordinates_max_project_round.min(axis=1)>=0,:]
                spots_z_slices = np.argmax(image_stack[coordinates_max_project_round[:,0],coordinates_max_project_round[:,1],:], axis=1)
                spots_z_coordinates = np.zeros((spots_z_slices.__len__(),1), dtype='float')
                
                for i in range(spots_z_slices.__len__()):

                    spots_z_coordinates[i] = np.asarray(images_pd_df.loc[images_pd_df['z_slice']== str(spots_z_slices[i]+1)]
                                                     ['z_coordinate'].iloc[0]).astype('float')
                if coordinates_max_project==[]:
                    coordinates_max_project=np.ones((0,2), dtype='float')
                xyz_coordinates = np.append(np.asarray(coordinates_max_project).astype('float'), spots_z_coordinates, 1)

        return xyz_coordinates, coordinates_stack, final_spots
    
   
    def RADIAL_DIST_CALC(self, xyz_round, spot_nuc_labels, radial_dist_df, dist_img):
        radial_dist=[]
        eps=0.000001
        for i in range(xyz_round.__len__()):
            
            sp_dist = dist_img[xyz_round[i,0], xyz_round[i,1]]
            spot_lbl =np.int(spot_nuc_labels[i])
            if spot_lbl>0:
                cell_max = radial_dist_df.loc[radial_dist_df['label']==spot_lbl]['max_intensity'].iloc[0]
                sp_radial_dist= (cell_max-sp_dist)/(cell_max-1+eps)
            else:
                sp_radial_dist = np.nan
            radial_dist.append(sp_radial_dist)
    
        return np.array(radial_dist).astype(float)

    
    def RUN_BTRACK(self,label_stack, gui_params):
        print(label_stack)
        obj_from_generator = btrack.utils.segmentation_to_objects(label_stack, properties = ('bbox','area',
                                                                                    'perimeter',
                                                                                      'major_axis_length','orientation',
                                                                                       'solidity','eccentricity'))
        # initialise a tracker session using a context manager
        with btrack.BayesianTracker() as tracker:
            tracker = btrack.BayesianTracker()
                # configure the tracker using a config file
            tracker.configure_from_file('./BayesianTracker/models/cell_config.json') #CHANGED
            tracker.max_search_radius = gui_params.NucSearchRadiusSpinbox_current_value
            print("step 1")
            # append the objects to be tracked
            tracker.append(obj_from_generator)
            print("step 2")
            # set the volume
            tracker.volume=((0, 1200), (0, 1200), (-1e5, 64.))
            
            # track them (in interactive mode)
            tracker.track_interactive(step_size=100)
            print("step 3")
            # generate hypotheses and run the global optimizer
            tracker.optimize()
            print("step 4")

#                             tracker.export(os.path.join('/data2/cell_tracking/','tracking.h5'), obj_type='obj_type_1')

#                             # get the tracks in a format for napari visualization
#                             data, properties, graph = tracker.to_napari(ndim=2)

            tracks = tracker.tracks
            print("step 5")
        tracks_pd = pd.DataFrame()
        for i in range(len(tracks)):
            tracks_pd = tracks_pd.append(pd.DataFrame(tracks[i].to_dict()))
        print(tracks_pd)
        tracks_pd = tracks_pd[tracks_pd['dummy']==False]
        return tracks_pd
    
    #removed gui_params from the args because it was not used at all
    def deepcell_tracking(self,t_stack_nuc,masks_stack):
    
        tracker = CellTracking()
        tracked_data = tracker.track(np.copy(np.expand_dims(t_stack_nuc, axis=-1)), np.copy(np.expand_dims(masks_stack, axis=-1)))

        tracks_pd = pd.DataFrame()
        for i in range(len(tracked_data['y_tracked'])):
            labeled_nuc, number_nuc = label(np.squeeze(tracked_data['y_tracked'][i]))
            dist_props = regionprops_table(labeled_nuc, properties=('label','centroid', 'bbox','area', 'perimeter', 'major_axis_length','orientation', 'solidity','eccentricity'))
            single_image_df = pd.DataFrame(dist_props)
            single_image_df['t']=i*np.ones(len(single_image_df),dtype=int)
            single_image_df.rename(columns={ "label":"ID"}, inplace = True)
            single_image_df.rename(columns={ "centroid-0":"y"}, inplace = True)
            single_image_df.rename(columns={ "centroid-1":"x"}, inplace = True)
            single_image_df['z']=i*np.zeros(len(single_image_df),dtype=int)
            single_image_df['generation']=i*np.zeros(len(single_image_df),dtype=int)
            single_image_df['parent']=single_image_df['ID']
            single_image_df['root']=single_image_df['ID']
            new_columns=['ID','t','x','y','z','parent','root','generation', 'area', 'major_axis_length', 'solidity', 'bbox-1', 'bbox-2', 'bbox-3', 'eccentricity', 'bbox-0', 'orientation', 'perimeter']
            new_df = single_image_df[new_columns]
            tracks_pd = pd.concat([tracks_pd,new_df], axis=0)
        
        return tracks_pd
    
    def mutual_information(self, hgram):
     
        pxy = hgram / float(np.sum(hgram))
        px = np.sum(pxy, axis=1) # marginal for x over y
        py = np.sum(pxy, axis=0) # marginal for y over x
        px_py = px[:, None] * py[None, :] # Broadcast to multiply marginals
        # Now we can do the calculation using the pxy, px_py 2D arrays
        nzs = pxy > 0 # Only non-zero pxy values contribute to the sum
        return np.sum(pxy[nzs] * np.log(pxy[nzs] / px_py[nzs]))
    
    def zero_pad_patch(self, input_image, desired_size=(256,256)):
    
        padded_image1 = np.zeros(desired_size, dtype=float)
        w,h= input_image.shape
        shift_w = int(desired_size[0]/2)
        shift_h = int(desired_size[1]/2)
        padded_image1[:w, :h] = input_image
        padded_image = np.roll(np.roll(padded_image1, int(shift_w-w/2),axis=0), int(shift_h-h/2), axis=1)    

        return padded_image
