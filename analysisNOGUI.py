#new
import parameters_class
# f = open("parameters_class.py", mode="r")
# print(f.readlines())
# f.close()
#old
import numpy as np
import cv2
from scipy.ndimage import label
from scipy import ndimage
from PIL import Image
from skimage.feature import peak_local_max
from skimage.draw import circle_perimeter
import math
from skimage.filters import median, gaussian, sobel
from skimage.morphology import disk, binary_closing, skeletonize, binary_opening, binary_erosion, white_tophat
from skimage.filters import threshold_li
from skimage.segmentation import watershed, find_boundaries
from skimage.exposure import equalize_adapthist
from skimage.measure import regionprops, regionprops_table
import pandas as pd
from cellpose import models, utils
import tensorflow as tf
from tensorflow.keras import backend as K
import pickle
from skimage.transform import rescale, resize, downscale_local_mean
from deepcell.applications import NuclearSegmentation
class ImageAnalyzer(object):
    
#     with open('/data2/test_images/for_training/trained_model1.pickle', 'rb') as handle:  ### Random forest model
#     with open('/data2/test_images/for_training/all_trained_models/Nearest_Neighbors.pickle', 'rb') as handle:  ### Nearest Neighbors model    
        
#     with open("/data2/test_images/for_training/all_trained_models/NeuralNet_deep.pickle", 'rb') as handle: ### neural netowork model (8,16,32,64,32,16,8)
#         model1 = pickle.load(handle)
    #edited
    #made all parameters optional, for gui init only use gui params, for nogui init use input_params dict
    #analysisgui and inoutresourcegui are instances of the classes that are instantiated in hitips/batchanalyzer
    def __init__(self, is_gui, input_params=None, analysisgui=None, inout_resource_gui=None):
        #all from the input excel file->dict
        if is_gui == False:
            self.input_params = input_params
            self.gui_params = parameters_class.Gui_Params(is_gui,input_params = input_params)
            self.spot_params_dict =self.INITIALIZE_SPOT_ANALYSIS_PARAMS()
        else:
            self.AnalysisGui = analysisgui
            self.inout_resource_gui = inout_resource_gui
            self.gui_params = parameters_class.Gui_Params(is_gui,analysisgui=analysisgui, inout_resource_gui=inout_resource_gui)
            self.spot_params_dict =self.INITIALIZE_SPOT_ANALYSIS_PARAMS()
        
    
    #edited
    #the only major edits here are how the params are called- all come from the gui_params class instantiation
    def nuclei_segmenter(self, input_img, pixpermic=None):
        
        #self.AnalysisGui.nuclei_detection_method.currentText()444
        if self.gui_params.NucDetectMethod == "Int.-based":
            
            first_thresh = self.gui_params.NucSeparationSlider*2.55 #separationslider
            second_thresh = 255-(self.gui_params.NucDetectionSlider*2.55) #detectionslider
            Cell_Size = self.gui_params.NucleiAreaSlider
            
            boundary, mask = self.segmenter_function(input_img, cell_size=Cell_Size, 
                                                     first_threshold=first_thresh, second_threshold=second_thresh)
            
        if self.gui_params.NucDetectMethod == "Marker Controlled":
            pixpermic=1
          
            Cell_Size = self.gui_params.NucleiAreaSlider
            max_range = np.sqrt(Cell_Size/3.14)*2/float(pixpermic)
            nuc_detect_sldr = self.gui_params.NucDetectionSlider
            first_thresh = np.ceil((1-(nuc_detect_sldr/100))*max_range).astype(int)
            
            second_thresh = self.gui_params.NucSeparationSlider
            
            print(Cell_Size, first_thresh, second_thresh)
            boundary, mask = self.watershed_scikit(input_img, cell_size=Cell_Size, 
                                                     first_threshold=first_thresh, second_threshold=second_thresh)

        if self.gui_params.NucDetectMethod == "CellPose-GPU":
            pixpermic = 0.1 #commented out earlier, why?
            Cell_Size = self.gui_params.NucleiAreaSlider
            cell_diameter = np.sqrt(Cell_Size/(float(pixpermic)*float(pixpermic)))*2/3.14
            
            boundary, mask = self.cellpose_segmenter(input_img, use_GPU=1, cell_dia=cell_diameter)
            
        if self.gui_params.NucDetectMethod == "CellPose-CPU":
            
            Cell_Size = self.gui_params.NucleiAreaSlider
            cell_diameter = np.sqrt(Cell_Size/(float(pixpermic)*float(pixpermic)))*2/3.14
            
            boundary, mask = self.cellpose_segmenter(input_img, use_GPU=0, cell_dia=cell_diameter)
            
        if self.gui_params.NucDetectMethod == "CellPose-Cyto":
            
            Cell_Size = self.gui_params.NucleiAreaSlider
            cell_diameter = np.sqrt(Cell_Size/(float(pixpermic)*float(pixpermic)))*2/3.14
            
            boundary, mask = self.cellpose_segmenter(input_img, use_GPU=1, cell_dia=cell_diameter)
                
        if self.gui_params.NucDetectMethod == "DeepCell":
            
            Cell_Size = self.gui_params.NucleiAreaSlider
            cell_diameter = Cell_Size/100
            
            boundary, mask = self.deepcell_segmenter(input_img, cell_dia=cell_diameter)
        
#         kernel = np.ones((7,7), np.uint8)
#         first_pass = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel)
#         boundary = find_boundaries(first_pass, connectivity=1, mode='thick', background=0)
#         mask= (255*first_pass).astype('uint8')
#         boundary= (255*boundary).astype('uint8')
        #ask for an output directory to return the file
        return boundary, mask

#edited
#I moved the deepcell import to the top
    def deepcell_segmenter(self, input_img, cell_dia=None):
        
        app = NuclearSegmentation()
        im = np.expand_dims(input_img, axis=-1)
        im = np.expand_dims(im, axis=0)
        
        masks1 = app.predict(im, image_mpp=cell_dia)
        masks = np.squeeze(masks1)
        boundary = find_boundaries(masks, connectivity=1, mode='thick', background=0)
        boundary_img = (255*boundary).astype('uint8')
        resized_bound = boundary_img
        filled1 = ndimage.binary_fill_holes(resized_bound)
        mask= (255*filled1).astype('uint8')-resized_bound
        boundary= resized_bound.astype('uint8')

        return boundary, mask
#edited
#only edits to gui_params
    def cellpose_segmenter(self, input_img, use_GPU, cell_dia=None):
        
        if self.gui_params.NucRemoveBoundaryCheckBox == True: #remove boundary
            img_uint8 = input_img
        else:
            img_uint8 = cv2.copyMakeBorder(input_img,5,5,5,5,cv2.BORDER_CONSTANT,value=0)

        if self.gui_params.NucDetectMethod == "CellPose-Cyto":
            model = models.Cellpose(gpu=use_GPU, model_type='cyto')
        else:
            model = models.Cellpose(gpu=use_GPU, model_type='nuclei')
        masks, flows, styles, diams = model.eval(img_uint8, diameter=cell_dia, flow_threshold=None)
                
        boundary = find_boundaries(masks, connectivity=1, mode='thick', background=0)

        if self.gui_params.NucRemoveBoundaryCheckBox == True:

            boundary_img = (255*boundary).astype('uint8')
            filled1 = ndimage.binary_fill_holes(boundary_img)
            mask1= (255*filled1).astype('uint8')-boundary_img
            kernel = np.ones((3,3), np.uint8)
            mask = cv2.morphologyEx(mask1.astype('uint8'), cv2.MORPH_OPEN, kernel)
#             kernel = np.ones((1,1),np.uint8)
#             erosion = cv2.erode(mask,kernel,iterations = 1)

            
            boundary_img = find_boundaries(mask, connectivity=1, mode='thick', background=0)
            resized_bound = cv2.resize((255*boundary_img).astype('uint8'),(input_img.shape[1],input_img.shape[0]))
            filled1 = ndimage.binary_fill_holes(resized_bound)
        else:
            boundary_img = (255*boundary[3:boundary.shape[0]-3,3:boundary.shape[1]-3]).astype('uint8')
            resized_bound = cv2.resize(boundary_img,(input_img.shape[1],input_img.shape[0]))
            filled1 = ndimage.binary_fill_holes(resized_bound)
        
        mask1= (255*filled1).astype('uint8')-resized_bound
        kernel = np.ones((11,11), np.uint8)
        mask = cv2.morphologyEx(mask1, cv2.MORPH_OPEN, kernel)
        mask = 255*(mask==255).astype('uint8')
        mask[mask>0]=255
        boundary= resized_bound.astype('uint8')

        return boundary, mask

    def segmenter_function(self, input_img, cell_size=None, first_threshold=None, second_threshold=None):
    
        img_uint8 = cv2.copyMakeBorder(input_img,5,5,5,5,cv2.BORDER_CONSTANT,value=0)
        
        
        ## First blurring round
        if (cell_size %2)==0:
            cell_size = cell_size + 1
        median_img = cv2.medianBlur(img_uint8,cell_size)
        gaussian_blurred = cv2.GaussianBlur(median_img,(5,5),0)
        ## Threhsolding and Binarizing
        ret, thresh = cv2.threshold(gaussian_blurred,0,255,cv2.THRESH_BINARY_INV+cv2.THRESH_OTSU)
        bin_img = (1-thresh/255).astype('bool')
        ## Binary image filling
        filled = ndimage.binary_fill_holes(bin_img)
        filled_int= (filled*255).astype('uint8')
        ## Gray2RGB to feed the watershed algorithm
        img_rgb  = cv2.cvtColor(img_uint8,cv2.COLOR_GRAY2RGB)
        boundary = img_rgb
        boundary = boundary - img_rgb
        ## Distance trasform and thresholing to set the watershed markers
        dt = cv2.distanceTransform(filled.astype(np.uint8), 2, 3)
        dt = ((dt - dt.min()) / (dt.max() - dt.min()) * 255).astype(np.uint8)
        _, dt = cv2.threshold(dt, first_threshold, 255, cv2.THRESH_BINARY)
        lbl, ncc = label(dt)
        lbl = lbl * (255 / (ncc + 1))
        lbl = lbl.astype(np.int32)
        ## First round of Watershed transform
        cv2.watershed(img_rgb, lbl)
        ## Correcting image boundaries
        boundary[lbl == -1] = [255,255,255]
        boundary[0,:] = 0
        boundary[-1,:] = 0
        boundary[:,0] = 0
        boundary[:, -1] = 0
        b_gray = cv2.cvtColor(boundary,cv2.COLOR_BGR2GRAY)
        diff = filled_int-b_gray

        kernel = np.ones((11,11), np.uint8)
        first_pass = cv2.morphologyEx(diff, cv2.MORPH_OPEN, kernel)

        ## Second round of marker generation and watershed 
        kernel = np.ones((5,5),np.uint8)
        aa = first_pass.astype('uint8')
        erosion = cv2.erode(aa,kernel,iterations = 1)
        kernel = np.ones((11,11), np.uint8)
        opening = cv2.morphologyEx(erosion, cv2.MORPH_OPEN, kernel)
        blur = cv2.GaussianBlur(opening,(11,11),50)
        ret2, thresh2 = cv2.threshold(blur,0,255,cv2.THRESH_BINARY_INV+cv2.THRESH_OTSU)
        dt = cv2.distanceTransform(255-thresh2, 2, 3)
        dt = ((dt - dt.min()) / (dt.max() - dt.min()) * 255).astype(np.uint8)
        _, dt = cv2.threshold(dt, second_threshold, 255, cv2.THRESH_BINARY)
        lbl, ncc = label(dt)
        lbl = lbl * (255 / (ncc + 1))
        lbl = lbl.astype(np.int32)
        cv2.watershed(img_rgb, lbl)
        ########
        boundary = img_rgb
        boundary = boundary - img_rgb

        boundary[lbl == -1] = [255,255,255]
        boundary_img = boundary[3:boundary.shape[0]-3,3:boundary.shape[1]-3]
        bound_gray = cv2.cvtColor(boundary_img,cv2.COLOR_BGR2GRAY)
        resized_bound = cv2.resize(bound_gray,(input_img.shape[1],input_img.shape[0]))

        kernel = np.ones((3,3),np.uint8)
        boundary = cv2.dilate(resized_bound,kernel,iterations = 1)
        filled1 = ndimage.binary_fill_holes(boundary)
        fin= 255*filled1-boundary
        mask = ndimage.binary_fill_holes(fin)
        mask = (255*mask).astype(np.uint8)

        return boundary, mask

    def watershed_scikit(self, input_img, cell_size=None, first_threshold=None, second_threshold=None):
        
        img_uint8 = cv2.copyMakeBorder(input_img,5,5,5,5,cv2.BORDER_CONSTANT,value=0)
        
        med_scikit = median(img_uint8, disk(1))
        thresh = threshold_li(med_scikit)
        binary = med_scikit > thresh
        filled = ndimage.binary_fill_holes(binary)
        filled_blurred = gaussian(filled, 1)
        filled_int= (filled_blurred*255).astype('uint8')
        
# #         edge_sobel = sobel(img_uint8)
# #         enhanced = 50*edge_sobel/edge_sobel.max() + img_uint8
# #         enhanced.astype('uint8')
# #         med_scikit = median(img_uint8, disk(5))
        thresh = threshold_li(filled_int)
        binary = filled_int > thresh
        filled = ndimage.binary_fill_holes(binary)
        filled_int = binary_opening(filled, disk(5))
        filled_int = ndimage.binary_fill_holes(filled_int)
#         filled_blurred = gaussian(openeed, 3)
        
#         thresh = threshold_li(filled_int)
#         binary = filled_int > thresh
        #binary = binary_erosion(filled_int, disk(5))
        distance = ndimage.distance_transform_edt(filled_int)
        binary1 = distance > first_threshold
        distance1 = ndimage.distance_transform_edt(binary1)
        binary2 = distance1 > second_threshold

        labeled_spots, num_features = label(binary2)
        spot_labels = np.unique(labeled_spots)    

        spot_locations = ndimage.measurements.center_of_mass(binary2, labeled_spots, spot_labels[spot_labels>0])

        mask = np.zeros(distance.shape, dtype=bool)
        if spot_locations:
            mask[np.ceil(np.array(spot_locations)[:,0]).astype(int), np.ceil(np.array(spot_locations)[:,1]).astype(int)] = True
        markers, _ = ndimage.label(mask)
        labels = watershed(-distance, markers, mask=binary, compactness=0.5, watershed_line=True)
        boundary = find_boundaries(labels, connectivity=1, mode='thick', background=0)
        boundary_img = (255*boundary[3:boundary.shape[0]-3,3:boundary.shape[1]-3]).astype('uint8')
        resized_bound = cv2.resize(boundary_img,(input_img.shape[1],input_img.shape[0]))
        filled1 = ndimage.binary_fill_holes(resized_bound)
        
        mask= (255*filled1).astype('uint8')-resized_bound
        boundary= resized_bound.astype('uint8')
        
        return boundary, mask
    def max_z_project(self, image_stack):
        
        z_imglist=[]
        
        for index, row in image_stack.iterrows():
            if row['ImageName']=="dask_array":
                im = row["Type"].compute()
            else: 
                im = Image.open(row['ImageName'])
            z_imglist.append( np.asarray(im))
        z_stack = np.stack(z_imglist, axis=2)
        max_project = z_stack.max(axis=2)
        
        return max_project
#edited
#changed the way it checks the methods to comparing to strings
#ex. used to be if str(params_to_pass[0]) == '0': and I changed to if str(params_to_pass[0]) == "LOG":
#this is because with the GUI, lists of the possible param values are created with the widget, but in headless
#we need to compare directly
    def SpotDetector(self, input_image_raw, nuclei_image, spot_channel):
        #no need to update spot analysis params in nogui ver, removed

        if self.gui_params.NucRemoveBoundaryCheckBox=='All':
            
            params_to_pass= self.spot_params_dict['Ch1']
        else:
            params_to_pass= self.spot_params_dict[spot_channel]
        
        
        uint8_max_val = 255
#         w,h = input_image.shape
#         dim = (2*w,2*h)
#         input_image1 = cv2.resize(input_image, dim, interpolation = cv2.INTER_AREA)
        
#         noise = np.random.normal(1*input_image_raw.mean(), 1*input_image_raw.std(), input_image_raw.shape)
#         input_image_raw = input_image_raw +noise
        
        input_image = cv2.normalize(input_image_raw, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)

        input_image1 = input_image
        ## First blurring round
        median_img = cv2.medianBlur(nuclei_image,11)
        ## Threhsolding and Binarizing
        ret, thresh = cv2.threshold(median_img,0,uint8_max_val,cv2.THRESH_BINARY_INV+cv2.THRESH_OTSU)
        bin_img = (1-thresh/uint8_max_val).astype('bool')
        ## Binary image filling
        filled = ndimage.binary_fill_holes(bin_img)
        struct = ndimage.generate_binary_structure(2, 2)
        filled = ndimage.binary_dilation(filled, structure=struct).astype(filled.dtype)
        filled = ndimage.binary_dilation(filled, structure=struct).astype(filled.dtype)
        #### this part is for removing bright junk in the image################
        labeled_nuc, num_features_nuc = label(filled)
        props = regionprops_table(labeled_nuc, input_image, properties=('label', 'area', 'max_intensity', 'mean_intensity'))
        props_df = pd.DataFrame(props)
        mean_intensity_ave=props_df['mean_intensity'].mean()
        max_intensity_max=props_df['max_intensity'].max()
        for ind, row in props_df.iterrows():
            
            if row['mean_intensity'] > 2*mean_intensity_ave:
                
                input_image1[labeled_nuc==row['label']]=0
        input_image1[input_image1>max_intensity_max]=0 
        if self.gui_params.SpotPerChSpinBox>1:
            input_image1 = rescale(input_image1.copy(), self.gui_params.SpotPerChSpinBox, anti_aliasing=False, preserve_range=True)

        input_image1 = cv2.normalize(input_image1, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)    
########################################
        sig=int(params_to_pass[3])
        if str(params_to_pass[0]) == "LOG": #0
            log_result = ndimage.gaussian_laplace(input_image1, sigma=sig)
            if self.gui_params.SpotPerChSpinBox>1:
                log_result =  resize(log_result, input_image_raw.shape, anti_aliasing=False, preserve_range=True)

            if str(params_to_pass[1]) == "Auto": 

                ret_log, thresh_log = cv2.threshold(log_result.astype("uint8"),0,255,cv2.THRESH_BINARY_INV+cv2.THRESH_OTSU)
                bin_img_log = (1-thresh_log/255).astype('bool')
                
                if self.gui_params.SpotPerChSpinBox>1:
                    struct = ndimage.generate_binary_structure(2, 2)
                    bin_img_log = ndimage.binary_dilation(bin_img_log, structure=struct).astype(filled.dtype)

            if str(params_to_pass[1]) == "Manual": 
                
                manual_threshold = np.ceil(params_to_pass[2]*2.55).astype(int)
                thresh_log = log_result.astype("uint8") > manual_threshold
                bin_img_log = thresh_log
                if self.gui_params.SpotPerChSpinBox>1:
                    struct = ndimage.generate_binary_structure(2, 2)
                    bin_img_log = ndimage.binary_dilation(bin_img_log, structure=struct).astype(filled.dtype)
            
            spots_img_log = (bin_img_log*255).astype('uint8')
            kernel = np.ones((3,3), np.uint8)
            spot_openned_log = cv2.morphologyEx(spots_img_log, cv2.MORPH_OPEN, kernel)
            final_spots = np.multiply(spot_openned_log,filled)
            spots_df, bin_img_log, labeled_spots = self.spots_information(final_spots, input_image_raw)   
            
        if str(params_to_pass[0]) == "Laplacian of Gaussian": #1
            
            result_gaussian = ndimage.gaussian_filter(input_image1, sigma=sig)
            
            if str(params_to_pass[1]) == "Auto":
                
                ret_log, thresh_log = cv2.threshold(result_gaussian,0,255,cv2.THRESH_BINARY_INV+cv2.THRESH_OTSU)
                bin_img_g = (1-thresh_log/255).astype('bool')
                
            if str(params_to_pass[1]) == "Manual":
                
                manual_threshold = np.ceil(params_to_pass[2]*2.55).astype(int)
                
                thresh_log = result_gaussian > manual_threshold
                bin_img_g = thresh_log
#CHANGED
            spots_img_g = ((bin_img_g>0)*255).astype('uint8') 
            kernel = np.ones((3,3), np.uint8)
            spot_openned_g = cv2.morphologyEx(spots_img_g, cv2.MORPH_OPEN, kernel)
            final_spots = np.multiply(spot_openned_g,filled)
            spots_df, bin_img_g, labeled_spots = self.spots_information(final_spots, input_image_raw)
#             final_spots = ((bin_img_g>0)*255).astype('uint8') 
        
        if str(params_to_pass[0]) == "IntensityThreshold":
            
            if str(params_to_pass[1]) == "Auto":
                
                ret_log, thresh_log = cv2.threshold(input_image,0,255,cv2.THRESH_BINARY_INV+cv2.THRESH_OTSU)
                bin_img_g = (1-thresh_log/255).astype('bool')
                
            if str(params_to_pass[1]) == "Manual":
                
                manual_threshold = np.ceil(params_to_pass[2]*2.55).astype(int)
                
                thresh_log = input_image > manual_threshold
                bin_img_g = thresh_log
            
            spots_img_g = (bin_img_g*255).astype('uint8')
            kernel = np.ones((3,3), np.uint8)
            spot_openned_g = cv2.morphologyEx(spots_img_g, cv2.MORPH_OPEN, kernel)
            final_spots = np.multiply(spot_openned_g,filled)
        ### center of mass calculation
        if str(self.gui_params.SpotLocationCbox_currentText) == 'CenterOfMass':
            
#             labeled_spots, num_features = label(final_spots)
#             spot_labels = np.unique(labeled_spots)
            
#             bin_img = (final_spots/uint8_max_val).astype('bool')
#             ## Binary image filling
#             masked_spots = np.multiply(input_image,bin_img)
            
#             spot_locations = ndimage.measurements.center_of_mass(masked_spots, labeled_spots, spot_labels[spot_labels>0])
            spot_locations = list(spots_df['center_of_mass_coords'].to_numpy())
            ###### Brightest spot calculation
        if str(self.gui_params.SpotLocationCbox_currentText) == "MaxIntensity":
            
#             labeled_spots, num_features = label(final_spots)
#             spot_labels = np.unique(labeled_spots)
#             bin_img = (final_spots/uint8_max_val).astype('bool')
#             masked_spots = np.multiply(input_image,bin_img)
#             spot_locations = peak_local_max(masked_spots, labels=labeled_spots, num_peaks_per_label=1)
            spot_locations = list(spots_df['max_intensity_coords'].to_numpy())
            ##### Centroid calculation
        if str(self.gui_params.SpotLocationCbox_currentText) == "Centroid":
            
            labeled_spots, num_features = label(final_spots)
            spot_labels = np.unique(labeled_spots)
            
            spot_locations = ndimage.measurements.center_of_mass(final_spots, labeled_spots, spot_labels[spot_labels>0])
                        
        return spot_locations, spots_df
    
    def spots_information(self, bin_img_log, max_project ):
        
        labeled_spots, num_features = label(bin_img_log)
        if num_features>0:
            props = regionprops_table(labeled_spots, max_project,  properties=( 'label',  'area', 'solidity','coords'))
            props_df = pd.DataFrame(props)
#             props_df = props_df.loc[(props_df["area"]<25) & (props_df["area"]>5) & (props_df["solidity"]>0.98)]
            props_df = props_df.loc[props_df['area']>4]
            new_bin_img_log = np.zeros(bin_img_log.shape)
            for spot_no in range(len(props_df)):
                spot_patch=max_project[props_df['coords'].iloc[spot_no][:,0].min():props_df['coords'].iloc[spot_no][:,0].max(),
                                       props_df['coords'].iloc[spot_no][:,1].min():props_df['coords'].iloc[spot_no][:,1].max()]
                if (spot_patch.shape[0]>0)&(spot_patch.shape[1]>0):
                    spot_fwhm = np.multiply(spot_patch>spot_patch.max()/2,spot_patch)
                    new_bin_img_log[props_df['coords'].iloc[spot_no][:,0].min():props_df['coords'].iloc[spot_no][:,0].max(),
                                    props_df['coords'].iloc[spot_no][:,1].min():props_df['coords'].iloc[spot_no][:,1].max()] = spot_fwhm>0
            
            kernel = np.ones((3,3), np.uint8)
            new_bin_img_log = cv2.morphologyEx(new_bin_img_log, cv2.MORPH_CLOSE, kernel)
            
            labeled_spots, num_features = label(new_bin_img_log)
            if num_features>0:

                props = regionprops_table(labeled_spots, max_project, properties=( 'area', 'max_intensity', 'min_intensity', 'mean_intensity',
                                                        'perimeter', 'solidity','coords', 'bbox'))
                props_df = pd.DataFrame(props)
                masked_spots = np.multiply(max_project,new_bin_img_log)
                spot_labels = np.unique(labeled_spots)            
                props_df['center_of_mass_coords'] = ndimage.measurements.center_of_mass(masked_spots, labeled_spots, spot_labels[spot_labels>0])
                props_df['max_intensity_coords'] = props_df.apply(lambda row : tuple(row['coords'][max_project[row['coords'][:,0], 
                                                                                             row['coords'][:,1]].argmax(),:]), axis = 1)
        
                nospot_max_project = np.multiply(max_project,~(new_bin_img_log>0))
                r=4
                for i in range(len(props_df)):

                    y0 = max(props_df.loc[i, "bbox-0"]-r,0)
                    y1 = min(props_df.loc[i, "bbox-2"]+r, max_project.shape[0])
                    x0 = max(props_df.loc[i, "bbox-1"]-r,0)
                    x1 = min(props_df.loc[i, "bbox-3"]+r, max_project.shape[1])
                    spot_area_patch = nospot_max_project[y0:y1,x0:x1]
                    nonzero_patch = spot_area_patch[np.nonzero(spot_area_patch)]

                    props_df.loc[i,"spot_area_mean"]= nonzero_patch.mean()
                    props_df.loc[i,"spot_area_std"]= nonzero_patch.std()
                    props_df.loc[i,"spot_area_median"]=np.median(nonzero_patch)
                    
                    props_df.loc[i,"spot_max_to_min"]=props_df.loc[i, "max_intensity"]/props_df.loc[i, "min_intensity"]
                    props_df.loc[i,"spot_to_area_mean"]=props_df.loc[i, "mean_intensity"]/props_df.loc[i, "spot_area_mean"]
                    props_df.loc[i,"spot_max_to_area_mean"]=props_df.loc[i, "max_intensity"]/props_df.loc[i, "spot_area_mean"]
                # props_df = props_df.loc[(props_df["spot_to_area_mean"]>1.2)& (props_df["spot_max_to_area_mean"]>1.2)]
                props_df = props_df.loc[(props_df["spot_max_to_area_mean"]>1.2)]
                props_df = props_df.reset_index(drop=True)
#                 print("before filtering...")
#                 print(props_df)
#                 true_ind = self.model1.predict(props_df[['max_intensity', 'min_intensity', 'mean_intensity','spot_area_mean', 'spot_area_std', 'spot_area_median']].to_numpy())
#                 props_df = props_df.loc[np.where(true_ind)[0]].reset_index(drop=True)
#                 print("after filtering...")
#                 print(props_df)
        else:
            
            props_df = pd.DataFrame(columns=['area', 'max_intensity', 'min_intensity', 'mean_intensity', 'center_of_mass_coords',
                                             'max_intensity_coords', 'perimeter', 'solidity','coords', 'bbox', 'spot_area_mean',
                                             'spot_area_std', 'spot_area_median'])
            new_bin_img_log = np.zeros(bin_img_log.shape)
        
        return props_df,new_bin_img_log, labeled_spots
        
    def COORDINATES_TO_CIRCLE(self, coordinates,ImageForSpots, circ_radius =5):
        
        circles = np.zeros((ImageForSpots.shape), dtype=np.uint8)
        
        if coordinates.any():
            
            for center_y, center_x in zip(coordinates[:,0], coordinates[:,1]):
                    circy, circx = circle_perimeter(center_y, center_x, circ_radius, shape=ImageForSpots.shape)
                    circles[circy, circx] = 255

        return circles
    
    def SPOTS_TO_BOUNDARY(self, final_spots):
        
        labeled_spots, num_features = label(final_spots)
        boundary = find_boundaries(labeled_spots, connectivity=1, mode='thick', background=0)
        spot_boundary = (255*boundary).astype('uint8')
        
        return spot_boundary
    
    #EDITED
    #changed to grab from gui_params again
    def INITIALIZE_SPOT_ANALYSIS_PARAMS(self):

        self.spot_params_dict={
            
            "Ch1": np.array([self.gui_params.spotanalysismethod,self.gui_params.thresholdmethod,
                    self.gui_params.ThresholdSlider, self.gui_params.SensitivitySpinBox,
                            self.gui_params.SpotsPerChSpinBox]),
            "Ch2": np.array([self.gui_params.spotanalysismethod,self.gui_params.thresholdmethod,
                    self.gui_params.ThresholdSlider, self.gui_params.SensitivitySpinBox,
                            self.gui_params.SpotsPerChSpinBox]),
            "Ch3": np.array([self.gui_params.spotanalysismethod,self.gui_params.thresholdmethod,
                    self.gui_params.ThresholdSlider, self.gui_params.SensitivitySpinBox,
                            self.gui_params.SpotsPerChSpinBox]),
            "Ch4": np.array([self.gui_params.spotanalysismethod,self.gui_params.thresholdmethod,
                    self.gui_params.ThresholdSlider, self.gui_params.SensitivitySpinBox,
                            self.gui_params.SpotsPerChSpinBox]),
            "Ch5": np.array([self.gui_params.spotanalysismethod,self.gui_params.thresholdmethod,
                    self.gui_params.ThresholdSlider, self.gui_params.SensitivitySpinBox,
                            self.gui_params.SpotsPerChSpinBox])
            }
        return self.spot_params_dict
    
    def UPDATE_SPOT_ANALYSIS_PARAMS(self):
        
        if self.AnalysisGui.spotchannelselect.currentText()=='All':
            
            self.spot_params_dict={
            
            "Ch1": np.array([self.AnalysisGui.spotanalysismethod.currentIndex(),self.AnalysisGui.thresholdmethod.currentIndex(),
                    self.AnalysisGui.ThresholdSlider.value(), self.AnalysisGui.SensitivitySpinBox.value(),
                            self.AnalysisGui.SpotPerChSpinBox.value()],dtype=int),
            "Ch2": np.array([self.AnalysisGui.spotanalysismethod.currentIndex(),self.AnalysisGui.thresholdmethod.currentIndex(),
                    self.AnalysisGui.ThresholdSlider.value(), self.AnalysisGui.SensitivitySpinBox.value(),
                            self.AnalysisGui.SpotPerChSpinBox.value()],dtype=int),
            "Ch3": np.array([self.AnalysisGui.spotanalysismethod.currentIndex(),self.AnalysisGui.thresholdmethod.currentIndex(),
                    self.AnalysisGui.ThresholdSlider.value(), self.AnalysisGui.SensitivitySpinBox.value(),
                            self.AnalysisGui.SpotPerChSpinBox.value()],dtype=int),
            "Ch4": np.array([self.AnalysisGui.spotanalysismethod.currentIndex(),self.AnalysisGui.thresholdmethod.currentIndex(),
                    self.AnalysisGui.ThresholdSlider.value(), self.AnalysisGui.SensitivitySpinBox.value(),
                            self.AnalysisGui.SpotPerChSpinBox.value()],dtype=int),
            "Ch5": np.array([self.AnalysisGui.spotanalysismethod.currentIndex(),self.AnalysisGui.thresholdmethod.currentIndex(),
                    self.AnalysisGui.ThresholdSlider.value(), self.AnalysisGui.SensitivitySpinBox.value(),
                            self.AnalysisGui.SpotPerChSpinBox.value()],dtype=int)
            }
        
        if self.AnalysisGui.spotchannelselect.currentText()=='Ch1':
            
            self.spot_params_dict["Ch1"] =    np.array([self.AnalysisGui.spotanalysismethod.currentIndex(),self.AnalysisGui.thresholdmethod.currentIndex(),
                    self.AnalysisGui.ThresholdSlider.value(), self.AnalysisGui.SensitivitySpinBox.value(),
                            self.AnalysisGui.SpotPerChSpinBox.value()],dtype=int)
            
        if self.AnalysisGui.spotchannelselect.currentText()=='Ch2':
            
            self.spot_params_dict["Ch2"] =    np.array([self.AnalysisGui.spotanalysismethod.currentIndex(),self.AnalysisGui.thresholdmethod.currentIndex(),
                    self.AnalysisGui.ThresholdSlider.value(), self.AnalysisGui.SensitivitySpinBox.value(),
                            self.AnalysisGui.SpotPerChSpinBox.value()],dtype=int)
        
        if self.AnalysisGui.spotchannelselect.currentText()=='Ch3':
            
            self.spot_params_dict["Ch3"] =    np.array([self.AnalysisGui.spotanalysismethod.currentIndex(),self.AnalysisGui.thresholdmethod.currentIndex(),
                    self.AnalysisGui.ThresholdSlider.value(), self.AnalysisGui.SensitivitySpinBox.value(),
                            self.AnalysisGui.SpotPerChSpinBox.value()],dtype=int)
            
        if self.AnalysisGui.spotchannelselect.currentText()=='Ch4':
            
            self.spot_params_dict["Ch4"] =    np.array([self.AnalysisGui.spotanalysismethod.currentIndex(),self.AnalysisGui.thresholdmethod.currentIndex(),
                    self.AnalysisGui.ThresholdSlider.value(), self.AnalysisGui.SensitivitySpinBox.value(),
                            self.AnalysisGui.SpotPerChSpinBox.value()],dtype=int)
            
        if self.AnalysisGui.spotchannelselect.currentText()=='Ch5':
            
            self.spot_params_dict["Ch5"] =    np.array([self.AnalysisGui.spotanalysismethod.currentIndex(),self.AnalysisGui.thresholdmethod.currentIndex(),
                    self.AnalysisGui.ThresholdSlider.value(), self.AnalysisGui.SensitivitySpinBox.value(),
                            self.AnalysisGui.SpotPerChSpinBox.value()],dtype=int)
    
    def UPDATE_SPOT_ANALYSIS_GUI_PARAMS(self):
        
        if self.AnalysisGui.spotchannelselect.currentText()=='Ch1':
            
            self.AnalysisGui.spotanalysismethod.setCurrentIndex(np.array(self.spot_params_dict["Ch1"][0]).astype(int))
            self.AnalysisGui.thresholdmethod.setCurrentIndex(np.array(self.spot_params_dict["Ch1"][1]).astype(int))
            self.AnalysisGui.ThresholdSlider.setValue(np.array(self.spot_params_dict["Ch1"][2]).astype(int))
            self.AnalysisGui.SensitivitySpinBox.setValue(np.array(self.spot_params_dict["Ch1"][3]).astype(int))
            self.AnalysisGui.SpotPerChSpinBox.setValue(np.array(self.spot_params_dict["Ch1"][4]).astype(int))
        
        if self.AnalysisGui.spotchannelselect.currentText()=='Ch2':
            
            self.AnalysisGui.spotanalysismethod.setCurrentIndex(np.array(self.spot_params_dict["Ch2"][0]).astype(int))
            self.AnalysisGui.thresholdmethod.setCurrentIndex(np.array(self.spot_params_dict["Ch2"][1]).astype(int))
            self.AnalysisGui.ThresholdSlider.setValue(np.array(self.spot_params_dict["Ch2"][2]).astype(int))
            self.AnalysisGui.SensitivitySpinBox.setValue(np.array(self.spot_params_dict["Ch2"][3]).astype(int))
            self.AnalysisGui.SpotPerChSpinBox.setValue(np.array(self.spot_params_dict["Ch2"][4]).astype(int))
        
        if self.AnalysisGui.spotchannelselect.currentText()=='Ch3':
            
            self.AnalysisGui.spotanalysismethod.setCurrentIndex(np.array(self.spot_params_dict["Ch3"][0]).astype(int))
            self.AnalysisGui.thresholdmethod.setCurrentIndex(np.array(self.spot_params_dict["Ch3"][1]).astype(int))
            self.AnalysisGui.ThresholdSlider.setValue(np.array(self.spot_params_dict["Ch3"][2]).astype(int))
            self.AnalysisGui.SensitivitySpinBox.setValue(np.array(self.spot_params_dict["Ch3"][3]).astype(int))
            self.AnalysisGui.SpotPerChSpinBox.setValue(np.array(self.spot_params_dict["Ch3"][4]).astype(int))
        
        if self.AnalysisGui.spotchannelselect.currentText()=='Ch4':
            
            self.AnalysisGui.spotanalysismethod.setCurrentIndex(np.array(self.spot_params_dict["Ch4"][0]).astype(int))
            self.AnalysisGui.thresholdmethod.setCurrentIndex(np.array(self.spot_params_dict["Ch4"][1]).astype(int))
            self.AnalysisGui.ThresholdSlider.setValue(np.array(self.spot_params_dict["Ch4"][2]).astype(int))
            self.AnalysisGui.SensitivitySpinBox.setValue(np.array(self.spot_params_dict["Ch4"][3]).astype(int))
            self.AnalysisGui.SpotPerChSpinBox.setValue(np.array(self.spot_params_dict["Ch4"][4]).astype(int))
        
        if self.AnalysisGui.spotchannelselect.currentText()=='Ch5':
            
            self.AnalysisGui.spotanalysismethod.setCurrentIndex(np.array(self.spot_params_dict["Ch5"][0]).astype(int))
            self.AnalysisGui.thresholdmethod.setCurrentIndex(np.array(self.spot_params_dict["Ch5"][1]).astype(int))
            self.AnalysisGui.ThresholdSlider.setValue(np.array(self.spot_params_dict["Ch5"][2]).astype(int))
            self.AnalysisGui.SensitivitySpinBox.setValue(np.array(self.spot_params_dict["Ch5"][3]).astype(int))
            self.AnalysisGui.SpotPerChSpinBox.setValue(np.array(self.spot_params_dict["Ch5"][4]).astype(int))
        
        
