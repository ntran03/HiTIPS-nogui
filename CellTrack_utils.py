from xml.dom import minidom
import os
import pandas as pd
import numpy as np
import imageio
import matplotlib.pyplot as plt
import cv2
import tensorflow as tf
from tensorflow.keras import backend as K
import matplotlib.image as mpimg
from scipy import ndimage, misc
from PIL import Image
from skimage.draw import circle_perimeter

# parse an xml file by name (reads the metadata file)
# data_path = '/data2/test_images/Christopher/201211-MYH9-SIR-MS2-9hr_20201211_182502/AssayPlate_MatricalBioscience_MGB096-1-2-LG-L'
class data_handler(object):
    def init(self):

        print('Initialized')
    #needs both the path to the metadata file and the name of the file
    def return_metadata_df(self, data_path, metafilename):
        metadatafilename = os.path.join(data_path, metafilename)
        print(metadatafilename)
        mydoc = minidom.parse(metadatafilename)
        PATH_TO_FILES = os.path.split(metadatafilename)[0]
        items = mydoc.getElementsByTagName('bts:MeasurementRecord')

        metadatafilename_mrf = os.path.join(data_path,'MeasurementDetail.mrf')
        mydoc_mrf = minidom.parse(metadatafilename_mrf)
        PATH_TO_FILES = os.path.split(metadatafilename_mrf)[0]
        items_mrf = mydoc_mrf.getElementsByTagName('bts:MeasurementChannel')
        df_cols = ["ImageName", "column", "row", "time_point", "field_index", "z_slice", "channel", 
                   "x_coordinates", "y_coordinates","z_coordinate", "action_index", "action", "Type", "Time", "PixPerMic"]
        rows = []

        for i in range(items.length):

            fullstring = items[i].firstChild.data
            substring = "Error"


            if fullstring.find(substring) == -1:
                if items[i].attributes['bts:Type'].value=='IMG':
                    rows.append({

                         "ImageName": os.path.join(PATH_TO_FILES, items[i].firstChild.data), 
                         "column": items[i].attributes['bts:Column'].value, 
                         "row": items[i].attributes['bts:Row'].value, 
                         "time_point": items[i].attributes['bts:TimePoint'].value, 
                         "field_index": items[i].attributes['bts:FieldIndex'].value, 
                         "z_slice": items[i].attributes['bts:ZIndex'].value, 
                         "channel": items[i].attributes['bts:Ch'].value,
                         "x_coordinates": items[i].attributes['bts:X'].value,
                         "y_coordinates": items[i].attributes['bts:Y'].value,
                         "z_coordinate": items[i].attributes['bts:Z'].value,
                         "action_index": items[i].attributes['bts:ActionIndex'].value,
                         "action": items[i].attributes['bts:Action'].value, 
                         "Type": items[i].attributes['bts:Type'].value, 
                         "Time": items[i].attributes['bts:Time'].value,
                         "PixPerMic": items_mrf[0].attributes['bts:HorizontalPixelDimension'].value
                    })

        out_df = pd.DataFrame(rows) #, columns = df_cols)

        return out_df
    
    def time_lapse_read(self, select_color_nuc, select_color_spot):
    
        from Analysis import ImageAnalyzer
        analyzer_fun=ImageAnalyzer()
        time_series_nuc=[]
        time_series_spot=[]
        time_series_vid=[]
        time_points = select_color_nuc["time_point"].unique().astype(int)
        z_slices = select_color_nuc["z_slice"].unique().astype(int)
        all_spots=[]
        for t in time_points:
            select_time_nuc=select_color_nuc.loc[select_color_nuc['time_point'] == str(t)]
            z_imglist_nuc=[]
            select_time_spot=select_color_spot.loc[select_color_spot['time_point'] == str(t)]
            z_imglist_spot=[]
            for z in z_slices:

                select_z_nuc = select_time_nuc.loc[select_time_nuc['z_slice'] == str(z)]
                select_z_spot = select_time_nuc.loc[select_time_nuc['z_slice'] == str(z)]
                im_nuc = Image.open(select_z_nuc['ImageName'].iloc[0])
                im_spot = Image.open(select_z_spot['ImageName'].iloc[0])
                z_imglist_nuc.append( np.asarray(im_nuc))
                z_imglist_spot.append( np.asarray(im_spot))
            z_stack_nuc = np.stack(z_imglist_nuc, axis=2)
            z_stack_spot = np.stack(z_imglist_spot, axis=2)
            max_project_nuc = z_stack_nuc.max(axis=2)
            max_project_spot = z_stack_spot.max(axis=2)

            img_uint8_nuc = cv2.normalize(max_project_nuc, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
            img_uint8_spot = cv2.normalize(max_project_spot, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)

            coordinates_max_project,_ =analyzer_fun.SpotDetector(input_image=img_uint8_spot,nuclei_image=img_uint8_nuc, 
                                                                 detect_method='LOG', threshold_method='auto', kernel_size=3)

            if coordinates_max_project.__len__()>0:

                coordinates_max_project_round = np.round(np.asarray(coordinates_max_project)).astype('int')

                coordinates_max_project = np.array(coordinates_max_project)[coordinates_max_project_round.min(axis=1)>=0,:].tolist()
                coordinates_max_project_round = coordinates_max_project_round[coordinates_max_project_round.min(axis=1)>=0,:]
                spots_z_slices = np.argmax(z_stack_spot[coordinates_max_project_round[:,0],coordinates_max_project_round[:,1],:], axis=1)
                spots_z_coordinates = np.zeros((spots_z_slices.__len__(),1), dtype='float')


                for i in range(spots_z_slices.__len__()):

                    spots_z_coordinates[i] = np.asarray(select_color_spot.loc[select_color_spot['z_slice']== str(spots_z_slices[i]+1)]
                                                     ['z_coordinate'].iloc[0]).astype('float')
                xyz_coordinates = np.append(np.asarray(coordinates_max_project).astype('float'), spots_z_coordinates, 1)

            all_spots.append(xyz_coordinates)
            time_series_nuc.append(np.asarray(img_uint8_nuc))
            time_series_spot.append(np.asarray(img_uint8_spot))
        t_stack_nuc = np.stack(time_series_nuc, axis=0)
        t_stack_spot = np.stack(time_series_nuc, axis=0)

        return t_stack_nuc, t_stack_spot, all_spots
    
    def COORDINATES_TO_CIRCLE(self, coordinates,ImageForSpots):
        
        circles = np.zeros((ImageForSpots.shape), dtype=np.uint8)
        
        if coordinates.any():
            
            for center_y, center_x in zip(coordinates[:,0], coordinates[:,1]):
                    circy, circx = circle_perimeter(center_y, center_x, 7, shape=ImageForSpots.shape)
                    circles[circy, circx] = 255

        return circles
