
class Gui_Params(object):
    
    def __init__(self, is_gui, input_params=None, analysisgui=None, inout_resource_gui=None):
        self.is_gui = is_gui
        if is_gui == False:
            self.output_dir = input_params["output_dir"]
            self.NucInfoChkBox_check_status = input_params["nuc_info"] #self.AnalysisGui.NucInfoChkBox.isChecked()
            self.SpotsLocation_check_status = input_params["spots_location"] #self.AnalysisGui.SpotsLocation.isChecked()
            self.SpotLocationCbox_currentText = input_params["spot_coordinates"] #self.AnalysisGui.SpotLocationCbox.currentText()
            self.SpotsDistance_check_status = input_params["spots_distance"] #self.AnalysisGui.SpotsDistance.isChecked()
            self.NucMaskCheckBox_status_check = input_params["nuc_mask"] #self.AnalysisGui.NucMaskCheckBox.isChecked()
            self.NucMaxZprojectCheckBox_status_check = input_params["nuclei_z_project"] #self.AnalysisGui.NucMaxZprojectCheckBox.isChecked()
            self.SpotMaxZProject_status_check = input_params["spot_z_project"] #self.AnalysisGui.SpotMaxZProject.isChecked()
            self.SpotCh1CheckBox_status_check = input_params["ch1_spot"] #self.AnalysisGui.SpotCh1CheckBox.isChecked()
            self.SpotCh2CheckBox_status_check = input_params["ch2_spot"] #self.AnalysisGui.SpotCh2CheckBox.isChecked()
            self.SpotCh3CheckBox_status_check = input_params["ch3_spot"] #self.AnalysisGui.SpotCh3CheckBox.isChecked()
            self.SpotCh4CheckBox_status_check = input_params["ch4_spot"] #self.AnalysisGui.SpotCh4CheckBox.isChecked()
            self.SpotCh5CheckBox_status_check = input_params["ch5_spot"] #self.AnalysisGui.SpotCh5CheckBox.isChecked()
            self.NucleiChannel_index = str(input_params["nuclei_channel"])[-1] #self.AnalysisGui.NucleiChannel.currentIndex()
            self.NumCPUsSpinBox_value = 8#self.inout_resource_gui.NumCPUsSpinBox.value()
            self.Cell_Tracking_check_status = input_params["cell_tracking"] #self.AnalysisGui.Cell_Tracking.isChecked()
            self.NucTrackingMethod_currentText = input_params["nuc_track_method"] #self.AnalysisGui.NucTrackMethod.currentText()
            self.NucSearchRadiusSpinbox_current_value = input_params["nuc_search_radius"] #self.AnalysisGui.NucSearchRadiusSpinbox.value()
            self.SpotSearchRadiusSpinbox_current_value = input_params["spot_search_radius"] #self.AnalysisGui.SpotSearchRadiusSpinbox.value()
            
            self.MintrackLengthSpinbox_current_value = input_params["min_track_length"] #self.AnalysisGui.MintrackLengthSpinbox.value()
            
            #for analysis, variable names same as in analysisGui
            self.NucDetectMethod = input_params["nuclei_detection_method"]
            self.NucDetectionSlider = input_params["nuclei_detection"]
            self.NucSeparationSlider = input_params["nuclei_separation"]
            self.NucleiAreaSlider = input_params["nuclei_area"]
            self.NucRemoveBoundaryCheckBox = input_params["remove_boundary_nuclei"]
            self.spotchannelselect = input_params["spot_channel_select"]
            self.SpotPerChSpinBox = input_params["spots_per_channel"]
            #spot_params_dict initialization (only need to do it onceas you can only pick one channel either way)
            self.spotanalysismethod = input_params["ch1_spot_detection_method"] #self.AnalysisGui.spotanalysismethod.currentIndex()
            self.thresholdmethod = input_params["ch1_spot_threshold_method"] #self.AnalysisGui.thresholdmethod.currentIndex()
            self.ThresholdSlider = input_params["ch1_spot_threshold_value"] #self.AnalysisGui.ThresholdSlider.value()
            self.SensitivitySpinBox = input_params["ch1_kernel_size"] #self.AnalysisGui.SensitivitySpinBox.value()
            self.SpotsPerChSpinBox = input_params["ch1_spots/ch"] #self.AnalysisGui.SpotPerChSpinBox.value()
        if is_gui == True:
            self.AnalysisGui = analysisgui
            self.inout_resource_gui = inout_resource_gui
            self.NucInfoChkBox_check_status = self.AnalysisGui.NucInfoChkBox.isChecked()
            self.SpotsLocation_check_status = self.AnalysisGui.SpotsLocation.isChecked()
            self.SpotLocationCbox_currentText = self.AnalysisGui.SpotLocationCbox.currentText()
            self.SpotsDistance_check_status = self.AnalysisGui.SpotsDistance.isChecked()
            self.NucMaskCheckBox_status_check = self.AnalysisGui.NucMaskCheckBox.isChecked()
            self.NucMaxZprojectCheckBox_status_check = self.AnalysisGui.NucMaxZprojectCheckBox.isChecked()
            self.SpotMaxZProject_status_check = self.AnalysisGui.SpotMaxZProject.isChecked()
            self.SpotCh1CheckBox_status_check = self.AnalysisGui.SpotCh1CheckBox.isChecked()
            self.SpotCh2CheckBox_status_check = self.AnalysisGui.SpotCh2CheckBox.isChecked()
            self.SpotCh3CheckBox_status_check = self.AnalysisGui.SpotCh3CheckBox.isChecked()
            self.SpotCh4CheckBox_status_check = self.AnalysisGui.SpotCh4CheckBox.isChecked()
            self.SpotCh5CheckBox_status_check = self.AnalysisGui.SpotCh5CheckBox.isChecked()
            self.NucleiChannel_index = self.AnalysisGui.NucleiChannel.currentIndex() + 1
            self.NumCPUsSpinBox_value = self.inout_resource_gui.NumCPUsSpinBox.value()
            self.Cell_Tracking_check_status = self.AnalysisGui.Cell_Tracking.isChecked()
            self.NucTrackingMethod_currentText = self.AnalysisGui.NucTrackMethod.currentText()
            self.NucSearchRadiusSpinbox_current_value = self.AnalysisGui.NucSearchRadiusSpinbox.value()
            self.SpotSearchRadiusSpinbox_current_value = self.AnalysisGui.SpotSearchRadiusSpinbox.value()
            
            self.MintrackLengthSpinbox_current_value = self.AnalysisGui.MintrackLengthSpinbox.value()
            
            self.NucDetectMethod = self.AnalysisGui.NucDetectMethod.currentText()
            self.NucDetectionSlider = self.AnalysisGui.NucDetectionSlider.value()
            self.NucSeparationSlider = self.AnalysisGui.NucSeparationSlider.value()
            self.NucleiAreaSlider = self.AnalysisGui.NucleiAreaSlider.value()
            self.NucRemoveBoundaryCheckBox = self.AnalysisGui.NucRemoveBoundaryCheckBox.isChecked()
            self.spotchannelselect = self.AnalysisGui.spotchannelselect.currentText()
            self.SpotPerChSpinBox = self.AnalysisGui.SpotPerChSpinBox.value()
            #spot params dict
            self.spotanalysismethod = self.AnalysisGui.spotanalysismethod.currentIndex()
            self.thresholdmethod = self.AnalysisGui.thresholdmethod.currentIndex()
            self.ThresholdSlider = self.AnalysisGui.ThresholdSlider.value()
            self.SensitivitySpinBox = self.AnalysisGui.SensitivitySpinBox.value()
            self.SpotsPerChSpinBox = self.AnalysisGui.SpotPerChSpinBox.value()