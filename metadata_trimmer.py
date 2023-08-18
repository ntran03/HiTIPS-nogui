from xml.dom import minidom
import os
import shutil

#if you are planning on reusing the resulting temp file after an interrupt, make sure to rename
#to something different before running or it will write to the same file again
class metadata_trimmer(object):
    #initializes a new temp file with the path to the original metadata file and the name of the file
    def __init__(self, metadata_path, metadataname):
        self.metadatafilename = os.path.join(metadata_path, metadataname)
        #temp file will always have same name- can change this later
        self.tempfilename = os.path.join(metadata_path,'temp_metadata_file.mlf')
        #make a copy of the metadata file called temp_metadata_file
        shutil.copyfile(self.metadatafilename, self.tempfilename)
        print("copied")
    
    #trims the lines of the metadata file by the col,row,time,and fov of the corresponding image
    #don't be silly and use time as a variable name q-q
    def trim_metadata(self, col, row, t, fov):
        #parse metadata file
        self.mydoc = minidom.parse(self.tempfilename)
        self.items = self.mydoc.getElementsByTagName('bts:MeasurementRecord')
        #go through metadata file, look for corresponding lines and delete them
        for i in range(self.items.length):
            #print(self.items[i].attributes['bts:Column'].value)
            if self.items[i].attributes['bts:Column'].value == str(col) and self.items[i].attributes['bts:Row'].value == str(row) and self.items[i].attributes['bts:TimePoint'].value == str(t) and self.items[i].attributes['bts:FieldIndex'].value == str(fov):
                #self.mydoc.remove(self.items[i])
                self.items[i].parentNode.removeChild(self.items[i])
        #write the edited version of the file to the same file
        with open(self.tempfilename, "w" ) as fs: 
            fs.write(self.mydoc.toxml() )
            fs.close()
    #delete the metadata file after a successful run of the nuclear masks
    def del_temp(self):
        os.remove(self.tempfilename)

#testing
#trimmer = metadata_trimmer("C:\\Users\\tranne\\Downloads\\HiTIPS-main")
##trimmer.trim_metadata(10, 6, 1, 1)
#print("done")