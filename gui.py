# -*- coding: utf-8 -*-
"""
@author: gmeier
DMS gui
"""
import tkinter as tk
import json
import os
import DMS_processing_multiprocessing



class DMS(tk.Tk):
    """
    class to start gui for a new dataanalysis
    Allows user to define input directory, output directory, DNA sequence and reading frame of target

    """
    def __init__(self):
        super().__init__()
        
        
        
        self.title('DMS data processing')
        self.geometry('800x200')
        self.overall_canvas=tk.Canvas(self,width=200,height=400)
        self.overall_canvas.grid(row=0)
        tk.Label(self.overall_canvas,text='input file directory').grid(row=0, column=0)
        tk.Label(self.overall_canvas,text='output file directory').grid(row=1,column=0)
        tk.Label(self.overall_canvas,text='job name').grid(row=2,column=0)
        tk.Label(self.overall_canvas,text='enter positions to be analysed (comaseparated, no spaces)').grid(row=3, column=0)
        tk.Label(self.overall_canvas,text='enter reference DNA sequence').grid(row=4, column=0)
        #defining input file directory and searich for all bam files in directory over Entry validation
        self.job_name=tk.Entry(self.overall_canvas)
        self.job_name.grid(row=2, column=2)          

        self.position_list=tk.Entry(self.overall_canvas)
        self.position_list.grid(row=3, column=2) 

        self.ref_seq=tk.Entry(self.overall_canvas)
        self.ref_seq.grid(row=4, column=2) 
             
        self.input_file_directory=tk.Entry(self.overall_canvas)
        self.input_file_directory.grid(row=0, column=2)
        
        self.output_file_directory=tk.Entry(self.overall_canvas)
        self.output_file_directory.grid(row=1, column=2)
        
        
        
        self.process_button=tk.Button(self, text='process', command=self.process_DMS_data)
        self.process_button.grid(row=10,column=0)
        
        self.var=tk.IntVar()
        self.var2=tk.IntVar()
        
        self.Json_file_available=tk.Checkbutton(self, text='do you wish to load configurations form a Json file?', command=self.old_Json,variable=self.var2)
        self.Json_file_available.grid(row=6,column=0)
        
        self.readingFrame_checkbutton= tk.Checkbutton(self,text='are you working with several reading frames?', command=self.reading_frame,variable=self.var)
        self.readingFrame_checkbutton.grid(row=2,column=0)
        self.close_button=tk.Button(self, text='close', command=self.destroy)
        self.close_button.grid(row=10, column=1)
        
        
   
        
        
    def reading_frame(self):
        """
        Creates entry point for reading frame and frameshift offset
        
        """
        
        if self.var.get() ==1:
            self.readingframecan=tk.Canvas(self,width=200, height=100)
            self.readingframecan.grid(row=3)
            self.readframe=tk.Frame(self.readingframecan)
            self.readframe.grid()
            tk.Label(self.readframe,text='position of last base of reading frame 1').grid(row=0,column=0)
            self.rf=tk.Entry(self.readframe)
            self.rf.grid(row=0, column=1)
            tk.Label(self.readframe,text='frameshift offset').grid(row=1,column=0)
            self.fo=tk.Entry(self.readframe)
            self.fo.grid(row=1, column=1)
            
        else:
            self.readingframecan.destroy()
            
    def old_Json(self):
        """
        entry point to load a existing Json file from an old project.
        Usefull if reanalysis has to be done
        """
        
        if self.var2.get()==1:
            self.oldjsoncan=tk.Canvas(self,width=200, height=100)
            self.oldjsoncan.grid(row=7)
            self.jsonframe=tk.Frame(self.oldjsoncan)
            self.jsonframe.grid()
            tk.Label(self.jsonframe, text='Location of Json file').grid(row=0, column=0)
            self.jsonentry=tk.Entry(self.jsonframe)
            self.jsonentry.grid(row=0, column=1)


            
            
    def create_Json_object(self):
        """Write data to a json object. either from new entry or if  loading of old jsonfile is chosen from old jsonfile
        
        """
        
        print('creating json object')
        jsonfile=Json_file()
        if self.var2.get()==1:
            jsonfile.load_old_json(self.jsonentry.get())
            jsonfile.data_files=jsonfile.find_input_files()
            with open(jsonfile.outputdir+'/Json_'+jsonfile.job_name+'.json','w') as write_file:
           
                json.dump(jsonfile.create_data(), write_file,sort_keys = True, indent = 4)   
            
            
        else:
            jsonfile.inputdir=(self.input_file_directory.get())
            jsonfile.outputdir=(self.output_file_directory.get())
            jsonfile.job_name=self.job_name.get()
            jsonfile.data_files=jsonfile.find_input_files()
            jsonfile.position_list=(self.position_list.get()).split(',')
            jsonfile.reference_sequence=(self.ref_seq.get())
            if self.var.get() ==1:
                jsonfile.readingframes = True
                jsonfile.frameshift_position=int(self.rf.get())
                jsonfile.frameshift_offset=int(self.fo.get())
#producing sequences for readingframes and save to json
                jsonfile.frame1=jsonfile.reference_sequence[:jsonfile.frameshift_position]
                jsonfile.frame2=jsonfile.reference_sequence[jsonfile.frameshift_position-jsonfile.frameshift_offset:]
                

            jsonfile.create_json()   
        return(jsonfile)

    def process_DMS_data(self):
        """
        Starts data processing from Json file after all neccessary input have been saved.
        
        """
        jsonfile=self.create_Json_object()
        
        #start of processing
        DMS_processing_multiprocessing.run_analysis(jsonfile.outputdir+'/Json_'+jsonfile.job_name+'.json')


      

#Creat and wirte Jsone file with input settings (used to retriev input data for dms scripts)

class Json_file(object):
    """
    Class to create an empty json file
    
    """
    
    def  __init__      (self):
        self.job_name=None
        self.inputdir=None
        self.outputdir=None
        self.readingframes=False
        self.frameshift_position=0
        self.frameshift_offset=0
        self.data_files=None
        self.position_list=None
        self.reference_sequence=None
        self.frame1=None
        self.frame2=None
        
   
    def create_data(self):
        """
        creates data and returns a data dict  
        
        """
        
        data_dict={}
        data_dict['job_name']=self.job_name
        data_dict['inputdir']=self.inputdir
        data_dict['outputdir']=self.outputdir
        data_dict['readingframes']=self.readingframes
        
        
        data_dict['frameshift_position']=int(self.frameshift_position)
        data_dict['frameshift_offset']=int(self.frameshift_offset)
        data_dict['data_files']=self.data_files
        data_dict['position_list']=self.position_list
        data_dict['reference_sequence']=self.reference_sequence
        data_dict['frame1']=self.frame1
        data_dict['frame2']=self.frame2
        return data_dict
        
#loads existing json file    
    def load_old_json(self,json_directory):
        """
        used to load a preexisting Json file 

        """
        
        
        with open(json_directory,'r') as load_file:
            data=json.load(load_file)
        
        self.inputdir=data['inputdir']
        self.outputdir=data['outputdir']
        self.job_name=data['job_name']
        self.readingframes=data['readingframes']
        self.frameshift_position=data['frameshift_position']
        self.frameshift_offset=data['frameshift_offset']
        self.position_list=data['position_list']
        self.reference_sequence=data['reference_sequence']
        self.frame1=data['frame1']
        self.frame2=data['frame2']
 
    def find_input_files(self):
        """
        This function looks for all input files in the input directory
        Files must have the ending 'sorted.bam'
        
        """
        bamfile_list=[]
        for file in os.listdir(self.inputdir):
            if file.endswith('sorted.bam'):
                bamfile_list.append(file)
                
        return bamfile_list        
        
#creates a json file and dump data_dict       
    def create_json(self):
        with open(self.outputdir+'/Json_'+self.job_name+'.json','w') as write_file:
           
            json.dump(self.create_data(), write_file,sort_keys = True, indent = 4)
            print ('writing json file')