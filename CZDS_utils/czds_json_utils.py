import regex as re
import os
from pathlib import Path
import sys
from numpy import diff
import numpy as np
import json
import math
from scipy.stats import linregress
from lmfit import Model
from sklearn.metrics import root_mean_squared_error
#from xcbgen.xtypes import true_values

def poly2(x, a, b, c):
        """Defining the quadratic function to perform the EXAFS normalization."""
        polinomio = (a*x + b)*x + c
        return polinomio
    
class XDI():
    xdi_dict   = {}
    __xdi_filename = ''
    
    def __init__(self, filename, normalize=False):
        
        self.__xdi_filename = filename
        
        lines = []
        metadata_lines = []
        data_lines = []
        
        xdifile = Path(filename)
        if xdifile.is_file():
            with open(xdifile) as f:
                lines = [line.rstrip() for line in f]
        else:
            raise ValueError(f"'{filename}' is not a valid file. Check the file name.")
    
       
        #regex patterns formalizing xdi syntax
        #semantics is fully ignored at this stage,
        #in order to keep the approach fully general
        numeric_data_pattern   = re.compile(r'^[^#]([\s*][+-]?(?=\.\d|\d)(?:\d+)?(?:\.?\d*))?.*$')
        version_field_pattern  = re.compile(r'^#[\s+]XDI/[\w\W\s]+$')
        metadata_field_pattern = re.compile(r'^#[\s+][\w]+\.[\w]+:[\s+][\w\W\s]+$')
    
        field_eol_pattern      = re.compile(r'^#[\s+][/]+$')
        header_eol_pattern     = re.compile(r'^#[\s+][-]+$')
        float_pattern          = re.compile(r'([\s+]?[+-]?(?=\.\d|\d)(?:\d+)?(?:\.?\d*))')
        
        
        
        for line in lines :
            if metadata_field_pattern.search(line) or version_field_pattern.search(line):
                metadata_lines.append(line)
            elif numeric_data_pattern.search(line):
                numbers = re.findall( float_pattern, line)
                data_lines.append( numbers)
        
        #the code above separates data from metadata
        
        version_pattern  = r'^#[\s+]XDI/(?P<version>[\w\W\s]+)$'
        metadata_pattern = r'^#[\s+](?P<namespace>[\w]+)\.(?P<attribute>[\w]+):[\s+](?P<value>[\w\W\s]+)$'
    
        self.xdi_dict = {"metadata":{}, "data":{}}
    
        for line in metadata_lines:
            match_fields = re.search(version_pattern, metadata_lines[0])
            if match_fields:
                self.xdi_dict['metadata']['Version']=match_fields.group('version')
            
            match_fields = re.search(metadata_pattern, line)
            if match_fields:
                if match_fields.group('namespace') in self.xdi_dict['metadata']:
                    self.xdi_dict['metadata'][match_fields.group('namespace')][match_fields.group('attribute')]=match_fields.group('value')
                else:
                    self.xdi_dict['metadata'][match_fields.group('namespace')]={}
                    self.xdi_dict['metadata'][match_fields.group('namespace')][match_fields.group('attribute')]=match_fields.group('value')
        
        ncolumns = len(self.xdi_dict['metadata']['Column'])
        for ncolumn in range(ncolumns): 
            ncolumn_str = str(ncolumn+1)
            values = []
            for x in data_lines:
                if len(x) == ncolumns:
                    values.append( float(x[ncolumn]) ) 
            self.xdi_dict['data'][ncolumn_str]=values
        
        #from here on some semantic analysis is performed 
        #to identify the type of experiment
        #and the type of measurement  
        self.xdi_dict['metadata']["Experiment"] = {}
        self.xdi_dict['metadata']["Experiment"]["Measurement"] = "unknown"
        for ncolumn in range(ncolumns):
           ncolumn_str = str(ncolumn+1)
           if ( "itrans" in self.xdi_dict['metadata']['Column'][ncolumn_str].lower() ):
               self.xdi_dict['metadata']["Experiment"]["Measurement"] = "transmission"
               break
           if "fluorescence" in self.xdi_dict['metadata']['Column'][ncolumn_str].lower():
               self.xdi_dict['metadata']["Experiment"]["Measurement"] = "fluorescence"
               break
           if ("transmission" in self.xdi_dict['metadata']['Column'][ncolumn_str].lower()) or ("raw" in self.xdi_dict['metadata']['Column'][ncolumn_str].lower()):
               self.xdi_dict['metadata']["Experiment"]["Measurement"] = "raw transmission"
               break
           if "normalized" in self.xdi_dict['metadata']['Column'][ncolumn_str].lower():
               self.xdi_dict['metadata']["Experiment"]["Measurement"] = "normalized transmission"
               break
        
        #classify type of XAS: XANES or EXAFS
        self.xdi_dict['metadata']["Experiment"]["Type"]="undefined"
        for ncolumn in range(ncolumns): 
            ncolumn_str = str(ncolumn+1)
            if "energy" in self.xdi_dict['metadata']['Column'][ncolumn_str].lower():
                energy_range = max( self.xdi_dict['data'][ncolumn_str] ) - min( self.xdi_dict['data'][ncolumn_str] )
                if energy_range < 400:
                    self.xdi_dict['metadata']["Experiment"]["Type"] = "XANES"
                else:
                    self.xdi_dict['metadata']["Experiment"]["Type"] = "EXAFS" 
                break
            
        if normalize: self.normalize(True)
    
    def __str__(self):
        return json.dumps( self.xdi_dict, indent=5, default=str)
    
    def get_dict(self):
        return self.xdi_dict

    def set_dict(self, value):
        self.xdi_dict = value
        
    def save(self, filename=''):
        json_filename = ''
        
        if filename=='':
            json_filename = Path(self.__xdi_filename).stem + ".json"
        else:
            json_filename = filename
            
        with open(json_filename, 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=4)
        
        
    def __find_E0(self, energies, absorption):
        """Finds the absorption edge E0

        Keyword arguments:
        energies -- the array of energies from the experiment.
        absorption -- the absorption data, ln(i0/itrans) for transmission or i0/ifluo for fluorescence.

        Returns:
        The absorption edge E0 location in x (E0x) and y (E0y) axis.
        """

        x = energies
        y =  absorption
        dydx = diff(y)/diff(x)
        x_new = x[:-1]
        dydx2 = diff(dydx)/diff(x_new)
        x_new2 = x_new[:-1]
    
        #Test for the largest derivative, only if the second derivative is positive.
        #Searching only up to half the graph. Should be to the first peak, but for now didn't implemented it. 
        #TODO, SEARCH ONLY UP TO THE FIRST PEAK
        range_search_E0_start=0
        range_search_E0_end=int(len(energies)/2)
        deriv_max = 0
        point_max = 0 
        for point in range(range_search_E0_start,range_search_E0_end):
            if dydx[point] > deriv_max:
                if dydx2[point] > 0:
                    deriv_max = dydx[point]
                    point_max = point

        E0x = x[point_max]
        E0y = y[point_max]
    
        return E0x, E0y
    
    def __fit_pre_edge(self, energies_array, absorption_array, start_pre_edge_x_index, end_pre_edge_x_index):
        #Define the background
        #print('indexes start finish',start_pre_edge_x_index,end_pre_edge_x_index)
        background_x = energies_array[start_pre_edge_x_index:end_pre_edge_x_index]
        background_y = absorption_array[start_pre_edge_x_index:end_pre_edge_x_index]
        #Linear model for the pre-edge.
        #Perform the linear fit.
        resultado_fit = linregress(background_x, background_y)
        
        #Extrapolate for the whole range.
        predicted_initial_range = np.array(resultado_fit.intercept + 
                                           resultado_fit.slope*energies_array)
    
        # Normalize the pre-edge
        #df['Absorption']=df['Absorption'].values - predicted_initial_range
        #normalized_pre_edge = absorption_array - predicted_initial_range
    
        #return normalized_pre_edge, predicted_initial_range
        return predicted_initial_range
    
    def __poly2(x, a, b, c):
        """Defining the quadratic function to perform the EXAFS normalization."""
        polinomio = (a*x + b)*x + c
        return polinomio
    
    def __EXAFS_normalization(self, energies_array, absorption_array):
        
        #Find E0
        E0x, E0y = self.__find_E0(energies_array, absorption_array)
        
        #Get E0 index
        E0x_index = np.where(energies_array==E0x)[0]
        
        #Get the ranges where the pre-edge is located
        start_to_E0 = E0x - energies_array[0]
        start_to_E0_point_spacing = E0x_index/start_to_E0
      
        #Getting the starting energy for the pre-edge fit and its index in the array
        start_pre_edge_x = energies_array[0]
        start_pre_edge_x_index = np.where(energies_array==start_pre_edge_x)[0][0]
        
        #Here I am defining the end of the pre-edge region, with 30 points before E0
        pre_edge_end_fit_energy_from_E0 = 30
        end_pre_edge_x = start_pre_edge_x + (start_to_E0 - pre_edge_end_fit_energy_from_E0)
        idx = np.searchsorted(energies_array, end_pre_edge_x, side="left")
        #I do not remember what this test does, need to check
        if idx > 0 and (idx == len(energies_array) or math.fabs(end_pre_edge_x - energies_array[idx-1]) < math.fabs(end_pre_edge_x - energies_array[idx])):
            end_pre_edge_x_index = idx
        else:
            end_pre_edge_x_index = idx
    
        #Normalize the pre-edge with the obtained ranges.
        linear_fit_pre_edge = self.__fit_pre_edge(energies_array, absorption_array, start_pre_edge_x_index, end_pre_edge_x_index)
    
        #Starting the EXAFS Normalization.
        #It will perform quadratic fits in a range of starting points from 100 points after E0 to 200 points after E0.
        #The lowest RMSE of all the fits will be chosen as the best fit.
        #Define initial point for the starting range
        exafs_start_shit_from_E0 = 100
        start_exafs = np.where(energies_array > E0x + exafs_start_shit_from_E0)[0][0]
    
        #Define ending point for the starting range
        exafs_end_shit_from_E0 = 200
        if (E0x + exafs_end_shit_from_E0) > energies_array[-1]:
            #print('DEBUG1', E0x + exafs_end_shit_from_E0, energies_array[-1])
            end_exafs = np.where(energies_array > energies_array[-1] - 50)[0][0]
        else:
            #print('DEBUG2', E0x + exafs_end_shit_from_E0, energies_array[-1])
            end_exafs = np.where(energies_array > E0x + exafs_end_shit_from_E0)[0][0]
        
        #Starting a very large RMSE
        rmse_min = 1000
    
        #Define final point for the EXAFS FIT, it will always be this one for all fits, only the starting point changes.
        #It is 30 points before the end of the energies_array.
        last_point_to_fit_exafs = np.where(energies_array > energies_array[-1] - 30)[0][0]
    
        for npt in range(start_exafs, end_exafs):
            np_init = npt #Initial point of the interval
            np_final = last_point_to_fit_exafs #End of the array
    
            #Slice the data
            data_x = energies_array[np_init:np_final]
            data_y = absorption_array[np_init:np_final]
            
            #Perform the fit
            model_poly2 = Model(poly2)
            params = model_poly2.make_params(a=1.0, b=1.0, c=1.0)
            fit_result = model_poly2.fit(data_y, params, x=data_x)
            rmse = root_mean_squared_error(data_y, fit_result.best_fit)
            
            #Test if the RMSE is smaller then the smallest found yet.
            if abs(rmse) < abs(rmse_min):
                rmse_min = rmse
                npt_min = npt
    
        #Defining the range from the lowest RMSE found.
        data_x = energies_array[npt_min:last_point_to_fit_exafs]
        data_y = absorption_array[npt_min:last_point_to_fit_exafs]
        
        #Perform the fit again, on the best range found.
        final_exafs_fit = model_poly2.fit(data_y, params, x=data_x)
    
        #Infer the values
        predicted_exafs = 0
        predicted_exafs = np.array(final_exafs_fit.best_values['a']*energies_array**2 + 
                                   final_exafs_fit.best_values['b']*energies_array +
                                   final_exafs_fit.best_values['c'])
        
        #Perform the normalization
        edge_step = predicted_exafs[E0x_index] - linear_fit_pre_edge[E0x_index]
        normalized_spectra = (absorption_array - linear_fit_pre_edge) / (edge_step)
    
        data_x = energies_array[npt_min:last_point_to_fit_exafs]
        data_y = normalized_spectra[npt_min:last_point_to_fit_exafs]
        final_exafs_fit = model_poly2.fit(data_y, params, x=data_x)
        predicted_exafs = 0
        predicted_exafs = np.array(final_exafs_fit.best_values['a']*energies_array**2 + 
                                   final_exafs_fit.best_values['b']*energies_array +
                                   final_exafs_fit.best_values['c'])
        
        #Perform the flattening, with subtraction mode. One can also choose to divide the curve, but it is not the standard procedure
        method = 'subtract'
        flattened_curve = []
        if method == 'divide':
            #print('Normalizing EXAFS regions with division method')
            for point in range(len(energies_array)):
                if energies_array[point] < E0x:
                    flattened_curve.append(normalized_spectra[point])
                else:
                    flatting = normalized_spectra[point]/predicted_exafs[point]
                    #print(normalized_spectra[point], predicted_exafs[point], predicted_exafs[E0x_index][0], - predicted_exafs[point] + predicted_exafs[E0x_index][0])
                    flattened_curve.append(flatting) #[0])
        elif method == 'subtract':
            #print('Normalizing EXAFS regions with subtraction method')
            for point in range(len(energies_array)):
                if energies_array[point] < E0x:
                    flattened_curve.append(normalized_spectra[point])
                else:
                    flatting = normalized_spectra[point] - predicted_exafs[point] + predicted_exafs[E0x_index][0]
                    #print(normalized_spectra[point], predicted_exafs[point], predicted_exafs[E0x_index][0], - predicted_exafs[point] + predicted_exafs[E0x_index][0])
                    flattened_curve.append(flatting) #[0])


        return energies_array, flattened_curve

    def __XANES_normalization(self, energies_array, absorption_array):
        """Recebe um dataset com dados de energia e absorbância e normaliza o dataset."""

        #Call the E0 calculation function and map it to the input spectra.
        E0x, E0y = self.__find_E0(energies_array, absorption_array)
        E0x_index = np.where(energies_array==E0x)[0]

        #Define the initial and final points of the linear fit of the background.
        ##This definition seemed good for the the tested spectra, might be changed after more tests are performed.
        start_pre_edge_x_index = int(E0x_index/15)
        end_pre_edge_x_index = int(E0x_index/3)

        #Linear model for the pre-edge.
        data_x = energies_array[start_pre_edge_x_index:end_pre_edge_x_index] #background.iloc[:, 0].values
        data_y = absorption_array[start_pre_edge_x_index:end_pre_edge_x_index] #background.iloc[:, 1].values

        #Perform the linear fit.
        linear_fit_pre_edge = self.__fit_pre_edge(energies_array, absorption_array, start_pre_edge_x_index, end_pre_edge_x_index)

        #Perform the normalization
        edge_step = absorption_array[E0x_index] - linear_fit_pre_edge[E0x_index]
       
        pre_edge_normalized_spectra = (absorption_array - linear_fit_pre_edge) / E0y #(edge_step)

        #Start the needed variables.
        slope_min = 100000

        np_start = E0x_index[0]
        np_end = len(pre_edge_normalized_spectra)

        for npt in range(np_start,len(pre_edge_normalized_spectra-5)): 
            np_init = npt #Initial point of the interval

            #Define the range for the linear fit.
            data_x = energies_array[np_init:np_end] 
            data_y = pre_edge_normalized_spectra[np_init:np_end] 

            #Perform the fit
            resultado_fit = linregress(data_x, data_y)
            predicted_resulto_fit = np.array(resultado_fit.intercept +
                                         resultado_fit.slope*data_x)

            #Test if the slope is smaller then the smallest found yet.
            if abs(resultado_fit.slope) < abs(slope_min):
                slope_min = resultado_fit.slope
                min_intercept = resultado_fit.intercept
                npt_min = npt

        #Defining the data in range
        data_x = energies_array[npt_min:np_end] #faixa_final.iloc[:, 0].values
        data_y = pre_edge_normalized_spectra[npt_min:np_end] #faixa_final.iloc[:, 1].values

        #Perform the linear fit again. TODO: No need, we could just get from the already found parameters, although I did not saved the intercepts.
        resultado_fit_final = linregress(data_x, data_y)

        #Extrapolate for the whole range.
        xwide = energies_array 
        predicted_faixa_final = np.array(resultado_fit_final.intercept +
                                         resultado_fit_final.slope*xwide)

        #Normalized spectra.
        fit_final = pre_edge_normalized_spectra/predicted_faixa_final

        return energies_array, fit_final

    def __get_column_with_content(self, name):
        ncolumns = len(self.xdi_dict['metadata']['Column'])
        ncolumn_str = ""
        column_index = -1
        for ncolumn in range(ncolumns): 
            ncolumn_str = str(ncolumn+1)
            if name.lower() in self.xdi_dict['metadata']['Column'][ncolumn_str].lower():
                column_index = ncolumn
                break
        
        if column_index == -1 :
            return  []
        else:        
            return self.xdi_dict['data'][ncolumn_str] 
    
    
    def normalize(self, status=False):
        
        #get the data. Data not found will be represented by empty arrays
        energy   = np.array( self.__get_column_with_content("energy") )
        itrans   = np.array( self.__get_column_with_content("itrans") )
        irefer   = np.array( self.__get_column_with_content("irefer") )
        i0       = np.array( self.__get_column_with_content("i0")     )
        raw_data = np.array( self.__get_column_with_content("raw")    )
        fluorescence_emission  = np.array( self.__get_column_with_content("fluorescence")     )
        normalized_adsorbance  = np.array( self.__get_column_with_content("normalized transmission")     )
       
        if status:
            print( "energy : ",  bool(energy.size), " itrans : ",  bool(itrans.size),
                   " irefer : ",  bool(irefer.size), " i0 : ",  bool(i0.size),
                 " raw_data : ",  bool(raw_data.size), " fluorescence_emission : ",  bool(fluorescence_emission.size), 
                 " normalized_adsorbance : ",  bool(normalized_adsorbance.size) )
             
            print("metadata.Experiment.Measurement :", self.xdi_dict['metadata']["Experiment"]["Measurement"] )
            print("metadata.Experiment.Type : ", self.xdi_dict['metadata']["Experiment"]["Type"] ) 
                
        enorm   = None
        mu_norm = None
        
        if self.xdi_dict['metadata']["Experiment"]["Measurement"] == "transmission":
            
            if self.xdi_dict['metadata']["Experiment"]["Type"] == "EXAFS":
                
                enorm, mu_norm = self.__EXAFS_normalization(energy, np.log(i0/itrans) )
                    
            elif self.xdi_dict['metadata']["Experiment"]["Type"] == "XANES":
                
                enorm, mu_norm = self.__XANES_normalization(energy, np.log(i0/itrans) )
                
        elif self.xdi_dict['metadata']["Experiment"]["Measurement"] == "raw transmission":
            
            if self.xdi_dict['metadata']["Experiment"]["Type"] == "EXAFS":
                
                enorm, mu_norm = self.__EXAFS_normalization(energy, raw_data )
                
            elif self.xdi_dict['metadata']["Experiment"]["Type"] == "XANES":
                
                enorm, mu_norm = self.__XANES_normalization(energy, raw_data )
            
        elif self.xdi_dict['metadata']["Experiment"]["Measurement"] == "normalized transmission":
            
            if self.xdi_dict['metadata']["Experiment"]["Type"] == "EXAFS":
                
                enorm, mu_norm = self.__EXAFS_normalization(energy, normalized_adsorbance )
                
            elif self.xdi_dict['metadata']["Experiment"]["Type"] == "XANES":      
                
                enorm, mu_norm = self.__XANES_normalization(energy, normalized_adsorbance )
        
        elif self.xdi_dict['metadata']["Experiment"]["Measurement"] == "fluorescence":
         
            if self.xdi_dict['metadata']["Experiment"]["Type"] == "EXAFS":
                
                enorm, mu_norm = self.__EXAFS_normalization(energy, fluorescence_emission )
                
            elif self.xdi_dict['metadata']["Experiment"]["Type"] == "XANES":      
                
                enorm, mu_norm = self.__XANES_normalization(energy, fluorescence_emission )
        
        if not ("normalized" in self.xdi_dict['data']):
            self.xdi_dict['data']["normalized"]={}
            
        self.xdi_dict['data']["normalized"]["mu"]     =  mu_norm
     
               