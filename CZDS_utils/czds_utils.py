import regex as re
import os
import pandas as pd
import numpy as np
from numpy import diff
from scipy.stats import linregress
import math
from lmfit import Model
from sklearn.metrics import root_mean_squared_error
import matplotlib.pyplot as plt #DEBUG
from sklearn.metrics import mean_squared_error

class read():
    def xdi(filepath):
        secoes = [
            "Element.symbol",
            "Element.edge",
            "Mono.d_spacing",
            "Mono.name",
            "Sample.formula",
            "Sample.name",
            "Sample.prep",
            "Sample.temperature",
            "Sample.reference",
            "Detector.I0",
            "Detector.I1",
            "Detector.I2",
            "Facility.Name",
            "Beamline.Name",
            "Beamline.name",
            "Facility.name",
            "Beamline.xray_source",
            "Beamline.Storage_Ring_Current",
            "Beamline.I0",
            "Beamline.I1",
            "Scan.start_time",
            "Scan.end_time",
            "ScanParameters.Start",
            "ScanParameters.ScanType",
            "ScanParameters.E0",
            "ScanParameters.Legend",
            "ScanParameters.Region1",
            "ScanParameters.Region2",
            "ScanParameters.Region3",
            "ScanParameters.End"
            ]
        regex = '|'.join(map(re.escape, secoes))
    
        with open(filepath, 'r') as texto:
            linhas = texto.read()
            matches = re.findall(f'({regex}):\\s*(.*)', linhas)
        dicio = {}
        for match in matches:
            secao, valor = match[0], match[1]
            secao_primaria, secao_secundaria = secao.split('.')
            if secao_primaria not in dicio:
                dicio[secao_primaria] = {}
            dicio[secao_primaria][secao_secundaria] = valor
            
        match = re.search(r'#---+', linhas, re.MULTILINE)
        if match:
            tabela_inicio = match.end()
            tabela_linhas = linhas[tabela_inicio:].strip().split('\n')
            valores_tabela = []
            for linha in tabela_linhas:
                if re.match(r'(\s+\d+\.\d+\s+){2,4}', linha):
                    valores_tabela.append([float(valor) for valor in linha.split()])
        else:
            # Se não encontrar o início da tabela, definir valores_tabela como None
            valores_tabela = None
    
        with open(filepath, "r") as arquivo:
            texto = arquivo.read()
    
        # Encontrar nomes das colunas
        nomes_colunas = re.findall(r"# Column\.\d+: (\w+)", texto)
        # Separar as linhas de dados e criar dicionário associando cada valor ao seu respectivo nome de coluna
        dados = {}
        linhas = texto.splitlines()
        match = re.search(r'#---+', texto, re.MULTILINE)
        #print('Match', match) #DEBUG
        for nome in nomes_colunas:
            dados[nome] = []
        if match:
            inicio = match.end()
            #print('inicio',inicio) #DEBUG
            tab_linhas = texto[inicio:].strip().split('\n')
            #print('tab_linhas', tab_linhas) #DEBUG
            for linha in tab_linhas:
                #print('linha',linha) #DEBUG
                if re.match(r'((\s+)?(-)?\d+\.\d+\s+){1,3}(-)?\d+\.\d+(\s+)?', linha):
                    valores = linha.split()
                    #print('valores', valores) #DEBUG
                    for nome,valor in zip(nomes_colunas,valores):
                        dados[nome].append(float(valor))

#TRY TO GET THE EXPERIMENT DATA AND TAG THE EXPERIMENT TYPE BASED ON COLUMN NAMES
        #print(dados)
        try:
            energy = np.array(dados["energy"])
        except KeyError:
            try:
                energy = np.array(dados["energy(ev)"])
            except:
                energy = None
        try:
            i0 = np.array(dados["i0"])
        except KeyError:
            i0 = None
        try:
            itrans = np.array(dados["itrans"])
            experiment_type = "Transmission"
        except KeyError:
            itrans = None
        try:
            irefer = np.array(dados["irefer"])
        except KeyError:
            irefer = None
        try:
            fluorescence_emission = np.array(dados["fluorescence_emission"])
            experiment_type = "Fluorescence"
        except KeyError:
            fluorescence_emission = None
        try:
            raw_data = np.array(dados["raw"])
            experiment_type = "Transmission Raw"
        except KeyError:
            try:
                raw_data = np.array(dados["transmission"])
                experiment_type = "Transmission Raw"
            except KeyError:
                raw_data = None

        try:
            normalized_adsorbance = np.array(dados["normalized_absorbance"])
            experiment_type = "Normalized Transmission"
        except KeyError:
            normalized_adsorbance = None

#CHECK EXPERIMENT TYPE
        #print('experiment_type', experiment_type)
        #print('keys', nomes_colunas)
        if experiment_type == "Transmission":
            return experiment_type, energy, i0, itrans, np.log(i0/itrans)
        elif experiment_type == "Transmission Raw":
            return experiment_type, energy, raw_data
        elif experiment_type == "Normalized Transmission":
            return experiment_type, energy, normalized_adsorbance
        elif experiment_type == "Fluorescence":
            return experiment_type, energy, fluorescence_emission
        else:
            print("ERROR, EXPERIMENT TYPE NOT DETECTED")

# OLD STUFF RETURN A DATAFRAME, DELETE AFTER TESTS
#        df = pd.DataFrame()
#        df['Energy'] = energy
#        df['i0'] = i0
#        df['itrans'] = itrans
#        df['Absorption'] = np.log(df['i0'].values/df['itrans'].values)
#        return df

class XASNormalization():
    def find_E0(energies, absorption):
        """Finds the absorption edge E0

        Keyword arguments:
        energy -- the array of energies from the experiment.
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
        for point in range(range_search_E0_start,range_search_E0_end):
            if dydx[point] > deriv_max:
                if dydx2[point] > 0:
                    deriv_max = dydx[point]
                    point_max = point

#DEBUG STUFF, SHOULD BE DELETED AT SOME POINT
        #print('E0 found at', point_max, 'with derivative', dydx[point_max], 'and second derivative', dydx2[point_max])            
        # plt.plot(x, y, label='Input Data', color='black') #DEBUG
        # plt.plot(x_new, dydx, color='red', label='Derivative') #DEBUG
        # plt.plot(x_new2, dydx2, color='blue', label='Second derivative') #DEBUG
        # plt.scatter(x[point_max],y[point_max], color='black', label='E0')
        # plt.xlim(x[range_search_E0_start],x[range_search_E0_end])
        # #plt.ylim(-0.05,0.05)
        # plt.legend()
        # plt.show()
        
        E0x = x[point_max]
        E0y = y[point_max]
    
        return E0x, E0y
    
    def fit_pre_edge(energies_array, absorption_array, start_pre_edge_x_index, end_pre_edge_x_index):
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

    def xas_type(energy):
        """Test if the energy range is from a XANES or EXAFS experiment.
        XANES range is up to 400 eV.

        Keyword arguments:
        energy -- the array of energies from the experiment

        Returns:
        A string with the name of experiment type, either XANES or EXAFS
        """
        maximum = max(energy)
        minimum = min(energy)
        energy_range = float(maximum) - float(minimum)
        if energy_range < 400:
            return "XANES"
        else:
            return "EXAFS"  

    def poly2(x, a, b, c):
        """Defining the quadratic function to perform the EXAFS normalization."""
        polinomio = a*x**2 + b*x + c
        return polinomio

    def EXAFS_normalization(energies_array, absorption_array, debug=False):
        
        #Find E0
        E0x, E0y = XASNormalization.find_E0(energies_array, absorption_array)
        
        #Get E0 index
        E0x_index = np.where(energies_array==E0x)[0]
        
        #Get the ranges where the pre-edge is located
        start_to_E0 = E0x - energies_array[0]
        start_to_E0_point_spacing = E0x_index/start_to_E0
        if debug == True:
            print('Data starts at', energies_array[0], 'eV.')
            print('Data ends at', energies_array[-1], 'eV.')
            print('E0 is at', E0x, 'eV')
            print('Interval from start to E0 is', start_to_E0, 'with spacing', start_to_E0_point_spacing, 'eV')
    
        #Getting the starting energy for the pre-edge fit and its index in the array
        start_pre_edge_x = energies_array[0]
        start_pre_edge_x_index = np.where(energies_array==start_pre_edge_x)[0][0]
        if debug == True:
            print('Starting point for the pre-edge fit is', start_pre_edge_x, 'eV, point number:', start_pre_edge_x_index)
    
        #Here I am defining the end of the pre-edge region, with 30 points before E0
        pre_edge_end_fit_energy_from_E0 = 30
        end_pre_edge_x = start_pre_edge_x + (start_to_E0 - pre_edge_end_fit_energy_from_E0)
        idx = np.searchsorted(energies_array, end_pre_edge_x, side="left")
        #I do not remember what this test does, need to check
        if idx > 0 and (idx == len(energies_array) or math.fabs(end_pre_edge_x - energies_array[idx-1]) < math.fabs(end_pre_edge_x - energies_array[idx])):
            end_pre_edge_x_index = idx
        else:
            end_pre_edge_x_index = idx
    
        if debug == True:
            print('Ending point for the pre-edge fit is', energies_array[end_pre_edge_x_index], 'eV, point number:', end_pre_edge_x_index)
            print('E0 to end point of fitting distance is,', E0x - energies_array[end_pre_edge_x_index], 'eV')
    
        #Normalize the pre-edge with the obtained ranges.
        linear_fit_pre_edge = XASNormalization.fit_pre_edge(energies_array, absorption_array, start_pre_edge_x_index, end_pre_edge_x_index)
    
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
        
        if debug == True:
            print('EXAFS FIT STARTING POINT WILL BE FROM', start_exafs, '(', energies_array[start_exafs] , 'eV )', 
                  'to point', end_exafs, '(', energies_array[end_exafs], 'eV )')
    
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
            model_poly2 = Model(XASNormalization.poly2)
            params = model_poly2.make_params(a=1, b=1, c=1)
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

        # #DEBUG Plots
        if debug == True:
            reference_line_x = range(int(energies_array[0]), int(energies_array[-1])) #Reference line
            reference_line_y = [1] * len(reference_line_x) #Reference line
            plt.plot(reference_line_x, reference_line_y, color='black', linestyle='dotted')
            plt.plot(energies_array, absorption_array, label='Input Spectra', color='black')
            plt.plot(energies_array, linear_fit_pre_edge, label='Fit Pre-Edge', color='red', linestyle='dashed')
            plt.scatter(start_pre_edge_x, absorption_array[start_pre_edge_x_index], color='black', label='Background start')
            plt.scatter(energies_array[end_pre_edge_x_index], absorption_array[end_pre_edge_x_index], color='red', label='Background end')
            plt.scatter(E0x, absorption_array[E0x_index], label='E0', color='blue')
            plt.scatter(energies_array[npt_min], absorption_array[npt_min], label='EXAFS start', color='green')
            plt.scatter(energies_array[last_point_to_fit_exafs], absorption_array[last_point_to_fit_exafs], label='EXAFS end', color='gray')
            plt.plot(energies_array, predicted_exafs, label='EXAFS Fit', color='blue', linestyle='dashed')
            plt.plot(energies_array, normalized_spectra, label='Normalized Spectra', color='blue', linestyle='solid')
            plt.plot(energies_array, flattened_curve, label='Flattened Spectra', color='gray', linestyle='solid')
            plt.xlim(int(energies_array[0]), int(energies_array[-1]))
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.show()

        return energies_array, flattened_curve

    def XANES_normalization(energies_array, absorption_array, debug=False):
        """Recebe um dataset com dados de energia e absorbância e normaliza o dataset."""

        #Call the E0 calculation function and map it to the input spectra.
        E0x, E0y = XASNormalization.find_E0(energies_array, absorption_array)
        E0x_index = np.where(energies_array==E0x)[0]

        #Define the initial and final points of the linear fit of the background.
        ##This definition seemed good for the the tested spectra, might be changed after more tests are performed.
        start_pre_edge_x_index = int(E0x_index/15)
        end_pre_edge_x_index = int(E0x_index/3)

        #Linear model for the pre-edge.
        data_x = energies_array[start_pre_edge_x_index:end_pre_edge_x_index] #background.iloc[:, 0].values
        data_y = absorption_array[start_pre_edge_x_index:end_pre_edge_x_index] #background.iloc[:, 1].values

        #Perform the linear fit.
        linear_fit_pre_edge = XASNormalization.fit_pre_edge(energies_array, absorption_array, start_pre_edge_x_index, end_pre_edge_x_index)

        #Perform the normalization
        edge_step = absorption_array[E0x_index] - linear_fit_pre_edge[E0x_index]
        if debug:
            print('energies_array[E0x_index]', energies_array[E0x_index])
            print('linear_fit_pre_edge[E0x_index]', linear_fit_pre_edge[E0x_index])
            print('edge_step', edge_step)
            print('E0x', E0x, 'E0y', E0y)
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

        if debug==True:
            #For visualizations and DEBUG.
            plt.plot(energies_array,fit_final,label='Normalized Spectra', color='red')
            plt.legend(fontsize=16)
            plt.xlabel("Energy (eV)")
            plt.ylabel("Absorption") #, fontproperties=font_manager.FontProperties(size=12))
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            plt.scatter(E0x,float(fit_final[E0x_index]), label='E0 value')
            plt.scatter(energies_array[npt_min],fit_final[npt_min], color='b', label='Optimal Point')
            plt.legend()
            plt.show()

        return energies_array, fit_final

