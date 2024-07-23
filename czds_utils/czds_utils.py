import re
import os
import pandas as pd
import numpy as np
from numpy import diff
from scipy.stats import linregress


def czd_read_xdi(filepath):
    filepath = "CeO2.xdi"
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
    for nome in nomes_colunas:
        dados[nome] = []
    if match:
        inicio = match.end()
        tab_linhas = texto[inicio:].strip().split('\n')
        for linha in tab_linhas:
            if re.match(r'(\s+\d+\.\d+\s+){2,4}', linha):
                valores = linha.split()
                for nome,valor in zip(nomes_colunas,valores):
                    dados[nome].append(float(valor))
    try:
        energy = dados["energy"]
    except KeyError:
        energy = None
    try:
        i0 = dados["i0"]
    except KeyError:
        i0 = None
    try:
        itrans = dados["itrans"]
    except KeyError:
        itrans = None
    try:
        irefer = dados["irefer"]
    except KeyError:
        irefer = None

    df = pd.DataFrame()
    df['Energy'] = energy
    df['i0'] = i0
    df['itrans'] = itrans
    df['Absorption'] = np.log(df['i0'].values/df['itrans'].values)
    return df

def find_E0(energies, absorption):
    x = energies
    y =  absorption
    dydx = diff(y)/diff(x)
    x_new = x[:-1]
    dydx2 = diff(dydx)/diff(x_new)
    x_new2 = x_new[:-1]

    #Test for the largest derivative, only if the second derivative is positive.
    #Searching only up to half the graph. Should be to the first peak, but for now didn't implemented it. 
    range_search_E0_start=0
    range_search_E0_end=int(len(energies)/2)
    deriv_max = 0
    for point in range(range_search_E0_start,range_search_E0_end):
        if dydx[point] > deriv_max:
            if dydx2[point] > 0:
                deriv_max = dydx[point]
                point_max = point

    print('E0 found at', point_max, 'with derivative', dydx[point_max], 'and second derivative', dydx2[point_max])            

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
    print('indexes start finish',start_pre_edge_x_index,end_pre_edge_x_index)
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
