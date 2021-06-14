import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import*
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pylab
from RQAnalysis import*
from Utility import*
from DynamicSystems import*
from itertools import permutations
from matplotlib import pyplot as plt


########################### Defining Functions for Symbollic Recurrence Plots ###########################################

def perm_symbols(embedding):
    return(list(permutations(range(embedding))))

def symbolization(embedding):
    l = perm_symbols(embedding)
    d = dict()
    for i in range(len(l)):
        d[l[i]] = i+1
    return(d)

def symbolic_function(x,embedding):
    sort_index = tuple(np.argsort(np.array(x)))
    dict_symbols = symbolization(embedding)
    symbols = perm_symbols(embedding)
    for i in symbols:
        if (i==sort_index):
            function_result = dict_symbols[i]
            # print(str(x) + "--->" + str(function_result))
    return(function_result)

def symbols_for_all_time_series(timeseries,embedding):
    n = len(timeseries)
    m = embedding
    d = dict()
    for i in range(n-m+1):
        # print(timeseries[i:i+m])
        d["x_"+str(i+1)]= symbolic_function(timeseries[i:i+m],embedding)
        # print(d)
    return(d)
def plot_matrix(d):
    values = list(d.values())
    matrix = np.zeros((len(d),len(d)))
    for i in range(len(values)):
        for j in range(len(values)):
            if values[i]==values[j]:
                matrix[i][j] = values[i]
    return matrix    
def plot_SRP(matrix):
    matrix[matrix ==0] = np.nan
    plt.imshow(matrix,cmap="cool",interpolation="none",vmin=0)
    plt.show()   
###############################################################################################################
def SRR(mat):
    total = 0
    N = len(mat)
    for i in range(N):
        for j in range(N):
            if(mat[i][j]>0):
                total = total+1
    SRRate = total/N**2
    return SRRate

################################################################################################################
def SRR_Pi(mat):
    SRR_Pi = dict()
    SRRate = SRR(mat) 
    symbols_list = list(np.unique(mat))
    N = len(mat)
    for symbol in symbols_list:
        total = 0
        for i in range(N):
            for j in range(N):
                if mat[i][j]==symbol:
                    total = total+1
        SRR_Pi[symbol] = round(total/((N**2)*SRRate),2)
    return SRR_Pi
###############################################################################################################

#diagonal function is the each diagonal.
"""------------------------------------DETERMINISM (DET)--------------------------------"""
def diagonal(matrix,offset =0):
    mat = np.array(matrix)
    ret = mat.diagonal(offset)
    return(list(ret))

def frequency_distribution_diagonal_lines(recurrence_matrix):
    N = len(recurrence_matrix)
    frequency_array = np.zeros(N+1)
    for i in range(-N+1,N):
        dia_element = diagonal(recurrence_matrix,i)
        int_dia_element = list(np.array(dia_element).astype("int"))
        # print(int_dia_element)
        str_dia_element = ''.join(map(str, int_dia_element))
        diagonal_len_list = list(map(len, str_dia_element.split("0")))
        for j in range(len(diagonal_len_list)):
            frequency_array[diagonal_len_list[j]] = frequency_array[diagonal_len_list[j]]+1
        freq_dictionary = freq_dict(frequency_array)
        # freq_dictionary = change_dict_by_one(freq_dictionary)
    return(freq_dictionary)
    
def determinism_SRP(recurrence_matrix,q0,q1):
    """ Determinism (DET). --> Percentage of Recurrence Points that form diagonal lines
        DET =  Sum(l=lmin-->N,lP(l))/Sum(l=1-->N,lP(l))"""
    P = frequency_distribution_diagonal_lines(recurrence_matrix)
    N = len(recurrence_matrix)
    num =0
    deno =0
    for l in range(q0,q1):
        try:
            num = num + l*P[l]
        except:
            continue
    for l in range(1,N+1):
        try:
            deno = deno + l*P[l]
        except:
            continue
    det = num/deno
    return det
####################################################################################################
#each column element extracting function 
"""-----------------------------------LAMINARITY ----> LAM----------------------------------------- """
import collections
def vertical(matrix,column=0):
    l = matrix[:,column]
    return list(l)

def frequency_distribution_vertical_lines(recurrence_matrix):
    N = len(recurrence_matrix)
    frequency_array = np.zeros(N+1)
    for i in range(-N+1,N):
        vert_element = vertical(recurrence_matrix,i)
        int_vert_element = list(np.array(vert_element).astype("int"))
        # print(int_dia_element)
        str_vert_element = ''.join(map(str, int_vert_element))
        vertical_len_list = list(map(len, str_vert_element.split("0")))
        for j in range(len(vertical_len_list)):
            frequency_array[vertical_len_list[j]] = frequency_array[vertical_len_list[j]]+1
        freq_dictionary = freq_dict(frequency_array)
        # freq_dictionary = change_dict_by_one(freq_dictionary)
    return(freq_dictionary)

def laminarity_SRP(recurrence_matrix,q0,q1):
    """ Laminarity (LAM). --> Percentage of Recurrence Points that form vertical lines
        LAM =  Sum(v=vmin-->N,vP(v))/Sum(v=1-->N,vP(v))"""
    P = frequency_distribution_vertical_lines(recurrence_matrix)
    N = len(recurrence_matrix)
    num =0
    deno =0
    for v in range(q0,q1):
        try:
            num = num + v*P[v]
        except:
            continue
    for v in range(1,N):
        try:
            deno = deno + v*P[v]
        except:
            continue
    lam = num/deno
    return lam
"""-----------------------------------------------------------------------------------------------"""
###################################################################################################
""" --------------------------------------AVERAGE DIAGONAL LINE ------------------------------------"""   
def average_diagonal_line(recurrence_matrix,q0,q1):
    """ Average diagonal line length (L).
        DET =  Sum(l=lmin-->N,lP(l))/Sum(l=1-->N,P(l))"""
    P = frequency_distribution_diagonal_lines(recurrence_matrix)
    num =0
    deno =0
    N = len(recurrence_matrix)
    for l in range(q0,q1):
        try:
            num = num + l*P[l]
        except:
            continue
    for l in range(q0,q1):
        try:
            deno = deno + P[l]
        except:
            continue
    try:
        L = num/deno
        return L
    except:
        print("ZeroDivisionError")
#####################################################################################################

"""-------------------------------AVERAGE VERTICAL TIME (TT)--------------------------------------------"""
        
def average_vertical_line(recurrence_matrix,q0,q1):
    """ Average Vertical Line (AVL)---> The average length of the vertical lines.
        DET =  Sum(v=vmin-->N,vP(v))/Sum(l=1-->N,P(v))"""
    P = frequency_distribution_vertical_lines(recurrence_matrix)
    num =0
    deno =0
    N = len(recurrence_matrix)
    for v in range(q0,q1):
        try:
            num = num + v*P[v]
        except:
            continue
    for v in range(q0,q1):
        try:
            deno = deno + P[v]
        except:
            continue
    try:
        TT = num/deno
        return TT
    except:
        print("ZeroDivisionError")


####################################################################################################
def Quantification_of_Recurrence_Plot_SRP(matrix):
    matrix_num = convert_nan_to_num(matrix)
    print("SRR -->",SRR(matrix_num))
    print("SRR_Pi-->",SRR_Pi(matrix_num))
    matrix_num_b = convert_into_binary(matrix_num)
    N = len(matrix_num_b)
    print("DET -->",determinism_SRP(matrix_num_b,q0=2,q1=N))
    print("LAM -->",laminarity_SRP(matrix_num_b,q0=2,q1=N))
    # print("RATIO -->",ratio_determinism_recurrence_rate(recurrence_matrix))
    # print("L-->",average_diagonal_line(recurrence_matrix,q0=2,q1=len(recurrence_matrix)))
    # print("TT -->",trapping_time(recurrence_matrix,q0=2,q1=len(recurrence_matrix)))
    # print("L_max-->",set_longest_diagonal_line(recurrence_matrix))
    # print("V_max-->",set_longest_vertical_line(recurrence_matrix))
    # # print("DIV -->",divergence(recurrence_matrix))
    # # print("ENTR -->", set_entropy_diagonal_lines(recurrence_matrix))
