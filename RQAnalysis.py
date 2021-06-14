# Import all the important libraries
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
sns.set()
from Utility import*


def rec_plot_binary(s, eps=0.10, steps=10):
    d = pdist(s)
    d = np.floor(d/eps)
    d[d>steps] = steps
    Z = squareform(d)
    for i in range(len(Z)):
        for j in range(len(Z)):
            if (Z[i][j]>0):
                Z[i][j] =1
            else:
                Z[i][j] =0
    return Z


# Definition of all functions of Quantification Analysis
"""-------------------------------Quantification of Recurrence Plots------------------------------"""
#######################################################################################################
"""-------------------------------RECURRENCE RATE(RR)---------------------------------------------"""

def total_number_of_recurrence(recurrence_matrix):
    total_number_of_recurrence_points = 0
    N = len(recurrence_matrix)
    for i in range(N):
        for j in range(N):
            if (recurrence_matrix[i][j]==1):
                total_number_of_recurrence_points = total_number_of_recurrence_points + 1
    print("Total number of recurrence points-->",total_number_of_recurrence_points)
    return(total_number_of_recurrence_points)
    
def recurrence_rate(recurrence_matrix):
    """ Recurrence rate (RR).--> 1/N^2*(sum (i->N,j->N) R(i,j)) 
    
    Corresponds to the correlation sum """
    total_number_of_recurrence_points = total_number_of_recurrence(recurrence_matrix)
    N = len(recurrence_matrix)
    RR = total_number_of_recurrence_points/(N**2)
    return (RR) 


###########################################################################################
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
    
def determinism(recurrence_matrix,l_min):
    """ Determinism (DET). --> Percentage of Recurrence Points that form diagonal lines
        DET =  Sum(l=lmin-->N,lP(l))/Sum(l=1-->N,lP(l))"""
    P = frequency_distribution_diagonal_lines(recurrence_matrix)
    N = len(recurrence_matrix)
    num =0
    deno =0
    for l in range(l_min,N+1):
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
    for i in range(N):
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

def laminarity(recurrence_matrix):
    """ Laminarity (LAM). --> Percentage of Recurrence Points that form vertical lines
        LAM =  Sum(v=vmin-->N,vP(v))/Sum(v=1-->N,vP(v))"""
    P = frequency_distribution_vertical_lines(recurrence_matrix)
    N = len(recurrence_matrix)
    num =0
    deno =0
    for v in range(2,N):
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
    
"""---------------------------------------RATIO(DET/RR)-------------------------------------------"""    
    
def ratio_determinism_recurrence_rate(recurrence_matrix):
        """ Ratio determinism / recurrence rate (DET/RR). """
        try:
            Ratio = recurrence_rate(recurrence_matrix)/determinism(recurrence_matrix)
            return(Ratio)
        except:
            print("ZeroDivisionError")
#######################################################################################################

""" --------------------------------------AVERAGE DIAGONAL LINE ------------------------------------"""   
def average_diagonal_line(recurrence_matrix):
    """ Average diagonal line length (L).
        DET =  Sum(l=lmin-->N,lP(l))/Sum(l=1-->N,P(l))"""
    P = frequency_distribution_diagonal_lines(recurrence_matrix)
    num =0
    deno =0
    N = len(recurrence_matrix)
    for l in range(2,N+1):
        try:
            num = num + l*P[l]
        except:
            continue
    for l in range(2,N+1):
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
"""-------------------------------TRAPPING TIME (TT)--------------------------------------------"""
        
def trapping_time(recurrence_matrix):
    """ Trapping Time (TT)---> The average length of the vertical lines.
        DET =  Sum(v=vmin-->N,vP(v))/Sum(l=1-->N,P(v))"""
    P = frequency_distribution_vertical_lines(recurrence_matrix)
    num =0
    deno =0
    N = len(recurrence_matrix)
    for v in range(2,N+1):
        try:
            num = num + v*P[v]
        except:
            continue
    for v in range(2,N+1):
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
"""--------------------------------LONGEST DIAGONAL LINE (L_max)-----------------------------"""
    
def set_longest_diagonal_line(recurrence_matrix):
    """ Set longest diagonal line length (L_max). """
    P = frequency_distribution_diagonal_lines(recurrence_matrix)
    # print(set(P.keys()))
    P_nonzero = remove_zero_from_dict(P)
    L_max = second_largest(list(P_nonzero.keys()))
    return(L_max)
####################################################################################################
"""--------------------------------LONGEST VERTICAL LINE (V_max)-----------------------------"""
    

def set_longest_vertical_line(recurrence_matrix):
    """ Set longest vertical line length (V_max). """
    P = frequency_distribution_vertical_lines(recurrence_matrix)
    # print(set(P.keys()))
    P_nonzero = remove_zero_from_dict(P)
    V_max = second_largest(list(P_nonzero.keys()))
    return(V_max)
####################################################################################################
"""--------------------------------DIVERGENCE (DIV)-----------------------------"""
            
def divergence(recurrence_matrix):
        """ Divergence (DIV). """
        L_max = set_longest_diagonal_line(recurrence_matrix)
        DIV = 1/L_max
        return DIV
####################################################################################################
"""-------------------------------- ENTROPY(ENTR)-----------------------------"""
    
def set_entropy_diagonal_lines(recurrence_matrix):
    """ Set entropy of diagonal lines (ENTR). """
    P = frequency_distribution_diagonal_lines(recurrence_matrix)
    ENTR = 0
    N = len(recurrence_matrix)
    for l in range(2,N+1):
        try:
            ENTR = ENTR - P[l]*np.log(P[l])
        except:
            continue
    return(ENTR)

def Quantification_of_Recurrence_Plot(recurrence_matrix):
    print("RR -->",recurrence_rate(recurrence_matrix))
    print("DET -->",determinism(recurrence_matrix,l_min=2))
    print("LAM -->",laminarity(recurrence_matrix))
    # print("RATIO -->",ratio_determinism_recurrence_rate(recurrence_matrix))
    # print("L-->",average_diagonal_line(recurrence_matrix))
    # print("TT -->",trapping_time(recurrence_matrix))
    print("L_max-->",set_longest_diagonal_line(recurrence_matrix))
    print("V_max-->",set_longest_vertical_line(recurrence_matrix))
    # print("DIV -->",divergence(recurrence_matrix))
    # print("ENTR -->", set_entropy_diagonal_lines(recurrence_matrix))

def generate_plot(recurrence_matrix):
    plt.figure(figsize = (6,20))
    plt.imshow(recurrence_matrix,cmap="gray")
    plt.show()
