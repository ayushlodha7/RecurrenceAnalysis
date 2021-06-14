import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pylab
from RQAnalysis import*

sns.set()

def second_largest(array):
    second_ele = list(set(np.sort(array)))[-2]
    return(second_ele)

def freq_dict(array):
    index_array = list(range(len(array)))
    d_ = dict(zip(index_array,array))
    return (d_)

def change_dict_by_one(d):
    new_d = dict()
    d_keys = list(d.keys())
    d_values = list(d.values())
    for i in range(len(d)):
        new_d[d_keys[i]-1] = d_values[i]
    return(new_d)


# Converting the NaN values in Matrix to Num Values

def convert_nan_to_num(matrix):
    matrix_num = np.nan_to_num(matrix)
    return (matrix_num)

# Converting Colorful Matrix to Zero and One Matrix.

def convert_into_binary(matrix):
    matrix_b = np.zeros([len(matrix),len(matrix)])
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if (matrix[i][j]>0):
                matrix_b[i][j] =1
            else:
                matrix_b[i][j] =0
    return matrix_b

# Remove zero from Dictionary
def remove_zero_from_dict(dic_in):
    dic_out = {x:y for x,y in dic_in.items() if y!=0}
    return(dic_out)
