# -*- coding: utf-8 -*-
"""Young_Modulus.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1tSMFjI5d-TaROU5CTcF9pxFgwkktVfNK

OBJECTIVES:
- Read .dat file
- Plot Stress x Strain
- Be able to make a linear numerical extrapolation with the data
- Based on that, extract the Young Modulus

# NODE DATA
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from google.colab import drive
drive.mount('/content/drive')

## First, define path to data
Stress_Strain_node = '/content/drive/MyDrive/Projetos_pessoais/Tese_Polimi/Stress_Strain_node.dat' #In this case, I'm usign a file from my Drive, but you may also use this by specifying any directory from your PC
Header_rows_node = 0 ## If Abaqus .dat file -> 0

## Reading file
data_node = pd.read_csv(Stress_Strain_node,sep='\s+',header=Header_rows_node)
data_node = pd.DataFrame(data_node)
print(data_node) #Check if data was read correctly

## Ploting results
Strain_node = data_node['X']
Stress_node = data_node['Stress_Strain']
plt.plot(Strain_node, Stress_node,'r')
plt.xlim([data_node['X'].min(), data_node['X'].max()])
plt.ylim([data_node['Stress_Strain'].min(), data_node['Stress_Strain'].max()])
plt.grid()
plt.title("Aluminium Coupon Stress Strain curve")
plt.xlabel("Strain")
plt.ylabel("Stress [GPa]")
plt.show()

## Now, let's try to extract the Young Modulus by interpolating the data into a linear function
linear_node = data_node.interpolate(method='linear')
rows_node = len(linear_node)
E_node = linear_node.iloc[rows_node-2,1]/linear_node.iloc[rows_node-2,0]
print('Young Modulus (Node) =', E_node, 'GPa')

"""# ELEMENT DATA"""

## First, define path to data
Stress_Strain_element = '/content/drive/MyDrive/Projetos_pessoais/Tese_Polimi/Stress_Strain_element.dat'
Header_rows_element = 0 ## If Abaqus .dat file -> 0

## Reading file
data_element = pd.read_csv(Stress_Strain_element,sep='\s+',header=Header_rows_element)
data_element = pd.DataFrame(data_element)
print(data_element) #Check if data was read correctly

## Ploting results
Strain_element = data_element['X']
Stress_element = data_element['Stress_Strain']
plt.plot(Strain_element, Stress_element,'r')
plt.xlim([data_element['X'].min(), data_element['X'].max()])
plt.ylim([data_element['Stress_Strain'].min(), data_element['Stress_Strain'].max()])
plt.grid()
plt.title("Aluminium Coupon Stress Strain curve")
plt.xlabel("Strain")
plt.ylabel("Stress [GPa]")
plt.show()

## Now, let's try to extract the Young Modulus by interpolating the data into a linear function
linear_element = data_element.interpolate(method='linear')
rows_element = len(linear_element)
E_element = linear_element.iloc[rows_element-2,1]/linear_element.iloc[rows_element-2,0]
print('Young Modulus (Element) =', E_element, 'GPa')

"""# COMPARISON BETWEEN RESULTS

The inserted material data was E = 71.7 GPa

Comparing with both results, we get:
"""

E = 71.7
# Comparing with Node data
error_node = abs(100*((E_node/E)-1))

# Comparing with Element data
error_element = abs(100*((E_element/E)-1))

print('The error of the Young Modulus obtained with data from node is of', error_node,'%')
print('The error of the Young Modulus obtained with data from element is of', error_element,'%')