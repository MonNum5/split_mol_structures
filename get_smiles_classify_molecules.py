# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 08:28:24 2020

@author: hall_ce
"""
import pandas as pd
import pubchempy as pcp
import numpy as np

import copy
import os

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw

xlsx_path = 'Thermodataengine_names.xlsx'

xl = pd.ExcelFile(xlsx_path)

# reduce molecule names to unique names
for i, fam in enumerate(list(xl.sheet_names)):  # see all sheet names

    df_i = xl.parse(fam)
    if i ==0:
        df = df_i
    df = df.append(df_i)
 
df_dum = df.groupby(1).first()
df_dum['name'] = df_dum.index
df_dum['formula'] = df_dum[0]
df_dum['alternative_names'] = df_dum[2]
df_dum['family'] = df_dum[4]
df_dum['smiles'] = np.empty((df_dum.shape[0],1))*np.nan
df_dum['smarts'] = np.empty((df_dum.shape[0],1))*np.nan
df_dum = df_dum.drop(columns = ['Unnamed: 0',3])
df = df_dum[['name','alternative_names','formula','family','smiles','smarts']]
df = df.reset_index(drop=True)

# get smiles key, draw picture, classifiy hydrocarbon family
     
for i in df.index:
            
    # reformat alternative names
    alt_names = copy.deepcopy(df['alternative_names'][i])
    if type(alt_names) != type(np.nan):
        alt_names = alt_names[1:-1]
        alt_names = alt_names.split('; ')
        df['alternative_names'][i] = alt_names
    
    # search compound by name on pubchempy
    name = df['name'][i]
    c = pcp.get_compounds(name, 'name')
    
    # if c is none check using alt names
    if not c and type(alt_names) != type(np.nan):
        for alt_name in df['alternative_names'][i]:
            c = pcp.get_compounds(alt_name, 'name')
            df['name'][i] = alt_name
            if c:          
                break
        if not c:
            print(f"No match found for {name}")
    
    if c:
        if len(c)>1:           
            smiles = [comp.isomeric_smiles for comp in c]
            print(i,name, 'Multiple machtes', smiles)
        else:
            smiles = c[0].isomeric_smiles
            print(i,name, smiles)
    else:
        smiles = np.nan
    
    df['smiles'][i] = smiles   
    
    # convert to smarts
    if type(smiles) == str:
        mol = Chem.MolFromSmiles(smiles)
        smarts = Chem.MolToSmarts(mol, isomericSmiles=True)
        df['smarts'][i] = smarts
 
df.to_excel('Thermodataengine_names_smiles.xlsx')
ddd   

# New xlsx file to write to
writer = pd.ExcelWriter('Thermodataengine_names_and_structures.xlsx', engine='xlsxwriter')

if not os.path.isdir('formula_images'):
    os.mkdir('formula_images')
    
if not os.path.isdir(os.path.join('formula_images',fam)):
    os.mkdir(os.path.join('formula_images',fam))
        
df.to_excel(writer, fam)

workbook  = writer.book
worksheet = writer.sheets[fam] 
worksheet.set_default_row(80)

for i in df.index:  
    
    smiles = df['smiles'][i]
    name = df['name'][i]
    #draw molecule
    
    if type(smiles) != type(np.nan):
        print(smiles)
        smiles = smiles
        mol = Chem.MolFromSmiles(smiles)
        mol
        Draw.MolToFile(mol,os.path.join('formula_images',fam,f'{name}.png'), size=(100,100))
        worksheet.insert_image(f'I{i+2}', os.path.join('formula_images',fam,f'{name}.png'))
        
writer.save()


'''
#to smarts
mol = Chem.MolFromSmiles('CCCCC1=CC=C(C=C1)C=C')
mol  
sma = Chem.MolToSmarts(mol, isomericSmiles=True)
sma = Chem.MolFromSmiles(sma)
sma
'''