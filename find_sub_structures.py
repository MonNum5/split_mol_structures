# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 08:58:06 2020

@author: hall_ce
"""
import os
import pandas as pd

# path to smart substructures
xlsx_path = 'smart_sub_structures.xlsx'

df = pd.read_excel(xlsx_path)
ssss
writer = pd.ExcelWriter('smart_sub_structures_with_images.xlsx', engine='xlsxwriter')
        
df.to_excel(writer,'substructures')

workbook  = writer.book
worksheet = writer.sheets['substructures'] 
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

