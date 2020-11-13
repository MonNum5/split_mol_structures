# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 09:41:44 2020

@author: hall_ce
"""

import pandas as pd
import pubchempy as pcp
import numpy as np
import requests
import re



def get_headings(dict_, key):
    '''
    returns list of values for key
    '''
    return [i[key] for i in dict_]



#check if pyhsical and chemical values are available

property_dict = {
                 'Refractive Index':{'name':'refractive_index', 'unit':'-'},
                 'Surface Tension':{'name':'surface_tension','unit':'mN/m'},
                 'Autoignition Temperature':{'name':'autoignition_temperature','unit':'C'},
                 'Viscosity':{'name':'viscosity_dynamic','unit':'mPa*s'},
                 'Heat of Combustion':{'name':'net_heat_of_combustion','unit':'MJ/kg'},
                 'Heat of Vaporization':{'name':'heat_of_vaporization', 'unit':'kJ/kg'},
                 'Density':{'name':'density', 'unit':'kg/m3'},
                 'Vapor Pressure':{'name':'vapor_pressure', 'unit':'hPa'},
                 'Boiling Point':{'name':'boiling_temperature', 'unit':'C'},
                 'Melting Point':{'name':'freezing_point', 'unit':'C'},
                 'Flash Point':{'name':'flash_point', 'unit':'C'}
                 }

unit_list = ['K', 'F','C','mPa*s','mN/m','kg/m3','MJ/kg','hPa','kJ/kg','-']
drop_list = ['at','X10','@','to']
replace_unit = {
        'mm Hg':'mmHg',
        'g/cu cm':'g/cm3',
        'Â°C/D':'C',
        'Â°C':'C',
        'Â°F':'F',   
        'mPa.s':'mPa*s',
        'g/ml':'g/mL',
        'g/cmÂ³':'g/cm3',
        'DYNES/CM':'dyne/cm',
        'CAL/G':'cal/g',
        'PSI':'psi',
        'G/ML':'g/mL',
        'dyn/cm':'dynes/cm',
        'AT':'atm'
        }

# units
def convert_to_standard_unit(prop, unit, c, molmass):
    conversion_dict={'kg/m3':{'g/mL':1e-3*c,'g/cm3':1e-3*c},
                     'g/mL':{'kg/m3':1e3*c,'g/cm3':1*c},
                     'g/cm3':{'kg/m3':1e3*c,'g/mL':1*c},
                     'mN/m':{'dyne/cm':1*c,'dynes/cm':1*c},
                     'N/m':{'mN/m':1e-3*c},
                     'dyne/cm':{'mN/m':1*c,'dynes/cm':1*c},
                     'dynes/cm':{'mN/m':1*c,'dyne/cm':1*c},
                     'kJ/kg.K':{'J/gC':1*c,'J/gK':1*c,'kJ/(kg*K)':1*c},
                     'J/gC':{'kJ/kg.K':1*c,'J/gK':1*c,'kJ/(kg*K)':1*c},
                     'J/gK':{'J/gC':1*c,'kJ/kg.K':1*c,'kJ/(kg*K)':1*c},
                     'kJ/(kg*K)':{'J/gC':1*c,'kJ/kg.K':1*c,'J/gK':1*c},
                     'MJ/kg':{'kJ/kg':1000*c,'BTU/lb':429.923*c,'kJ/mol':c*molmass},
                     'kJ/mol':{'MJ/kg':abs(c/molmass), 'kJ/kg':abs(c/molmass*1000)},
                     'kJ/kg':{'MJ/kg':1e-3*c,'BTU/lb':0.429923*c},
                     'J/g':{'MJ/kg':1e-3*c,'BTU/lb':0.429923*c,'kJ/kg':c},
                     'BTU/lb':{'kJ/kg':2.326*c,'MJ/kg':0.002326*c},
                     'cal/g':{'MJ/kg':0.004184,'kJ/kg':4.184},
                     'cSt':{'mm2/s':1*c},
                     'cP':{'mPa*s':1*c},
                     'mm2/s':{'cSt':1*c},
                     'C':{'F':(c*9/5)+32,'K':c+273.15},
                     'F':{'C':(c-32)*5/9,'K':(c+459.67)*5/9},
                     'K':{'C':c-273.15,'F':(c*9/5)-459.67},
                     'mmHg':{'hPa':1.33322*c,'bar':0.00133322*c},
                     'hPa':{'mmHg':0.750062*c,'bar':1e-3*c},
                     'atm':{'bar':1.01325,'hPa':101.325},
                     'psi':{'hPa':68.9476*c,'bar':0.0689476*c},
					 }
    # check if unit not standard unit
    standard_unit = property_dict[prop]['unit']
    if unit != standard_unit:
        if unit in conversion_dict:
            val = conversion_dict[unit][standard_unit]
            new_unit = standard_unit
        else:
            print(f"{unit} for {prop} not conversion dict")
            val = c
            new_unit = unit
    else:
        val = c
        new_unit = unit
    return val, new_unit
    

df = pd.DataFrame([])
k = 0

# open xlsx data
df_names = pd.read_excel('Thermodataengine_names_smiles.xlsx')

for o, name in enumerate(df_names['name']):
    print(o/len(df_names['name'])*100)
    c = pcp.get_compounds(name, 'name')
    if len(c)==1:
    
        cid = c[0].cid
        response = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?response_type=display")
        dum = response.json()
        if 'Record' in dum:
            dum=dum['Record']
            section = dum['Section']
            headings = get_headings(section, 'TOCHeading')
            
            if 'Chemical and Physical Properties' in headings:
                ind = headings.index('Chemical and Physical Properties')
                chem_phy_sec = section[ind]['Section']
                headings = get_headings(chem_phy_sec, 'TOCHeading')
                
                #get comp data
                ind = headings.index('Computed Properties')
                comp_data = exp_data = chem_phy_sec[ind]['Section']
                headings_comp = get_headings(comp_data, 'TOCHeading')
                #get molar mass
                ind = headings_comp.index('Exact Mass')
                molar_mass = comp_data[ind]['Information'][0]['Value']['Number'][0]
                molar_mass_unit = comp_data[ind]['Information'][0]['Value']['Unit']
                # check if experimental values are available
                if 'Experimental Properties' in headings:
                    ind = headings.index('Experimental Properties')
                    exp_data = chem_phy_sec[ind]['Section']
                    headings = get_headings(exp_data, 'TOCHeading')
                    
                    # iterate through properties
                    for i in exp_data:
                
                        if i['TOCHeading'] in property_dict.keys():
                            for j in i['Information']:
                                #get reference
                                try:
                                    df.at[k, 'source'] = j['Reference'][0]
                                except:
                                    df.at[k, 'source'] = f"Reference Number {j['ReferenceNumber']}"
                                    
                                
                                #get value
                                if 'StringWithMarkup' in j['Value']:
                                    val_string = j['Value']['StringWithMarkup'][0]['String']
                                    
                                    #replace unwanted str tokens
                                    drop_list = ['at','X10','@','to']
                                    for drop_token in drop_list:
                                        if drop_token in val_string:
                                            val_string = val_string.replace(drop_token,'')
                                    
                                    # check if unit values have to be replaced
                                    for repl_unit in replace_unit:
                                        if repl_unit in val_string:
                                            val_string = val_string.replace(repl_unit,replace_unit[repl_unit])
                                            
                                    # insert space after digit
                                    digits = re.findall('(\d+(?:\.\d+)?)', val_string)
                                    for dig in digits:
                                        dig_ind = val_string.index(dig)
                                        if dig_ind+len(dig)+1<len(val_string) and val_string[dig_ind+len(dig)+1] !=' ':
                                            val_string = val_string[:dig_ind+len(dig)]+' '+val_string[dig_ind+len(dig):]
                    
                                    val_string_split = val_string.split(' ')
                                    val_string_split = [val_i for val_i in val_string_split if val_i not in [' ','']]
                                    
                                    val_ = False
                                    prop_name = i['TOCHeading']
                                    name_new = property_dict[prop_name]['name']
                                    df.at[k,'name'] = name
                                    df.at[k,'property'] = name_new
                                    df.at[k,'CID'] = dum['RecordNumber']
                                    df.at[k,'name_cid'] = dum['RecordTitle']                 
                                    for p in range(len(val_string_split)):
                                        try:
                                            if float(val_string_split[p]):
                                                val = val_string_split[p]
                                                unit = val_string_split[p+1]
                                                
                                                try:
                                                    val, unit = convert_to_standard_unit(prop_name, unit, float(val), molar_mass)
                                                except:
                                                    print('did not work for',prop_name, unit, val, molar_mass)
                                                    
                                                if name == 'Refractive Index':
                                                    unit = '-'
                                                elif name == 'Density' and 0 < val < 1.5 and type(unit) == type(np.nan):
                                                    val = val*1000
                                                    unit = 'kg/m3'
                                                if not val_:
                                                    df.at[k, 'value'] = val
                                                    if unit in unit_list: 
                                                        df.at[k, 'unit'] = unit
                                                    else:
                                                        df.at[k, 'untracked_unit'] = unit
                                                    val_ = True
                                                else:
                                                    df.at[k, 'dep'] = val
                                                    if unit in unit_list: 
                                                        df.at[k, 'dep_unit'] = unit
                                                    else:
                                                        df.at[k, 'untracked_dep_unit'] = unit
                                        except:
                                            print(f"Error {name}")
                                            pass
                                else:
                                    print(f"Error {name} no 'StringWithMarkup'")
                                    pass
                                        
                                k +=1
        else:
            print(f"Error {name} no 'Record'")
                        



