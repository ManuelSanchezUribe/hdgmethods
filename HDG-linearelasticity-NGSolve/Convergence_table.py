#from pandas import DataFrame
import pandas as pd
import numpy as np
from math import log10 as log

class Convergence_table:
    def __init__(self, h, Ndof, Err_list):
        self.Err_list= Err_list
        self.h = h
        self.Ndof = Ndof
    def create(self, errors):
        self.errors = errors
        self.orders = self.Compute_orders(errors)
        self.table_dict = {}
        if self.h is not None:
            self.table_dict['h'] = self.h
        if self.Ndof is not None:
            self.table_dict['Ndof'] = self.Ndof
        for i, key in enumerate(self.Err_list):
            self.table_dict['Error '+key] = self.errors[i]
            self.table_dict['Order '+key] = self.orders[i]
    def Compute_orders(self, errors):
        # Convert
        self.errors = [[item[i] for item in errors] for i in range(len(errors[0]))]
        orders = []
        if self.h is not None:
            for j in range( len(self.Err_list)):
                col = [] 
                i = 0
                while i < len(self.errors[j]):
                    if i>0:
                        if self.errors[j][i] == 0:
                            col.append(0)
                        else:
                            col.append(log(self.errors[j][i-1]/self.errors[j][i])/log(self.h[i-1]/self.h[i]))
                        #
                    else:
                        col.append(0)
                    i += 1
                orders.append(col)
        return orders

    def Show_in_Jupyter(self):
        pd.set_option('display.width', 120)
        pd.options.display.float_format = '{:.4e}'.format
        df = pd.DataFrame(self.table_dict)
        convert_dict = {}
        if self.h is not None:
            convert_dict['h'] = float
        if self.Ndof is not None:
            convert_dict['Ndof'] = int
        for j in range(len(self.Err_list)):
            convert_dict['Error '+self.Err_list[j]] = float
            #convert_dict['order '+Err_list[j]] = float
            df['Order '+self.Err_list[j]] = df['Order '+self.Err_list[j]].apply(lambda x: '{:.4f}'.format(x))
        df = df.astype(convert_dict)
        return df
