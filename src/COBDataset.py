#!/usr/bin/python

from COBDatabase import *

class COBDataset(COBDatabase):
    def __init__(self,name,description):
        super(COBDataset,self).__init__()
        self.name = name
        self.description = description
        self.expr_vals = pd.DataFrame()
        self.id = None
        self.FPKM = None

    def from_csv(self,filename,FPKM=True,sep="\t"):
        ''' Returns a COBDataset from a CSV '''
        self.expr_vals = pd.read_table(filename,sep=sep)
        try:
            self.expr_vals[self.expr_vals.columns] = self.expr_vals[self.expr_vals.columns].convert_objects(convert_numeric = True)
        except e:
            exit("csv expression values must be numbers")

    def save(self):
        self.id = self.query("SELECT MAX(id) as MID FROM datasets;").iloc[0]['MID'] + 1
        # add the entry into the database
        self.execute("INSERT INTO datasets (id, name, description) VALUES ({}, '{}', '{}')",self.id,self.name,self.description)
