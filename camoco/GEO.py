import pandas as pd
import urllib
import os
from camoco.Tools import log

class Soft(object):
    def __init__(self, name, type):
        self.name = name
        self.type = type
        self.info = dict()
        self.headers = []
        self.tbl = pd.DataFrame()
    def update_info(self,key,val):
        if key not in self.info:
            self.info[key] = val
        elif val == self.info[key]:
            self.info[key] = val
        elif key in self.info and not isinstance(self.info[key],list):
            self.info[key] = [self.info[key], val]
        else:
            self.info[key].append(val)

    def __add__(self,other):
        ''' combines two Soft instances into a single instance '''
        # complain if types do not match up
        if self.type != other.type:
            raise Exception("Can only combine similar types of Soft")
        #combine names
        if self.name != other.name:
            self.name = "{}-{}".format(self.name,other.name)
        # combine info
        for key,val in other.info.items():
            self.update_info(key,val)
            #combine headers
            self.headers.extend(other.headers)
        if not self.tbl.empty and not other.tbl.empty:
            self.tbl = pd.concat([self.tbl,other.tbl]).drop_duplicates()
        return self
    
      
    @staticmethod
    def wget(id,force=False):
        ''' Downloads the GEO series from the internets into PWD'''
        if os.path.exists("{}_family.soft.gz".format(id)) and force == False:
            log("{} already exists",id)
            return
        try:
            log("Fetching {}",id)
            gse = urllib.request.urlretrieve(
                "ftp://ftp.ncbi.nlm.nih.gov/geo/series/{}nnn/{}/soft/{}_family.soft.gz".format(id[0:len(id)-3],id,id),
                "{}_family.soft.gz".format(id)
            )
        except Exception as e:
            log("Could not download {}",id)
        
    
class Series(Soft):
    def __init__(self,name):
        super().__init__(name)

class Family(object):
    def __init__(self):
        self.database = None
        self.series = None
        self.platform = None
        self.samples = []
    
    @classmethod
    def from_file(cls,filename):
        self = cls() 
        with open(filename,'r') as IN:
            in_data_table = False
            cur_soft = None
            cur_data = list()
            for i,line in enumerate(IN):
                line = line.strip()
                if line.startswith('^'):
                    if cur_soft:
                        # Add the filled SOFT to Family
                        if cur_soft.type == 'Sample':
                            self.samples.append(cur_soft)
                        else:
                            setattr(self,cur_soft.type.lower(),cur_soft)
                    # WE have a new SOFT
                    type,name = line.replace('^','').replace(' = ','=').split('=',1)
                    cur_soft = Soft(name,type=type.lower().capitalize())
                    cur_data = list()
                elif line.startswith('!') and 'table_begin' in line:
                    in_data_table = True
                elif line.startswith('!') and 'table_end' in line:
                    in_data_table = False
                    # Create DataFrame and append to SOFT
                    cur_headers = cur_data.pop(0)
                    cur_soft.tbl = pd.DataFrame.from_records(data=cur_data,columns=cur_headers)
                    cur_soft.tbl.index = cur_soft.tbl.icol(0)
                    cur_data = list()
                elif line.startswith("!"):
                    # add info to 
                    key,val = map(str.strip,line.replace('!'+cur_soft.type+'_','').split('=',1))
                    cur_soft.update_info(key,val)
                elif line.startswith('#'):
                    # Columns descriptions
                    cur_soft.headers.append(line)
                elif in_data_table:
                    cur_data.append(line.replace('"','').split('\t'))
        # If multiple samples, create SeriesMatrix
        if len(self.samples) > 0:
            self.series_matrix = pd.DataFrame([x.tbl.VALUE for x in self.samples]).transpose()
            # Label the columns
            self.series_matrix.columns = [x.name for x in self.samples]
            # Relabel to orfs
            new_index = []
            self.series_matrix.index = self.platform.tbl['ORF'][self.series_matrix.index].values
        return self


    def __add__(self,family):
        ''' combines families into a single instance '''
        self.database += family.database
        self.series += family.series
        self.platform +=  family.platform
        self.samples.extend(family.samples)
        return self
        
