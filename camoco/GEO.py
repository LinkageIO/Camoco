import pandas as pd
import numpy as np
import urllib
import itertools
import os
from camoco.Tools import log
import warnings

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
        soft = Soft("{}-{}".format(self.name,other.name),type=self.type)
        if self.type != other.type:
            raise Exception("Can only combine similar types of Soft")
        if self.type == 'Series' and self.info['platform_id'] != other.info['platform_id']:
            warnings.warn("About to combine {} with the different platforms: {} != {}".format(
                self.type,
                self.info['platform_id'],
                other.info['platform_id']
            ))
        # combine info
        for key,val in other.info.items():
            soft.update_info(key,val)
            #combine headers
            soft.headers = self.headers + other.headers
        if not self.tbl.empty and not other.tbl.empty:
            soft.tbl = pd.concat([self.tbl,other.tbl]).drop_duplicates()
        return soft

    
      
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
        super().__init__(name,'Series')

class Sample(Soft):
    def __init__(self,name):
        super().__init__(name,'Sample')

    def is_raw(self):
        ''' guesses whether or not the samples values have been previously transformed
            into log space '''
        if max(list(map(float,self.tbl.VALUE.values))) > 200:
            # It is very unlikely a log transformed value would exceed 200
            # If you can come up with a better solution, by all means ...
            return True
        else:
            return False

    def is_empty(self):
        return self.tbl.empty

    def transform(self,func=np.log2):
        ''' Transform the samples table VALUES in place'''
        if not self.is_raw():
            log.warn("Attempting to perform transormation of apparently non raw data")
        self.tbl.VALUE = func(list(map(float,self.tbl.VALUE.values)))
        self.name = self.name
        # make sure we didnt introduce and -Inf values
        self.tbl.loc[self.tbl.VALUE == float('-Inf'),'VALUE'] = np.nan
        return True

    def __add__(self,other):
        if self.is_raw() and not other.is_raw():
            log.warn('WARNING! attempting to combine {} and {} which are not both normalized')
        super().__add__(self,other)

class Platform(Soft):
    def __init__(self,name):
        super().__init__(name,'Platform')
    
    def __repr__(self):
        return '''
            Platform: 
        '''

class Family(object):
    def __init__(self):
        self.database = None
        self.series = None
        self.platform = None
        self.samples = []

    def to_keepfile(self,filename,group_max_r2=0.99,keep_hint=None):
        ''' Creates a tsv file for each sample containing sample information. Included are two columns 
            used for further filtering of the aggregated samples. The 'Keep' columns is a hard filter
            used to remove experiments (columns) completely. The "Group" column is for biological replicates.

            See Also: GEO.Family.filter_from_keepfile()'''
        # Extract info from samples
        if os.path.exists(filename):
            self.log('{} exists ... skipping',filename)
            return
        info = pd.DataFrame([sample.info for sample in self.samples])[['geo_accession','title','description','characteristics_ch1','series_id']]
        info.sort("title") # Replicates most likely to have similar names
        if keep_hint:
            # This first join joins all the words together
            info.insert(1,'Keep',[keep_hint in " ".join(
                    # This second join is because when the comprehension comes to the string
                    # it thinks its a list of chars. Just join everything together.
                    ["".join(map(str.lower,str(x))) if isinstance(x,list) else str(x).lower() for x in row]
                ) for i,row in info[['title','description','characteristics_ch1']].iterrows()
            ])
        else:
            info.insert(1,'Keep',False)
        info.insert(1,'Group',self._guess_groups(self.series_matrix(),group_max_r2))
        info.to_csv(filename,sep='\t',index=False)


    @classmethod
    def from_file(cls,filename,normalize=True):
        self = cls() 
        with open(filename,'r') as IN:
            in_data_table = False
            cur_soft = None
            cur_data = list()
            for i,line in enumerate(IN):
                line = line.strip()
                if line.startswith('^'):
                    if cur_soft: # Add the filled SOFT to Family
                        if cur_soft.type == 'Sample':
                            if cur_soft.is_raw() and normalize:
                                log("Normalizing {}",cur_soft.name)
                                cur_soft.transform()
                            self.samples.append(cur_soft)
                        else:
                            setattr(self,cur_soft.type.lower(),cur_soft)
                    # WE have a new SOFT
                    type,name = line.replace('^','').replace(' = ','=').split('=',1)
                    type = type.lower().capitalize()
                    if type == 'Series':
                        cur_soft = Series(name)
                    elif type == 'Sample':
                        cur_soft = Sample(name)
                    elif type == 'Platform':
                        cur_soft = Platform(name)
                    else:
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
                    # Turn -Inf into NaNs
                    cur_soft.tbl[cur_soft.tbl == float('-Inf')]  = np.nan
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
            return self

    def probe2gene(self,probe_list):
        ''' returns a mapping from probes to genes '''
        mapping = []
        for probe in probe_list:
            if probe in self.platform.tbl.index:
                orf = self.platform.tbl.ix[probe]['ORF']
                if orf == '':
                    mapping.append(probe)
                else:
                    mapping.append(orf)
            else:
                mapping.append(probe)
        return mapping

    @staticmethod
    def _uniqify_columns(df_columns):
        ''' Given a list of column names, returns a similar sized list with 
            unique names '''
        seen = set()
        new = list()
        for item in df_columns:
            fudge = 1
            newitem = item
            while newitem in seen: # keep increasing fudge until not seen 
                fudge += 1
                newitem = "{}_{}".format(item, fudge)
            new.append(newitem)
            seen.add(newitem)
        return new

    @staticmethod
    def _guess_groups(dataframe,max_r2=0.99,max_namediff=0.8):
        ''' Given a data frame, this method checks to see that each column has a correlation
            below the max_r2. If it is above, a new column is created using the mean.'''
        # Calculate correlation
        cors = dataframe.corr()
        # Each column starts in its own group
        column_groups = list(range(0,len(cors)))
        # Iterate over upper triangular, this guarantees that 
        # we dont overwrite lower numbered groups
        #
        #   #OOOOOOOOOOO   <- group
        #   i              <- row number (along d)
        #  -------------   <- matrix
        #  |d   x    x     <- row (x means r2 > max)
        #  | d
        #  |  d
        #  |   d
        #  |    d
        for i,row in enumerate(np.triu(cors.as_matrix(),k=1)):
            if any(row > max_r2):
                # higher numbered columns are highly correlated with current column
                which = np.where(row > max_r2)[0]
                from difflib import SequenceMatcher
                for match in which:
                    diffratio = SequenceMatcher(None,dataframe.columns[i],dataframe.columns[match]).ratio()
                    # if column group is already assigned, keep assignment
                    if column_groups[i] != i and diffratio > max_namediff:
                        # Here, samples have high expression correlation AND similar names (i.e. rep1 vs rep2)
                        group = column_groups[i]
                        log("{} is {} correlated with {}",
                            dataframe.columns[i],
                            diffratio,
                            dataframe.columns[match]
                        )
                    else:
                        # Otherwise start your own group
                        group = i
                    for x in which:
                        column_groups[x] = group
        return column_groups

    def series_matrix(self,keepfile=None):
        # If multiple samples, create SeriesMatrix
        organisms = set([sample.info['organism_ch1'] for sample in self.samples])
        if len(organisms) > 1:
            raise ValueError("Cannot combine series from different organisms: {}".format(",".join()))
        if len(self.samples) > 0:
            series_matrix = pd.DataFrame(
                data=[x.tbl.VALUE for x in self.samples],
                index=[x.info['title'] for x in self.samples]
            ).transpose()
            # Label the rows
            old_index = series_matrix.index
            series_matrix.index = self.probe2gene(series_matrix.index.values)
            # Uniqify the columns
            series_matrix.columns = self._uniqify_columns(series_matrix.columns)
            # Collpase down highly correlated columns (biological replicates)
            series_matrix = series_matrix.astype('float')
            #series_matrix = self._collapse_columns(series_matrix)
            # replace NANs or empty strings with old labels
            if keepfile is not None:
                # Read in the keepfile
                keepfile = pd.read_table(keepfile,sep='\t')
                keepfile.title = self._uniqify_columns(keepfile.title)
                # Reorder the columns based on the keepfile
                series_matrix = series_matrix.reindex_axis(keepfile.title,axis=1)
                # Make sure the above worked
                assert all(series_matrix.columns == keepfile.title)
                # Filter based on the keepfile
                groups = series_matrix[series_matrix.columns[keepfile.Keep]].groupby(
                    keepfile[keepfile.Keep].Group.values, axis=1
                )
                series_matrix = groups.apply(lambda group: group.mean(skipna=True,axis=1))
                # Relabel the columns
                series_matrix.columns = [group[0] if len(group) == 1 else "(BioRep) "+";".join(group) for group in groups.groups.values()]
            return series_matrix

        else:
            raise ValueError("Not enough samples to build matrix")


    def __add__(self,family):
        ''' combines families into a single instance '''
        # make a new one
        fam = Family()
        fam.database = self.database + family.database
        fam.series = self.series + family.series
        fam.platform = self.platform + family.platform
        fam.samples = self.samples + family.samples
        return fam
    def __radd__(self,something):
        if isinstance(something,int):
            return self
        else:
            return something + self

    def __str__(self):
        return '''
            Family: {},
            Summary: {},
            Title: {},
            Series Platform: {}
            Num Samples: {}
        '''.format(
                self.series.name,
                self.series.info['summary'],
                self.series.info['title'],
                self.series.info['platform_id'],
                len(self.samples)
        )

    def __repr__(self):
        return str(self)

