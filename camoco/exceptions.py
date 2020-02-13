class CamocoError(Exception):
    def __init__(self,msg=''):
        self.message=msg

class CoexExistsError(CamocoError):
    pass

class CoexEmptyError(CamocoError):
    pass
