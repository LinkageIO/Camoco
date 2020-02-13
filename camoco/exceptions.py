class CamocoError(Exception):
    def __init__(self,msg=''):
        self.message=msg

class CoexNetEmptyError(CamocoError):
    pass
