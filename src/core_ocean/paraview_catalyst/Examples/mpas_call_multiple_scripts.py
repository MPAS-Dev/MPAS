# Substitute operation_* with your own script name, e.g. mpas_T_KE
import operation_a
import operation_b

def RequestDataDescription(datadescription):
    operation_a.RequestDataDescription(datadescription)
    operation_b.RequestDataDescription(datadescription)

def DoCoProcessing(datadescription):
    operation_a.DoCoProcessing(datadescription)
    operation_b.DoCoProcessing(datadescription)


