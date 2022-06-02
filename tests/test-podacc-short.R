

library("podacc")

data = amiodarone_tggates_29day_BMD

head(data)

list_bmd_values = data[,"BMD"]

results = runPODAccMethod(list_bmd_values)

plotPODAccResults(list_bmd_values,results)

plotPODAccResults(list_bmd_values,results,titleplot="DEMO Accumulation Plot Results",xlab_text = "BMD (mg/kg/day)")


