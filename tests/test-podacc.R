

library("podacc")


data = amiodarone_tggates_29day_BMD

head(data)

list_bmd_values = data[,"BMD"]

aux = returnAntimode(list_bmd_values,printHist=FALSE)
firstmode= aux[1]
antimode = aux[2]

# Get accumulation plot values
output = generateAccumulationValuesFromListBMDValues(list_bmd_values)
bmd_Values_acc = output[,1]
accumulated = output[,2]


plot(bmd_Values_acc,accumulated,log="x")

# Calculate method for the whole curve first - Needs package scam
elbows_all_curve = getMaximumCurvaturePoints(bmd_Values_acc,accumulated,bmd_Values_acc[length(bmd_Values_acc)])

# Calculate using first antimode as the end point of the curve
elbows_anti_mode = getMaximumCurvaturePoints(bmd_Values_acc,accumulated,10^antimode)

all_elbows = setdiff(union(elbows_all_curve,elbows_anti_mode),NA)

if(!is.na(all_elbows[1])){
  # Print plot - choose closest to first mode
  difffirstmode = abs(all_elbows-10^firstmode)
  elbow = all_elbows[which.min(difffirstmode)]
}else{
  elbow = NA
}

if(is.na(elbow)){
  elbow = median(bmdValues[1:2])
}

titleplot="DEMO Accumulation Plot Results"

xlab_text = "BMD (mg/kg/day)"

plot(bmd_Values_acc,accumulated,log="x",main=titleplot,xlab=xlab_text,ylab="Accumulation")

abline(v=(10^firstmode),col="red", lwd=3)
abline(v=(10^antimode),col="blue", lwd=3)

for(e in 1:length(all_elbows)){
  abline(v=all_elbows[e],col="green", lwd=3,lty=2)
}
abline(v=elbow,col="green", lwd=3)

legend(x="bottomright",legend=c("First Mode","Antimode",expression("POD"[Accum])),fill=c("red","blue","green"))

