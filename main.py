import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import pandas as pd
import re
import os.path
import os

narrowSlitInstrumentalLineWidth = .0779634





def process(runtime, fileindex, parentFolder, magfield, wavelength, Element, peakLocation, off = 0, maxlim=120):
    deltaLambda = []
    shotnumber = []
    notCold = []
    for z in range(runtime):
        filenm = str(fileindex)
        if os.path.isfile(parentFolder+"/"+filenm+".asc"):
            shotnumber.append(filenm)
            #Takes the ascii file and will extract the information per row

            index1 = int(35+off)
            index2 = int(42+off)
            with open(f"{parentFolder}/{filenm}.asc") as fl:
                lines = fl.readlines()
                gateDelay = re.findall(r'\d+', lines[index1])
                print(lines[index1])
                gateDelay = gateDelay[0]+ " ns"
                lines = lines[index2:]

            #Processes the file to remove metadata and other unneccesary entries
            with open(f"{parentFolder}/txt/{filenm}.txt", "w") as file:
            
                for line in lines:
                    index = line.index(",") #the first entry is just the column number, not needed
                    line = line[index+1:-2] #removes the column number and the ,\n at the end
                    file.write(f"{line}\n")


            col=list(range(len(lines))) #not strictly neccessary
            df = pd.read_csv(f"{parentFolder}/txt/{filenm}.txt", delimiter=',', names=col)
            df_calculation = df
            #df_calculation = df_calculation/np.max(df_calculation) #normalization
            df_calculation['mean'] = df_calculation.mean(axis=1) #final y-axis
            dfList = df_calculation['mean'].tolist()


            #raw ccd (rotated 90 degrees)
            plt.contourf(df.T, cmap="gist_heat")
            plt.title(Element)
            plt.xlabel("pixels")
            plt.ylabel("pixels")
            plt.figtext(.15,.5,f"Gate Delay: {gateDelay}\nGuide Field: {magfield}\nWavelength: {wavelength} nm", color='w')
            plt.savefig(f"{parentFolder}/ccd/{filenm}.png", bbox_inches='tight')
            plt.clf()

            #scaling col to be in nm knowing the peak's wavelength
            #do this part :)


            col = np.array(col)
            dfList = np.array(dfList)


            #filters the intensity vs pixel plot
            filterer = sc.signal.lfilter
            n = 20
            b = [1/n]*n
            a = 1
            filtered = filterer(b, a, dfList)
            for i in range(20):
                filtered[i] = filtered[21] #bug

            mini = np.min(filtered)
            #then finds the min of that to subtract from the original
            #this gets most values very close to 0
            for i in range(len(col)):
               dfList[i] = np.abs(dfList[i] - mini)
               #if dfList[i] < .0015:
               # dfList[i] = 0

            #simple lorentzian, replace for more complex graphs
            def Lorentz(x,amp1, cen1, wid1):
                return amp1/(1+((x-cen1)/wid1)**2)
            
            def Gauss(x, cen1, std1):
                return np.min(dfList)+(np.max(dfList)*np.exp(-((x-cen1)**2)/(2*std1**2)))

            par, cor = sc.optimize.curve_fit(Gauss, col, dfList, p0=[peakLocation, 12])
            error = np.sqrt(np.diag(cor))

            HgCalibrationSlope = (366.32840-365.01580)/(570-404) #for narrow slit
            offset = (HgCalibrationSlope*1024)*((par[0]/1024)) #takes the resolution of the camera and multiplies it by how far along the center of the peak is

            Calibrated = []
            for i in col:
                Calibrated.append(i*HgCalibrationSlope+(wavelength-offset)) #linear regression, wavelength = pixel*slope + offset


            parNM, corNM = sc.optimize.curve_fit(Gauss, Calibrated, dfList, p0=[wavelength, .1])
            errorNM = np.sqrt(np.diag(corNM))

            lineWidth = np.abs(parNM[1])
            lineWidthError = np.abs(errorNM[1])

            lorenz = []
            for i in range(len(col)):
                lorenz.append(Gauss(Calibrated[i],parNM[0], lineWidth))

            print(parNM)
            print(errorNM)

            if lineWidth**2-narrowSlitInstrumentalLineWidth**2 > 0:
                truelineWidth = np.sqrt(lineWidth**2-narrowSlitInstrumentalLineWidth**2)
                notCold.append([filenm,truelineWidth])
            else:
                truelineWidth = 0

            #intensity vs pixels
            plt.plot(Calibrated, dfList)
            plt.plot(Calibrated, lorenz)
            plt.title(Element)
            plt.xlabel("wavelength (nm)")
            plt.ylabel("intensity (arb. units)")
            plt.ylim([0,maxlim])
            plt.legend(["Blue: filtered peak", "Orange: fitted peak"])
            plt.figtext(.15,.5,f"Gate Delay: {gateDelay}\nGuide Field: {magfield}\nWavelength: {wavelength} nm\nLine Width: {lineWidth:.5f} nm\nsqrt(Δλ_measured^2-Δλ_instrumental^2): {truelineWidth: .5f}")
            plt.savefig(f"{parentFolder}/intensity/{filenm}.png", bbox_inches='tight')
            plt.clf()

            deltaLambda.append(lineWidth)

            #df2 = pd.read_csv(f"{parentFolder}/txt/{filenm}.txt", delimiter=',', names=Calibrated)
            #plt.contourf(list(df2.columns.values),col,df2.T, cmap="gist_heat")
            #plt.title(Element)
            #plt.xlabel("wavelength (nm)")
            #plt.ylabel("pixels")
            #plt.figtext(.15,.5,f"Gate Delay: {gateDelay}\nGuide Field: {magfield}\nWavelength: {wavelength} nm", color='w')
            #plt.savefig(f"{parentFolder}/ccd/{filenm}.png", bbox_inches='tight')
            #plt.clf()

            #cleans up space by removing the txt's after itself
            os.remove(f"{parentFolder}/txt/{filenm}.txt")

        fileindex += 1
    return deltaLambda, shotnumber, notCold, lorenz #only the last lorentzian will get returned

#####HELIUM#########
parentHe = "./Helium/ion"

fileHe = 138250
runHe = 6
magHe = "225 G"
waveHe = 657.50154
EleHe = "He"

peakHe = 420
###################
#####Mercury#########
parentHg = "./Mercury"

fileHg = 100001
runHg = 1
magHg = "225 G"
waveHg = 365.0158
EleHg = "Hg I"

peakHg = 410
###################
#####Argon#########
parentAr = "./Argon/Ion225G656nm"

fileAr = 138263
runAr = 6
magAr = "225 G"
waveAr = 657.50154
EleAr = "Ar"

peakAr = 420
###################


deltaAr, shotAr, ColdAr, LorAr = process(runAr, fileAr, parentAr, magAr, waveAr, EleAr, peakAr, off=1, maxlim=500)
deltaHe, shotHe, ColdHe, LorHe = process(runHe, fileHe, parentHe, magHe, waveHe, EleHe, peakHe, off=1)
deltaHg, shotHg, ColdHg, LorHg = process(runHg, fileHg, parentHg, magHg, waveHg, EleHg, peakHg, maxlim=40000)


for i in range(len(deltaHe)):
    deltaHe[i] = deltaHe[i]/waveHe

for i in range(len(deltaHg)):
    deltaHg[i] = deltaHg[i]/waveHe


col = []
for i in range(1024):
    col.append(i)
plt.plot(col,LorHe, color="r")
plt.plot(col,LorHg, color="b")
plt.plot(col,LorAr, color="g")
plt.legend(["Helium", "Mercury", "Argon"])
plt.show()
plt.clf()


plt.plot(shotHe, deltaHe, 'ro')
plt.plot(shotHg, deltaHg, 'o')


plt.title("Δλ/λhe")
plt.legend(["Helium", "Mercury", "Argon"])
plt.xlabel("shotnumber")
plt.xticks(rotation=90)
plt.ylabel("Δλ/λhe")
plt.show()
plt.clf()

for i in range(len(deltaAr)):
    deltaAr[i] = deltaAr[i]/waveAr
for i in range(len(deltaHg)):
    deltaHg[i] = deltaHg[i]/waveAr

plt.plot(shotAr, deltaAr, 'go')
plt.plot(shotHg, deltaHg, 'o')
plt.xlabel("shotnumber")
plt.xticks(rotation=90)
plt.ylabel("Δλ/λar")
plt.legend(["Argon", "Mercury"])
plt.show()

