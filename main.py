import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

line = "beta100"
peaks = [656.279, 486.135, 434.0472, 365.0158] #alpha, beta, gamma, mercury
n=4

def DoubleGauss(x, mini, con1, cen1, std1, con2, cen2, std2):
    return mini+con1*np.exp(-((x-cen1)**2)/(2*std1**2)) + con2*np.exp(-((x-cen2)**2)/(2*std2**2))
def TripleGauss(x, mini, con1, cen1, std1, con2, cen2, std2, con3, cen3, std3):
    return mini+con1*np.exp(-((x-cen1)**2)/(2*std1**2)) + con2*np.exp(-((x-cen2)**2)/(2*std2**2))+ con3*np.exp(-((x-cen3)**2)/(2*std3**2))

average = []
pixel = []
with open(f"./{line}shotsWITHOUTX.txt") as f:
    lines = f.readlines() #new lines are seperate x values, commas are seperate shots, values are the y value
    for i in range(len(lines)):
        row = lines[i].split(",") 
        row[-1] = row[-1].strip() #removes the newline from the last element
        row = list(map(int, row)) #converts string to int
        average.append((sum(row))/len(row)) #averages every shot
        pixel.append(i+1)


data, error = sc.optimize.curve_fit(DoubleGauss, pixel, average, p0=[np.min(average),np.max(average),1580, 15, np.max(average), 1650, 15])
mini, con1, cen1, std1, con2, cen2, std2 = data

HgCalibrationSlope = 0.00323738 #for narrow slit, peaks at 1625, 1773, 2031 pixels
offset = (HgCalibrationSlope*3352)*((cen2/3352)) #takes the resolution of the camera and multiplies it by how far along the center of the peak is

Calibrated = []
for i in pixel:
    Calibrated.append(i*HgCalibrationSlope+(peaks[n-3]-offset)) #linear regression, wavelength = pixel*slope + offset

data, error = sc.optimize.curve_fit(DoubleGauss, Calibrated, average, p0=[np.min(average),np.max(average),Calibrated[int(cen1)], .1, np.max(average), Calibrated[int(cen2)], .1])
#data, error = sc.optimize.curve_fit(TripleGauss, Calibrated, pixel, p0=[np.min(average),np.max(average),1625, 15, np.max(average)/5.5, 1773, 15, np.max(average)/7, 2031, 15])
mini, con1, cen1, std1, con2, cen2, std2 = data
#mini, con1, cen1, std1, con2, cen2, std2, con3, cen3, std3 = data
var = np.diag(error)
uncert = np.sqrt(var)
fit = DoubleGauss(Calibrated, mini, con1, cen1, std1, con2, cen2, std2)
#fit = TripleGauss(pixel, mini, con1, cen1, std1, con2, cen2, std2, con3, cen3, std3)

#Rydberg
indexOfRefraction = 1.00029
cen2I = cen2*indexOfRefraction
cen1I = cen1*indexOfRefraction
unc1I = (std1/np.sqrt(2*np.log(2)))*indexOfRefraction
unc2I = (std2/np.sqrt(2*np.log(2)))*indexOfRefraction

Rh = (1/cen2I) * ((4*n**2)/(n**2-4))*(10)**9 #hydrogen rydberg
Rd = (1/cen1I) * ((4*n**2)/(n**2-4))*(10)**9 #deuterium rydberg
aRh = "x"#uncertainty in Rh
aRd = "x"#uncertainty in Rd

#mass ratio
idealized = 1/.01097373156816 * ((4*n**2)/(n**2-4))
mdmh = (cen2I-idealized)/((cen2I-idealized)-(cen2I-cen1I)) #mass ratio
amdmh = np.sqrt((unc2I**2/(cen1I-idealized)**2)+(((cen2I-idealized)**2*(unc1I)**2)/(cen1I-idealized)**4)) #mass ratio uncertainty

print(f"center deuterium: {cen1}\ncenter hydrogen: {cen2}\nD error: {std1/np.sqrt(2*np.log(2))}\nH error: {std2/np.sqrt(2*np.log(2))}\nHydrogen Rydberg: {Rh}\nDeuterium Rydberg: {Rd}\nα_Rh: {aRh}\nα_Rd: {aRd}\nMass Ratio: {mdmh}\nα_md/mh: {amdmh}")
#print(f"peak 1: {cen1}\npeak 2: {cen2}\npeak 3: {cen3}\n")
with open(f"./output/{line}shots.txt", "w", encoding='utf-8') as f:
    f.write(f"center deuterium: {cen1}\ncenter hydrogen: {cen2}\nD error: {std1/np.sqrt(2*np.log(2))}\nH error: {std2/np.sqrt(2*np.log(2))}\nHydrogen Rydberg: {Rh}\nDeuterium Rydberg: {Rd}\nα_Rh: {aRh}\nα_Rd: {aRd}\nAverage_Rh: x\nAverage_Rd: x\nMass Ratio: {mdmh}\nα_md/mh: {amdmh}\n x\n x")
    #f.write(f"peak 1: {cen1}\npeak 2: {cen2}\npeak 3: {cen3}\n")


plt.plot(Calibrated, average)
plt.plot(Calibrated, fit)
plt.legend(["Raw Data","Fitted Gaussian"])
plt.xlabel("wavelength (nm)")
plt.ylabel("counts (relative intensity)")
plt.savefig(f"./output/{line}shotsFit.jpg", bbox_inches='tight')

plt.show()
