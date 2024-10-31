import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

line = "gamma50"

def DoubleGauss(x, mini, con1, cen1, std1, con2, cen2, std2):
    return mini+con1*np.exp(-((x-cen1)**2)/(2*std1**2)) + con2*np.exp(-((x-cen2)**2)/(2*std2**2))

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


data, error = sc.optimize.curve_fit(DoubleGauss, pixel, average, p0=[np.min(average),np.max(average),1580, 43/3, np.max(average), 1640, 43/3])
mini, con1, cen1, std1, con2, cen2, std2 = data
var = np.diag(error)
uncert = np.sqrt(var)
fit = DoubleGauss(pixel, mini, con1, cen1, std1, con2, cen2, std2)

print(f"center deuterium: {cen1}\ncenter hydrogen: {cen2}\nD error: {std1/np.sqrt(2*np.log(2))}\nH error: {std2/np.sqrt(2*np.log(2))}")
with open(f"./output/{line}shots.txt", "w") as f:
    f.write(f"center deuterium: {cen1}\ncenter hydrogen: {cen2}\nD error: {std1/np.sqrt(2*np.log(2))}\nH error: {std2/np.sqrt(2*np.log(2))}")


plt.plot(pixel, average)
plt.plot(pixel, fit)
plt.legend(["Raw Data","Fitted Gaussian"])
plt.xlabel("pixels")
plt.ylabel("counts (relative intensity)")
plt.savefig(f"./output/{line}shotsFit.jpg", bbox_inches='tight')

plt.show()
