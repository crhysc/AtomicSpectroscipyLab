import numpy as np
filenm = ["alpha100shots.txt", "beta100shots.txt", "gamma50shots.txt"]
ns = [((4*25)/(25-4)),((4*16)/(16-4)),((4*9)/(9-4))]
Rh = []
Rd = []
ah = []
ad = []
mdmh = []
amdmh = []
for i in range(len(filenm)):
    with open(f"./output/{filenm[i]}") as f:
        lines = f.readlines()
        Rh.append(lines[4][18:-1])
        Rd.append(lines[5][19:-1])
        ah.append(lines[3][8:-1])
        ad.append(lines[2][8:-1])
        mdmh.append(lines[10][11:-1])
        amdmh.append(lines[11][9:-1])
Rh = list(map(float, Rh))
Rd = list(map(float, Rd))
ah = list(map(float, ah))
ad = list(map(float, ad))
mdmh = list(map(float, mdmh))
amdmh = list(map(float, amdmh))

weightsh = []
weightsd = []
oneoverh = 0
oneoverd = 0
for i in range(len(ah)):
    oneoverh = oneoverh + (1/ah[i])
for i in range(len(ad)):
    oneoverd = oneoverd + (1/ad[i])
for i in range(len(filenm)):
    weightsh.append(3/(((ah[i])**2)*oneoverh**2))
    weightsd.append(3/(((ad[i])**2)*oneoverd**2))

fh = 0
sh = 0
fd = 0
sd = 0
for i in range(len(weightsh)):
    fh = fh+(weightsh[i] * (ns[i]**2))
    sh = sh + (weightsh[i]*ns[i])**2
    fd = fd+(weightsd[i] * (ns[i]**2))
    sd = sd + (weightsd[i]*ns[i])**2
deltah = 3*fh-sh
deltad = 3*fd-sd

amh = np.sqrt(sum(weightsh)/deltah)
amd = np.sqrt(sum(weightsd)/deltad)

Rhtrue = sum(Rh)/len(Rh)
Rdtrue = sum(Rd)/len(Rd)

aRh = (amh*(Rhtrue/1000000000)**2) * 1000000000
aRd = (amd*(Rdtrue/1000000000)**2) * 1000000000

mdmhtrue = sum(mdmh)/len(mdmh)
amdmhtrue = np.sqrt(amdmh[0]**2+amdmh[1]**2+amdmh[2]**2)/len(mdmh)

print(f"{aRh} {aRd}")
for i in range(len(filenm)):
    with open(f"./output/{filenm[i]}", "r", encoding='utf-8') as f:
        lines = f.readlines()
        lines[6] = f"α_Rh: {aRh}\n"
        lines[7] = f"α_Rd: {aRd}\n"
        lines[8] = f"Average Rh: {Rhtrue}\n"
        lines[9] = f"Average Rd: {Rdtrue}\n"
        lines[12] = f"Average md/mh: {mdmhtrue}\n"
        lines[13] = f"α_Average md/mh: {amdmhtrue}\n"
    with open(f"./output/{filenm[i]}", "w", encoding='utf-8') as f:
        f.writelines(lines)