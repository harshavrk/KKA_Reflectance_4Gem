## Author: Sanat Kumar Gogoi 29th March 2026
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rcParams["figure.constrained_layout.use"] = True


rfs = np.loadtxt('Diamond_specR.dat', float)
(d, d1) = np.shape(rfs)

l = rfs[:,0]
en = 1239.9219/l
R = rfs[:,1]/100

IGgI = np.zeros([d,d])
igm = np.zeros([d,d])
IGhI = np.zeros([d,d])
dnf = np.zeros([d,d])
n = np.zeros(d)


l0 = l[0]-1
ds = 1
a = 199
b = 1001
g = 950
h = 980


def calc_m(matrix,lamda):
    m,n = np.shape(matrix)
    for i in range(m):
        for j in range(m):
            l0jds = l0 + j*ds       #ds = 1, step size
            if l0jds ==lamda:
                matrix[i,j] = 0
            else:
                matrix[i,j] = 1e9*(1/np.pi)* (np.log(R[j]))*l[i]/(lamda**2 -l0jds**2)
    sum_mat = np.sum(matrix, axis=1)  
    return matrix, sum_mat

IGgI, sg = calc_m(IGgI, g)
IGhI, sh = calc_m(IGhI, h)

def calc_PP(matrix,sum_mat):
    PP = (sum_mat - (matrix[:, 0] + matrix[:, -1]) / 2) * 1e-9
    return PP

PPg = calc_PP(IGgI, sg)
PPh = calc_PP(IGhI, sh)
Pg = PPg[750]
Ph = PPh[780]


ag = np.log(abs((a+g)/(a-g)))
bg = np.log(abs((b+g)/(b-g)))
ah = np.log(abs((a+h)/(a-h)))
bh = np.log(abs((b+h)/(b-h)))
A  = (Ph*bg-Pg*bh)/(ag*bh-ah*bg)       #Evaluating the coefficient A 
B  = (Pg*ah-ag*Ph)/(ag*bh-ah*bg)       #Evaluating the coefficient B



for i in range(d):
    for j in range(d):
        l0jds = l0 + j*ds
        if l[i]==l0jds:
            igm[i,j] = 0
        else:
            igm[i,j] = 1e9*(1/np.pi)* (np.log(R[j]))*l[i]/((l[i])**2 -l0jds**2)

s = np.sum(igm, axis=1)


thetam = calc_PP(igm, s)
Ac = np.log(np.abs((a+l)/(a-l)))
Bc = np.log(np.abs((b+l)/(b-l)))
thetafin = A * Ac + B * Bc + thetam

denom = 1 + R - 2*np.cos(thetafin) * np.sqrt(R)
u = (1 -R)/denom
v = (-2 * np.sqrt(R)* np.sin(thetafin))/denom

n2 = 2.3673
n1 = 2.5981
wl2 = 800.99
wl1 = 208.0

kl = (v * l * (1e-9)) / (4 * np.pi * 2 * 3.5*1e-3)
for i in range(d):
    for j in range(d):
        l0jds = l0 + j*ds
        if (l[i] == l0jds or l0jds == wl1 or l0jds == wl2):
            dnf[i, j] = 0
        else:
            denom1 = (l[i]**2-l0jds**2)*(wl2**2-l0jds**2)*(wl1**2-l0jds**2)
            dnf[i, j] = (kl[j]/(l[i]**2))*(l0jds**3/denom1)

sf = np.sum(dnf, axis=1)

wl21 = wl2**2 - wl1**2
wl12 = wl1**2 - wl2**2
wl11 = wl1**2
wl22 = wl2**2
step = l[-1] - l[-2]
two_pi = 2 / np.pi

l2 = l**2
term1 = (n1 - 1) * (wl11 / l2) * ((wl22 - l2) / wl21)
term2 = (n2 - 1) * (wl22 / l2) * ((wl11 - l2) / wl12)
cor = sf - (dnf[:, 1] + dnf[:, -2]) / 2
term3 = two_pi * (l2 - wl22) * (l2 - wl12) * step * cor
n = 1 + term1 + term2 + term3



nfin = np.zeros(d)
for i in range(d):
    if 200 < i < 500:
        nfin[i] = n[i]
    else:
        nfin[i] = 0

navg = np.sum(nfin) / 300  
print(f"Mean refractive index (400-700 nm): {navg}")

# Plotting
fig, ax1 = plt.subplots(figsize=(10, 8))
ax1.plot(l, thetafin, 'b')
ax1.set_xlabel('wavelength', fontsize=16)
ax1.set_ylabel('Phase (radian)', color='b', fontsize=16)
ax1.tick_params(axis='y', labelcolor='b')

ax2 = ax1.twinx()
ax2.plot(l, R, 'r')
ax2.set_ylabel('Specular Reflectance', color='r', fontsize=16)
ax2.tick_params(axis='y', labelcolor='r')

plt.title('Diamond (SCD)')
plt.grid(True)
plt.show()

# Refractive index plot
plt.figure()
plt.plot(l, n, 'b')
plt.xlabel('wavelength')
plt.ylabel('Refractive Index')
plt.title('DSKK estimate for Diamond (SCD)')
plt.grid(True)
plt.show()

