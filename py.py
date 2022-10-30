from cgitb import grey
from cmath import sqrt
import math
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from skimage import io, color

def CIEDE2000(Lab_1, Lab_2) :
    C_25_7 = 6103515625 # 25**7
    
    L1, a1, b1 = Lab_1[0], Lab_1[1], Lab_1[2]
    L2, a2, b2 = Lab_2[0], Lab_2[1], Lab_2[2]
    C1 = math.sqrt(a1**2 + b1**2)
    C2 = math.sqrt(a2**2 + b2**2)
    C_ave = (C1 + C2) / 2
    G = 0.5 * (1 - math.sqrt(C_ave**7 / (C_ave**7 + C_25_7)))
    
    L1_, L2_ = L1, L2
    a1_, a2_ = (1 + G) * a1, (1 + G) * a2
    b1_, b2_ = b1, b2
    
    C1_ = math.sqrt(a1_**2 + b1_**2)
    C2_ = math.sqrt(a2_**2 + b2_**2)
    
    if b1_ == 0 and a1_ == 0: h1_ = 0
    elif a1_ >= 0: h1_ = math.atan2(b1_, a1_)
    else: h1_ = math.atan2(b1_, a1_) + 2 * math.pi
    
    if b2_ == 0 and a2_ == 0: h2_ = 0
    elif a2_ >= 0: h2_ = math.atan2(b2_, a2_)
    else: h2_ = math.atan2(b2_, a2_) + 2 * math.pi

    dL_ = L2_ - L1_
    dC_ = C2_ - C1_    
    dh_ = h2_ - h1_
    if C1_ * C2_ == 0: dh_ = 0
    elif dh_ > math.pi: dh_ -= 2 * math.pi
    elif dh_ < -math.pi: dh_ += 2 * math.pi        
    dH_ = 2 * math.sqrt(C1_ * C2_) * math.sin(dh_ / 2)
    
    L_ave = (L1_ + L2_) / 2
    C_ave = (C1_ + C2_) / 2
    
    _dh = abs(h1_ - h2_)
    _sh = h1_ + h2_
    C1C2 = C1_ * C2_
    
    if _dh <= math.pi and C1C2 != 0: h_ave = (h1_ + h2_) / 2
    elif _dh  > math.pi and _sh < 2 * math.pi and C1C2 != 0: h_ave = (h1_ + h2_) / 2 + math.pi
    elif _dh  > math.pi and _sh >= 2 * math.pi and C1C2 != 0: h_ave = (h1_ + h2_) / 2 - math.pi 
    else: h_ave = h1_ + h2_
    
    T = 1 - 0.17 * math.cos(h_ave - math.pi / 6) + 0.24 * math.cos(2 * h_ave) + 0.32 * math.cos(3 * h_ave + math.pi / 30) - 0.2 * math.cos(4 * h_ave - 63 * math.pi / 180)
    
    h_ave_deg = h_ave * 180 / math.pi
    if h_ave_deg < 0: h_ave_deg += 360
    elif h_ave_deg > 360: h_ave_deg -= 360
    dTheta = 30 * math.exp(-(((h_ave_deg - 275) / 25)**2))
    
    R_C = 2 * math.sqrt(C_ave**7 / (C_ave**7 + C_25_7))  
    S_C = 1 + 0.045 * C_ave
    S_H = 1 + 0.015 * C_ave * T
    
    Lm50s = (L_ave - 50)**2
    S_L = 1 + 0.015 * Lm50s / math.sqrt(20 + Lm50s)
    R_T = -math.sin(dTheta * math.pi / 90) * R_C

    k_L, k_C, k_H = 1, 1, 1
    
    f_L = dL_ / k_L / S_L
    f_C = dC_ / k_C / S_C
    f_H = dH_ / k_H / S_H
    
    dE_00 = math.sqrt(f_L**2 + f_C**2 + f_H**2 + R_T * f_C * f_H)
    return dE_00
def colorDiff(px1 ,px2) :
    res = CIEDE2000(px1,px2)
    return res
def histoSens(histo) :
    dist = 0
    avrg = np.average(histo)
    for i in histo :
        dist += (abs(i-avrg))
    return dist
def histoNorm(histo) :
    return [i/sum(histo) for i in histo]
def getSens(img ,sensScale) :
    histoR = np.zeros(16)
    histoG = np.zeros(16)
    histoB = np.zeros(16)

    for i in range(img.shape[0]) :
        for j in range(img.shape[1]) :
            histoR[img[i][j][0] // 16] += 1
            histoG[img[i][j][1] // 16] += 1 
            histoB[img[i][j][2] // 16] += 1
        loading = round(((i+1)/img.shape[0]) * 100 , 2)
        print (f'Calculating sens: {loading}%', end="\r") 

    histoR = histoNorm(histoR)
    histoG = histoNorm(histoG)
    histoB = histoNorm(histoB)

    return (histoSens(histoR)*0.299 + histoSens(histoG)*0.587 + histoSens(histoB)*0.144) * sensScale 
def borders(img, dist):
    input('')
    temp = []
    for i in range(1, img.shape[0] - 1):
        for j in range(1, img.shape[1] - 1):
            up = False
            left = False
            #up = colorDiff(img[i][j], img[i + 1][j]) > dist
            down = colorDiff(img[i][j], img[i - 1][j]) > dist
            right = colorDiff(img[i][j], img[i][j + 1]) > dist
            #left = colorDiff(img[i][j], img[i][j - 1]) > dist
            if (up or down or left or right):
                temp.append({'i': i, 
                            'j': j,
                            'col' : (max(colorDiff(img[i][j], img[i][j + 1]),colorDiff(img[i][j], img[i - 1][j])) - dist) * 16 })
        loading = round(((i+1)/img.shape[0]) * 100 , 2)
        print (f'Extracting boders: {loading}%', end="\r")
    return temp

#importing image
img = io.imread('face.jpg')

#calculating colors sensitivity based on RGB valeus
sens  = getSens(img , 8)
print(f'sense for this image is {sens}')

#converting image to CIELAB
img = color.rgb2lab(img)

#extracting borders pixels
border = borders(img ,sens)

#creating white image same size as imported image
plain = np.full((img.shape[0], img.shape[1], 3), 255, dtype=np.uint8)

#draw the borders on tha plain image
for i in border:
    grey = 255-min(255,int(i['col'])) 
    plain[i['i']][i['j']] = [grey, grey, grey]

#show the result
plt.imshow(plain)
plt.show()
print(f'image borders sucsesfully extracted !')


