
"""
V 1.3
1.这个脚本主要的目标是计数300的瞬时辐射的观测
2.挑出BAT和XRT的数据，去掉绿色，显示BAT推导到的数据与XRT的数据
3. 具有不同仪器谱指数区分功能：BAT的谱指数用-2，XRT的谱指数用-2
4. 计数MDP时，若有XRT的数据，则忽略BAT的数据

"""

import numpy as np
import matplotlib as plt
import os
import re
from scipy import integrate
import matplotlib.pyplot as plt
from fpdf import FPDF
pdf = FPDF()

t_start = 1.0
t_end = 300
miu = 0.3     #调制因子
eff = 0.1     # 探测效率
area = 300    # 有效探测面积
bat_photon_index = 2.0  #BAT平均谱指数，假设是powerlaw谱
xrt_photon_index = 2.0  #XRT平均谱指数，假设是powerlaw谱

def prompt_MDP(dataname):
    # ---------------数据读取
    dat=[]
    bat = []
    xrt = []
    path_png = '../GRB_prompt_MDP(3)_300s/'
    num = 0
    with open(dataname, 'r') as f:
        bat_start=bat_end=xrt_start=xrt_end=0
        for line in f.readlines():
            num += 1
            if 'batSNR5flux' in line:
                bat_start = num
            if 'batSNR5gamma' in line:
                bat_end = num-2
            if 'xrtwtflux' in line:
                xrt_start = num
            if 'xrtwtgamma' in line:
                xrt_end = num-2
    with open(dataname, 'r') as f:
        for line in f.readlines()[bat_start:bat_end]:
            bat.append(re.split(r'\s+', line))
        bat = np.array(bat)
        if len(bat) == 0 :
            pass
        else:
            bat = bat[:, :-1]
            bat = bat.astype(np.float)
    with open(dataname, 'r') as f:
        for line in f.readlines()[xrt_start:xrt_end]:
            xrt.append(re.split(r'\s+', line))
        xrt = np.array(xrt)
        xrt = xrt[:, :-1]
        xrt = xrt.astype(np.float)

    #---------------------- 转换光子数
    def N_count(flux, index=xrt_photon_index):
        N = flux * integrate.quad(lambda E: E ** (-index), 2, 10.0)[0] / \
            integrate.quad(lambda E: E * E ** (-index), 0.3, 10)[0] / 1.6e-9
        return N

    #print('%3.3f count/cm2/s'%N_count(2.4e-8))
    if len(bat) == 0 :
            pass
    else:
        flux = np.array(bat[:, 3])
        flux_err = np.array(bat[:, 4])
        bat_count = [N_count(flux[i],index=bat_photon_index) for i in range(len(flux))]
        bat_count_err = [N_count(flux_err[i],index=bat_photon_index) for i in range(len(flux_err))]
        bat = np.column_stack((bat, bat_count, bat_count_err))
    
    xrt_flux = np.array(xrt[:, 3])
    xrt_flux_err = np.array(xrt[:, 4])
    xrt_count = [N_count(xrt_flux[i]) for i in range(len(xrt_flux))]
    xrt_count_err = [N_count(xrt_flux_err[i]) for i in range(len(xrt_flux_err))]
    xrt = np.column_stack((xrt, xrt_count, xrt_count_err))
    #print(bat_count)
    # ---------------------------画图
    fig, ax = plt.subplots()
    if len(bat) == 0 :
            pass
    else:
        x = bat[:, 0]
        xerr = bat[:, 1]
        xerr_ = bat[:, 2]
        y = bat[:, -2]
        yerr = bat[:, -1]
        plt.errorbar(x, y, yerr=yerr, xerr=xerr, fmt='o',label='BAT(flux to count)')
    
    xrt_x = xrt[:, 0]
    xrt_xerr = xrt[:, 1]
    xrt_xerr_ = xrt[:, 2]
    xrt_y = xrt[:, -2]
    xrt_yerr = xrt[:, -1]
    plt.errorbar(xrt_x, xrt_y, yerr=xrt_yerr, xerr=xrt_xerr, fmt='o',color='red',label='XRT(flux to count)')
    plt.xlabel('Time since BAT trigger (s)')
    plt.ylabel(r'2-10 keV (Count/cm$^2$/s)')
    plt.title('Swift BAT-XRT data of %s'%dataname)
    plt.loglog()
    
    #-----------------------------------合并BAT与XRT反推光子的数据
    xrt_x0= xrt[0, 0]
    print (xrt_x0)
    if len(bat) == 0 :
            bat = xrt
    else:
        bat = bat[bat[:,0] < xrt_x0, :]
        bat = np.row_stack((bat,xrt))
    
    x = bat[:, 0]
    xerr = bat[:, 1]
    xerr_ = bat[:, 2]
    y = bat[:, -2]
    yerr = bat[:, -1]
  #  plt.errorbar(x, y, yerr=yerr, xerr=xerr, fmt='o',label='BAT')

    # -----------------------------------计算总MPD
    t = xerr-xerr_
    N_cm = sum(y*t)
    N_total = N_cm*eff*area
  #  print(N_cm, N_total)
    MDP = 4.29/(miu*np.sqrt(N_total))*100
  #  print(MDP)
    # ----------------------------t_start 之后观测到的数据点及画图
    bat2 = bat[bat[:, 0] > t_start, :]
    bat2 = bat2[bat2[:, 0] < t_end, :]
    #print(bat2)
    x2 = bat2[:, 0]
    xerr2 = bat2[:, 1]
    xerr2_ = bat2[:, 2]
    y2 = bat2[:, 6]
    yerr2 = bat2[:, 7]
   # plt.errorbar(x2, y2, yerr=yerr2, xerr=xerr2, fmt='*')
    plt.axvline(t_start, label='t=%s s'%t_start, color='green')
    plt.axvline(t_end, label='t=%s s'%t_end, color='green')
    # ----------------------------t_start 之后观测到的MDP
    t2 = xerr2-xerr2_
    N_cm2 = sum(y2*t2)
    N_total2 = N_cm2*eff*area
   # print(N_cm2, N_total2)
    MDP2 = 4.29/(miu*np.sqrt(N_total2))*100
  #  print(MDP2)
    fig.text(0.2, 0.2, 'MDP = %2.2f %%' % MDP2,color='red', fontsize=12, fontweight='bold')
    plt.legend(loc = 'upper left')
    #plt.legend()
    plt.savefig(path_png+dataname+' %2.2f%%'%MDP2+'.png')
    plt.show()
    return MDP2


def compresspdf(path):
    imagelist = [x for x in os.listdir(path) if os.path.splitext(x)[1] == '.png']
    os.chdir(path)
    tex = 'SXP observed time is %3.2f to %3.2f s;\n' \
          'SXP polarization modulation factor is %2.2f;\n' \
          'SXP detection efficiency is %2.2f\n' \
          'SXP effective area is %3.2f cm^2\n' %(t_start,t_end, miu,eff,area)
    for image in imagelist:
        print(image)
        pdf.add_page()
        pdf.image(image, 0, 0)
        pdf.set_xy(10,180)
        pdf.set_font('Arial','',18)
        pdf.write(10,tex)
    pdf.output('GRB_prompt_MPD.pdf', 'F')
    pdf.close()
    os.chdir('../')

def MDP_infile(path):
    batalist = [x for x in os.listdir(path) if os.path.splitext(x)[1] == '.txt']
    os.chdir(path)
    MDP_sta = []
    for i in range(len(batalist)):
        print(batalist[i])
        MDP_sta.append(prompt_MDP(batalist[i]))
    os.chdir('../')
    return MDP_sta


print('start process')
out_file = 'GRB_prompt_MDP(3)_300s'
if (not(os.path.exists(out_file))) :
    os.mkdir(out_file)
bata_path = 'GRB'
MDP_all = MDP_infile(bata_path)
MDP_png = out_file
compresspdf(MDP_png)

print(MDP_all)
np.savetxt('GRB_prompt_MDP(3)_300s.txt', MDP_all)
