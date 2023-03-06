#********************************************************************************************************
#------------This code is based on the parametrized trajectory----------------------------------
#------------method given by T. Su., J. Chem. Phys. 100, 4703 (1994).---------------------------
#------------Proton transfer reaction rates can be computed provided the VOCs have--------------
#------------proton affinity higher than H3O or NH4 and for electron transfer rates if----------
#------------the ionization energy less than NO and O2 molecules, respectively.-----------------
#------------Please leave your comments at: quantumsimm@gmail.com-------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------Instructions to run the code--------------------------------------
#--------------------put CSV file in the same folder as the code and run: python rate-file.py---
#********************************************************************************************************
import math
from math import*
import pandas as pd

from tkinter import *
window = Tk()
window.title("Chemical Ionization Mass Spectrometry Rate Calculation")
window.geometry("900x540")
#window.config(bg='#9A9ACD')

label5 = Label(window, text = "Rate Coefficient Calculation",fg='#000080',font=("Times New Roman",18)).place(x = 35,y = 10)
label7 = Label(window, text = "Results",fg='#000080',font=("Times New Roman",18)).place(x = 600,y = 10)
#label5.grid(padx=5,pady=10)
#------------------------------------------text1-------------------------------
label1 = Label(window, text = "Enter CAS Number", fg='blue',font=("Times New Roman",14)).place(x = 15,y = 50)
#label1.grid(padx=5,pady=35)
label8 = Label(window, text = "Ex. 87-40-1", fg='blue',font=("Times New Roman",14)).place(x = 380,y = 50)
label8 = Label(window, text = "Td", fg='blue',font=("Times New Roman",14)).place(x = 380,y = 90)
label8 = Label(window, text = "K", fg='blue',font=("Times New Roman",14)).place(x = 380,y = 130)

label2 = Label(window, text = "Reduced Field (E/N)", fg='blue',font=("Times New Roman",14)).place(x = 15,y = 90)
#label2.grid(padx=5,pady=20)

label3 = Label(window, text = "Enter Temperature (T)", fg='blue',font=("Times New Roman",14)).place(x = 15,y = 130)
#label3.grid(padx=5,pady=5)


data=StringVar()    #IntVar(), DoubleVar()
data1=StringVar()
data2=StringVar() 
#-------------------------------------radio button--------------
label4 = Label(window, text = "Select reagent ion:",fg='blue',font=("Times New Roman",14)).place(x = 15,y = 180)
def viewSelected():
    global res
    choice = var.get()
    if choice == 1:
        output = 19.02
    elif choice == 2:
        output = 18.04
    elif choice == 3:
        output = 30.01
    elif choice == 4:
        output = 31.99
    else:
        output = "Invalid selection"
    res = output
    return Label(window).place(x = 430,y = 150)

var=IntVar()
Radiobutton(window, text="h3o+",variable=var,value=1,command=viewSelected,font=("Times New Roman",14)).place(x = 70,y = 220)
Radiobutton(window, text="nh4+",variable=var,value=2,command=viewSelected,font=("Times New Roman",14)).place(x = 70,y = 250)
Radiobutton(window, text="no+",variable=var,value=3,command=viewSelected,font=("Times New Roman",14)).place(x = 70,y = 280)
Radiobutton(window, text="o2+",variable=var,value=4,command=viewSelected,font=("Times New Roman",14)).place(x = 70,y = 310)
#--------------------------------------------------------------

def myfunction():
    cas = data.get()
    field = data1.get()
    temp = data2.get()
    field1 = int(field)
    temp1 = int(temp) 
#    emptylabel.config(text="CAS number of Molecule :"+cas)
#    emptylabel1.config(text="Electric field (E/N) : {field1}")
#    emptylabel2.config(text="Temperature (K) : {temp1}")
    df = pd.read_csv('data.csv')
    dataframe = df.to_html(justify='center')
    #print(dataframe)
    #print(df)
    rslt_df = df[df['CAS'] == cas]
    rslt_df1 = rslt_df.to_string(index = False)
    blankIndex=[''] * len(rslt_df)
    rslt_df.index=blankIndex
    #print(rslt_df)
    m = rslt_df['Mass']
    #print(m)
    dipole_moment1 = rslt_df['Dipole Moment']
    dipole_moment = float(dipole_moment1)
    #print(dipole_moment)
    #dipole_moment = "{:.2f}".format(dipole_moment)
    plz1 = rslt_df['Polarizability']
    plz = float(plz1)
    #print(plz)
    ebyn = field1
    temp = temp1
    value12 = res
#    value12 = 18.03                              #mass of ion
    #print(value12)    
    alpha = plz*10**(-30)
    DM = dipole_moment*3.336*10**(-30)
    #==============================velocity from E/N====v = K_0*N_0*E/N=========
    K_o = 2.81 #cm^2 V^-1 S^-1
    N_o = 2.687*10**(19) # 10 19 cm −3
    #EoverN = float(input("\n Enter the reduced electric field (in Td): "))  #1 Td=10^(-17) V cm^2
    redu_field = ebyn*10**(-17)
    v = K_o*N_o*redu_field*(10**(-2))         #v in ms^(-1)
    #print ("drift velocity is :", v, "ms^-1")
    #================================extra conversion factors===================
    mmol = m
    mion = value12
    T = temp
    vel = v*v
    avg_n_o = 28.8
    mass_conv = 1.67*10**(-27)
    k  = 1.38*10**-23
    #=================================KE-ion=====================================
    m1 = mion*mass_conv
    M1 = avg_n_o*mass_conv
    E1 = 0.5*(m1*vel)
    E2 = 0.5*(M1*vel)
    E3 = 1.5*(k*T)
    E3_ev = E3*6.24*10**(18)            #conversion to eV
    KE_ion = E1+E2+E3
    KE_ion_ev = KE_ion*6.24*10**(18)
    KE_ion_ev1 = "{:.4f}".format(KE_ion_ev)
    #print(KE_ion_ev1)
    #================================KE-COM======================================
    ms = (mmol)/(mmol+mion)
    ke = (KE_ion_ev-E3_ev)
    KE_cm = E3_ev+(ms*ke)
    KE_cm2 = float(KE_cm)
    KE_cm1 = "{:.4f}".format(KE_cm2)
    #print(KE_cm1)
    #==================================k-Lang====================================
    mu = (mion*mmol)/(mion+mmol)
    mu1 = mu*1.67*10**(-27)
    epsilon = 8.854*10**(-12)     # J^-1.C^2.m^-1
    q = 1.602*10**(-19)           # C
    #mu = 2.38*10**(-26)           # reduced mass
    part1 = math.pi*alpha*q**2
    part2 = mu1*epsilon
    K_L = math.sqrt(part1/part2)
    K_L1 =10**(6)*K_L

    
    #res  = "{:.2f}".format(v)
    kinetic_energy_ion = "{:.2f}".format(KE_ion_ev)
    kinetic_energy_com = KE_cm
    #rate_Lang = "{:.2e}".format(K_L1)
    rate_Lang = K_L1
    #print(rate_Lang)
    #===============================T_effective==============================
    second_term = vel/(3*k)                                                 #average of N2 & O2
    mass = (mmol)*(mion+avg_n_o)/(mion+mmol)
    mass_1 = mass*mass_conv
    T_eff1 = T+second_term*mass_1
    T_eff = int(T_eff1)
    #print(T_eff)
    #==============================values-neede===========================
    c1 = 0.727143
    c2 = 3.71823
    c3 = 0.586920
    c4 = 4.97894
    tau = dipole_moment/math.sqrt(plz*T)
    eps = dipole_moment/math.sqrt(plz*KE_cm)
    #print(eps)
    theta = c3*(c4+math.log(tau))
    #==============================k_cap-final================================
    if eps > 1.5:
        S = math.exp(-2*(eps-1.5))
        a = tau**(0.4)
        b = eps*eps
        K_c1 = c1*a*b*S
        K_c2 = c2*(1-S)
        K_c3 = math.sin(theta)
        K_c4 = tau**(0.6)
        K_c5 = math.sqrt(eps-0.5)
        K_c6 = K_c2*K_c3*K_c4*K_c5
        K_c = 1+K_c1+K_c6
        K_cap = K_L1*K_c
#        print ("..............>>>>>\n      The Capture rate constant k_cap is :  ", K_cap, "10^(-9)cm^3 s^-1", "@", i, "m/s")
    elif eps <= 1.5:
        S = 1
        a = tau**(0.4)
        b = eps*eps
        K_c1 = c1*a*b*S
        K_c2 = c2*(1-S)*math.sin(theta)
        K_c = 1+K_c1+K_c2
        K_cap = K_L1*K_c
#        print ("..............>>>>>\n      The Capture rate constant k_cap is :  ", K_cap, "10^(-9)cm^3 s^-1", "@", i, "m/s")

    #print(K_cap)
    #K_cap = "{:.2e}".format(K_cap)
    #print(K_cap)
    emptylabel3.config(text="Dipole Moment of the Molecule: %.2f"" Debye\r\n" % dipole_moment)
    emptylabel4.config(text="Polarizability of the Molecule: %.2f"" Ang^3\r\n" % plz)
    emptylabel5.config(text="Mass of Ion : %.2f"" au\r\n" % value12)
    emptylabel6.config(text="Mass of the Molecule : %.2f"" au\r\n" % mmol)
    emptylabel7.config(text="KE of the Ion : %.4f"" eV\r\n" % KE_ion_ev)
    emptylabel8.config(text="Center-of-mass KE : %.4f"" eV\r\n" % KE_cm2)
    emptylabel9.config(text="Effective Temperature : %.0f"" K\r\n" % T_eff)
    emptylabel10.config(text="Langevin Rate Coefficient : %.2e"" cm^3 s^-1\r\n" % rate_Lang)
    emptylabel11.config(text="Parametrized Rate Coefficient : %.2e"" cm^3 s^-1\r\n" % K_cap)
    
def cleartext():
    textbox1.delete(0, "end")
    textbox2.delete(0, "end")
    textbox3.delete(0, "end")

def close():
    window.destroy()

button1=Button(window,command=myfunction,text="Submit",fg='blue',font=("Times New Roman",14)).place(x = 40,y = 360)
#button3=Button(window,command=cleartext,text="Clear",fg='blue',font=("Times New Roman",14)).place(x = 180,y = 360)
button4=Button(window,command=close,text="Exit",fg='blue',font=("Times New Roman",14)).place(x = 320,y = 360)
#button1.grid(row=9,column=1,sticky=W)
#---------------------------------------text2---------------------------------
textbox1=Entry(window,textvariable=data,fg='black',font=("Times New Roman",14)).place(x = 185,y = 50)
#textbox1.grid(padx=5,pady=35)

textbox2=Entry(window,textvariable=data1,fg='black',font=("Times New Roman",14)).place(x = 185,y = 90)
#textbox2.grid(padx=5,pady=20)

textbox3=Entry(window,textvariable=data2,fg='black',font=("Times New Roman",14)).place(x = 185,y = 130)
#textbox3.grid(padx=5,pady=10)
#============================================================================================
'''emptylabel=Label(window,fg='green',font=("Arial",14))
#emptylabel.grid(row=5,column=1,sticky=W,pady=10)
emptylabel.place(x = 430,y = 50)

emptylabel1=Label(window,fg='green',font=("Arial",14))
#emptylabel1.grid(row=6,column=1,sticky=W,pady=10)
emptylabel1.place(x = 430,y = 80)

emptylabel2=Label(window,fg='green',font=("Arial",14))
#emptylabel2.grid(row=8,column=1,sticky=W,pady=10)
emptylabel2.place(x = 430,y = 110)'''

emptylabel3=Label(window,fg='green',font=("Times New Roman",14))
#emptylabel2.grid(row=8,column=1,sticky=W,pady=10)
emptylabel3.place(x = 490,y = 50)

emptylabel4=Label(window,fg='green',font=("Times New Roman",14))
#emptylabel2.grid(row=8,column=1,sticky=W,pady=10)
emptylabel4.place(x = 490,y = 80)

emptylabel5=Label(window,fg='green',font=("Times New Roman",14))
#emptylabel2.grid(row=8,column=1,sticky=W,pady=10)
emptylabel5.place(x = 490,y = 110)


emptylabel6=Label(window,fg='green',font=("Times New Roman",14))
#emptylabel2.grid(row=8,column=1,sticky=W,pady=10)
emptylabel6.place(x = 490,y = 140)

emptylabel7=Label(window,fg='green',font=("Times New Roman",14))
#emptylabel2.grid(row=8,column=1,sticky=W,pady=10)
emptylabel7.place(x = 490,y = 170)

emptylabel8=Label(window,fg='green',font=("Times New Roman",14))
#emptylabel2.grid(row=8,column=1,sticky=W,pady=10)
emptylabel8.place(x = 490,y = 200)

emptylabel9=Label(window,fg='green',font=("Times New Roman",14))
#emptylabel2.grid(row=8,column=1,sticky=W,pady=10)
emptylabel9.place(x = 490,y = 230)

emptylabel10=Label(window,fg='green',font=("Times New Roman",14))
#emptylabel2.grid(row=8,column=1,sticky=W,pady=10)
emptylabel10.place(x = 490,y = 260)

emptylabel11=Label(window,fg='green',font=("Times New Roman",14))
#emptylabel2.grid(row=8,column=1,sticky=W,pady=10)
emptylabel11.place(x = 490,y = 290)


label9 = Label(window, text = "---------------------------------------------------------------------------------------------------------------------------------------------", fg='#2F2F4F',font=("Times New Roman",14)).place(x = 25,y = 410)
label6 = Label(window, text = "Currently available molecules with data can be found in .csv file attached. User can update .csv file with \n dipole moment and polarizability values and get the required parameters.",fg='#2F2F4F',font=("Times New Roman",14)).place(x = 50,y = 440)
label11 = Label(window, text = "*Dipole moment and polarizability values are computed using DFT calculations \n for molecules with no experimental values available.",fg='#2F2F4F',font=("Times New Roman",10)).place(x = 450,y = 360)
#button2=Button(window,text="click here",fg='blue',font=("Times New Roman",14)).place(x = 30,y = 470)
label10 = Label(window, text = "------------------------------Parametrized trajectory method, T. Su., J. Chem. Phys. 100, 4703 (1994)---------------------------------", fg='#4A4A70',font=("Times New Roman",12)).place(x = 50,y = 500)

#=============================================result box=============================================
#display = Text(window,fg='green',font=("Arial",18)).place(x = 450, y = 50, width=404, height=401)
#display.pack(padx=21,pady=34)
#Creating a text box widget
#display_box=Text(window,height=21, width=34).place(x = 430,y = 50)
#==================================================================================================
window.mainloop()
