# -*- coding: utf-8 -*-

import numpy as np
import math 
import sys # for debug
import xlrd
import matplotlib.pyplot as plt

print("Program Start !!")
np.seterr(divide='raise') #debug code

freq = 3000
dt = 1.0/freq


# 未入力データチェック及び計算制御
flagThrust = False
flagChamberPressure = False
flagTankPressure = False
flagSupplyPressure = False

#--------------------------------
# 初期値                         -
#--------------------------------
input_book = xlrd.open_workbook('input_file.xlsx')
sheet = input_book.sheet_by_index(0)
# print (sheet.cell(13,1).value)

cartridge_type = sheet.cell(1,1).value             # カートリッジタイプ？？
fuel_type = sheet.cell(3,1).value                   # 燃料タイプ？？


if(cartridge_type.count('L') == 1):
  cartridge_type = 'L'
elif(cartridge_type.count('K') == 1):
  cartridge_type = 'K'
elif(cartridge_type.count('J') == 1):
  cartridge_type = 'J'
else:
  print("Error: incorrect cartridge type.")
  sys.exit()

if(fuel_type.count("3") == 1):
  fuel_type = 3
elif(fuel_type.count('5') == 1):
  fuel_type = 5
elif(fuel_type.count('7') == 1):
  fuel_type = 7
else:
  print("Error: incorrect fuel type.")
  sys.exit()



injector_diamiter = sheet.cell(15,1).value         # インジェクタ径[mm]
cartridge_length = sheet.cell(16,1).value         # インジェクタ長さ[mm]
insulator_inner_diameter = sheet.cell(17,1).value  # インシュレータ内径[mm]
insulator_length = sheet.cell(18,1).value          # インシュレータ長さ[mm]
bp_inner_diameter = sheet.cell(19,1).value         # B.P内径[mm]
bp_depth = sheet.cell(20,1).value                  # B.P深さ[mm]
bp_length = sheet.cell(21,1).value                 # B.P全長[mm]
bp_plate_thickness = sheet.cell(22,1).value        # プレート部厚さ[mm]
central_hole_diameter = sheet.cell(23,1).value     # 中心孔径[mm]
surrouding_hole_diameter = sheet.cell(24,1).value # 周囲孔径[mm]
surrouding_hole_number = sheet.cell(25,1).value       # 周囲孔数[個]
acc_inner_diameter = sheet.cell(26,1).value        # ACC内径[mm]
acc_length = sheet.cell(27,1).value               # ACC長さ[mm]
engaging_length = sheet.cell(28,1).value            # 噛み合わせ部長さ[mm]
throat_diameter = sheet.cell(29,1).value           # スロート径[mm]
oxidizer_volume = sheet.cell(7,4).value          # 酸化剤搭載容積[cc]
fuel_outer_diameter = sheet.cell(9,4).value*10**(-3)   # 燃料外径[m]
initial_port_diameter = sheet.cell(10,4).value*10**(-3) # 初期ポート径[m]
fuel_length = sheet.cell(11,4).value*10**(-3)          # 燃料長さ[m]
fuel_mass = sheet.cell(12,4).value                # 燃料質量[g]
a = sheet.cell(13,4).value                       # a係数
n = sheet.cell(14,4).value                        # n係数
Before_All_Mass = sheet.cell(15,4).value*10**(-3)
After_All_Mass = sheet.cell(16,4).value*10**(-3)
Before_FGC_Mass = sheet.cell(17,4).value*10**(-3)
After_FCG_Mass = sheet.cell(18,4).value*10**(-3)
Before_ACC_Mass = sheet.cell(19,4).value*10**(-3)
After_ACC_Mass = sheet.cell(20,4).value*10**(-3)
freq = sheet.cell(22,4).value
Pa = sheet.cell(23,4).value                      # 大気圧
T = sheet.cell(24,4).value
kappa = sheet.cell(26,4).value                   # κappa
Pe_Pe = sheet.cell(27,4).value          # Pe/Pe


#↓これは後で自動化したい↓
actuation_start_time = sheet.cell(7,1).value          # 作動開始 
combustion_start_time = sheet.cell(8,1).value         # 燃焼開始
combustion_stable_start_time = sheet.cell(9,1).value # 燃焼安定開始
combustion_end_time = sheet.cell(10,1).value          # 燃焼終了
actuation_end_time = sheet.cell(11,1).value           # 作動終了
calibration_start_time = sheet.cell(12,1).value         # 較正開始
calibration_end_time = sheet.cell(13,1).value            # 較正終了

actuation_start_index = math.ceil(actuation_start_time * freq)
combustion_start_index = math.ceil(combustion_start_time * freq)
combustion_stable_start_index = math.ceil(combustion_stable_start_time * freq)
combustion_end_index = math.ceil(combustion_end_time * freq)
actuation_end_index = math.ceil(actuation_end_time * freq)
calibration_start_index = math.ceil(calibration_start_time * freq)
calibration_end_index = math.ceil(calibration_end_time * freq)





#--------------------------------
# データロード                    -
#--------------------------------
print("\nLoading Data...\n")
# url = input('Data >> ')
url = 'data.csv'
Log = np.loadtxt(url, delimiter=',', skiprows=1)

t = Log[:,0]
Thrust = Log[:,1]
ChamberPressure = Log[:,2]
TankPressure = Log[:,3]
SupplyPressure = Log[:,4]
# TankMass = Log[:,5]

length = len(t)

#--------------------------------
# 未入力データチェック             -
#--------------------------------
for i in range(length):
  if(not Thrust[i] == 0):
    flagThrust = True

  if(not ChamberPressure[i] == 0):
    flagChamberPressure = True

  if(not TankPressure[i] == 0):
    pass# flagTankPressure = True #未対応のためpass

  if(not SupplyPressure[i] == 0):
    flagSupplyPressure = True

print("-------------------------------")
print("- CheckData                   -")
print("-------------------------------")
print("Thrust:          ",flagThrust)
print("ChamberPressure: ",flagChamberPressure)
print("TankPressure:    ",flagTankPressure)
print("SupplyPressure:  ",flagSupplyPressure)
print("-------------------------------\n")


#--------------------------------
# データ較正                     -
#--------------------------------
Thrust_calibration = 0
ChamberPressure_calibration = 0
TankPressure_calibration = 0
SupplyPressure_calibration = 0
# TankMass_calibration = 0

calibration_index = calibration_start_index
while(calibration_index < calibration_end_index):
  Thrust_calibration += Thrust[calibration_index]
  ChamberPressure_calibration += ChamberPressure[calibration_index]
  TankPressure_calibration += TankPressure[calibration_index]
  SupplyPressure_calibration += SupplyPressure[calibration_index]
  # TankMass_calibration += TankMass[calibration_index]
  
  calibration_index += 1

Thrust_calibration /= (calibration_end_index - calibration_start_index)
ChamberPressure_calibration /= (calibration_end_index - calibration_start_index)
TankPressure_calibration /= (calibration_end_index - calibration_start_index)
SupplyPressure_calibration /= (calibration_end_index - calibration_start_index)
# TankMass_calibration /= (calibration_end_index - calibration_start_index)

for i in range(length):
  Thrust[i] -= Thrust_calibration
  ChamberPressure[i] += -ChamberPressure_calibration + Pa
  # TankPressure[i] += -TankPressure_calibration + Pa
  SupplyPressure[i] += -SupplyPressure_calibration + Pa
  # TankMass[i] -= TankMass_calibration
# 較正完了--------------------------------------------------------


#- 作動時間までのデータを削除 -----------------------------------------------
t = t[actuation_start_index:actuation_end_index]
Thrust = Thrust[actuation_start_index:actuation_end_index]
ChamberPressure = ChamberPressure[actuation_start_index:actuation_end_index]
TankPressure = TankPressure[actuation_start_index:actuation_end_index]
SupplyPressure = SupplyPressure[actuation_start_index:actuation_end_index]
# TankMass = TankMass[actuation_start_index:actuation_end_index]
length = len(t)
for i in range(length): #時間を作動時間=0にセット
  t[i] = i*dt
# 各時間再計算
actuation_start_time = 0.0       # 作動開始 
combustion_start_time = (combustion_start_index - actuation_start_index)/freq      # 燃焼開始
combustion_stable_start_time = (combustion_stable_start_index - actuation_start_index)/freq
combustion_end_time = combustion_start_time + (combustion_end_index - combustion_start_index)/freq      # 燃焼終了
actuation_end_time = length/freq  # 作動終了
# 各インデックス再計算
actuation_start_index = int(0)
combustion_start_index = int(combustion_start_time*freq)
combustion_stable_start_index = int(combustion_stable_start_time*freq)
combustion_end_index = int(combustion_end_time*freq)
actuation_end_index = int(actuation_end_time*freq)
#------------------------------------------------------------------------

#--------------------------------
# 初期計算                       -
#--------------------------------
injector_area = np.pi * injector_diamiter**2 / 4.0 * 10**(-6)
if cartridge_type == 'J': 
  cartridge_inner_diameter = 44.0
elif cartridge_type == 'K':
  cartridge_inner_diameter = 52.0
elif cartridge_type == 'L':
  cartridge_inner_diameter = 69.0
elif cartridge_type == 'M':
  cartridge_inner_diameter = 92.0
else:
  print('Error: incorrect cartridge type.')
  sys.exit()
mean_supply_pressure = np.average(SupplyPressure[combustion_start_index:combustion_end_index]) # 平均酸化剤供給圧
cartridge_volume = np.pi * cartridge_inner_diameter**2 / 4.0 * cartridge_length
insulator_volume = np.pi * cartridge_inner_diameter**2 / 4.0 - np.pi * insulator_inner_diameter**2 / 4.0 * insulator_length
fuel_volume = (2/3*np.pi*(1/math.tan(np.deg2rad(8)))*(21.875**2+((fuel_outer_diameter*10**3/2)*((fuel_outer_diameter*10**3-2)/2))+((fuel_outer_diameter*10**3-2)/2)**2)+((np.pi*(fuel_outer_diameter*10**3)**2/4)*(fuel_length*10**3-(2/math.tan(np.deg2rad(8))))))-((np.pi*(initial_port_diameter*10**3)**2/4)*fuel_length*10**3)
hole_volume = np.pi * central_hole_diameter**2 / 4.0 + np.pi * surrouding_hole_diameter**2 / 4.0 * bp_plate_thickness * surrouding_hole_number
bp_volume = np.pi * bp_inner_diameter**2 / 4.0 * bp_plate_thickness + hole_volume
acc_volume = np.pi * acc_inner_diameter**2 / 4.0 * acc_length
fuel_bp = np.pi * cartridge_inner_diameter**2 / 4.0 * engaging_length
bp_acc_nozzle = np.pi * acc_inner_diameter**2 / 4.0 *engaging_length * 2
chamber_volume = (cartridge_volume + bp_volume + acc_volume - insulator_volume - fuel_volume - fuel_bp - bp_acc_nozzle) * 10**(-9)
throat_area = (np.pi * throat_diameter**2 / 4.0) * 10**(-6)
L_astarisk = chamber_volume / throat_area
oxidizer_initial_density = 0.031*mean_supply_pressure**6 - 0.8787*mean_supply_pressure**5 + 9.4257*mean_supply_pressure**4 - 50.329*mean_supply_pressure**3 + 144.29*mean_supply_pressure**2 - 280.32*mean_supply_pressure + 1244.5
fuel_density = fuel_mass / fuel_volume * 10**6 #!!fuel_volumeは別途計算が必要かも
epsilon = ((2/(kappa+1))**(1/(kappa-1)))/((Pe_Pe**(1/kappa))*math.sqrt(((kappa+1)/(kappa-1))*(1-(Pe_Pe**((kappa-1)/kappa)))))

flow_history = np.zeros(length)
oxidizer_mass_flow_rate = np.zeros(length)
Average_Cf = np.zeros(length)
oxidizer_mass_flow_flux = np.zeros(length)
port_diameter = np.zeros(length)
fuel_recession_speed = np.zeros(length)
fuel_mass_flow_rate = np.zeros(length)
propellant_mass_flow_rate = np.zeros(length)
of_rate = np.zeros(length)
c_astarisk = np.zeros(length)
theory_c_astarisk = np.zeros(length)
c_astarisk_rate = np.zeros(length)
specific_impulse_C = np.zeros(length)
specific_impalse = np.zeros(length)
Pe = np.zeros(length)
Cf_root = np.zeros(length)
Cf_pressure = np.zeros(length)
theory_Cf_history = np.zeros(length)
eta_Cf = np.zeros(length)


#--------------------------------
# main calculation              -
#--------------------------------

print("Calulating...\n")



# 流量履歴
for i in range(length):
  if (combustion_start_time < t[i]):
    flow_history[i] = injector_area * math.sqrt(2.0*oxidizer_initial_density*abs(SupplyPressure[i]-ChamberPressure[i])*10**6)
  else:
    flow_history[i] = 0
    
# 流量係数
flow_coefficient = (oxidizer_volume*oxidizer_initial_density/10**6)/(flow_history[combustion_start_index:combustion_end_index].sum()/freq)

for i in range(length):
  # 酸化剤質量流量
  oxidizer_mass_flow_rate[i] = flow_coefficient * flow_history[i]
  
  # Average Cf
  if (combustion_start_time < t[i]):
    Average_Cf[i] = Thrust[i]/(throat_area*ChamberPressure[i]*10**6)
  else:
    Average_Cf[i] = 0
    
  if (combustion_start_time < t[i]):
    # 酸化剤流量流束 !!ポート直径(n-1)
    oxidizer_mass_flow_flux[i] = oxidizer_mass_flow_rate[i] / (np.pi*(0.5*port_diameter[i-1])**2)
    # 燃料後退速度 !!酸化剤流量速束
    fuel_recession_speed[i] = (a*oxidizer_mass_flow_flux[i]**n)/1000
    # print(fuel_recession_speed[i])
    # print(a,n)
  else:
    oxidizer_mass_flow_flux[i] = 0
    fuel_recession_speed[i] = 0

  if (combustion_start_time <= t[i]):
    # ポート直径
    port_diameter[i] = 2*(fuel_recession_speed[i]*t[i] + initial_port_diameter/2)
  else:
    port_diameter[i] = 0  

  # 燃料質量流量      
  fuel_mass_flow_rate[i] = 2 * np.pi * fuel_length * fuel_density * fuel_recession_speed[i] * 0.5*port_diameter[i]

  # 推進剤質量流量
  propellant_mass_flow_rate[i] = oxidizer_mass_flow_rate[i] + fuel_mass_flow_rate[i]

  # O/F
  if(combustion_start_time < t[i]):
    of_rate[i] = oxidizer_mass_flow_rate[i] / fuel_mass_flow_rate[i]
  else:
    of_rate[i] = 0

  # 特性排気速度
  if(combustion_start_time < t[i]):
    c_astarisk[i] = (ChamberPressure[i] * throat_area * 10.0**6) / propellant_mass_flow_rate[i]

  # 理論特性排気速度
  if(combustion_start_time < t[i]):
    if (fuel_type == 3):
      theory_c_astarisk[i] = -0.05041*of_rate[i]**6 + 1.446427*of_rate[i]**5 - 15.7928*of_rate[i]**4 + 81.54317*of_rate[i]**3 -215.172*of_rate[i]**2 + 404.7733*of_rate[i] + 892.1382
    elif (fuel_type == 5):
      theory_c_astarisk[i] = -0.04795*of_rate[i]**6 + 1.386782*of_rate[i]**5 - 15.2965*of_rate[i]**4 + 80.03966*of_rate[i]**3 - 214.696*of_rate[i]**2 + 407.7471*of_rate[i] + 887.4279 
    elif (fuel_type == 7):
      theory_c_astarisk[i] = -0.04668*of_rate[i]**6 + 1.355494*of_rate[i]**5 - 15.029*of_rate[i]**4 + 79.19622*of_rate[i]**3 - 214.555*of_rate[i]**2 + 410.3958*of_rate[i] + 886.1227
    else:
      print('Error: incorrect fuel type')
    
  # 特性排気速度効率
  if (combustion_start_time < t[i]):
    c_astarisk_rate[i] = c_astarisk[i] / theory_c_astarisk[i]

  # 比推力C*,Cf
  specific_impulse_C[i] = c_astarisk[i] * c_astarisk_rate[i] * Average_Cf[i]

  # 比推力 
  if (combustion_start_time < t[i]):
    specific_impalse[i] = Thrust[i] / (9.8*propellant_mass_flow_rate[i])
  else:
    specific_impalse[i] = 0
    
  # Pe
  Pe[i] = ChamberPressure[i]*Pe_Pe
  
  # Cf_ルート項
  Cf_root[i] = math.sqrt(((2*(kappa**2)/(kappa-1))*((2/(kappa+1))**((kappa+1)/(kappa-1)))*(1-(Pe[i]/ChamberPressure[i])**((kappa-1)/kappa))))
  
  # Cf_圧力項
  Cf_pressure[i] = ((Pe[i]-Pa)/ChamberPressure[i])*epsilon
  
  # 理論Cf時間履歴
  theory_Cf_history[i] = Cf_root[i] + Cf_pressure[i]
  
  # ηCf
  eta_Cf[i] = Average_Cf[i] / theory_Cf_history[i]



mean_Thrust = np.average(Thrust[combustion_start_index:combustion_end_index])
total_impulse = np.sum(Thrust[combustion_start_index:actuation_end_index])/freq
mean_combustion_pressure = np.average(ChamberPressure[combustion_start_index:combustion_end_index])

print("-------------------------------")
print("- Result                      -")
print("-------------------------------")

print("インジェクタ径:  ",injector_diamiter)
print("スロート径: 　　 ",throat_diameter)
print("燃焼時間: 　　　 ",combustion_end_time-combustion_start_time)
print("作動時間: 　　　 ",actuation_end_time-actuation_start_time)
print("平均推力:        ",mean_Thrust)
print("全力積:          ",total_impulse)
print("酸化剤供給圧:　  ",mean_supply_pressure)
print("燃焼圧:          ",mean_combustion_pressure)
print("流量係数: 　　　 ",flow_coefficient)
print("酸化剤質量流量:  ",np.average(oxidizer_mass_flow_rate[combustion_start_index:combustion_end_index]))
print("燃料質量流量:  　",np.average(fuel_mass_flow_rate[combustion_start_index:combustion_end_index]))
print("推進剤質量流量:  ",np.average(propellant_mass_flow_rate[combustion_start_index:combustion_end_index]))
print("Average O/F:   　",np.average(of_rate[combustion_start_index:combustion_end_index]))
print("酸化剤質量流束:  ",np.average(oxidizer_mass_flow_flux[combustion_start_index:combustion_end_index]))
print("燃料後退速度: 　 ",np.average(fuel_recession_speed[combustion_start_index:combustion_end_index])*1000)
print("スライバ率:　　  ",(1-(((port_diameter[combustion_end_index])**2-initial_port_diameter**2)/(fuel_outer_diameter**2-initial_port_diameter**2)))*100)
print("推力係数:        ",np.average(Average_Cf[combustion_stable_start_index:combustion_end_index]))
print("推力係数効率:　  ",np.average(eta_Cf[combustion_stable_start_index:combustion_end_index])*100)
print("比推力:        　",total_impulse/((Before_All_Mass-After_All_Mass)+(oxidizer_volume/10**6*oxidizer_initial_density))/9.8)
print("比推力効率:　　  ",np.average(eta_Cf[combustion_stable_start_index:combustion_end_index])*np.average(c_astarisk[combustion_start_index:combustion_end_index])/np.average(theory_c_astarisk[combustion_start_index:combustion_end_index])*100)
print("特性排気速度:　  ",np.average(c_astarisk[combustion_start_index:combustion_end_index]))
print("特性排気速度効率:",np.average(c_astarisk[combustion_start_index:combustion_end_index])/np.average(theory_c_astarisk[combustion_start_index:combustion_end_index])*100)
print("--------------------------------\n")

output_file = open('output.csv','w',encoding="shift_jis")
output_file.write("燃焼時間,"+str(combustion_end_time-combustion_start_time)+"\n")
output_file.write("作動時間,"+str(actuation_end_time-actuation_start_time)+"\n")
output_file.write("酸化剤質量流量,"+str(np.average(oxidizer_mass_flow_rate[combustion_start_index:combustion_end_index]))+"\n")
output_file.write("燃料質量流量,"+str(np.average(fuel_mass_flow_rate[combustion_start_index:combustion_end_index]))+"\n")
output_file.write("推進剤質量流量,"+str(np.average(propellant_mass_flow_rate[combustion_start_index:combustion_end_index]))+"\n")
output_file.write("Average O/F,"+str(np.average(of_rate[combustion_start_index:combustion_end_index]))+"\n")
output_file.write("燃料後退速度,"+str(np.average(fuel_recession_speed[combustion_start_index:combustion_end_index])*1000)+"\n")
output_file.write("平均推力,"+str(mean_Thrust)+"\n")
output_file.write("全力積,"+str(total_impulse)+"\n")
output_file.write("平均比推力,"+str(total_impulse/((Before_All_Mass-After_All_Mass)+(oxidizer_volume/10**6*oxidizer_initial_density))/9.8)+"\n")
output_file.close()


print("Calculation Compleate !!\n")

fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(12,5))
plt.subplots_adjust(wspace=0.4, hspace=0.4)


axes[0,0].plot(t,Thrust, color="red", linewidth=0.5,label='Thrust')
axes[0,0].set_xlabel('Time[sec]')
axes[0,0].set_ylabel('Thrust[N]')
axes[0,0].legend(loc=1,fontsize=8) 
axes[0,0].grid(True)

axes[0,1].plot(t, ChamberPressure, linewidth=0.5,label='ChamberPressure')
axes[0,1].plot(t, SupplyPressure, linewidth=0.5,label='SupplyPressure')
axes[0,1].plot(t, TankPressure, linewidth=0.5,label='TankPressure')
axes[0,1].set_xlabel('Time[sec]')
axes[0,1].set_ylabel('Pressure[MPa]')
axes[0,1].legend(loc=1,fontsize=8)
axes[0,1].grid(True)

axes[0,2].plot(t, c_astarisk, linewidth=0.5,label='c*')
axes[0,2].set_xlabel('Time[sec]')
axes[0,2].set_ylabel('C*')
axes[0,2].set_ylim(0,2000)
axes[0,2].legend(loc=1,fontsize=8)
axes[0,2].grid(True)

axes[1,0].plot(t, specific_impalse, linewidth=0.5,label='Isp')
axes[1,0].set_xlabel('Time[sec]')
axes[1,0].set_ylabel('Isp')
axes[1,0].set_ylim(0,250)
axes[1,0].legend(loc=1,fontsize=8)
axes[1,0].grid(True)

axes[1,1].plot(t, of_rate, linewidth=0.5,label='O/F')
axes[1,1].set_xlabel('Time[sec]')
axes[1,1].set_ylabel('O/F')
axes[1,1].set_ylim(0,10)
axes[1,1].legend(loc=1,fontsize=8)
axes[1,1].grid(True)

axes[1,2].plot(t, oxidizer_mass_flow_rate, linewidth=0.5,label='Oxidizer Mass Flow rate')
axes[1,2].plot(t, fuel_mass_flow_rate, linewidth=0.5,label='Fuel Mass Flow rate')
axes[1,2].set_xlabel('Time[sec]')
axes[1,2].set_ylabel('Mass Flow rate[kg/sec]')
axes[1,2].set_ylim(0,0.5)
axes[1,2].legend(loc=1,fontsize=8)
axes[1,2].grid(True)


plt.show()