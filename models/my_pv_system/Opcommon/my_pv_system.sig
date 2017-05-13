[AcqGr1]
1=my_pv_system/SM_master/port1(1)|signal1.Irradiance(W/m^2)|1|1|6
10=my_pv_system/SM_master/port2(4)|signal1.Vgrid|2|4|11
11=my_pv_system/SM_master/port2(5)|signal1.Ig|2|5|11
12=my_pv_system/SM_master/port2(6)|signal1.P_PV|2|6|11
13=my_pv_system/SM_master/port2(7)|signal1.P_PV_MEAN|2|7|11
14=my_pv_system/SM_master/port2(8)|signal1.Gate(1)|2|8|11
15=my_pv_system/SM_master/port2(9)|signal1.Gate(2)|2|9|11
16=my_pv_system/SM_master/port2(10)|signal1.Gate(3)|2|10|11
17=my_pv_system/SM_master/port2(11)|signal1.Gate(4)|2|11|11
2=my_pv_system/SM_master/port1(2)|signal1.I_PV|1|2|6
3=my_pv_system/SM_master/port1(3)|signal1.V_PV|1|3|6
4=my_pv_system/SM_master/port1(4)|signal1.I_Diode|1|4|6
5=my_pv_system/SM_master/port1(5)|signal1.Temperature|1|5|6
6=my_pv_system/SM_master/port1(6)|signal1.Vph|1|6|6
7=my_pv_system/SM_master/port2(1)|signal1.Vinv|2|1|11
8=my_pv_system/SM_master/port2(2)|signal1.Vdc|2|2|11
9=my_pv_system/SM_master/port2(3)|signal1.Igrid|2|3|11
nbsignals=17
[AcqGr2]
1=my_pv_system/SM_master/port3(1)|signal1.Computation time|1|1|4
2=my_pv_system/SM_master/port3(2)|signal1.Real step size|1|2|4
3=my_pv_system/SM_master/port3(3)|signal1.Idle time|1|3|4
4=my_pv_system/SM_master/port3(4)|signal1.Overrun times|1|4|4
5=my_pv_system/SS_controller/port2(1)|signal1.Computation time|2|1|4
6=my_pv_system/SS_controller/port2(2)|signal1.Real step size|2|2|4
7=my_pv_system/SS_controller/port2(3)|signal1.Idle time|2|3|4
8=my_pv_system/SS_controller/port2(4)|signal1.Overrun times|2|4|4
nbsignals=8
[SM_master]
1=my_pv_system/SC_user_interface/port1(1)|signal1.Irradiance_con|1|1|2
2=my_pv_system/SC_user_interface/port1(2)|signal1.Temperature_con|1|2|2
nbsignals=2
