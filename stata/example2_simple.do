
cd "N:\_COMMITMENTS\LC-REVIEW\data"
*cd "S:\_COMMITMENTS\LC-REVIEW\data"
adopath ++ "N:\_GFORMULA\version 1.14 beta"
*adopath ++ "S:\_GFORMULA\version 1.14 beta"

cap log close
clear all
set more off
log using example2_simple_200421.log,replace

use example,replace

*************************************************************************
*EXAMPLE 2 :  birth weight, all BMI, and BE - for Table 3
*************************************************************************

***********************
*TCE of bw- to check mediation results
***********************
capture program drop tce
program define tce, rclass
   preserve
   expand 2, generate(interv)
   expand 2 if interv == 0, generate(interv2)
   replace interv = -1  if interv2 ==1
   drop interv2 
   tab interv
   replace y = . if interv != -1
   replace std_bw = 0 if interv == 0
   replace std_bw= 1 if interv == 1
   replace std_bw2 = 0 if interv == 0
   replace std_bw2= 1 if interv == 1
 foreach n of numlist 1/4{
	replace inter`n'_std_bw = 0 if interv == 0
   replace inter`n'_std_bw= c`n' if interv == 1
 }
   reg y  std_bw std_bw2 c1 c2 c3 c4 matBMI inter2_-inter4_ inter1_std_bw inter1_matBMI
   predict predY, xb
 
 su predY if(interv == -1)
  scalar obs=r(mean)
  return scalar OBS=r(mean)
  
 su predY if(interv == 0)
  scalar y0=r(mean)
  return scalar Y0=r(mean)
 
 su predY if(interv == 1)
 scalar y1=r(mean)
 return scalar Y1=r(mean)
 
 scalar tce=y1-y0
 return scalar TCE=tce
 scalar list
 restore
 end

tce
bootstrap obs=r(OBS) y0=r(Y0) y1=r(Y1) tce=r(TCE) , reps(1000): tce


***********************************
*only bmi 12 as M, on log scale, non-linear
#delimit;
gformula y c1 c2 c3 c4  matBMI std_bw std_bw2 l_bmi12 l_bmi2_12 inter1_std_bw inter2_std_bw inter3_std_bw
    inter4_std_bw inter1_matBMI inter1_l_bmi12,
mediation outcome(y) exposure(std_bw) mediator(l_bmi12)
base_confs(c1 c2 c3 c4 matBMI) 
 control(l_bmi12:2) 
commands(y:regress,  l_bmi12:regress)
equations(
	y:           std_bw std_bw2 l_bmi12 l_bmi2_12 c1 c2 c3 c4 matBMI inter1_std_bw inter2_std_bw inter3_std_bw
                        inter4_std_bw inter1_matBMI inter1_l_bmi12,
	l_bmi12:     std_bw std_bw2                  c1 c2 c3 c4 matBMI inter1_std_bw inter2_std_bw inter3_std_bw
                        inter4_std_bw inter1_matBMI
	)
linexp	
derived(std_bw2 l_bmi2_12 inter1_std_bw inter2_std_bw inter3_std_bw inter4_std_bw inter1_matBMI inter1_l_bmi12) 
derrules(std_bw2: std_bw^2, l_bmi2_12: l_bmi12^2, inter1_std_bw:c1*std_bw, inter2_std_bw:c2*std_bw, 	  
     inter3_std_bw: c3*std_bw, inter4_std_bw: c4*std_bw, inter1_matBMI: c1*matBMI, inter1_l_bmi12: c1*l_bmi12) 
minsim samples(1000) moreMC simulations(10000) replace seed(2803);
#delimit cr
***********************************
*only bmi 12 as M, on log scale, linear
#delimit;
gformula y c1 c2 c3 c4  matBMI std_bw  l_bmi12 ,
mediation outcome(y) exposure(std_bw) mediator(l_bmi12)
base_confs(c1 c2 c3 c4 matBMI) 
 control(l_bmi12:2) 
commands(y:regress,  l_bmi12:regress)
equations(
	y:           std_bw         l_bmi12          c1 c2 c3 c4 matBMI ,
	l_bmi12:     std_bw                          c1 c2 c3 c4 matBMI	)
linexp	
minsim samples(1000) moreMC simulations(10000) replace seed(2803);
#delimit cr

ex
