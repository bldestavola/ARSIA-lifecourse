cd "N:\_COMMITMENTS\LC-REVIEW\data"
adopath ++ "N:\_GFORMULA\version 1.14 beta"
*cd "S:\_COMMITMENTS\LC-REVIEW\data"
*adopath ++ "S:\_GFORMULA\version 1.14 beta"
cap log close
clear all
set more off
set seed 210804
log using example1_200421.log,replace

use example,replace

*************************************************************************
*EXAMPLE 1 :  birth weight,  BMI AT 12, and BE
*************************************************************************
su l_bmi12 l_bmi2_12
su ht*
egen std_l_bmi7=std(l_bmi7)
egen std_ht7=std(ht7)
su std*

egen std_l_bmi12=std(l_bmi12)
gen std_l_bmi12_2=std_l_bmi12^2
gen std_inter=std_bw*std_l_bmi12
gen inter1_std_l_bmi12=c1* std_l_bmi12
gen inter=std_bw* std_l_bmi12

su inter*



*************************************************************************
*DIRECT EFFECTS OF BW (with MX interaction)

*for CDE_1 at m=0
#delimit;
gformula y c1 c2 c3 c4 std_bw std_bw2 matBMI std_l_bmi12 inter,
mediation outcome(y) exposure(std_bw) mediator(std_l_bmi12)
base_confs(c1 c2 c3 c4 std_bw matBMI)
 control(std_l_bmi12:0) 
commands(y:regress, std_l_bmi12:regress)
equations(
	y:           std_bw std_bw2 std_l_bmi12 c1 c2 c3 c4 matBMI inter ,
	std_l_bmi12: std_bw std_bw2             c1 c2 c3 c4 matBMI  
	)
linexp	
derived(std_bw2 inter) derrules(std_bw2:std_bw^2, inter:(std_bw)*(std_l_bmi12))
minsim samples(1000) moreMC simulations(10000) ;
#delimit cr

*for CDE_1 at +l
#delimit;
gformula y c1 c2 c3 c4 std_bw std_bw2 matBMI std_l_bmi12 inter,
mediation outcome(y) exposure(std_bw) mediator(std_l_bmi12)
base_confs(c1 c2 c3 c4 std_bw matBMI)
 control(std_l_bmi12:1) 
commands(y:regress, std_l_bmi12:regress)
equations(
	y:           std_bw std_bw2 std_l_bmi12 c1 c2 c3 c4 matBMI inter ,
	std_l_bmi12: std_bw std_bw2             c1 c2 c3 c4 matBMI  
	)
linexp	
derived(std_bw2 inter) derrules(std_bw2:std_bw^2, inter:(std_bw)*(std_l_bmi12))
minsim samples(1000) moreMC simulations(10000) ;
#delimit cr

*for CDE_1 at -1
#delimit;
gformula y c1 c2 c3 c4 std_bw std_bw2 matBMI std_l_bmi12 inter,
mediation outcome(y) exposure(std_bw) mediator(std_l_bmi12)
base_confs(c1 c2 c3 c4 std_bw matBMI)
 control(std_l_bmi12:-1) 
commands(y:regress, std_l_bmi12:regress)
equations(
	y:           std_bw std_bw2 std_l_bmi12 c1 c2 c3 c4 matBMI inter ,
	std_l_bmi12: std_bw std_bw2             c1 c2 c3 c4 matBMI  
	)
linexp	
derived(std_bw2 inter) derrules(std_bw2:std_bw^2, inter:(std_bw)*(std_l_bmi12))
minsim samples(1000) moreMC simulations(10000) ;
#delimit crex
*/



*************************************************************************
*DIRECT EFFECTS OF BW (with MX interaction) and additional control for  bmi7

*for CDE_1 at m=0
#delimit;
gformula y c1 c2 c3 c4 std_bw std_bw2 matBMI std_l_bmi12 inter std_l_bmi7  ,
mediation outcome(y) exposure(std_bw) mediator(std_l_bmi12)
base_confs(c1 c2 c3 c4 std_bw matBMI) post_confs(std_l_bmi7 )
 control(std_l_bmi12:0) 
commands(y:regress, std_l_bmi12:regress, std_l_bmi7:regress,   )
equations(
	y:           std_bw std_bw2 std_l_bmi12 c1 c2 c3 c4 matBMI inter std_l_bmi7   ,
	std_l_bmi12: std_bw std_bw2             c1 c2 c3 c4 matBMI       std_l_bmi7  ,
	std_l_bmi7 : std_bw std_bw2             c1 c2 c3 c4 matBMI                   ,
	)
linexp	
derived(std_bw2 inter) derrules(std_bw2:std_bw^2, inter:(std_bw)*(std_l_bmi12))
minsim samples(1000) moreMC simulations(10000) ;
#delimit cr



*for CDE_1 at m=1
#delimit;
gformula y c1 c2 c3 c4 std_bw std_bw2 matBMI std_l_bmi12 inter std_l_bmi7  ,
mediation outcome(y) exposure(std_bw) mediator(std_l_bmi12)
base_confs(c1 c2 c3 c4 std_bw matBMI) post_confs(std_l_bmi7 )
 control(std_l_bmi12:1) 
commands(y:regress, std_l_bmi12:regress, std_l_bmi7:regress,   )
equations(
	y:           std_bw std_bw2 std_l_bmi12 c1 c2 c3 c4 matBMI inter std_l_bmi7   ,
	std_l_bmi12: std_bw std_bw2             c1 c2 c3 c4 matBMI       std_l_bmi7  ,
	std_l_bmi7 : std_bw std_bw2             c1 c2 c3 c4 matBMI                   ,
	)
linexp	
derived(std_bw2 inter) derrules(std_bw2:std_bw^2, inter:(std_bw)*(std_l_bmi12))
minsim samples(1000) moreMC simulations(10000) ;
#delimit cr


*for CDE_1 at m=-1
#delimit;
gformula y c1 c2 c3 c4 std_bw std_bw2 matBMI std_l_bmi12 inter std_l_bmi7  ,
mediation outcome(y) exposure(std_bw) mediator(std_l_bmi12)
base_confs(c1 c2 c3 c4 std_bw matBMI) post_confs(std_l_bmi7 )
 control(std_l_bmi12:-1) 
commands(y:regress, std_l_bmi12:regress, std_l_bmi7:regress,   )
equations(
	y:           std_bw std_bw2 std_l_bmi12 c1 c2 c3 c4 matBMI inter std_l_bmi7   ,
	std_l_bmi12: std_bw std_bw2             c1 c2 c3 c4 matBMI       std_l_bmi7  ,
	std_l_bmi7 : std_bw std_bw2             c1 c2 c3 c4 matBMI                   ,
	)
linexp	
derived(std_bw2 inter) derrules(std_bw2:std_bw^2, inter:(std_bw)*(std_l_bmi12))
minsim samples(1000) moreMC simulations(10000) ;
#delimit cr
ex
*for CDE_1 at m=1

/*************************************************************************
*DIRECT EFFECTS OF BMI12 (with MX interaction)
*for CDE_2: TCE for l_bmi12
*/

capture program drop tce
program define tce, rclass
   preserve
   expand 2, generate(interv)
   expand 2 if interv == 0, generate(interv2)
   replace interv = -1  if interv2 ==1
   drop interv2 
   tab interv
   replace y = . if interv != -1
   replace std_l_bmi12 = 0 if interv == 0
   replace std_l_bmi12 = 1 if interv == 1
   replace std_l_bmi12_2 = 0 if interv == 0
   replace std_l_bmi12_2 = 1 if interv == 1
   replace inter1_std_l_bmi12 = 0 if interv == 0
   replace inter1_std_l_bmi12 = 1 if interv == c1
   replace inter = 0 if interv == 0
   replace inter = 1 if interv == std_bw


reg y std_l_bmi12 std_l_bmi12_2 std_bw std_bw2 c1 c2 c3 c4 matBMI inter2_-inter4_ inter1_std_bw inter1_matBMI inter1_std_l_bmi12 inter
   predict predY, xb
   tabstat predY,by(interv)
  
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


*tce
bootstrap obs=r(OBS) y0=r(Y0) y1=r(Y1) tce=r(TCE) , reps(1000): tce
ex

*create a new variable that takes value 0 at the mean and 1 at 1SD
gen X2=(l_bmi12-3)/.2
scatter X2 l_bmi12, yline(0 1)
medeff (regress y X2 c1 c2 c3 c4 std_bw matBMI) (regress GLUC60 SYSBT40 BMIbin BMISYS male),
mediate(SYSBT40) treat(BMIbin) interact(BMISYS) sims(1000)
			



