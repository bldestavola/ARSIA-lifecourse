cd "N:\_COMMITMENTS\LC-REVIEW\data"
adopath ++ "N:\_GFORMULA\version 1.14 beta"
*cd "S:\_COMMITMENTS\LC-REVIEW\data"
*adopath ++ "S:\_GFORMULA\version 1.14 beta"

cap log close
clear all
set more off
set seed 210404
log using example2_growth_interventional_200421.log,replace

use example, clear


*INTERVENTIONAL EFFECTS FOR AVAILABLE BMI MEASURES
*NON LINEAR MODELS

************************
*G-formula- MONTE CARLO STEPS
***********************
capture program drop gform_growth_interv
program define gform_growth_interv, rclass

*MONTE CARLO STEP
preserve
expand 1000
sort id
qui by id:gen original=_n==1
						
****************************** models for M1**************************************************
reg	l_bmi7 std_bw std_bw2 c1 c2 c3 c4 matBMI inter1_std_bw inter2_std_bw inter3_std_bw inter4_std_bw inter1_matBMI ///
    if original==1
	reg,coeflegend
	
cap drop U7
qui gen U7=rnormal()

foreach v of numlist 0/1{		  
#delimit ;
	gen bmi7_`v'`v'=_b[_cons] +_b[std_bw]*`v'+_b[std_bw2]*`v'^2  
              +_b[c1]*c1+_b[c2]*c2+_b[c3]*c3+_b[c4]*c4+_b[matBMI]*matBMI
			  +_b[inter1_std_bw]*c1*`v'+_b[inter2_std_bw]*c2*`v'+_b[inter3_std_bw]*c3*`v'
			  +_b[inter4_std_bw]*c4*`v'+_b[inter1_matBMI]*matBMI*c1
			  +U7*e(rmse);
#delimit cr
}
su bmi7_*

reg	l_bmi9 l_bmi7 std_bw std_bw2 c1 c2 c3 c4 matBMI inter1_std_bw inter2_std_bw inter3_std_bw inter4_std_bw inter1_matBMI ///
    if original==1
cap drop U9
qui gen U9=rnormal()
foreach v of numlist 0/1{		  
#delimit ;
	gen bmi9_`v'`v'=_b[_cons] +_b[std_bw]*`v'+_b[std_bw2]*`v'^2  
	           +_b[l_bmi7]*bmi7_`v'`v' 
              +_b[c1]*c1+_b[c2]*c2+_b[c3]*c3+_b[c4]*c4+_b[matBMI]*matBMI
			  +_b[inter1_std_bw]*c1*`v'+_b[inter2_std_bw]*c2*`v'+_b[inter3_std_bw]*c3*`v'
			  +_b[inter4_std_bw]*c4*`v'+_b[inter1_matBMI]*matBMI*c1
			  +U9*e(rmse);
#delimit cr
}

su U*
su bmi*_*


*************M2************************************
*conditional models
reg	l_bmi10 l_bmi9 l_bmi7 std_bw std_bw2 c1 c2 c3 c4 matBMI inter1_std_bw inter2_std_bw inter3_std_bw inter4_std_bw inter1_matBMI ///
    if original==1
cap drop U10
qui gen U10=rnormal()
foreach v of numlist 0/1{	
#delimit ;
	gen bmi10_`v'`v'_cond=_b[_cons] +_b[std_bw]*`v'+_b[std_bw2]*`v'^2  
	          +_b[l_bmi7]*bmi7_`v'`v' +_b[l_bmi9]*bmi9_`v'`v' 
              +_b[c1]*c1+_b[c2]*c2+_b[c3]*c3+_b[c4]*c4+_b[matBMI]*matBMI
			  +_b[inter1_std_bw]*c1*`v'+_b[inter2_std_bw]*c2*`v'+_b[inter3_std_bw]*c3*`v'
			  +_b[inter4_std_bw]*c4*`v'+_b[inter1_matBMI]*matBMI*c1
			  +U10*e(rmse);
#delimit cr
}

reg	l_bmi11 l_bmi10 l_bmi9 l_bmi7 std_bw std_bw2 c1 c2 c3 c4 matBMI inter1_std_bw inter2_std_bw inter3_std_bw inter4_std_bw inter1_matBMI ///
    if original==1
cap drop U11
qui gen U11=rnormal()
foreach v of numlist 0/1{	
#delimit ;
	gen bmi11_`v'`v'_cond=_b[_cons] +_b[std_bw]*`v'+_b[std_bw2]*`v'^2  
	           +_b[l_bmi7]*bmi7_`v'`v' +_b[l_bmi9]*bmi9_`v'`v' 
			   +_b[l_bmi10]*bmi10_`v'`v'_cond 
               +_b[c1]*c1+_b[c2]*c2+_b[c3]*c3+_b[c4]*c4+_b[matBMI]*matBMI
			   +_b[inter1_std_bw]*c1*`v'+_b[inter2_std_bw]*c2*`v'+_b[inter3_std_bw]*c3*`v'
			   +_b[inter4_std_bw]*c4*`v'+_b[inter1_matBMI]*matBMI*c1
			   +U11*e(rmse);
#delimit cr
}


reg	l_bmi12 l_bmi11 l_bmi10 l_bmi9 l_bmi7 std_bw std_bw2 c1 c2 c3 c4 matBMI inter1_std_bw inter2_std_bw inter3_std_bw inter4_std_bw inter1_matBMI ///
    if original==1
cap drop U12
qui gen U12=rnormal()
foreach v of numlist 0/1{	
#delimit ;
	gen bmi12_`v'`v'_cond=_b[_cons] +_b[std_bw]*`v'+_b[std_bw2]*`v'^2  
	           +_b[l_bmi7]*bmi7_`v'`v' +_b[l_bmi9]*bmi9_`v'`v' 
			   +_b[l_bmi10]*bmi10_`v'`v'_cond+_b[l_bmi11]*bmi11_`v'`v'_cond 
              +_b[c1]*c1+_b[c2]*c2+_b[c3]*c3+_b[c4]*c4+_b[matBMI]*matBMI
			  +_b[inter1_std_bw]*c1*`v'+_b[inter2_std_bw]*c2*`v'+_b[inter3_std_bw]*c3*`v'
			  +_b[inter4_std_bw]*c4*`v'+_b[inter1_matBMI]*matBMI*c1
			  +U12*e(rmse);
#delimit cr
}


*************M2' ************************************
*marginal models
reg	l_bmi10 std_bw std_bw2 c1 c2 c3 c4 matBMI inter1_std_bw inter2_std_bw inter3_std_bw inter4_std_bw inter1_matBMI ///
    if original==1
cap drop U10
qui gen U10=rnormal()
foreach v of numlist 0/1{	
#delimit ;
	gen bmi10_`v'`v'_marg=_b[_cons] +_b[std_bw]*`v'+_b[std_bw2]*`v'^2  
              +_b[c1]*c1+_b[c2]*c2+_b[c3]*c3+_b[c4]*c4+_b[matBMI]*matBMI
			  +_b[inter1_std_bw]*c1*`v'+_b[inter2_std_bw]*c2*`v'+_b[inter3_std_bw]*c3*`v'
			  +_b[inter4_std_bw]*c4*`v'+_b[inter1_matBMI]*matBMI*c1
			  +U10*e(rmse);
#delimit cr
}

reg	l_bmi11 l_bmi10  std_bw std_bw2 c1 c2 c3 c4 matBMI inter1_std_bw inter2_std_bw inter3_std_bw inter4_std_bw inter1_matBMI ///
    if original==1
cap drop U11
qui gen U11=rnormal()
foreach v of numlist 0/1{	
#delimit ;
	gen bmi11_`v'`v'_marg=_b[_cons] +_b[std_bw]*`v'+_b[std_bw2]*`v'^2  
	          +_b[l_bmi10]*bmi10_`v'`v'_marg 
              +_b[c1]*c1+_b[c2]*c2+_b[c3]*c3+_b[c4]*c4+_b[matBMI]*matBMI
			  +_b[inter1_std_bw]*c1*`v'+_b[inter2_std_bw]*c2*`v'+_b[inter3_std_bw]*c3*`v'
			  +_b[inter4_std_bw]*c4*`v'+_b[inter1_matBMI]*matBMI*c1
			  +U11*e(rmse);
#delimit cr
}


reg	l_bmi12 l_bmi11 l_bmi10 std_bw std_bw2 c1 c2 c3 c4 matBMI inter1_std_bw inter2_std_bw inter3_std_bw inter4_std_bw inter1_matBMI ///
    if original==1
cap drop U12
qui gen U12=rnormal()
foreach v of numlist 0/1{	
#delimit ;
	gen bmi12_`v'`v'_marg=_b[_cons] +_b[std_bw]*`v'+_b[std_bw2]*`v'^2  
	          +_b[l_bmi10]*bmi10_`v'`v'_marg +_b[l_bmi11]*bmi11_`v'`v'_marg  
              +_b[c1]*c1+_b[c2]*c2+_b[c3]*c3+_b[c4]*c4+_b[matBMI]*matBMI
			  +_b[inter1_std_bw]*c1*`v'+_b[inter2_std_bw]*c2*`v'+_b[inter3_std_bw]*c3*`v'
			  +_b[inter4_std_bw]*c4*`v'+_b[inter1_matBMI]*matBMI*c1
			  +U12*e(rmse);
#delimit cr
}

su bmi*_*_*


******************************MODEL for Y ***************************************************
#delimit ;
reg	y std_bw std_bw2 l_bmi7 l_bmi9 l_bmi10 l_bmi11 l_bmi12 l_bmi2_12
       c1 c2 c3 c4 matBMI inter1_std_bw inter2_std_bw inter3_std_bw inter4_std inter1_matBMI 
	   inter1_l_bmi12 if original==1;
#delimit cr		   

cap drop UY
qui gen UY=rnormal()

cap drop y_00* y_01* y_10* y_11*

foreach m in cond marg{
  foreach a of numlist 0,1 {     /* index for X */
   foreach b of numlist 0,1 {   /* index for M_1   */
    foreach c of numlist 0,1 {  /* index for M_2 */

   #delimit;
	qui gen y_`a'`b'`c'_`m'
			 =_b[_cons] +_b[std_bw]*`a'+_b[std_bw2]*`a'^2      
	           +_b[l_bmi7]*bmi7_`b'`b' +_b[l_bmi9]*bmi9_`b'`b' 
			   +_b[l_bmi10]*bmi10_`c'`c'_`m'+_b[l_bmi11]*bmi11_`c'`c'_`m' 
			   +_b[l_bmi12]*bmi12_`c'`c'_`m'+_b[l_bmi2_12]*(bmi12_`c'`c'_`m')^2   
			   +_b[c1]*c1+_b[c2]*c2+_b[c3]*c3+_b[c4]*c4 +_b[matBMI]*matBMI    
	           +_b[inter1_std_bw]*c1*`a'+_b[inter2_std_bw]*c2*`a'+_b[inter3_std_bw]*c3*`a'
			   +_b[inter4_std_bw]*c4*`a'+_b[inter1_matBMI]*matBMI*c1
			   +_b[inter1_l_bmi12]*c1*bmi12_`c'`c'_`m'
			   +e(rmse)*UY;
	#delimit cr	
   }
  }
 }
}	

su y_*

********************************RESULTS**********************************************************************

*DE: only X set to exposed versus none, conditional draws
	qui gen d7=y_100_cond-y_000_cond

*IE_1: X and M1 set to exposed versus only X set to exposed, marginal draws
	qui gen d8=y_110_marg-y_100_marg
	
*IE_2: X and M2 set to exposed versus X only set to exposed, marginal draws
	qui gen d9=y_101_marg-y_100_marg
	
*MD: the extra bit
	qui gen d10=(y_111_cond-y_111_marg)-(y_100_cond-y_100_marg)

*TOTAL MEDIATED
	qui gen nie=y_111_cond-y_100_cond

*TOTAL
	qui gen tce=y_111_cond-y_000_cond

 su tce* nie* d*	


*POST THE RESULTS

 qui summ tce
	return scalar tce=r(mean)

	qui summ nie
	return scalar nie=r(mean)

	qui summ d7
	return scalar d7=r(mean)

	qui summ d8
	return scalar d8=r(mean)

	qui summ d9
	return scalar d9=r(mean)

	qui summ d10
	return scalar d10=r(mean)

restore

	
end

gform_growth_interv 

bootstrap tce=r(tce) nie=r(nie) de=r(d7) nie1=r(d8) nie2=r(d9)  rem=r(d10), reps(2): gform_growth_interv

ex

