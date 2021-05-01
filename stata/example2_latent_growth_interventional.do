cd "N:\_COMMITMENTS\LC-REVIEW\data"
adopath ++ "N:\_GFORMULA\version 1.14 beta"
*cd "S:\_COMMITMENTS\LC-REVIEW\data"
*adopath ++ "S:\_GFORMULA\version 1.14 beta"

cap log close
clear all
set more off
set seed 210404
log using example2_latent_growth_interventional_200421.log,replace

use example, clear


*INTERVENTIONAL EFFECTS FOR INTERCEPT AND SLOPE USING ALL AVAILABLE BMI MEASURES TO DEFINE ALPHA AND BETA
*NON LINEAR MODELS
************************
*G-formula- MONTE CARLO STEPS
***********************
capture program drop gform_interv
program define gform_interv, rclass

*MONTE CARLO STEP
preserve
expand 1000
sort id
qui by id:gen original=_n==1
qui by id:gen double k=_n
gen double newid=k*1000000+id
su n* id
sort newid
qui by newid: gen counter=_N
ta  counter 
drop counter 

****************************** 1-  model for M: latent growth**************************************************
*Mixed effects model for log(BMI)
reshape long l_bmi l_bmi2_ ,i(newid) j(age)
gen c_age=(age-10)
*gen c_age2=(age-10)^2
gen age_bw=c_age*std_bw
gen age_matBMI=matBMI*(age-10)
*gen inter1_matBMI_age=c1*matBMI*(age-10)
*gen inter1_std_bw_age=c1*std_bw*(age-10)

mixed l_bmi c_age std_bw std_bw2 c1 c2 c3 c4 matBMI age_bw age_matBMI ||newid: c_age if original==1,  cov(unstr) nolog stddev

predict  slope  intercept , reffects
predict  se_slope  se_intercept , reses

sort id k
qui by id: replace slope=slope[1]
qui by id: replace se_slope=se_slope[1]
qui by id: replace intercept=intercept[1]
qui by id: replace se_intercept=se_intercept[1]

matrix S=e(b)
matrix list S

scalar sd_e=exp(S[1,15])

*draw values from the level 2 error distribution
reshape wide l_bmi l_bmi2_ c_age age_bw age_matBMI, i(newid) j(age)
cap drop i s
gen i=rnormal(0,se_intercept)
gen s=rnormal(0,se_slope)
su i s
cap drop alpha
cap drop alpha_*
cap drop beta_*

#delimit ;
	gen alpha=_b[_cons] +_b[std_bw]*std_bw+_b[std_bw2]*std_bw^2  
              +_b[c1]*c1+_b[c2]*c2+_b[c3]*c3+_b[c4]*c4 +_b[matBMI]*matBMI
			  + (intercept+i);
#delimit cr					  
	gen beta =_b[c_age]  +_b[age_bw]*std_bw +_b[age_matBMI]*matBMI+(slope+s)
su alpha beta
corr alpha beta

*************M1************************************
foreach v of numlist 0/1{		  
#delimit ;
	gen alpha_`v'`v'=_b[_cons] +_b[std_bw]*`v'+_b[std_bw2]*`v'^2  
              +_b[c1]*c1+_b[c2]*c2+_b[c3]*c3+_b[c4]*c4+_b[matBMI]*matBMI
			  +(intercept+i);
#delimit cr
}
su alpha_*

*************M2************************************
*conditional model
reg beta alpha std_bw std_bw2 c1 c2 c3 c4 inter1_std_bw inter2_std_bw inter3_std_bw inter4_std_bw inter1_matBMI
cap drop UB
qui gen UB=rnormal()

foreach v of numlist 0/1{		  
	foreach w of numlist 0/1{		  
#delimit ;
	gen beta_`v'`w'_cond=_b[_cons] +_b[alpha]*alpha_`w'`w'
	          +_b[std_bw]*`v'+_b[std_bw2]*`v'^2
              +_b[c1]*c1+_b[c2]*c2+_b[c3]*c3+_b[c4]*c4 
			  +_b[inter1_std_bw]*c1*`v'+_b[inter2_std_bw]*c2*`v'+_b[inter3_std_bw]*c3*`v'
			  +_b[inter4_std_bw]*c4*`v' 
			  +_b[inter1_matBMI]*c1*matBMI
			  +e(rmse)*UB;
#delimit cr
	}
}

*marginal model
reg beta std_bw std_bw2 c1 c2 c3 c4 inter1_std_bw inter2_std_bw inter3_std_bw inter4_std_bw inter1_matBMI
cap drop UB
qui gen UB=rnormal()

foreach v of numlist 0/1{		  
#delimit ;
	gen beta_`v'`v'_marg=_b[_cons] +_b[std_bw]*`v'+_b[std_bw2]*`v'^2  
					+_b[c1]*c1+_b[c2]*c2+_b[c3]*c3+_b[c4]*c4 	
		        	  +_b[inter1_std_bw]*c1*`v'+_b[inter2_std_bw]*c2*`v'
					  +_b[inter3_std_bw]*c3*`v'+_b[inter4_std_bw]*c4*`v' 
					+_b[inter1_matBMI]*c1*matBMI
					+e(rmse)*UB;
#delimit cr
}
su alpha_* beta_*
gen inter1_i=c1*(alpha)
gen inter1_s=c1*(beta)


******************************MODEL for Y ***************************************************
*with estimated intercepts and slopes
reg	y std_bw std_bw2 alpha beta c1 c2 c3 c4 matBMI inter1_std_bw inter2_std_bw inter3_std_bw inter4_std_bw  ///
	       inter1_matBMI inter1_i inter1_s if original==1
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
					  +_b[alpha]*alpha_`b'`b'+ _b[beta]*beta_`c'`c'_`m'
					  +_b[c1]*c1+_b[c2]*c2+_b[c3]*c3+_b[c4]*c4 +_b[matBMI]*matBMI    
					  +_b[inter1_std_bw]*c1*`a'+_b[inter2_std_bw]*c2*`a'+_b[inter3_std_bw]*c3*`a'
					  +_b[inter4_std_bw]*c4*`a'+_b[inter1_matBMI]*matBMI*c1
					  +_b[inter1_i]*alpha_`b'`b'*c1
					  +_b[inter1_s]*beta_`c'`c'_`m'*c1				  
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

gform_interv 

bootstrap tce=r(tce) nie=r(nie) de=r(d7) nie1=r(d8) nie2=r(d9)  rem=r(d10), reps(2): gform_interv

ex
