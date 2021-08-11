/******************************************************************************/
/* Title: How much does route to diagnosis explain comorbidity 		      */
/*		inequalities in survival of patients with Diffuse Large       */
/*		B-cell lymphoma or Follicular Lymphoma                        */
/* Authors: Matthew J. Smith, Miguel Angel Luque Fernandez                    */
/* Date: 29-04-2021				                              */
/******************************************************************************/

/* Preliminaries */
    set more off
    set autotabgraphs on
    cd "V:\Manuscript"
    use "V:\Manuscript\cancerdata.dta", clear


/* Data preparation */
    gen bcmcat = cmcat
    recode bcmcat (2=1)
    fre bcmcat

    forvalues i = 1/15 {
      tab CM_HES_6_24_cm`i'
    }


* Generate outcome variables
	* Death before 6 months
	  tab dead6m, m

	* Death within 6 and 12 months
    gen dead1y6m = 1 if timesurv>=0.5 & timesurv<1
    recode dead1y6m (.=0) if timesurv>=0.5
    tab dead1y6m, m

	* Death within 12 and 36 months
    tab dead3y1y, m
    sum timesurv if dead3y1y==1
	
	* Death within 36 and 60 months
	  tab dead5y3y, m


/* Descriptive statistics */

  * DLBCL
      preserve
      drop if agediag==.
      keep if subtype==1
        * Age
        sum agediag
        * Gender
        fre sex 
        * Deprivation
        fre dep
        * Ethnicity
        fre bethnic
        * Comorbidity
        fre bcmcat
        * Route
        fre broute
        * Death within 6m
        fre dead6m, m
        * Death within 1y given 6m
        fre dead1y6m, m
        * Death within 3y given 1y
        fre dead3y1y, m
        * Death within 5y given 3y
        fre dead5y3y, m
      restore


  * FL
      preserve
      drop if agediag==.
      keep if subtype==2
        * Age
        sum agediag
        * Gender
        fre sex 
        * Deprivation
        fre dep
        * Ethnicity
        fre bethnic
        * Comorbidity
        fre bcmcat
        * Route
        fre broute
        * Death within 6m
        fre dead6m, m
        * Death within 1y given 6m
        fre dead1y6m, m
        * Death within 3y given 1y
        fre dead3y1y, m
        * Death within 5y given 3y
        fre dead5y3y, m
      restore


/* Probability of emergency route to diagnosis */
    
    * Plot the proportion of cases of diagnostic route
    catplot dep if subtype==1, percent(broute) recast(bar)

    catplot dep broute if subtype==1,  over(bcmcat) subtitle("DLBCL", ring(1) size(vlarge)) ///
      percent(broute bcmcat) recast(bar) /* blabel(bar, format(%4.1f)) */ ///
      b1title("", size(vsmall))  ylabel(0(5)27) ///
      bar(1, fcolor(black*0.2) lcolor(black*0.2)) ///
      bar(2, fcolor(black*0.3) lcolor(black*0.3)) ///
      bar(3, fcolor(black*0.5) lcolor(black*0.5)) ///
      bar(4, fcolor(black*0.7) lcolor(black*0.7)) ///
      bar(5, fcolor(black*0.9) lcolor(black*0.9)) ///
      legend(order(1 "Least deprived" 2 3 4 5 "Most deprived") ring(0) pos(12) ///
            cols(5)	title("", size(vsmall)) size(vsmall) ///
            region(lstyle(none)) nobox symy(3) symx(3)) ///
      scheme(s2mono) graphregion(fcolor(white)) ///
      name(DLBCLroute, replace) saving(DLBCLroute, replace)
	
    catplot dep broute if subtype==2,  over(bcmcat) subtitle("FL", ring(1) size(vlarge)) ///
      percent(broute bcmcat) recast(bar) /* blabel(bar, format(%4.1f)) */ ///
      b1title("", size(vsmall)) ylabel(0(5)27) ///
      bar(1, fcolor(black*0.2) lcolor(black*0.2)) ///
      bar(2, fcolor(black*0.3) lcolor(black*0.3)) ///
      bar(3, fcolor(black*0.5) lcolor(black*0.5)) ///
      bar(4, fcolor(black*0.7) lcolor(black*0.7)) ///
      bar(5, fcolor(black*0.9) lcolor(black*0.9)) ///
      legend(off) ///
      scheme(s2mono) graphregion(fcolor(white)) ///
      name(FLroute, replace) saving(FLroute, replace)
	
  * Combine graphs of DLBCL and FL 
    graph combine DLBCLroute.gph FLroute.gph, c(1)  ///
      xcom com scheme(s2mono) graphregion(fcolor(white)) ///
      name(DLBCLFLroute, replace) saving(DLBCLFLroute, replace)
    graph export DLBCLFLroute.pdf, replace


/* Net Survival estimates */

	// Sex must be coded as 1 and 2 (it will not accept anything else)

    *** First, set the data to be survival time
      // Use the stns command 
        stset timesurv, failure(dead==1) id(id3) 

      // Describe the dataset 
        stdescribe
        *stvary
        stsum, by(subtype)

      // Confidence intervals
        stci 
        stci, p(25)
        stci, rmean by(subtype)

      // Kaplan-Meier by subtype
        sts list
        sts graph, by(subtype) gwood name(NHLKM, replace)

      // Test equality of survivor functions
        sts test subtype


    *** For DLBCL, use the Cohort approach to estimate crude 5-year net survival	
      ** Year of diagnosis 2005-2010 - Cohort approach

      // General net survival at 1-, 3-, and 5-years after diagnosis
        stns list using Lifetable_2013_stns.dta 	///
          if subtype==1 /*& ydiag<=2010*/, ///
          age(agediagindays=_age) /// 
          period(dodiag=yearindays) strata(sex2 dep) rate(rate_day) /// 
          at(1 3 5,  unit(year)) 

      // Cohort survival estimates by comorbidity
        forvalues i = 1/2 {
          forvalues j = 0/1 {
            stns list using Lifetable_2013_stns.dta 	///
              if subtype==`i' & bcmcat==`j' & (dep==1 | dep==5), ///
              age(agediagindays=_age) period(dodiag=yearindays) rate(rate_day) ///
              strata(sex2 dep) by(dep) at(0.5 1 3 5, unit(year))
          }
        }

      // Cohort survival estimates by route
        forvalues i = 1/2 {
          forvalues j = 0/1 {
            stns list using Lifetable_2013_stns.dta 	///
              if subtype==`i' & broute==`j' & (dep==1 | dep==5), ///
              age(agediagindays=_age) period(dodiag=yearindays) rate(rate_day) ///
              strata(sex2 dep) by(dep) at(0.5 1 3 5, unit(year))
          }
        }

      // Graph the net survival for COMORBIDITIES, stratified by deprivation
          * DLBCL
            * For no comorbidities
            stns graph using Lifetable_2013_stns.dta 	///
                if subtype==1 & bcmcat==0 & (dep==1 | dep==5), ///
                age(agediagindays=_age) period(dodiag=yearindays) rate(rate_day) ///
                strata(sex2 dep) by(dep) ///
                ci tma(5) name(DLBCLcm0ns, replace) saving(DLBCLcm0ns, replace) ///
                title("(a) DLBCL") subtitle("No comorbidities", ring(0)) ///
                yti("Net survival probability")	xti("", size(vsmall)) ///
                xlabel(,nolabel) ///
                scheme(s2mono) graphregion(fcolor(white)) ///					
                legend(off)

                *risktable(0 1 3 5, order(1 "Least deprived   " 2 "2.   ") title(" ", size(vsmall)) ) ///

            * For at least one comorbidity
            stns graph using Lifetable_2013_stns.dta 	///
                if subtype==1 & bcmcat==1 & (dep==1 | dep==5), ///
                age(agediagindays=_age) period(dodiag=yearindays) rate(rate_day) ///
                strata(sex2 dep) by(dep) ///
                ci tma(5) name(DLBCLcm1ns, replace) saving(DLBCLcm1ns, replace) ///
                title("", size(vsmall)) subtitle("At least one comorbidity", ring(0)) ///
                yti("Net survival probability")	xti("", size(vsmall)) ///
                xlabel(,nolabel) ///
                scheme(s2mono) graphregion(fcolor(white)) ///					
                legend(off) 

          * FL
            * For no comorbidities
            stns graph using Lifetable_2013_stns.dta 	///
                if subtype==2 & bcmcat==0 & (dep==1 | dep==5), ///
                age(agediagindays=_age) period(dodiag=yearindays) rate(rate_day) ///
                strata(sex2 dep) by(dep) ///
                ci tma(5) name(FLcm0ns, replace) saving(FLcm0ns, replace) ///
                title("(b) FL") subtitle("No comorbidities", ring(0)) ///
                yti("", size(vsmall)) xti("", size(vsmall)) ///
                xlabel(,nolabel) ///
                ylabel(,nolabel) ///
                scheme(s2mono) graphregion(fcolor(white)) ///					
                legend(off) 

                *risktable(0 1 3 5, order(1 "Least deprived   " 2 "2.   ") title(" ", size(vsmall)) ) ///

            * For at least one comorbidity
            stns graph using Lifetable_2013_stns.dta 	///
                if subtype==2 & bcmcat==1 & (dep==1 | dep==5), ///
                age(agediagindays=_age) period(dodiag=yearindays) rate(rate_day) ///
                strata(sex2 dep) by(dep) ///
                ci tma(5) name(FLcm1ns, replace) saving(FLcm1ns, replace) ///
                title("", size(vsmall)) subtitle("At least one comorbidity", ring(0)) ///
                yti("", size(vsmall)) xti("", size(vsmall)) ///
                ylabel(,nolabel) ///	
                xlabel(,nolabel) ///
                scheme(s2mono) graphregion(fcolor(white)) ///					
                legend(off) 

          * DLBCL 
            * For elective 
            stns graph using Lifetable_2013_stns.dta 	///
                if subtype==1 & broute==0 & (dep==1 | dep==5), ///
                age(agediagindays=_age) period(dodiag=yearindays) rate(rate_day) ///
                strata(sex2 dep) by(dep) ///
                ci tma(5) name(DLBCLelns, replace) saving(DLBCLelns, replace) ///
                title("", size(vsmall)) subtitle("Elective route", ring(0)) ///
                yti("Net survival probability")	xti("", size(vsmall)) ///
                xlabel(,nolabel) ///
                scheme(s2mono) graphregion(fcolor(white)) ///					
                legend(off)

                *risktable(0 1 3 5, order(1 "Least deprived   " 2 "2.   ") title(" ", size(vsmall)) ) ///

            * For emergency
            stns graph using Lifetable_2013_stns.dta 	///
                if subtype==1 & broute==1 & (dep==1 | dep==5), ///
                age(agediagindays=_age) period(dodiag=yearindays) rate(rate_day) ///
                strata(sex2 dep) by(dep) ///
                ci tma(5) name(DLBCLemns, replace) saving(DLBCLemns, replace) ///
                title("", size(vsmall)) subtitle("Emergency route", ring(0)) ///
                yti("Net survival probability")	xti("Time (years)") ///
                scheme(s2mono) graphregion(fcolor(white)) ///					
                legend(off) 

          * FL
            * For elective
            stns graph using Lifetable_2013_stns.dta 	///
                if subtype==2 & broute==0 & (dep==1 | dep==5), ///
                age(agediagindays=_age) period(dodiag=yearindays) rate(rate_day) ///
                strata(sex2 dep) by(dep) ///
                ci tma(5) name(FLelns, replace) saving(FLelns, replace) ///
                title("", size(vsmall)) subtitle("Elective route", ring(0)) ///
                yti("", size(vsmall)) xti("", size(vsmall)) ///
                xlabel(,nolabel) ///
                ylabel(,nolabel) ///
                scheme(s2mono) graphregion(fcolor(white)) ///					
                legend(order(5 "Least deprived" 6 "Most deprived" 1 "95% CI") ///
                    size(vsmall) c(1) r(3) subtitle("", size(vsmall)) ///
                    ring(0) position(8) bmargin(large) nobox ///
                    lc(none) fc(white) ///
                    region(lstyle(none)) nobox symy(2) symx(7)) 

                *risktable(0 1 3 5, order(1 "Least deprived   " 2 "2.   ") title(" ", size(vsmall)) ) ///

            * For emergency
            stns graph using Lifetable_2013_stns.dta 	///
                if subtype==2 & broute==1 & (dep==1 | dep==5), ///
                age(agediagindays=_age) period(dodiag=yearindays) rate(rate_day) ///
                strata(sex2 dep) by(dep) ///
                ci tma(5) name(FLcm1ns, replace) saving(FLemns, replace) ///
                title("", size(vsmall)) subtitle("Emergency route", ring(0)) ///
                yti("", size(vsmall)) xti("Time (years)") ///
                ylabel(,nolabel) ///						
                scheme(s2mono) graphregion(fcolor(white)) ///					
                legend(off) 

            * Combine the graphs
            graph combine DLBCLcm0ns.gph DLBCLcm1ns.gph DLBCLelns.gph DLBCLemns.gph ///
                    FLcm0ns.gph FLcm1ns.gph FLelns.gph FLemns.gph, ///
                colf r(2) c(2) ycomm xcomm com ///
                scheme(s2mono) graphregion(fcolor(white)) ysize(7)

            graph export survival.pdf, replace





/******************************************************************************/
/* Mediation analysis */
/******************************************************************************/


/* Firstly for DLBCL */
	
	* Remove missing data for now...
	keep if agediag!=.
	
	* Drop the variables if necessary
	drop ind dpml dpme dpmh dniel dniee dnieh dndel dndee dndeh dtcel dtcee dtceh fpml fpme fpmh fniel fniee fnieh fndel fndee fndeh ftcel ftcee ftceh indte indie indd indf indte indie indd indf indte2 indie2 indd2 indf2 indte3 indie3 indd3 indf3
	
	* Generate the indicator variable
	range ind 1 3 3
	
	* Generate variable for PM estimates
	gen dpml = . 
	gen dpme = . 
	gen dpmh = . 
	
	* Generate variable for NIE estimates
	gen dniel = .
	gen dniee = .
	gen dnieh = .
	
	* Generate variable for NDE estimates
	gen dndel = .
	gen dndee = .
	gen dndeh = .
	
	* Generate variable for TCE estimates
	gen dtcel = .
	gen dtcee = .
	gen dtceh = .
	
	* Generate variable for PM estimates
	gen fpml = . 
	gen fpme = . 
	gen fpmh = . 
	
	* Generate variable for NIE estimates
	gen fniel = .
	gen fniee = .
	gen fnieh = .
	
	* Generate variable for NDE estimates
	gen fndel = .
	gen fndee = .
	gen fndeh = .
	
	* Generate variable for TCE estimates
	gen ftcel = .
	gen ftcee = .
	gen ftceh = .
	
	* Offsets
	local offset 0.1
	gen indte = ind - `offset'
	gen indie = ind + `offset'
	gen indd  = ind - `offset'
	gen indf  = ind + `offset'
			
		
	* Death within 1 year 
	gformula dead1y broute bcmcat sex agediag dep subtype if subtype==1, ///
		mediation out(dead1y) mediator(broute) ex(bcmcat) ///
		samples(1000) minsim moreMC base_confs(sex agediag dep)  ///
		commands(dead1y: logit, broute: logit) ///
		eq(dead1y: broute bcmcat agediag sex i.dep, ///
			broute: bcmcat agediag sex i.dep) ///
		impute(broute) imp_cmd(logit) ///
		imp_eq(broute: dead1y bcmcat agediag sex dep) imp_cycles(10) 	///
		seed(123) logOR obe control(broute:0) 

		replace dpml = r(pm) - 1.96*r(se_pm)  if ind==1
		replace dpme = r(pm) 				          if ind==1
		replace dpmh = r(pm) + 1.96*r(se_pm)  if ind==1
		
		replace dniel = exp(r(nie) - 1.96*r(se_nie))  if ind==1
		replace dniee = exp(r(nie)) 				 	        if ind==1
		replace dnieh = exp(r(nie) + 1.96*r(se_nie))  if ind==1
		
		replace dndel = exp(r(nde) - 1.96*r(se_nde))  if ind==1
		replace dndee = exp(r(nde)) 				 	        if ind==1
		replace dndeh = exp(r(nde) + 1.96*r(se_nde))  if ind==1
		
		replace dtcel = exp(r(tce) - 1.96*r(se_tce))  if ind==1
		replace dtcee = exp(r(tce)) 				 	        if ind==1
		replace dtceh = exp(r(tce) + 1.96*r(se_tce))  if ind==1

	* Death within 3 years given 1 year survival
	gformula dead3y1y broute bcmcat sex agediag dep subtype if subtype==1, ///
		mediation out(dead3y1y) mediator(broute) ex(bcmcat) ///
		samples(1000) minsim moreMC base_confs(sex agediag dep) ///
		commands(dead3y1y: logit, broute: logit) ///
		eq(dead3y1y: broute bcmcat agediag sex i.dep, ///
			broute: bcmcat agediag sex i.dep)	///
		impute(broute) imp_cmd(logit) ///
		imp_eq(broute: dead3y1y bcmcat agediag sex dep) imp_cycles(10) 	///
		seed(123) logOR obe control(broute:0) 		

		replace dpml = r(pm) - 1.96*r(se_pm)  if ind==2
		replace dpme = r(pm) 				          if ind==2
		replace dpmh = r(pm) + 1.96*r(se_pm)  if ind==2

		replace dniel = exp(r(nie) - 1.96*r(se_nie))  if ind==2
		replace dniee = exp(r(nie)) 				 	        if ind==2
		replace dnieh = exp(r(nie) + 1.96*r(se_nie))  if ind==2
		
		replace dndel = exp(r(nde) - 1.96*r(se_nde))  if ind==2
		replace dndee = exp(r(nde)) 				 	        if ind==2
		replace dndeh = exp(r(nde) + 1.96*r(se_nde))  if ind==2
		
		replace dtcel = exp(r(tce) - 1.96*r(se_tce))  if ind==2
		replace dtcee = exp(r(tce)) 				 	        if ind==2
		replace dtceh = exp(r(tce) + 1.96*r(se_tce))  if ind==2
		
	* Death within 5 years given 3 year survival
	gformula dead5y3y broute bcmcat sex agediag dep subtype if subtype==1, ///
		mediation out(dead5y3y) mediator(broute) ex(bcmcat) ///
		samples(1000) minsim moreMC base_confs(sex agediag dep) ///
		commands(dead5y3y: logit, broute: logit) ///
		eq(dead5y3y: broute bcmcat agediag sex i.dep, ///
			broute: bcmcat agediag sex i.dep)	///
		impute(broute) imp_cmd(logit) ///
		imp_eq(broute: dead5y3y bcmcat agediag sex dep) imp_cycles(10) 	///
		seed(123) logOR obe control(broute:0) 

		replace dpml = r(pm) - 1.96*r(se_pm)  if ind==3
		replace dpme = r(pm) 				          if ind==3
		replace dpmh = r(pm) + 1.96*r(se_pm)  if ind==3
		
		replace dniel = exp(r(nie) - 1.96*r(se_nie))  if ind==3
		replace dniee = exp(r(nie)) 				 	        if ind==3
		replace dnieh = exp(r(nie) + 1.96*r(se_nie))  if ind==3
		
		replace dndel = exp(r(nde) - 1.96*r(se_nde))  if ind==3
		replace dndee = exp(r(nde)) 				 	        if ind==3
		replace dndeh = exp(r(nde) + 1.96*r(se_nde))  if ind==3
		
		replace dtcel = exp(r(tce) - 1.96*r(se_tce))  if ind==3
		replace dtcee = exp(r(tce)) 				 	        if ind==3
		replace dtceh = exp(r(tce) + 1.96*r(se_tce))  if ind==3
		
	
		
	/* Plot of natural effects */
	twoway rcap dtcel dtceh indte, lc(red)   || scatter dtcee indte, mc(red)   ///
		|| rcap dndel dndeh ind,   lc(blue)  || scatter dndee ind,   mc(blue)  ///
		|| rcap dniel dnieh indie, lc(green) || scatter dniee indie, mc(green) ///
			xlabel(1 "12" 2 "36" 3 "60")  ///
			xtitle("Conditional months of follow up") ///
			title("DLBCL", size(m)) ///
			legend(order(2 "Total" 4 "Direct" 6 "Indirect") ring(0) pos(11) ///
					cols(1) region(lstyle(none)) nobox ///
					title("Natural effects", size(msmall))) ///
			ytitle("Odds Ratio") ylabel(1(0.2)2.05) yline(0) ///
			scheme(s2mono) graphregion(fcolor(white)) ///
			name(DLBCLallv3, replace)	saving(DLBCLallv3, replace)
			
	graph export DLBCLallv3.pdf, replace
		


/* Now for Follicular lymphomas */
	

	* Death within 1 year
	gformula dead1y broute bcmcat sex agediag dep subtype if subtype==2, ///
		mediation out(dead1y) mediator(broute) ex(bcmcat) ///
		samples(1000) minsim moreMC base_confs(sex agediag dep)  ///
		commands(dead1y: logit, broute: logit) ///
		eq(dead1y: broute bcmcat agediag sex i.dep, ///
			broute: bcmcat agediag sex i.dep) ///
		impute(broute) imp_cmd(logit) ///
		imp_eq(broute: dead1y bcmcat agediag sex dep) imp_cycles(10) 	///
		seed(123) logOR obe control(broute:0) 

		replace fpml = r(pm) - 1.96*r(se_pm)  if ind==1 
		replace fpme = r(pm) 			            if ind==1 
		replace fpmh = r(pm) + 1.96*r(se_pm)  if ind==1 
		
		replace fniel = exp(r(nie) - 1.96*r(se_nie))  if ind==1 
		replace fniee = exp(r(nie)) 				 	        if ind==1 
		replace fnieh = exp(r(nie) + 1.96*r(se_nie))  if ind==1 
		
		replace fndel = exp(r(nde) - 1.96*r(se_nde))  if ind==1
		replace fndee = exp(r(nde)) 				 	        if ind==1 
		replace fndeh = exp(r(nde) + 1.96*r(se_nde))  if ind==1 
		
		replace ftcel = exp(r(tce) - 1.96*r(se_tce))  if ind==1 
		replace ftcee = exp(r(tce)) 				 	        if ind==1 
		replace ftceh = exp(r(tce) + 1.96*r(se_tce))  if ind==1


	* Death within 3 years given 1 year survival
	gformula dead3y1y broute bcmcat sex agediag dep subtype if subtype==2, ///
		mediation out(dead3y1y) mediator(broute) ex(bcmcat) ///
		samples(1000) minsim moreMC base_confs(sex agediag dep) ///
		commands(dead3y1y: logit, broute: logit) ///
		eq(dead3y1y: broute bcmcat agediag sex i.dep, ///
			broute: bcmcat agediag sex i.dep)	///
		impute(broute) imp_cmd(logit) ///
		imp_eq(broute: dead3y1y bcmcat agediag sex dep) imp_cycles(10) 	///
		seed(123) logOR obe control(broute:0)		
		
		replace fpml = r(pm) - 1.96*r(se_pm)  if ind==2
		replace fpme = r(pm) 				          if ind==2 
		replace fpmh = r(pm) + 1.96*r(se_pm)  if ind==2 

		replace fniel = exp(r(nie) - 1.96*r(se_nie))  if ind==2 
		replace fniee = exp(r(nie)) 				 	        if ind==2
		replace fnieh = exp(r(nie) + 1.96*r(se_nie))  if ind==2
		
		replace fndel = exp(r(nde) - 1.96*r(se_nde))  if ind==2 
		replace fndee = exp(r(nde)) 				 	        if ind==2
		replace fndeh = exp(r(nde) + 1.96*r(se_nde))  if ind==2
		
		replace ftcel = exp(r(tce) - 1.96*r(se_tce))  if ind==2 
		replace ftcee = exp(r(tce)) 				 	        if ind==2
		replace ftceh = exp(r(tce) + 1.96*r(se_tce))  if ind==2
		
	* Death within 5 years given 3 year survival
	gformula dead5y3y broute bcmcat sex agediag dep subtype if subtype==2, ///
		mediation out(dead5y3y) mediator(broute) ex(bcmcat) ///
		samples(1000) minsim moreMC base_confs(sex agediag dep) ///
		commands(dead5y3y: logit, broute: logit) ///
		eq(dead5y3y: broute bcmcat agediag sex i.dep, ///
			broute: bcmcat agediag sex i.dep)	///
		impute(broute) imp_cmd(logit) ///
		imp_eq(broute: dead5y3y bcmcat agediag sex dep) imp_cycles(10) ///
		seed(123) logOR obe control(broute:0) 

		replace fpml = r(pm) - 1.96*r(se_pm)  if ind==3
		replace fpme = r(pm) 				          if ind==3
		replace fpmh = r(pm) + 1.96*r(se_pm)  if ind==3
		
		replace fniel = exp(r(nie) - 1.96*r(se_nie))  if ind==3
		replace fniee = exp(r(nie)) 				          if ind==3
		replace fnieh = exp(r(nie) + 1.96*r(se_nie))  if ind==3
		
		replace fndel = exp(r(nde) - 1.96*r(se_nde))  if ind==3
		replace fndee = exp(r(nde)) 				          if ind==3
		replace fndeh = exp(r(nde) + 1.96*r(se_nde))  if ind==3
		
		replace ftcel = exp(r(tce) - 1.96*r(se_tce))  if ind==3
		replace ftcee = exp(r(tce)) 				          if ind==3
		replace ftceh = exp(r(tce) + 1.96*r(se_tce))  if ind==3
		

		
	/* Summarise the estimates */
	list dtcee dtcel dtceh dndee dndel dndeh dniee dniel dnieh dpme dpml dpmh in 1/3
	list ftcee ftcel ftceh fndee fndel fndeh fniee fniel fnieh fpme fpml fpmh in 1/3

	
	/* Plots for the natural effects */
	twoway rcap ftcel ftceh indte , lc(red)   || scatter ftcee indte, mc(red)   ///
		|| rcap fndel fndeh ind ,   lc(blue)  || scatter fndee ind,   mc(blue)  ///
		|| rcap fniel fnieh indie , lc(green) || scatter fniee indie, mc(green) ///
			xlabel( 1 "12" 2 "36" 3 "60")  ///
			xtitle("Conditional months of follow up") ///
			title("FL", size(medium)) ///
			legend(off) ///
			ytitle("", size(vsmall)) ylabel(1(0.2)2.05) yline(0) ///
			scheme(s2mono) graphregion(fcolor(white)) ///
			name(FLallv3, replace) saving(FLallv3, replace)
	
	graph export FLallv3.pdf, replace
	
	
	/* Combine graphs of DLBCL and FL */
	graph combine DLBCLallv3.gph FLallv3.gph, ///
		ycom com scheme(s2mono) graphregion(fcolor(white)) ///
		name(DLBCLFLallv3, replace) saving(DLBCLFLallv3, replace)
		
	graph export DLBCLFLallv3.pdf, replace
	
	
/* Extract estimate of Proportion mediated for both lymphomas */
	
	twoway rcap dpml dpmh indd || scatter dpme indd || ///
			rcap fpml fpmh indf || scatter fpme indf, ///
			xlabel(1 "12" 2 "36" 3 "60")  ///
			xtitle("Conditional months of follow up") ///
			ytitle("Proportion mediated") ylabel(0(0.1)0.4) ///
			yline(0, lpattern(dash)) ///
			legend(order(2 "DLBCL" 4 "FL") ring(0) pos(2) ///
					cols(2)	title("Lymphoma", size(m)) ///
					region(lstyle(none)) nobox) ///
			scheme(s2mono) graphregion(fcolor(white)) ///
			name(DLBCLFLpmv3, replace) saving(DLBCLFLpmv3, replace)
			
	graph export DLBCLFLpmv3.pdf, replace
	


