#include <mpfr_cdf.h>

#include <count_temp_start.h>

void mpf_gamma_pdf(mpf_ptr pdf, mpf_ptr x, mpf_ptr a, mpf_ptr b) {
	size_t ntemp = NTEMP(mpf_gamma_pdf);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr pdfr = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr xr = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr ar = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr br = GET_NEXT_TEMP_PTR(temp);
	MPFR_FUN(set_f, xr, x);
	MPFR_FUN(set_f, ar, a);
	MPFR_FUN(set_f, br, b);
	mpfr_gamma_pdf(pdfr, xr, ar, br);
	MPFR_FUN(get_f,pdf, pdfr);
	mpfr_array_clear_free(temp,ntemp);
	return;
}
#define FUN_NAME mpf_gamma_pdf
#include <count_temp.h>
#define FUN_NAME1 mpf_gamma_pdf

void mpfr_gamma_pdf(mpfr_ptr pdf, mpfr_ptr x, mpfr_ptr a, mpfr_ptr b) {
	if (mpfr_sgn(x)<0 || mpfr_inf_p(x)) {
		MPFR_FUN(set_ui, pdf, 0);
		return;
	}
	size_t ntemp = NTEMP(mpfr_gamma_pdf);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr val = GET_NEXT_TEMP_PTR(temp);
	MPFR_FUN(sub_ui, tmp, a, 1); // a-1
	MPFR_FUN(log, val, x); // ln(x)
	MPFR_FUN(mul, val, val, tmp); // (a-1)ln(x)
	MPFR_FUN(div, tmp, x, b); // x/b
	MPFR_FUN(sub, val, val, tmp); // (a-1)ln(x)-x/b
	MPFR_FUN(log, tmp, b); // ln(b)
	MPFR_FUN(mul, tmp, tmp, a); //a*ln(b)
	MPFR_FUN(sub, val, val, tmp); // -a*ln(b)+(a-1)ln(x)-x/b
	MPFR_FUN(lngamma, tmp, a); // ln(Gamma(a))
	MPFR_FUN(sub, val, val, tmp); // -ln(Gamma(a))-a*ln(b)+(a-1)ln(x)-x/b
	MPFR_FUN(exp, val, val); // 1/Gamma(a)/b^a*x^(a-1)*exp(-x/b)
	MPFR_FUN(set, pdf, val);
	mpfr_array_clear_free(temp, ntemp);
	return;
}
#define FUN_NAME mpfr_gamma_pdf
#include <count_temp.h>
#define FUN_NAME1 mpfr_gamma_pdf

void mpf_gamma_cdf_Q(mpf_t Q, mpf_t x, mpf_t a, mpf_t b) {
	size_t ntemp = NTEMP(mpf_gamma_cdf_Q);
	mpfr_t *temp = mpfr_array_alloc(ntemp);
	mpfr_array_init(temp,ntemp);
	mpfr_t *Qr = GET_NEXT_TEMP(temp);
	mpfr_t *xr = GET_NEXT_TEMP(temp);
	mpfr_t *ar = GET_NEXT_TEMP(temp);
	mpfr_t *br = GET_NEXT_TEMP(temp);
	MPFR_FUN(set_f, *xr, x);
	MPFR_FUN(set_f, *ar, a);
	MPFR_FUN(set_f, *br, b);
	mpfr_gamma_cdf_Q(*Qr, *xr, *ar, *br);
	MPFR_FUN(get_f,Q,*Qr);
	mpfr_array_clear(temp,ntemp);
	free(temp);
	return;
}
#define FUN_NAME mpf_gamma_cdf_Q
#include <count_temp.h>
#define FUN_NAME1 mpf_gamma_cdf_Q

void mpfr_gamma_cdf_Q(mpfr_t Q, mpfr_t x, mpfr_t a, mpfr_t b) {
	#undef NPASS
	#define NPASS MAX(NTEMP(mpfr_gamma_cdf_Q_routine),NTEMP(mpfr_array_get_epsilon_routine))
	size_t ntemp = NTEMP(mpfr_gamma_cdf_Q);
	mpfr_t *temp = mpfr_array_alloc(ntemp);
	mpfr_array_init(temp,ntemp);
	mpfr_t *y = GET_NEXT_TEMP(temp);
	MPFR_FUN(div,*y,x,b);
	mpfr_t *eps = GET_NEXT_TEMP(temp);
	mpfr_t *pass_temp = GET_NEXT_TEMP(temp);
	mpfr_array_get_epsilon_routine(*eps, pass_temp);
	mpfr_gamma_cdf_Q_routine(Q, *y, a, *eps, pass_temp);
	mpfr_array_clear(temp,ntemp);
	free(temp);
	return;
}
#define FUN_NAME mpfr_gamma_cdf_Q
#include <count_temp.h>
#define FUN_NAME1 mpfr_gamma_cdf_Q

void mpfr_pdf_ugaussian(mpfr_ptr pdf, mpfr_ptr x) {
	if (mpfr_inf_p(x)) {
		MPFR_FUN(set_ui, pdf, 0);
		return;
	}
	size_t ntemp = NTEMP(mpfr_pdf_ugaussian);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr tmp1 = GET_NEXT_TEMP_PTR(temp);
	MPFR_FUN(pow_ui, tmp, x, 2);
	MPFR_FUN(div_ui, tmp, tmp, 2);
	MPFR_FUN(neg, tmp, tmp);
	MPFR_FUN(exp, tmp, tmp);
	MPFR_FUN(const_pi, tmp1);
	MPFR_FUN(mul_ui, tmp1, tmp1, 2);
	MPFR_FUN(sqrt, tmp1, tmp1);
	MPFR_FUN(div, tmp, tmp, tmp1);
	MPFR_FUN(set, pdf, tmp);
	mpfr_array_clear_free(temp, ntemp);
	return;
}
#define FUN_NAME mpfr_pdf_ugaussian
#include <count_temp.h>
#define FUN_NAME1 mpfr_pdf_ugaussian

void mpfr_pdf_gaussian(mpfr_ptr pdf, mpfr_ptr x, mpfr_ptr mean, mpfr_ptr stdev) {
	size_t ntemp = NTEMP(mpfr_pdf_gaussian);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
	MPFR_FUN(sub, tmp, x, mean);
	MPFR_FUN(div, tmp, tmp, stdev);
	mpfr_pdf_ugaussian(tmp, tmp);
	MPFR_FUN(div, tmp, tmp, stdev);
	MPFR_FUN(set, pdf, tmp);
	mpfr_array_clear_free(temp, ntemp);
	return;
}
#define FUN_NAME mpfr_pdf_gaussian
#include <count_temp.h>
#define FUN_NAME1 mpfr_pdf_gaussian

void mpfr_cdf_ugaussian_P(mpfr_ptr P, mpfr_ptr x) {
	if (mpfr_inf_p(x)) {
		if (mpfr_sgn(x)>0) {
			MPFR_FUN(set_ui, P, 1);
		} else {
			MPFR_FUN(set_ui, P, 0);
		}
	}
	size_t ntemp = NTEMP(mpfr_cdf_ugaussian_P);
	mpfr_ptr *temp = mpfr_ptr_array_alloc_init(ntemp);
	mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
	MPFR_FUN(set_ui, tmp, 2);
	MPFR_FUN(sqrt, tmp, tmp);
	MPFR_FUN(div, tmp, x, tmp); // x/sqrt(2)
	#if 0
	/* cdf_ugaussian_P(x) = 1/2 + 1/2*erf(x/sqrt(2)) */
	MPFR_FUN(erf, tmp, tmp);
	MPFR_FUN(add_ui, tmp, tmp, 1);
	MPFR_FUN(div_ui, tmp, tmp, 2);
	#endif /* 0 */
	/* cdf_ugaussian_P(x) = 1/2*erf(-x/sqrt(2)) */
	MPFR_FUN(neg, tmp, tmp);
	MPFR_FUN(erfc, tmp, tmp);
	MPFR_FUN(div_ui, tmp, tmp, 2);
	MPFR_FUN(set, P, tmp);
	mpfr_ptr_array_clear_free(temp, ntemp);
	return;
}
#define FUN_NAME mpfr_cdf_ugaussian_P
#include <count_temp.h>
#define FUN_NAME1 mpfr_cdf_ugaussian_P

void mpfr_cdf_ugaussian_Q(mpfr_ptr Q, mpfr_ptr x) {
	if (mpfr_inf_p(x)) {
		if (mpfr_sgn(x)>0) {
			MPFR_FUN(set_ui, Q, 0);
		} else {
			MPFR_FUN(set_ui, Q, 1);
		}
	}
	size_t ntemp = NTEMP(mpfr_cdf_ugaussian_Q);
	mpfr_ptr *temp = mpfr_ptr_array_alloc_init(ntemp);
	mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
	/* cdf_ugaussian_Q(x) = 1/2*erfc(x/sqrt(2)) */
	MPFR_FUN(set_ui, tmp, 2);
	MPFR_FUN(sqrt, tmp, tmp);
	MPFR_FUN(div, tmp, x, tmp);
	MPFR_FUN(erfc, tmp, tmp);
	MPFR_FUN(div_ui, tmp, tmp, 2);
	MPFR_FUN(set, Q, tmp);
	mpfr_ptr_array_clear_free(temp, ntemp);
	return;
}
#define FUN_NAME mpfr_cdf_ugaussian_Q
#include <count_temp.h>
#define FUN_NAME1 mpfr_cdf_ugaussian_Q

void mpfr_cdf_gaussian_P(mpfr_ptr P, mpfr_ptr x, mpfr_ptr mean, mpfr_ptr stdev) {
	size_t ntemp = NTEMP(mpfr_cdf_gaussian_P);
	mpfr_ptr *temp = mpfr_ptr_array_alloc_init(ntemp);
	mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
	MPFR_FUN(sub, tmp, x, mean);
	MPFR_FUN(div, tmp, tmp, stdev);
	mpfr_cdf_ugaussian_P(P, tmp);
	mpfr_ptr_array_clear_free(temp, ntemp);
	return;
}
#define FUN_NAME mpfr_cdf_gaussian_P
#include <count_temp.h>
#define FUN_NAME1 mpfr_cdf_gaussian_P

void mpfr_gamma_cdf_Q_routine(mpfr_t Q, mpfr_t x, mpfr_t a, mpfr_t eps, mpfr_t *temp) {
	#undef NPASS
	#define NPASS MAX(MAX(NTEMP(mpfr_gamma_cdf_Q_CF),NTEMP(mpfr_gamma_cdf_Q_series)),MAX(NTEMP(mpfr_gamma_cdf_P_series),NTEMP(mpfr_gamma_cdf_Q_asymp_unif)))
	if (mpfr_sgn(a)<=0)
		return;
	if (mpfr_sgn(x)<=0) {
		MPFR_FUN(set_ui,Q,1);
		return;
	}
	mpfr_t *pass_temp = GET_NEXT_TEMP(temp);
	mpfr_t *tmp = pass_temp; // the same memory
	MPFR_FUN(div_ui,*tmp, a, 2);
	if (mpfr_cmp(x,*tmp)<=0) {
		mpfr_gamma_cdf_P_series(Q, x, a, eps, pass_temp);
		MPFR_FUN(ui_sub, Q, 1, Q);
		return;
	}
	MPFR_FUN(sub, *tmp, x, a);
	MPFR_FUN(pow_ui, *tmp, *tmp, 2);
	if (mpfr_cmp_d(a,1.0e+06)>=0 && mpfr_cmp(*tmp, a)<0) {
		mpfr_gamma_cdf_Q_asymp_unif(Q, x, a, eps, pass_temp);
		return;
	}
	if (mpfr_cmp_d(a,0.2)<0 && mpfr_cmp_d(x,5.0)<0) {
		mpfr_gamma_cdf_Q_series(Q, x, a, eps, pass_temp);
		return;
	}
	MPFR_FUN(sqrt, *tmp, a);
	MPFR_FUN(sub, *tmp, a, *tmp);
	if (mpfr_cmp(x, *tmp)>0) {
		mpfr_gamma_cdf_Q_CF(Q, x, a, eps, pass_temp);
		return;
	}
	mpfr_gamma_cdf_P_series(Q, x, a, eps, pass_temp);
	MPFR_FUN(ui_sub, Q, 1, Q);
	return;
}

#define FUN_NAME mpfr_gamma_cdf_Q_routine
#include <count_temp.h>
#define FUN_NAME1 mpfr_gamma_cdf_Q_routine

void mpfr_gamma_cdf_Q_CF(mpfr_t Q_CF, mpfr_t x, mpfr_t a, mpfr_t eps, mpfr_t *temp) {
	#undef NPASS
	#define NPASS MAX(NTEMP(mpfr_gamma_cdf_D),NTEMP(mpfr_gamma_cdf_F_CF))
	mpfr_t* D = GET_NEXT_TEMP(temp);
	mpfr_t* pass_temp = GET_NEXT_TEMP(temp);
	mpfr_gamma_cdf_D(*D, x, a, pass_temp);
	mpfr_gamma_cdf_F_CF(Q_CF, x, a, eps, pass_temp);
	MPFR_FUN(mul,Q_CF,Q_CF,*D);
	MPFR_FUN(mul,Q_CF,Q_CF,a);
	MPFR_FUN(div,Q_CF,Q_CF,x);
	return;
}
#define FUN_NAME mpfr_gamma_cdf_Q_CF
#include <count_temp.h>
#define FUN_NAME1 mpfr_gamma_cdf_Q_CF

void mpfr_gamma_cdf_F_CF(mpfr_t F_CF, mpfr_t x, mpfr_t a, mpfr_t eps, mpfr_t *temp) {
	const size_t nmax = MPFR_CDF_MAX_ITER;
	mpfr_t *small = GET_NEXT_TEMP(temp);
	MPFR_FUN(pow_ui, *small, eps,3);
	mpfr_t *Cn = GET_NEXT_TEMP(temp);
	mpfr_t *Dn = GET_NEXT_TEMP(temp);
	mpfr_set_inf(*Cn,1);
	MPFR_FUN(set_d,*Dn,1.0);
	MPFR_FUN(set_d,F_CF,1.0);
	size_t n;
	for (n=2; n<nmax; n++) {
		mpfr_t *an = GET_NEXT_TEMP(temp);
		if (n%2) /* is odd */ { 
			MPFR_FUN(set_ui,*an, n-1);
			MPFR_FUN(div_ui,*an, *an, 2);
			MPFR_FUN(div,*an, *an, x);
		} else {
			MPFR_FUN(set_ui,*an, n);
			MPFR_FUN(div_ui,*an, *an, 2);
			MPFR_FUN(sub,*an, *an, a);
			MPFR_FUN(div,*an, *an, x);
		}
		 /* Dn = 1.0 + an * Dn; */
		MPFR_FUN(mul,*Dn,*Dn,*an);
		MPFR_FUN(add_ui,*Dn,*Dn,1);
		if (mpfr_cmpabs(*Dn,*small)<0) {
			MPFR_FUN(set, *Dn, *small);
		}
		/* Dn = 1.0 / Dn */
		MPFR_FUN(ui_div,*Dn,1,*Dn);
		/* Cn = 1.0 + an/Cn; */
		MPFR_FUN(div,*Cn, *an, *Cn);
		MPFR_FUN(add_ui,*Cn, *Cn, 1);
		if (mpfr_cmpabs(*Cn,*small)<0) {
			MPFR_FUN(set, *Cn, *small);
		}
		MPFR_FUN(mul,*an, *Cn, *Dn);
		MPFR_FUN(mul,F_CF, F_CF, *an);
		/* if (abs(an-1)<eps) break; */
		MPFR_FUN(sub_ui,*an,*an,1);
		if (mpfr_cmpabs(*an, eps)<0) break;
	}
	//printf("%s: iterations: %d\n", __func__, n);
	return;
}
#define FUN_NAME mpfr_gamma_cdf_F_CF
#include <count_temp.h>
#define FUN_NAME1 mpfr_gamma_cdf_F_CF

void mpfr_gamma_cdf_P_series(mpfr_t P, mpfr_t x, mpfr_t a, mpfr_t eps, mpfr_t *temp) {
	#undef NPASS
	#define NPASS MAX(NTEMP(mpfr_gamma_cdf_D),NTEMP(mpfr_exprel_n_CF))
	const size_t nmax = MPFR_CDF_MAX_ITER;
	mpfr_t *D = GET_NEXT_TEMP(temp);
	mpfr_t *tmp = GET_NEXT_TEMP(temp);
	mpfr_t *pass_temp = GET_NEXT_TEMP(temp);
	mpfr_gamma_cdf_D(*D,x,a,pass_temp);
	MPFR_FUN(mul_d,*tmp,a,0.995);
	if (mpfr_cmp(x,*tmp)>0 && mpfr_cmp_d(a,1.0e+5)>0) {
		mpfr_exprel_n_CF(P, x, a, eps, pass_temp);
		MPFR_FUN(mul, P, P, *D);
	}
	
	size_t nlow;
	if (mpfr_cmp(x,a)>0) {
		MPFR_FUN(sub, *tmp, x, a);
		nlow=MPFR_FUN(get_ui, *tmp);
	} else {
		nlow = 0;
	}
	mpfr_t *term = pass_temp; // overrides first pass_temp
	MPFR_FUN(set_ui, P, 1);
	MPFR_FUN(set_ui, *term, 1);
	size_t n;
	for (n=1; n<nlow && n<nmax; n++) {
		MPFR_FUN(add_ui, *tmp, a, n);
		MPFR_FUN(mul, *term, *term, x);
		MPFR_FUN(div, *term, *term, *tmp);
		MPFR_FUN(add, P, P, *term);
	}
	
	for ( ; n<nmax; n++) {
		MPFR_FUN(add_ui, *tmp, a, n);
		MPFR_FUN(mul, *term, *term, x);
		MPFR_FUN(div, *term, *term, *tmp);
		MPFR_FUN(add, P, P, *term);
		
		MPFR_FUN(div, *tmp, *term, P);
		if (mpfr_cmpabs(*tmp, eps)<0) break;
	}
	MPFR_FUN(mul, P, P, *D);
	//printf("%s: iterations: %d\n", __func__, n);
	return;
}
#define FUN_NAME mpfr_gamma_cdf_P_series
#include <count_temp.h>
#define FUN_NAME1 mpfr_gamma_cdf_P_series

void mpfr_gamma_cdf_Q_series(mpfr_t Q, mpfr_t x, mpfr_t a, mpfr_t eps, mpfr_t *temp) {
	const size_t nmax = MPFR_CDF_MAX_ITER;
	mpfr_t *term1 = GET_NEXT_TEMP(temp);
	mpfr_t *tmp = GET_NEXT_TEMP(temp);
	MPFR_FUN(add_ui, *tmp, a, 1); // a+1
	MPFR_FUN(lngamma, *tmp, *tmp); // ln(Gamma(a+1))
	MPFR_FUN(log, *term1, x); // ln(x)
	MPFR_FUN(mul, *term1, *term1, a); // a*ln(x)
	MPFR_FUN(sub, *term1, *term1, *tmp); // a*ln(x)-ln(Gamma(a+1))
	MPFR_FUN(expm1, *term1, *term1); // x^a/Gamma(a+1) - 1
	MPFR_FUN(neg, *term1, *term1); //1 - x^a/Gamma(a+1)
	/*
	MPFR_FUN(gamma, *tmp, *tmp); // Gamma(a+1)
	MPFR_FUN(pow, *term1, x, a); // x^a
	MPFR_FUN(div, *term1, *term1, *tmp); // x^a/Gamma(a+1)
	MPFR_FUN(ui_sub, *term1, 1, *term1); // 1-x^a/Gamma(a+1)
	*/
	mpfr_t *sum = GET_NEXT_TEMP(temp);
	MPFR_FUN(set_ui,*sum,1);
	MPFR_FUN(set_ui,*tmp,1);
	size_t n;
	mpfr_t *tmp2 = GET_NEXT_TEMP(temp);
	mpfr_t *tmp3 = GET_NEXT_TEMP(temp);
	for (n=1; n<nmax; n++) {
		MPFR_FUN(mul,*tmp,*tmp,x);
		MPFR_FUN(div_si,*tmp,*tmp,-(n+1)); // t=-t*x/(n+1)
		MPFR_FUN(add_ui, *tmp2, a, 1); // a+1
		MPFR_FUN(add_ui, *tmp3, a, n+1); // a+n+1
		MPFR_FUN(div, *tmp2, *tmp2, *tmp3); // (a+1)/(a+n+1)
		MPFR_FUN(mul, *tmp2, *tmp2, *tmp); // (a+1)/(a+n+1)*t
		MPFR_FUN(add, *sum, *sum, *tmp2); // sum = sum + (a+1)/(a+n+1)*t
		// if (abs(t/sum)<eps) break;
		MPFR_FUN(div,*tmp2, *tmp,*sum);
		if (mpfr_cmpabs(*tmp2,eps)<0) break;
	}
	MPFR_FUN(ui_sub,Q,1,*term1);
	MPFR_FUN(mul,Q,Q,a);
	MPFR_FUN(add_ui,*tmp,a,1);
	MPFR_FUN(div, Q, Q, *tmp);
	MPFR_FUN(mul, Q, Q, x);
	MPFR_FUN(mul, Q, Q, *sum); // (1 - term1) * a/(a+1) * x * sum
	MPFR_FUN(add, Q, Q, *term1); // term1 + ...
	//printf("%s: iterations: %d\n", __func__, n);
	return;
	/*
	const double pg21_d = -2.404113806319188570799476;  // PolyGamma[2,1]
	mpfr_t *pg21 = GET_NEXT_TEMP(temp);
	MPFR_FUN(set_d,*pg21,pg21_d);
	//const double lnx  = log(x);
	mpfr_t *lnx = GET_NEXT_TEMP(temp);
	MPFR_FUN(log,*lnx, x);
    //const double el   = M_EULER+lnx;
	mpfr_t *el = GET_NEXT_TEMP(temp);
	MPFR_FUN(const_euler,*el);
	mpfr_t *temp1 = GET_NEXT_TEMP(temp);
	#undef NPASS
	#define NPASS 10
	mpfr_t *c = GET_NEXT_TEMP(temp);
	//const double c1 = -el;
	MPFR_FUN(neg,c[0],*el);
	//const double c2 = M_PI*M_PI/12.0 - 0.5*el*el;
	MPFR_FUN(pow_ui,c[1],*el,2);
	MPFR_FUN(div_ui,c[1],c[1],2); // euler^2/2
	MPFR_FUN(const_pi,*temp1);
	MPFR_FUN(pow_ui,*temp1,*temp1,2);
	MPFR_FUN(div_ui,*temp1,*temp1,12); // pi^2/12
	MPFR_FUN(sub,c[1],*temp1,c[1]);
    //const double c3 = el*(M_PI*M_PI/12.0 - el*el/6.0) + pg21/6.0;
	MPFR_FUN(pow_ui,c[2],*el,2);
	MPFR_FUN(div_ui,c[2],c[2],6); // euler^2/6
	MPFR_FUN(sub,c[2],*temp1,c[2]);
	MPFR_FUN(mul,c[2],c[2],*el);
	MPFR_FUN(div_ui, *pg21, *pg21,6);
	MPFR_FUN(add,c[2],c[2],*pg21);
	*/
}
#define FUN_NAME mpfr_gamma_cdf_Q_series
#include <count_temp.h>
#define FUN_NAME1 mpfr_gamma_cdf_Q_series

void mpfr_gamma_cdf_D(mpfr_ptr D, mpfr_ptr x, mpfr_ptr a, mpfr_t *temp) {
	mpfr_ptr temp1 = GET_NEXT_TEMP_PTR(temp);
	MPFR_FUN(add_ui,temp1, a, 1);
	MPFR_FUN(lngamma,temp1, temp1);
	MPFR_FUN(log,D, x);
	MPFR_FUN(mul,D, D, a);
	MPFR_FUN(sub,D, D, x);
	MPFR_FUN(sub,D, D, temp1);
	MPFR_FUN(exp, D, D);
	return;
}
#define FUN_NAME mpfr_gamma_cdf_D
#include <count_temp.h>
#define FUN_NAME1 mpfr_gamma_cdf_D

void mpfr_gamma_cdf_Q_asymp_unif(mpfr_t Q, mpfr_t x, mpfr_t a, mpfr_t eps, mpfr_t *temp) {
	mpfr_set_nan(Q); // do not need it in my application
	return;
}
#define FUN_NAME mpfr_gamma_cdf_Q_asymp_unif
#include <count_temp.h>
#define FUN_NAME1 mpfr_gamma_cdf_Q_asymp_unif

void mpfr_exprel_n_CF(mpfr_t CF, mpfr_t x, mpfr_t N, mpfr_t eps, mpfr_t *temp) {
	size_t nmax = MPFR_CDF_MAX_ITER;
	mpfr_t *recur_big = GET_NEXT_TEMP(temp);
	MPFR_FUN(set_d, *recur_big, GSL_SQRT_DBL_MAX);
	size_t n;
	mpfr_t *eps2 = GET_NEXT_TEMP(temp);
	MPFR_FUN(mul_ui, *eps2, eps, 2);
	mpfr_t *An = GET_NEXT_TEMP(temp);
	mpfr_t *Bn = GET_NEXT_TEMP(temp);
	mpfr_t *Anm1 = GET_NEXT_TEMP(temp);
	mpfr_t *Anm2 = GET_NEXT_TEMP(temp);
	mpfr_t *Bnm1 = GET_NEXT_TEMP(temp);
	mpfr_t *Bnm2 = GET_NEXT_TEMP(temp);
	mpfr_t *an = GET_NEXT_TEMP(temp);
	mpfr_t *bn = GET_NEXT_TEMP(temp);
	MPFR_FUN(add_ui, *An, N, 1); // An = N+1
	MPFR_FUN(sub, *Bn, *An, x); // Bn = N+1-x
	MPFR_FUN(set_ui, *Anm1, 1);
	MPFR_FUN(set_ui, *Anm2, 0);
	MPFR_FUN(set_ui, *Bnm1, 1);
	MPFR_FUN(set_ui, *Bnm1, 1);
	MPFR_FUN(div, CF, *An, *Bn); // fn = An/Bn

	for (n=3; n<nmax; n++) {
		mpfr_t *tmp = GET_NEXT_TEMP(temp);
		MPFR_FUN(set, *Anm2, *Anm1); // Anm2 = Anm1
		MPFR_FUN(set, *Bnm2, *Bnm1); // Bnm2 = Bnm1
		MPFR_FUN(set, *Anm1, *An); // Anm1 = An
		MPFR_FUN(set, *Bnm1, *Bn); // Bnm1 = Bn
		/* an = ( GSL_IS_ODD(n) ? ((n-1)/2)*x : -(N+(n/2)-1)*x ) */
		if (n%2) { /* is odd */
			MPFR_FUN(mul_ui, *an, x, (n-1)/2);
		} else {
			MPFR_FUN(add_ui, *an, N, n/2-1);
			MPFR_FUN(mul, *an, *an, x);
			MPFR_FUN(neg, *an, *an);
		}
		MPFR_FUN(add_ui, *bn, N, n-1); // bn = N + n - 1
		MPFR_FUN(mul, *tmp, *bn, *Anm1); 
		MPFR_FUN(mul, *An, *an, *Anm2);
		MPFR_FUN(add, *An, *An, *tmp); // An = bn*Anm1 + an*Anm2
		MPFR_FUN(mul, *tmp, *bn, *Bnm1);
		MPFR_FUN(mul, *Bn, *an, *Bnm2); 
		MPFR_FUN(add, *Bn, *Bn, *tmp); // Bn = bn*Bnm1 + an*Bnm2
		if (mpfr_cmpabs(*An, *recur_big) || mpfr_cmpabs(*Bn, *recur_big)) {
			MPFR_FUN(div, *An, *An, *recur_big);
			MPFR_FUN(div, *Bn, *Bn, *recur_big);
			MPFR_FUN(div, *Anm1, *Anm1, *recur_big);
			MPFR_FUN(div, *Anm2, *Anm2, *recur_big);
			MPFR_FUN(div, *Bnm1, *Bnm1, *recur_big);
			MPFR_FUN(div, *Bnm2, *Bnm2, *recur_big);
		}
		MPFR_FUN(set, *tmp, CF);
		MPFR_FUN(div, CF, *An, *Bn);
		MPFR_FUN(div, *tmp, *tmp, CF);
		MPFR_FUN(ui_sub, *tmp, 1, *tmp);
		if (mpfr_cmpabs(*tmp, *eps2)<0) break;
	}
	//printf("%s: iterations: %d\n", __func__, n);
	return;
}
#define FUN_NAME mpfr_exprel_n_CF
#include <count_temp.h>
#define FUN_NAME1 mpfr_exprel_n_CF

void mpfr_erf_approx(mpfr_ptr erf, mpfr_ptr x) {
	size_t ntemp = NTEMP(mpfr_erf_approx);
	mpfr_ptr *temp = mpfr_ptr_array_alloc_init(ntemp);
	mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr tmp1 = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr a = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr pi = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr x2 = GET_NEXT_TEMP_PTR(temp);
	MPFR_FUN(const_pi, pi);
	MPFR_FUN(sub_ui, a, pi, 3);
	MPFR_FUN(mul_ui, a, a, 8);
	MPFR_FUN(ui_sub, tmp, 4, pi);
	MPFR_FUN(mul, tmp, tmp, pi);
	MPFR_FUN(mul_ui, tmp, tmp, 3);
	MPFR_FUN(div, a, a, tmp);
	
	MPFR_FUN(ui_div, tmp, 4, pi);
	MPFR_FUN(pow_ui, x2, x, 2);
	MPFR_FUN(mul, tmp1, x2, a);
	MPFR_FUN(add, tmp, tmp, tmp1); // 4/pi+a*x^2
	MPFR_FUN(add_ui, tmp1, tmp1, 1); // 1+a*x^2
	MPFR_FUN(div, tmp, tmp, tmp1);
	MPFR_FUN(mul, tmp, tmp, x2);
	MPFR_FUN(neg, tmp, tmp);
	
	MPFR_FUN(expm1, tmp, tmp);
	MPFR_FUN(neg, tmp, tmp);
	MPFR_FUN(sqrt, tmp, tmp);
	if (mpfr_sgn(x)<0) {
		MPFR_FUN(neg, tmp, tmp);
	}
	MPFR_FUN(set, erf, tmp);
	mpfr_ptr_array_clear_free(temp, ntemp);
	return;
}
#define FUN_NAME mpfr_erf_approx
#include <count_temp.h>
#define FUN_NAME1 mpfr_erf_approx

#include <count_temp_end.h>