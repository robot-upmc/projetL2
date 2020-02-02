#include <math.h>
#include <stdio.h>
#define EPS 0.0001

double f0(double x){
	return -2*exp(-(x+1)*(x+1));
}

double f1(double x){
	return ((x*x*x*x)/4)-2*(x*x)+x;
}
double f2(double x){
	return ((x*x)/2)-3*x*(log(x)-1);
}

double f0derivative(double x){
	return (4*exp(-(x+1)*(x+1)))*(x+1);
}

double f1derivative(double x){
	return (x*x*x)-4*x+1;
}

double f2derivative(double x){
	return x-3*log(x);
}

double f0secondderivative(double x){
	return 4*(-2*exp(-(x+1)*(x+1))*(x*x)-4*exp(-(x+1)*(x+1))*x-exp(-(x+1)*(x+1)));
}

double f1secondderivative(double x){
	return (3*x*x)-4;
}

double f2secondderivative(double x){
	return 1-(3/x);
}

double grabs(double x){
	if(x>0){
		return x;
	}
	else{
		return -x;
	}
}
double methodenewtonf0(double initial){
	double xk=initial;
	double xkp1=xk-(f0derivative(xk)/f0secondderivative(xk));
	int cpt=0;
	while(grabs(xkp1-xk)>EPS){
		cpt++;
		xk=xkp1;
		xkp1=xk-(f0derivative(xk)/f0secondderivative(xk));
	}
	printf("%d\n",cpt);
	return xkp1;
}

double methodenewtonf1(double initial){
	double xk=initial;
	double xkp1=xk-(f1derivative(xk)/f1secondderivative(xk));
	int cpt=0;
	while(grabs(xkp1-xk)>EPS){
		cpt++;
		xk=xkp1;
		xkp1=xk-(f1derivative(xk)/f1secondderivative(xk));
	}
	printf("%d\n",cpt);
	return xkp1;
}
double methodenewtonf2(double initial){
	double xk=initial;
	double xkp1=xk-(f2derivative(xk)/f2secondderivative(xk));
	int cpt=0;
	while(grabs(xkp1-xk)>EPS){
		cpt++;
		xk=xkp1;
		xkp1=xk-(f2derivative(xk)/f2secondderivative(xk));
	}
	printf("%d\n",cpt);
	return xkp1;
}
int main(){
	printf("x : %f\n",methodenewtonf1(-6));
	printf("x : %f\n",methodenewtonf1(-0.1));
	printf("x : %f\n",methodenewtonf1(7));
	return 0;
}
