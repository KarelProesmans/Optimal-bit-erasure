#include "MersenneTwister.h"
#include <math.h>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string.h>
#include <sstream>

MTRand mt;
double ran(void);//Gives a random number between 0 and 1.
double dV(double x);//V'(x) (i.e. minus force)
void eustep(double* g,double* v, double x,double dx, double c);//Calculates Gamma(x+dx) from the Euler-Lagrange equation. g=Gamma(x),v=Gamma'(x),c=2kT/(D*tau).
double* boltzdis(int n, double xmin, double xmax);//calculates initial Boltzman distribution. n is number of lattice points, xmin and xmax are the boundaries of the domain.
double* pt(double tfrac, double* path, double xmin, double xmax, int n);//calculates p(x,t). tfrac is t/tau, path is Gamma(x).
double*vleq(double tfrac, double xmin, double xmax, int n, double eps, double c);// Calculates V(x,t) for the protocol that ends in local equilibrium. eps is the erasure error.
double*vcalc(double tfrac, double*path, double xmin, double xmax, int n,double c);//calculates V(x,t) for the optimal protocol
double lowbound(double xmin, double xmax, int n, double eps);//calculates the lower bound for the amount of dissipation.
double upbound(double xmin, double xmax, int n, double eps);//calculates the upper bound for the amount of dissipation
double*gamleq(double tfrac, double xmin, double xmax, int n, double eps, double c);//calculates



int main()
{
	double var1,var2,var3;
	double*v=&var2;
	double*g=&var1;
	double*g2=&var3,xtar0;
	double xmin=-1.73,xmax=1.73;// boundaries of domain
	*g=xmin;
	*v=0;
	int n=4000,i;//c is 2kT/Dt, n the number of lattice points, dl the space between two lattice points, g=Gamma(x), v=Gamma'(x),errmin is the maximal error in the shooting algorithm.
	double c=4.0,err=1.0,errmin=0.000001,errmin0=errmin,l,dl=(xmax-xmin)/((double)n),xtar,v0=(*v),dv=0.01,test1=1.0,test2,eps=0.25,eps2=1-eps;
	FILE*filv1=fopen("optimal.txt","w");//fist column: t, second column: position, third column: V(x,t) for optimal protocol, fourth column: p(x,t) for optimal protocol
	FILE *fil2=fopen("leq.txt","w");//fist column: t, second column: position, third column: V(x,t) for optimal local-equilibrium protocol, fourth column: p(x,t) for optimal local-equilibrium protocol
	double*bol=boltzdis(n,xmin,xmax);
	double* path=(double*)malloc(n*sizeof(double));//path is gamma(x)
	int jb=-1,j;
	while(eps2>0)
	{
		++jb;
		eps2-=dl*bol[jb];
	}//jb is the lattice point that corresponds to x=f_0^{-1}(1-\epsilon).
//	while(c<10.1)//it is possible to do simulations for several tau in one run, by changing the constraint in the while loop.
	{
//	c=c*2.0;
	err=1.0;//initially set error to some value higher than errmin
	*g=xmin;
	v0=0.05;//initial guess for Gamma'(xmin)
	dv=0.00001;//After each step in the shooting algorithm, increase/decrease v0 by an amount maximally equal to dv.
	xtar=0.0;//xtar is the requested value of path[jb] or gamma(f_0^{-1}(1-\epsilon)).
	while(err>errmin)//in this routine we calculate gamma(x) for x<f_0^{-1}(1-\epsilon).
	{
		test2=test1;
		if((*g)-xtar>0)//If gamma(f_0^{-1}(1-\epsilon))>xtar, we need to make the initial value of gamma'(xmin) smaller.
		{
			v0-=ran()*dv;
			test1=-1.0;
		}else
		{
			v0+=ran()*dv;
			test1=1.0;
		}
		if(test2*test1<0)//if the signs of test2 and test1 differ, than v0 was between the two last estimates, which means that the algorithm is converging, so we can make dv smaller.
		{
			dv*=0.9999;
		}
		*g=xmin;*v=v0;
		for(i=0;i<jb;++i)
		{
			l=i*dl+xmin;
			eustep(g,v,l,dl,c);
			path[i]=(*g);
		}
		err=pow((*g)-xtar,2);//the error is defined as the squared distance between the target value of gamma(f_0^{-1}(1-\epsilon)) and the value of gamma(f_0^{-1}(1-\epsilon)) found if gamma'(xmin) is given by its current value.
	}
	xtar0=xtar;
	*g=xmax;//We now repeat the above procedure for x>f_0^{-1}(1-\epsilon).
	dv=0.0001;
	v0=0.05;
	err=1.0;
	errmin=errmin0;
	xtar=0.0;
	int k=0;
	while(err>errmin)
	{
		test2=test1;
		++k;
		if((*g)-xtar>0)
		{
			v0+=ran()*dv;
			test1=1.0;
		}else
		{
			v0-=ran()*dv;
			test1=-1.0;
		}
		if(test2*test1<0&&k>2)
		{
			dv*=0.9999;
		}
		*g=xmax;*v=v0;
		for(i=jb;i<n;++i)
		{
			l=(n-1+jb-i)*dl+xmin;
			for(j=0;j<10;++j)
			{
				eustep(g,v,l,-dl/10.0,c);
			}
			path[n-1+jb-i]=(*g);
		}
		err=pow((*g)-xtar,2);
		*g2=*g;
	}
	double w=0;//w is the total amount of work needed to complete the process.
	for(i=0;i<n-1;++i)
	{
		w+=dl*bol[i]*(log(dl*bol[i]/sqrt(pow((path[i+1]-path[i]),2))/bol[(int)((path[i]-xmin)/dl)])+c/2.0*pow((path[i]-(i*dl+xmin)),2));
	}
	printf("%f\t%f\t%f\t%f\n",2.0/c,w-(log(2.0)+eps*log(eps)+(1-eps)*log(1-eps)),upbound(xmin,xmax,n,eps)*c/2.0,lowbound(xmin,xmax,n,eps)*c/2.0);
	}
	double*V;
	double*V2;
	for(i=0;i<51;++i){//we calculate V(x,t/tau) and p(x,t/tau) for the optimal and the local-equilibrium protocol, for t/tau=0,0.02,...,0.98,1.
	V=vcalc(0.02*i,path,xmin,xmax,n,c/2.0);//V is the optimal potential.
	V2=vcalc(0.02*i,gamleq(0.01+49.0/50.0*0.02*i,xmin,xmax,n,eps,c/2.0),xmin,xmax,n,c/2.0);//V2 is the optimal potential for local equilibrium protocols.
	double* pres=pt(0.02*i,path,xmin,xmax,n);
	double* pleq=pt(0.02*i,gamleq(0.00001+0.02*i,xmin,xmax,n,eps,c/2.0),xmin,xmax,n);
	printf("%d\n",i);
	for(j=1;j<n;++j)
	{
		if(j%10==0)
		{
			fprintf(filv1,"%f\t%f\t%f\t%f\n",2.0/c*(00.02*i),xmin+j*dl,V[j],pres[j]);
			fprintf(fil2,"%f\t%f\t%f\t%f\n",2.0/c*(00.02*i),xmin+j*dl,V2[j],pleq[j]);
		}
	}
    }
    return 0;
}


double* boltzdis(int n, double xmin, double xmax)
{
	double*res=(double*)malloc(n*sizeof(double));
	int i;
	double nor=0;
	double V=0,dx=(xmax-xmin)/n;
	for(i=0;i<n;++i)
	{
		V+=dV(xmin+i*dx)*dx;
		res[i]=exp(-V);
		nor+=res[i]*dx;
	}
	for(i=0;i<n;++i)
	{
		res[i]/=nor;
	}
	return res;
}

void eustep(double* g,double* v, double t,double dt, double c)
{
	double a=c*pow((*v),2)*((*g)-t)+(*v)*((*v)*dV((*g))-dV(t));
	*v=(*v)+a*dt;
	*g=(*g)+(*v)*dt;
}

double* pt(double tfrac, double* path, double xmin, double xmax, int n)
{
	double* p=(double*)malloc(n*sizeof(double));
	double* pe=(double*)malloc(n*sizeof(double));
	double* f=(double*)malloc(n*sizeof(double));
	double* f0=(double*)malloc(n*sizeof(double));
	double* zpos=(double*)malloc(n*sizeof(double));
	double* ypos=(double*)malloc(n*sizeof(double));
	int i,j;
	for(i=0;i<n;++i)
	{
		p[i]=0;pe[i]=0;f[i]=0;f0[i]=0;zpos[i]=0;ypos[i]=0;
	}
	double dx=(xmax-xmin)/n,x,lam;
	p=boltzdis(n,xmin,xmax);
	f[0]=p[0]*dx;
	for(i=1;i<n;++i)
	{
		f[i]=f[i-1]+dx*p[i];
	}
	for(i=0;i<n;++i)
	{
		x=i*dx+xmin;
		zpos[i]=x+tfrac*(path[i]-x);
	}
	ypos[0]=xmin;
	pe[0]=f[0]*dx;
	for(i=1;i<n;++i)
	{
		ypos[i]=(zpos[i]+zpos[i-1])/2.0;
		pe[i]=(f[i]-f[i-1])/(zpos[i]-zpos[i-1]);
	}
	p[0]=pe[0];j=0;
	for(i=1;i<n;++i)
	{
		x=xmin+i*dx;
		while(x>ypos[j]&&j<n-1)
		{
			++j;
		}
		lam=(x-ypos[j-1])/(ypos[j]-ypos[j-1]);
		p[i]=lam*pe[j]+(1-lam)*pe[j-1];
		if(p[i]<0)
		{
			p[i]=0;
		}
	}
	free(f);free(f0);free(zpos);free(ypos);
	return p;
}

double*gamleq(double tfrac, double xmin, double xmax, int n, double eps, double c)
{
	int i,j;
	double sol,lam;
	double dx=(xmax-xmin)/n;
	double* p=boltzdis(n,xmin,xmax);
	double*p2=(double*)malloc(n*sizeof(double));
	double* f=boltzdis(n,xmin,xmax);
	double*f2=(double*)malloc(n*sizeof(double));
	double*gam=(double*)malloc(n*sizeof(double));
	for(i=0;i<n/2;++i)
	{
		p2[i]=2*(1-eps)*p[i];
	}
	for(i=n/2;i<n;++i)
	{
		p2[i]=2*eps*p[i];
	}
	f[0]=p[0]*dx;
	f2[0]=p2[0]*dx;
	for(i=1;i<n;++i)
	{
		f[i]=f[i-1]+p[i]*dx;
		f2[i]=f2[i-1]+p2[i]*dx;
	}
	j=0;
	for(i=0;i<n;++i)
	{
		while(f[i]>f2[j] &&j<n-1)
		{
			++j;
		}
		lam=(f[i]-f2[j-1])/(f2[j]-f2[j-1]+0.000000001);
		gam[i]=lam*(j)+(1-lam)*(j-1);
	}
	for(i=0;i<n;++i)
	{
		gam[i]=gam[i]*dx+xmin;
	}
	sol=0;
	return gam;
}


double*vcalc(double tfrac, double*path, double xmin, double xmax, int n,double c)
{
	int i,j,k;
	double*v=(double*)malloc(n*sizeof(double));
	double*p=pt(tfrac,path,xmin,xmax,n);
	double dt=0.00001,dx=(xmax-xmin)/n,sol,sol2,s1=0,s2=0;
	double*dpt=pt(tfrac-dt,path,xmin,xmax,n);
	double*f=(double*)malloc(n*sizeof(double));
	double*df=(double*)malloc(n*sizeof(double));
	for(i=0;i<n;++i)
	{
		v[i]=0;
		s1+=p[i]*dx;
		s2+=dpt[i]*dx;
	}
	for(i=0;i<n;++i)
	{
		v[i]=0;
		p[i]/=s1;
		dpt[i]/=s2;
	}
	f[0]=p[0]*dx;
	df[0]=dpt[0]*dx;
	for(i=1;i<n;++i)
	{
		f[i]=f[i-1]+p[i]*dx;
		df[i]=df[i-1]+dpt[i]*dx;
	}
	for(i=0;i<n;++i)
	{
		sol=0;
		for(j=0;j<i;++j)
		{
			sol2=(f[j]-df[j])/dt;
			sol+=sol2/(p[j]+0.001)*dx*c;
		}
		sol-=log(p[i]+0.00000001);
		v[i]=sol;
	}
	free(p);free(dpt);
	return v;
}


double upbound(double xmin, double xmax, int n, double eps)
{
	int i,j;
	double sol,lam;
	double dx=(xmax-xmin)/n;
	double* p=boltzdis(n,xmin,xmax);
	double*p2=(double*)malloc(n*sizeof(double));
	double* f=boltzdis(n,xmin,xmax);
	double*f2=(double*)malloc(n*sizeof(double));
	double*gam=(double*)malloc(n*sizeof(double));
	for(i=0;i<n/2;++i)
	{
		p2[i]=2*(1-eps)*p[i];
	}
	for(i=n/2;i<n;++i)
	{
		p2[i]=2*eps*p[i];
	}
	f[0]=p[0]*dx;
	f2[0]=p2[0]*dx;
	for(i=1;i<n;++i)
	{
		f[i]=f[i-1]+p[i]*dx;
		f2[i]=f2[i-1]+p2[i]*dx;
	}
	j=0;
	for(i=0;i<n;++i)
	{
		while(f[i]>f2[j] &&j<n-1)
		{
			++j;
		}
		lam=(f2[j]-f[i])/(f2[j]-f2[j-1]+0.000000001);
		gam[i]=lam*(j)+(1-lam)*(j-1);
	}
	sol=0;
	for(i=0;i<n;++i)
	{
		sol+=dx*p[i]*pow((gam[i]-(i))*dx,2);
	}
	free(p);free(p2);free(f);free(f2);free(gam);
	return sol;
}


double lowbound(double xmin, double xmax, int n, double eps)
{
	int i;
	double dx=(xmax-xmin)/n,sol=0,xmax2;
	double*p=boltzdis(n,xmin,xmax);
	xmax2=0;
	i=n/2;
	while(xmax2<0.5-eps&&i<n)//
	{
		xmax2+=dx*p[i];
		++i;
	}
	xmax2=i;
	for(i=n/2;i<xmax2;++i)
	{
		sol+=pow(xmin+dx*i,2)*dx*p[i];
	}
	free(p);
	return sol;
}


double ran(void)
{
	double x=mt.rand();
	while(x<=0||x>=1)
	{x=mt.rand();}
	return mt.rand();
}



double dV(double x)//V'(x)
{
	double eb=4.0;
	double xm=1.043;
	return eb*4*(pow(x/xm,3)/(xm)-x/(xm*xm));
}