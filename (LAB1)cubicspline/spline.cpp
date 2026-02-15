#include <stdio.h>
#include <math.h>
double f(double x)
{
return sin(x);
}
void progonka(double * y,double *h,const int N,double *c)
{
int i = 1;
double alfa[N+1], beta[N+1], hamma[N+1], delta[N+1], A[N+1], B[N+1];
alfa[1]=hamma[1]=delta[1]=0.0;
beta[1]=1.0;
for(i=2;i<=N;i++)
{
alfa[i]=h[i-1];
beta[i]=2*(h[i-1]+h[i]);
hamma[i]=h[i];
delta[i]=3*(((y[i]-y[i-1])/h[i])-((y[i-1]-y[i-2])/h[i-1]));
}
hamma[N]=0.0;
A[1]=-hamma[1]/beta[1];
B[1]=delta[1]/beta[1];
for(i=2;i<=N-1;i++)
{
A[i]=-hamma[i]/(alfa[i]*A[i-1]+beta[i]);
B[i]=(delta[i]-alfa[i]*B[i-1])/(alfa[i]*A[i-1]+beta[i]);
}
c[N]=(delta[N]-alfa[N]*B[N-1])/(alfa[N]*A[N-1]+beta[N]);
for(i=N;i>1;i--)
{
c[i-1]=A[i-1]*c[i]+B[i-1];
}
}
int main()
{
FILE *fdata = fopen("input.txt","r");
int N = 0, i = 0, j = 0;
double x0=0.0;
char one_char;
while((one_char = fgetc(fdata)) != EOF)
if (one_char == '\n') ++N;
N=N-1;
double x[N+1], y[N+1], xm[20*N+1], ym[20*N+1], h[N+1], a[N+1], b[N+1],
c[N+1], d[N+1];
rewind(fdata);
for(i=0;i<=N;i++)
{
fscanf(fdata,"%i\t%le\t%le\t%le\n",&j,&x[i],&y[i],&h[i]);
};
x0=x[0];
double hh = h[0];progonka(y,h,N,c);
for(i=1;i<N;i++)
{
a[i]=y[i-1];
b[i]=(y[i]-y[i-1])/h[i]-(h[i]/3)*(c[i+1]+2*c[i]);
d[i]=(c[i+1]-c[i])/(3*h[i]);
};
a[N]=y[N-1];
b[N]=(y[N]-y[N-1])/h[N]-(2.0/3.0)*h[N]*c[N];
d[N]=-c[N]/(3*h[N]);
FILE *foutput = fopen("output.txt","w");
hh /= 20.0;
for (i = 0; i <=20*N; i++)
{
xm[i] = x0+i*hh;
ym[i] = f(xm[i]);
}
double s = 0, eps = 0;
j=1;
for (i = 0; i <= 20*(N-1); i++)
{
s=a[j]+b[j]*(xm[i]-x[j-1])+c[j]*(xm[i]-x[j-1])*(xm[i]-x[j-1])+
d[j]*(xm[i]-x[j-1])*(xm[i]-x[j-1])*(xm[i]-x[j-1]);
eps = fabs(s-ym[i]);
fprintf(foutput,"[%i,%i]\t%le\t%le\t%le\t%le\n",i,j,xm[i],ym[i],s,eps);
if ((i!=0)&&((i)%20)==0) {j++;}
}
fclose(foutput);
return 0;
}