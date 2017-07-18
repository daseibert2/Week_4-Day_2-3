#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double a(double v)
{
    return v;
}

double fa(double R, double newR)
{
    double x;
    x=((1.0/pow(R,3))-1-(3.0/2)*(newR*newR))/R;
    return x;
}

double dtChange(double dt, double RNew1, double vNew1, double RNew2, double vNew2, double R, double v)
{
    double deltaR, deltaV, e = 1, deltadesired, phi, phi2;

    deltaR = fabs(RNew2-RNew1);
    deltaV = fabs(vNew2-vNew1);

    if(deltaV>deltaR)
        deltaR=deltaV;

    deltadesired = e*(fabs(R)+dt*fabs(v));
    phi = deltadesired/deltaR;

    if(deltaR<deltadesired)
        dt=dt*pow(phi,.2);

    if(deltaR>deltadesired)
        dt=dt*pow(phi,.25);

    return(dt);
}

double RK(double max,double r)
{
    double Vstar,Vstarstar,Vstarstarstar,Rstar,Rstarstar,Rstarstarstar,newR,newV,v2=0.0,t=0.0,t1=0.0;
    double newR1,newV1;
    double R=1+r, DT=0.01, DT1=DT,v=0,oldR=R, tol=0.01;

    printf("Cavitation: R(%lf): %lf\t\tV(%lf): %lf\t\tActual V: %lf\t\tdt: %lf\n\n",t,R,t,v,v2,DT);

    for(t=DT;t<=max;t+=DT1)
    {
        Rstar=R+DT/2*a(v);
        Vstar=v+DT/2*fa(R,v);
        Rstarstar=R+DT/2*a(Vstar);
        Vstarstar=v+DT/2*fa(Rstar,Vstar);
        Rstarstarstar=R+DT*a(Vstarstar);
        Vstarstarstar=v+DT*fa(Rstarstar,Vstarstar);
        newR=R+DT/6*(a(v)+2*a(Vstar)+2*a(Vstarstarstar)+a(Vstarstarstar));
        newV=v+DT/6*(fa(R,v)+2*fa(Rstar,Vstar)+2*fa(Rstarstar,Vstarstar)+fa(Rstarstarstar,Vstarstarstar));

        double R1=R;
        double v1=v;

        for(t1=DT1;t1<(DT1*2);t1+=(DT1/2))
        {
            R1=newR;
            v1=newV;
            Rstar=R1+DT1/2*a(v1);
            Vstar=v1+DT1/2*fa(R1,v1);
            Rstarstar=R1+DT1/2*a(Vstar);
            Vstarstar=v1+DT1/2*fa(Rstar,Vstar);
            Rstarstarstar=R1+DT1*a(Vstarstar);
            Vstarstarstar=v1+DT1*fa(Rstarstar,Vstarstar);
            newR1=R1+DT1/6*(a(v1)+2*a(Vstar)+2*a(Vstarstarstar)+a(Vstarstarstar));
            newV1=v1+DT1/6*(fa(R1,v1)+2*fa(Rstar,Vstar)+2*fa(Rstarstar,Vstarstar)+fa(Rstarstarstar,Vstarstarstar));
        }

        R=newR;
        v=newV;

        /*if(fabs(R-oldR)>5*tol)
            DT/=2;
        else if(fabs(R-oldR)<tol)
            DT*=2;*/

        v2=sqrt((2/pow(R,3))*log(R/(1+r))-(2.0/3)*(1-pow((1+r)/R,3)));
        printf("Cavitation: R(%lf): %lf\t\tV(%lf): %lf\t\tActual V: %lf\t\tdt: %lf\n\n",t,R,t,v,v2,DT1);
        DT1=dtChange(DT,newR,newV,newR1,newV1,R,v);

        //oldR=R;
    }

}



double b(double x,double y,double O)
{
    double w;
    w=O*(y-x);
    return w;
}

double fb(double x,double y,double z,double P)
{
    double w;
    w=x*(P-z)-y;
    return w;
}

double ffb(double x,double y,double z,double B)
{
    double w;
    w=x*y-(B*z);
    return w;
}

double dtChange1(double dt, double RNew1, double vNew1, double zNew1, double RNew2, double vNew2, double zNew2, double R, double v, double z, double O)
{
    double deltaR, deltaV, deltaZ,e = 0.001, deltadesired, phi, phi2;

    deltaR = fabs(RNew2-RNew1);
    deltaV = fabs(vNew2-vNew1);
    deltaZ=fabs(zNew2=zNew1);

    if(deltaV>deltaR){
        if(deltaV>deltaZ)
            deltaR=deltaV;
    }

    if(deltaZ>deltaR){
        if(deltaZ>deltaV)
            deltaR=deltaZ;
    }

    deltadesired = e*(fabs(R)+dt*fabs(b(R,v,O)));
    phi = deltadesired/deltaR;

    if(deltaR<deltadesired)
        dt=dt*pow(phi,.2);

    if(deltaR>deltadesired)
        dt=dt*pow(phi,.25);

    return(dt);
}

double Lorenz(double max, FILE *fptr)
{
    double Ystar,Ystarstar,Ystarstarstar,Xstar,Xstarstar,Xstarstarstar,Zstar,Zstarstar,Zstarstarstar,newZ,newY,newX,x=1.0,y=1.0,z=1.0,t=0.0;
    double O=10, DT=0.01, DT1=DT, t1=0.0,newX1,newY1,newZ1,B=8.0/3,P=28,oldx=x, tol=0.001;

    printf("Lorenz: X(%lf): %lf\t\tY(%lf): %lf\t\tZ(%lf): %lf\t\tdt:%lf\n\n",t,x,t,y,t,z,DT);

    for(t=DT;t<=max;t+=DT1)
    {
        Xstar=x+DT/2*b(x,y,O);
        Ystar=y+DT/2*fb(x,y,z,P);
        Zstar=z+DT/2*ffb(x,y,z,B);
        Xstarstar=x+DT/2*b(Xstar,Ystar,O);
        Ystarstar=y+DT/2*fb(Xstar,Ystar,Zstar,P);
        Zstarstar=z+DT/2*ffb(Xstar,Ystar,Zstar,B);
        Xstarstarstar=x+DT*b(Xstarstar,Ystarstar,O);
        Ystarstarstar=y+DT*fb(Xstarstar,Ystarstar,Zstarstar,P);
        Zstarstarstar=z+DT*ffb(Xstarstar,Ystarstar,Zstarstar,B);
        newX=x+DT/6*(b(x,y,O)+2*b(Xstar,Ystar,O)+2*b(Xstarstar,Ystarstar,O)+b(Xstarstarstar,Ystarstarstar,O));
        newY=y+DT/6*(fb(x,y,z,P)+2*fb(Xstar,Ystar,Zstar,P)+2*fb(Xstarstar,Ystarstar,Zstarstar,P)+fb(Xstarstarstar,Ystarstarstar,Zstarstarstar,P));
        newZ=z+DT/6*(ffb(x,y,z,B)+2*ffb(Xstar,Ystar,Zstar,B)+2*ffb(Xstarstar,Ystarstar,Zstarstar,B)+ffb(Xstarstarstar,Ystarstarstar,Zstarstarstar,B));

        double x1=x;
        double y1=y;
        double z1=z;

        for(t1=DT1;t1<(DT1*2);t1+=(DT1/2))
        {
            x1=newX;
            y1=newY;
            z1=newZ;

            Xstar=x1+DT/2*b(x1,y1,O);
            Ystar=y1+DT/2*fb(x1,y1,z1,P);
            Zstar=z1+DT/2*ffb(x1,y1,z1,B);
            Xstarstar=x1+DT/2*b(Xstar,Ystar,O);
            Ystarstar=y1+DT/2*fb(Xstar,Ystar,Zstar,P);
            Zstarstar=z1+DT/2*ffb(Xstar,Ystar,Zstar,B);
            Xstarstarstar=x1+DT*b(Xstarstar,Ystarstar,O);
            Ystarstarstar=y1+DT*fb(Xstarstar,Ystarstar,Zstarstar,P);
            Zstarstarstar=z1+DT*ffb(Xstarstar,Ystarstar,Zstarstar,B);
            newX1=x1+DT/6*(b(x1,y1,O)+2*b(Xstar,Ystar,O)+2*b(Xstarstar,Ystarstar,O)+b(Xstarstarstar,Ystarstarstar,O));
            newY1=y1+DT/6*(fb(x1,y1,z1,P)+2*fb(Xstar,Ystar,Zstar,P)+2*fb(Xstarstar,Ystarstar,Zstarstar,P)+fb(Xstarstarstar,Ystarstarstar,Zstarstarstar,P));
            newZ1=z1+DT/6*(ffb(x1,y1,z1,B)+2*ffb(Xstar,Ystar,Zstar,B)+2*ffb(Xstarstar,Ystarstar,Zstarstar,B)+ffb(Xstarstarstar,Ystarstarstar,Zstarstarstar,B));
        }

        x=newX;
        y=newY;
        z=newZ;

        fprintf(fptr,"%lf, %lf, %lf, %lf, %lf\n",t,x,y,z,DT1);
        printf("Lorenz: X(%lf): %lf\t\tY(%lf): %lf\t\tZ(%lf): %lf\t\tdt:%lf\n\n",t,x,t,y,t,z,DT1);
        DT1=dtChange1(DT,newX,newY,newZ,newX1,newY1,newZ1,x,y,z,O);

    }
}

int main()
{
    FILE *fptr;
    remove("Lorenz_Data.csv");
    fptr=fopen("Lorenz_Data.csv","w");
    fprintf(fptr,"t, X, Y, Z, DT\n");
    printf("Problem 7.1.1\n\n");
    RK(10,.5);
    printf("\n\nProblem 7.1.2\n\n");
    Lorenz(70,fptr);

    fclose(fptr);

    return 0;
}
