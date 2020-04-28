#include<iostream>
#include<cmath>
#include<fstream>


using namespace std;

double alpha=0,a=0,b=0,eps=0.0001,x1=9,x2=10,new_x1,new_x2,p_x1,p_x2;


double func(double x,double y){
	return (-x*x-13*x-y*y-4*y-37);
}

double f(double alp,double p1,double p2){
	return (-(x1+alp*p1)*(x1+alp*p1)-13*(x1+alp*p1)-(x2+alp*p2)*(x2+alp*p2)-4*(x2+alp*p2)-37);
}

void swann(double alphaa,double p1,double p2){
double h=5,delta,al=0,aln=0;
int k=0;

if (f(alphaa-h,p1,p2) >= f(alphaa,p1,p2) && f(alphaa+h,p1,p2) >= f(alphaa,p1,p2))
            {
                a = alphaa-h; b = alphaa+h;
            }
            else if (f(alphaa,p1,p2) >= f(alphaa-h,p1,p2) && f(alphaa,p1,p2) >= f(alphaa+h,p1,p2))
            {
                alphaa=2*(alphaa+1);    
            }
        

	
if((f(alphaa-h,p1,p2)<=f(alphaa,p1,p2))&&(f(alphaa+h,p1,p2)>=f(alphaa,p1,p2))){
alphaa=alphaa+h;
al=alphaa+pow(2,k)*h;	
}
else{
if((f(alphaa-h,p1,p2)>=f(alphaa,p1,p2))&&(f(alphaa+h,p1,p2)<=f(alphaa,p1,p2))){
alphaa=alphaa-h;
h=-h;
al=alphaa+pow(2,k)*h;	
}
else{
	
	a=alphaa-h; b=alphaa+h;
}
}


while(f(al,p1,p2)>=f(alphaa,p1,p2)){
k++;
alphaa=al;
al=alphaa+pow(2,k)*h;
}

if(h>0){
	a=alphaa;
	b=al;
}else{
	a=al;
	b=alphaa;
}


}


double kvad_interpol(double alphaa,double p1,double p2){
double h=0.5;
double alpha1=(a+b)/2;
double alpha2=alpha1+h;
double alpha3=0,q=0;


if(f(alpha1,p1,p2)>f(alpha2,p1,p2)){
alpha3=alpha1+2*h;
}else{
alpha3=alpha1-h;
}

while(abs(f(alpha1,p1,p2)-f(q,p1,p2))>=eps){

double A=(alpha2-alpha3)*(alpha3-alpha1)*(alpha1-alpha2);

double a=((alpha3-alpha2)*f(alpha1,p1,p2)+(alpha1-alpha3)*f(alpha2,p1,p2)+(alpha2-alpha1)*f(alpha3,p1,p2))/A;
double b=((alpha2*alpha2-alpha3*alpha3)*f(alpha1,p1,p2)+(alpha3*alpha3-alpha1*alpha1)*f(alpha2,p1,p2)+(alpha1*alpha1-alpha2*alpha2)*f(alpha3,p1,p2))/A;
q=-b/(2*a);

if ((f(alpha1,p1,p2)<f(alpha2,p1,p2)) && ((alpha3<q)&&(q<alpha1))){
	alpha2=alpha1;
	alpha1=q;
}
else{

	if((f(alpha1,p1,p2)<f(alpha2,p1,p2))&&((alpha1<q)&&(q<alpha2))){
	alpha3=alpha1; 	alpha1=q;
	}
	else{
		if((f(alpha1,p1,p2)>f(alpha2,p1,p2))&&((alpha2<q)&&(q<alpha3))){	
		alpha3=alpha2; alpha2=q;	
		}
		else{
			 alpha1=alpha2; alpha2=q;
		}
	}	
}


}
alphaa=q;
return alphaa;

}
void Spusk(ofstream &fout){
	
	
    p_x1=(-2*x1-13)/ sqrt((-2*x1-13)*(-2*x1-13)+(-2*x2-4)*(-2*x2-4));
    p_x2=(-2*x2-4)/ sqrt((-2*x1-13)*(-2*x1-13)+(-2*x2-4)*(-2*x2-4));
   
    alpha=0;
    
    swann(alpha,p_x1,p_x2);
	alpha=kvad_interpol(alpha,p_x1,p_x2);
   
    new_x1=x1+alpha*p_x1;
    new_x2=x2+alpha*p_x2;
     
     
     if(abs(func(new_x1,new_x2)-func(x1,x2))>=eps){
     	x1=new_x1;
     	x2=new_x2;
     	fout<<x1<<" "<<x2<<" "<<func(x1,x2)<<" "<<endl;
     	Spusk(fout);
	 }else{
	 	x1=new_x1;
     	x2=new_x2;
	 }
}


int main(){
	ofstream fout("data.txt",ios_base::app);

	fout<<x1<<" "<<x2<<" "<<func(x1,x2)<<" "<<endl;
	Spusk(fout);
	cout<<"Spusk: x1= "<<x1<<" "<<"x2= "<<x2<<" "<<"f(x1,x2)= "<<func(x1,x2)<<endl;
	system("python yx.py");
	return 0;;
}
