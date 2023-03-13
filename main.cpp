#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#define S 50.0 //size of diffraction grating slit
#define O 240.0 //size of image on screen
#define LAMBDA 1.0 //wave length
#define AMPLITUDE 10.0 //wave amplitude
#define D 120.0 //distance between slit and screen
#define B 5.0 //distance between slits

#define N 1000 // Number of data points
#define INTEGRATION_N 100
#define No_slits 1 // Number of slits  currently only 1 & 2 work
#define File_Name "Data1.txt" //name of the output file


double Real(double x, double y)
{
    double k = 2*M_PI/LAMBDA;
    double r = sqrt(D*D+(y-x)*(y-x));
    double amplitude = AMPLITUDE * cos(k*r)/sqrt(r);

    return amplitude;
}

double Imaginary(double x, double y)
{
    double k = 2*M_PI/LAMBDA;
    double r = sqrt(D*D+(y-x)*(y-x));
    double amplitude = AMPLITUDE * sin(k*r)/sqrt(r);

    return amplitude;
}

double Integration(double xmin, double xmax, double y, double n, double (*Function)(double, double))
{
    double h = (xmax-xmin)/(2.0*n);
    double sum = 0;
    double x = xmin+h;
    
    for(int i=0; i<n; i++)
    {
        sum += h*(Function(x+h,y)+4.0*Function(x,y)+Function(x-h,y))/3.0;
        x += 2.0*h;
    }
    

    return sum;
}

void save_to_file(double* values,string name)
{
    ofstream file;
    file.open(name);
    for(int i=0;i<2*N;i++)
    {
        file<<values[i]<<"\t";
        if(i%2==1){
            file<<endl;
        }
    }
}

double Calculation(double y,double (*Function)(double, double)) // calculates either real or imaginary intensity values and squares sum of them all
{
    double dy = double(B)/double(No_slits);
    double safe = 1;
    if(No_slits == 2){
        safe = 2;
    }
    double sum = 0;
    
    for(int i = 0;i<No_slits;i++)
    {
        if(i%No_slits == 0)
        {
            sum+=Integration((-S/2.0), (S/2.0), (y) , INTEGRATION_N, Function);
        }
        else{
            sum+=Integration((-S/2.0), (S/2.0), (y+safe*i*dy) , INTEGRATION_N, Function);
        }
        
    }
    

    return pow(sum,2);
}

double max_value(double* values){
    double max = -100000000000000;
    for(int i=0;i<N;i++){
        if(values[i]>max){
            max = values[i];
        }
    }
    return max;
}


int main()
{

    double *data = new double[N];
    double *final_data = new double[2*N];
    double *y = new double[N];

    double hy = O/(N-1);

    for(int i=0; i<N; i++)
    {
        y[i] = i*hy + (-O/2.0);
        if(No_slits>1){
            y[i]-=double(B)/double(No_slits);
        }
    }


    for(int i = 0;i<N;i++)
    {
        data[i]=Calculation(y[i],Real)+Calculation(y[i],Imaginary);
    }

    double normal = max_value(data);

    for(int i = 0;i<2*N;i++){
        if(i%2 == 0){
            if(No_slits>1){
                y[i/2]+=double(B)/double(No_slits);
            }
            final_data[i] = y[i/2];
        }
        else{
            final_data[i] = data[int(floor(i/2))]/normal;

        }

    }

    save_to_file(final_data,File_Name);

    return 0;
}