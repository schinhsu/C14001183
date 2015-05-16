#include <cstdlib>
#include <iostream>
#include <math.h>
#include <time.h>
#include <string>
#include <cstring>

using namespace std;

void sft(double *xR,double *xI,double *yR,double *yI, int N){
	int i,j;
	double wR,wI;
	for (i=0;i<N;i++){
		yR[i] = 0.0;
		yI[i] = 0.0;
		for (j=0;j<N;j++){
			wR = cos((-1)*i*j*2*M_PI/N);
			wI = sin((-1)*i*j*2*M_PI/N);
			yR[i] += (xR[j]*wR - xI[j]*wI); 
			yI[i] += (xR[j]*wI + xI[j]*wR);
		}
	}
}

void index_transform(double *xR,double *xI,double *yR,double *yI,int N,int p,int q,int r){
     int i,j,k,m,control,tmp;
     control = 0;
     tmp = N;
     while (r>0){
           tmp /= 5;
           if (control%2 == 0){
              for (k=0;k<N;k+=5*tmp){
                  for (i=0;i<5;i++){
                      for (j=0;j<tmp;j++){
                          yR[tmp*i+j+k] = xR[(5*j+i)+k];
                          yI[tmp*i+j+k] = xI[(5*j+i)+k];
                      }
                  }
              }
           }
           else {
                for (k=0;k<N;k+=5*tmp){
                    for (i=0;i<5;i++){
                        for (j=0;j<tmp;j++){
                            xR[tmp*i+j+k] = yR[(5*j+i)+k];
                            xI[tmp*i+j+k] = yI[(5*j+i)+k];
                        }
                    }
                }
           }
           r --;
           control ++;
     }
     while (q>0){
           tmp /= 3;
           if (control%2 == 0){
              for (k=0;k<N;k+=tmp*3){
                  for (i=0;i<3;i++){
                      for (j=0;j<tmp;j++){
                          yR[tmp*i+j+k] = xR[(3*j+i)+k];
                          yI[tmp*i+j+k] = xI[(3*j+i)+k];
                      }
                  }
              }
           }
           else {
                for (k=0;k<N;k+=tmp*3){
                    for (i=0;i<3;i++){
                        for (j=0;j<tmp;j++){
                            xR[tmp*i+j+k] = yR[(3*j+i)+k];
                            xI[tmp*i+j+k] = yI[(3*j+i)+k];
                        }
                    }
                }
           }
           q --;
           control ++;
     }
     while (p>0){
           tmp /= 2;
           if (control%2 == 0){
              for (k=0;k<N;k+=tmp*2){
                  for (i=0;i<2;i++){
                      for (j=0;j<tmp;j++){
                          yR[tmp*i+j+k] = xR[(2*j+i)+k];
                          yI[tmp*i+j+k] = xI[(2*j+i)+k];
                      }
                  }
              }
           }
           else {
                for (k=0;k<N;k+=tmp*2){
                    for (i=0;i<2;i++){
                        for (j=0;j<tmp;j++){
                            xR[tmp*i+j+k] = yR[(2*j+i)+k];
                            xI[tmp*i+j+k] = yI[(2*j+i)+k];
                        }
                    }
                }
           }
           p --;
           control ++;
     }
     if (control %2 == 1){
        for (i=0;i<N;i++){
            yR[i] = xR[i];
            yI[i] = xI[i];
        }
     }
}

void fft(double *yR,double *yI,int N,int p,int q,int r){
     int i,j,k,tmp=1;
     double theta,wR,wI,tmpR,tmpI;
     while (p>0){
           for (i=0;i<tmp;i++){
               theta = -2.0*i*M_PI/(tmp*2);
               wR = cos(theta);
               wI = sin(theta);
               for (j=i;j<N;j+=tmp*2){
                   k = j+tmp;
                   tmpR = wR*yR[k] - wI*yI[k];
                   tmpI = wR*yI[k] + wI*yR[k];
                   yR[k] = yR[j] - tmpR;
                   yI[k] = yI[j] - tmpI;
                   yR[j] += tmpR;
                   yI[j] += tmpI;
               }
           }
           tmp *=2;
           p --;
     }
     double const1,const2,const1R,const1I,const2R,const2I,const3R,const3I,const4R,const4I,tmp1R,tmp1I,tmp2R,tmp2I,tmp3R,tmp3I,tmp4R,tmp4I,w1R,w1I,w2R,w2I,w3R,w3I,w4R,w4I,theta1,theta2,theta3,theta4;
     int k1,k2,k3,k4;
     const1 = -2.0*M_PI/3;
     const1R = cos(const1);
     const1I = sin(const1);
     const2R = const1R;
     const2I = -const1I;
     while (q>0){
           for (i=0;i<tmp;i++){
               theta1 = -2.0*i*M_PI/(tmp*3);
               theta2 = theta1*2;
               w1R = cos(theta1);
               w1I = sin(theta1);
               w2R = cos(theta2);
               w2I = sin(theta2);
               for (j=i;j<N;j+=(tmp*3)){
                   k1 = j+tmp;
                   k2 = j+tmp*2;
                   tmp1R = w1R*yR[k1] - w1I*yI[k1];
                   tmp1I = w1R*yI[k1] + w1I*yR[k1];
                   tmp2R = w2R*yR[k2] - w2I*yI[k2];
                   tmp2I = w2R*yI[k2] + w2I*yR[k2];
                   yR[k2] = yR[j] + (tmp1R*const2R-tmp1I*const2I) + (tmp2R*const1R-tmp2I*const1I);
                   yI[k2] = yI[j] + (tmp1R*const2I+tmp1I*const2R) + (tmp2R*const1I+tmp2I*const1R);
                   yR[k1] = yR[j] + (tmp1R*const1R-tmp1I*const1I) + (tmp2R*const2R-tmp2I*const2I);
                   yI[k1] = yI[j] + (tmp1R*const1I+tmp1I*const1R) + (tmp2R*const2I+tmp2I*const2R);
                   yR[j] += (tmp1R + tmp2R);
                   yI[j] += (tmp1I + tmp2I);
               }
           }
           tmp *=3;
           q --;
     }
     const1 = -2.0*M_PI/5;
     const2 = const1*2;
     const1R = cos(const1);
     const1I = sin(const1);
     const2R = cos(const2);
     const2I = sin(const2);
     const3R = const2R;
     const3I = -const2I;
     const4R = const1R;
     const4I = -const1I;
     while (r>0){
           for (i=0;i<tmp;i++){
               theta1 = -2.0*i*M_PI/(tmp*5);
               theta2 = theta1*2;
               theta3 = theta1*3;
               theta4 = theta1*4;
               w1R = cos(theta1);
               w1I = sin(theta1);
               w2R = cos(theta2);
               w2I = sin(theta2);
               w3R = cos(theta3);
               w3I = sin(theta3);
               w4R = cos(theta4);
               w4I = sin(theta4);
               for (j=i;j<N;j+=(tmp*5)){
                   k1 = j+tmp;
                   k2 = j+tmp*2;
                   k3 = j+tmp*3;
                   k4 = j+tmp*4;                   
                   tmp1R = w1R*yR[k1] - w1I*yI[k1];
                   tmp1I = w1R*yI[k1] + w1I*yR[k1];
                   tmp2R = w2R*yR[k2] - w2I*yI[k2];
                   tmp2I = w2R*yI[k2] + w2I*yR[k2];
                   tmp3R = w3R*yR[k3] - w3I*yI[k3];
                   tmp3I = w3R*yI[k3] + w3I*yR[k3];
                   tmp4R = w4R*yR[k4] - w4I*yI[k4];
                   tmp4I = w4R*yI[k4] + w4I*yR[k4];
                   yR[k4] = yR[j] + (tmp1R*const4R-tmp1I*const4I) + (tmp2R*const3R-tmp2I*const3I) + (tmp3R*const2R-tmp3I*const2I) + (tmp4R*const1R-tmp4I*const1I);
                   yI[k4] = yI[j] + (tmp1R*const4I+tmp1I*const4R) + (tmp2R*const3I+tmp2I*const3R) + (tmp3R*const2I+tmp3I*const2R) + (tmp4R*const1I+tmp4I*const1R);
                   yR[k3] = yR[j] + (tmp1R*const3R-tmp1I*const3I) + (tmp2R*const1R-tmp2I*const1I) + (tmp3R*const4R-tmp3I*const4I) + (tmp4R*const2R-tmp4I*const2I);
                   yI[k3] = yI[j] + (tmp1R*const3I+tmp1I*const3R) + (tmp2R*const1I+tmp2I*const1R) + (tmp3R*const4I+tmp3I*const4R) + (tmp4R*const2I+tmp4I*const2R);
                   yR[k2] = yR[j] + (tmp1R*const2R-tmp1I*const2I) + (tmp2R*const4R-tmp2I*const4I) + (tmp3R*const1R-tmp3I*const1I) + (tmp4R*const3R-tmp4I*const3I);
                   yI[k2] = yI[j] + (tmp1R*const2I+tmp1I*const2R) + (tmp2R*const4I+tmp2I*const4R) + (tmp3R*const1I+tmp3I*const1R) + (tmp4R*const3I+tmp4I*const3R);
                   yR[k1] = yR[j] + (tmp1R*const1R-tmp1I*const1I) + (tmp2R*const2R-tmp2I*const2I) + (tmp3R*const3R-tmp3I*const3I) + (tmp4R*const4R-tmp4I*const4I);
                   yI[k1] = yI[j] + (tmp1R*const1I+tmp1I*const1R) + (tmp2R*const2I+tmp2I*const2R) + (tmp3R*const3I+tmp3I*const3R) + (tmp4R*const4I+tmp4I*const4R);
                   yR[j] += (tmp1R + tmp2R + tmp3R + tmp4R);
                   yI[j] += (tmp1I + tmp2I + tmp3I + tmp4I);
               }
           }
           tmp *=5;
           r --;
     }
}

int main(int argc, char *argv[])
{
    clock_t tic,toc;
    int i,j,p,q,r;
    long long int N=1;
    cin >> p >> q >> r;
    for (i=0;i<p;i++)
        N*=2;
    for (i=0;i<q;i++)
        N*=3;
    for (i=0;i<r;i++)
        N*=5;
//    cout << "N = " << N << endl;
    double *xR,*xI,*yR,*yI,*rR,*rI;
    xR = (double*) malloc(N*sizeof(double));
    xI = (double*) malloc(N*sizeof(double));
    yR = (double*) malloc(N*sizeof(double));
    yI = (double*) malloc(N*sizeof(double));
    rR = (double*) malloc(N*sizeof(double));
    rI = (double*) malloc(N*sizeof(double));
    for (i=0;i!=N;i++){
        xR[i] = i;
        xI[i] = 0.0;
    }
//    sft(xR,xI,rR,rI,N);
    index_transform(xR,xI,yR,yI,N,p,q,r);
    tic = clock();
    fft(yR,yI,N,p,q,r);
    toc = clock();
/*    for (i=0;i<N;i++){
//        cout << i << " " << rR[i] << " " << rI[i] << endl;
//        cout << i << " " << yR[i] << " " << yI[i] << endl;
          double compare1 = (rR[i]>yR[i]) ? rR[i]-yR[i] : yR[i]-rR[i];
          double compare2 = (rI[i]>yI[i]) ? rI[i]-yI[i] : yI[i]-rI[i];
          if (compare1>pow(10,6) || compare2>pow(10,6)) cout << "Wrong " << i << endl;
    }
*/
    cout << "Time: " << 1.0 * (toc-tic)/CLOCKS_PER_SEC << endl;
    system("PAUSE");
    return EXIT_SUCCESS;
}
