/* Optional likelihood accelerator for the documented linear Gaussian model.
   The Python filter is the reference and is always available. */
#include <math.h>
#include <stdlib.h>

double revision_filter(int n, int k, int ni, const double *time,
    const int *image, const int *instrument, const double *value,
    const double *error, const double *variance, const double *rate,
    const double *diffusion, const double *jitter, double zero,
    double reference, double year) {
    double *p = calloc((size_t)k*k, sizeof(double));
    double *m = calloc((size_t)k, sizeof(double));
    double *a = malloc((size_t)k*sizeof(double));
    double *v = malloc((size_t)k*sizeof(double));
    if (!p || !m || !a || !v) {free(p);free(m);free(a);free(v);return NAN;}
    for (int j=0;j<k;j++) p[j*k+j]=variance[j];
    double ll=0.;
    for(int z=0;z<n;z++) {
        double dt=z ? time[z]-time[z-1] : 0.;
        for(int j=0;j<k;j++) {a[j]=exp(-rate[j]*dt);m[j]*=a[j];}
        for(int j=0;j<k;j++) {
            for(int l=0;l<k;l++) p[j*k+l]*=a[j]*a[l];
            p[j*k+j]+=variance[j]*(-expm1(-2*rate[j]*dt))+diffusion[j]*dt;
        }
        int im=image[z];
        int indices[5]={0, 1+im, 1+ni+im, 1+2*ni+im,
                       1+3*ni+instrument[z]*ni+im};
        double h[5]={1.,1.,(time[z]-reference)/year,1.,1.};
        double innovation=value[z]-zero;
        for(int j=0;j<5;j++) innovation-=h[j]*m[indices[j]];
        for(int j=0;j<k;j++) {
            v[j]=0.;
            for(int l=0;l<5;l++) v[j]+=p[j*k+indices[l]]*h[l];
        }
        double s=error[z]*error[z]+jitter[im]*jitter[im];
        for(int j=0;j<5;j++) s+=h[j]*v[indices[j]];
        if (!(s>0.) || !isfinite(s)) {ll=NAN;break;}
        ll-=.5*(1.8378770664093453+log(s)+innovation*innovation/s);
        for(int j=0;j<k;j++) m[j]+=v[j]*innovation/s;
        for(int j=0;j<k;j++) for(int l=0;l<k;l++) p[j*k+l]-=v[j]*v[l]/s;
        for(int j=0;j<k;j++) for(int l=0;l<j;l++) {
            double symmetric=.5*(p[j*k+l]+p[l*k+j]);
            p[j*k+l]=p[l*k+j]=symmetric;
        }
    }
    free(p);free(m);free(a);free(v);
    return ll;
}
