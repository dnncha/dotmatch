/* MAGeCK-compatible EM loop. BSD-3-Clause; see LICENSE.
 * No fast-math, NumPy ABI, or linked BLAS dependency. SciPy special-function
 * capsules are checked before use. Input buffers stay alive while the GIL is
 * released; all mutable working storage belongs to this call.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>

typedef double (*special1)(double, int);
typedef double (*special2)(double, double, int);
static special1 loggamma_fn;
static special2 xlog1py_fn;

static int initialize_special(void) {
    PyObject *module = PyImport_ImportModule("scipy.special.cython_special");
    if (!module) return -1;
    PyObject *api = PyObject_GetAttrString(module, "__pyx_capi__");
    Py_DECREF(module);
    if (!api) return -1;
    PyObject *gamma = PyDict_GetItemString(api, "gammaln");
    PyObject *xlog = PyDict_GetItemString(api, "__pyx_fuse_1xlog1py");
    if (!gamma || !xlog) {
        Py_DECREF(api);
        PyErr_SetString(PyExc_ImportError, "Unsupported SciPy special-function C API");
        return -1;
    }
    loggamma_fn = (special1)PyCapsule_GetPointer(gamma, "double (double, int __pyx_skip_dispatch)");
    xlog1py_fn = (special2)PyCapsule_GetPointer(xlog, "double (double, double, int __pyx_skip_dispatch)");
    Py_DECREF(api);
    if (!loggamma_fn || !xlog1py_fn) {
        PyErr_Clear();
        PyErr_SetString(PyExc_ImportError, "Unsupported SciPy special-function capsule signature");
        return -1;
    }
    return 0;
}

/* Partial-pivot LU for the small ridge system. No extra regularization or
 * pseudoinverse fallback: a singular system is reported to the caller. */
static int solve(double *a, double *b, int p) {
    for (int k=0; k<p; k++) {
        int pivot=k;
        for (int i=k+1; i<p; i++)
            if (fabs(a[i*p+k]) > fabs(a[pivot*p+k])) pivot=i;
        if (a[pivot*p+k] == 0.0) return -1;
        if (pivot != k) {
            for (int j=0; j<p; j++) {
                double temp=a[k*p+j]; a[k*p+j]=a[pivot*p+j]; a[pivot*p+j]=temp;
            }
            double temp=b[k]; b[k]=b[pivot]; b[pivot]=temp;
        }
        for (int i=k+1; i<p; i++) {
            double factor=a[i*p+k]/a[k*p+k];
            for (int j=k+1; j<p; j++) a[i*p+j] -= factor*a[k*p+j];
            b[i] -= factor*b[k];
        }
    }
    for (int i=p-1; i>=0; i--) {
        double value=b[i];
        for (int j=i+1; j<p; j++) value -= a[i*p+j]*b[j];
        b[i]=value/a[i*p+i];
    }
    return 0;
}

static double likelihood(double k, double factorial, double mu, double alpha) {
    double square=mu*mu, variance=mu+alpha*square;
    double p=mu/variance, r=square/(variance-mu);
    if (!(r>0 && p>0 && p<=1) || isnan(k)) return 0;
    if (!(k>=0 && k<=INFINITY && k==floor(k))) return -INFINITY;
    double value=loggamma_fn(r+k,1)-factorial-loggamma_fn(r,1)
                 + r*log(p)+xlog1py_fn(k,-p,1);
    return isnan(value) ? 0 : value;
}

static PyObject *loop(PyObject *self, PyObject *args) {
    PyObject *objects[6]; Py_buffer buffers[6]={{0}};
    int n,c,other,estimate,update;
    double ridge;
    PyObject *result=NULL;
    double *storage=NULL;
    if (!PyArg_ParseTuple(args,"OOOOOOiiidpp",&objects[0],&objects[1],&objects[2],
                          &objects[3],&objects[4],&objects[5],&n,&c,&other,&ridge,&estimate,&update)) return NULL;
    if (n<1 || c<1 || other<1 || n>10000 || c>10000 || other>10000) {
        PyErr_SetString(PyExc_ValueError,"Invalid native model dimensions"); return NULL;
    }
    Py_ssize_t p=(Py_ssize_t)n+c, m=(Py_ssize_t)n*(2*other+1);
    if (p>10000 || m>1000000 || m*p>100000000) {
        PyErr_SetString(PyExc_ValueError,"Native model exceeds bounded working dimensions"); return NULL;
    }
    Py_ssize_t lengths[6]={m*p,m,m,p,n,m};
    for (int i=0;i<6;i++) {
        if (PyObject_GetBuffer(objects[i],&buffers[i],PyBUF_FORMAT|PyBUF_C_CONTIGUOUS)<0) goto cleanup;
        if (!buffers[i].format || strcmp(buffers[i].format,"d") || buffers[i].itemsize!=sizeof(double)
            || buffers[i].len!=lengths[i]*(Py_ssize_t)sizeof(double)
            || (uintptr_t)buffers[i].buf % _Alignof(double)) {
            PyErr_SetString(PyExc_ValueError,"Expected contiguous native float64 buffers of matching size"); goto cleanup;
        }
    }
    /* Seven m-vectors, two p-vectors, one n-vector, and three p*p matrices. */
    Py_ssize_t length=7*m+2*p+n+3*p*p;
    storage=PyMem_Calloc((size_t)length,sizeof(double));
    if (!storage) {PyErr_NoMemory(); goto cleanup;}
    const double *x=buffers[0].buf,*counts=buffers[1].buf,*sizes=buffers[2].buf,*alpha=buffers[5].buf;
    double *beta=storage, *eff=beta+p, *mu=eff+n, *res=mu+m, *weights=res+m;
    double *response=weights+m,*rounded=response+m,*factorial=rounded+m,*logmean=factorial+m;
    double *proposal=logmean+m,*gram=proposal+p,*regularized=gram+p*p,*factor=regularized+p*p;
    memcpy(beta,buffers[3].buf,p*sizeof(double)); memcpy(eff,buffers[4].buf,n*sizeof(double));
    int status=0, iteration=1, singular=0;
    Py_BEGIN_ALLOW_THREADS
    for (Py_ssize_t i=0;i<m;i++) {rounded[i]=nearbyint(counts[i]); factorial[i]=loggamma_fn(rounded[i]+1,1);}
    while (1) {
        for (Py_ssize_t i=0;i<m;i++) {
            double value=0;
            for (Py_ssize_t j=0;j<p;j++) value+=x[i*p+j]*beta[j];
            logmean[i]=value; mu[i]=sizes[i]*exp(value); res[i]=counts[i]-mu[i];
            double coefficient=1;
            if (estimate && i>=n) coefficient=i<n*(other+1) ? eff[i%n] : 1-eff[i%n];
            response[i]=(estimate ? res[i]*coefficient : res[i])/mu[i]+value;
        }
        if (estimate && update) {
            for (int guide=0;guide<n;guide++) {
                double minimum=INFINITY;
                for (int sample=1;sample<=other;sample++) {
                    Py_ssize_t i=(Py_ssize_t)sample*n+guide;
                    double value=0;
                    for (int j=0;j<n;j++) value+=x[i*p+j]*beta[j];
                    double baseline=sizes[i]*exp(value);
                    double difference=likelihood(rounded[i],factorial[i],baseline,alpha[i])
                                     -likelihood(rounded[i],factorial[i],mu[i],alpha[i]);
                    if (difference>100) difference=100;
                    if (isnan(difference) || difference<minimum) minimum=difference;
                }
                eff[guide]=1/(1+exp(minimum));
            }
        }
        double maximum=-INFINITY;
        int missing_weight=0;
        for (Py_ssize_t i=0;i<m;i++) {
            double w=1/(1/mu[i]+alpha[i]);
            if (estimate && i>=n) {
                double coefficient=i<n*(other+1) ? eff[i%n] : 1-eff[i%n];
                if (coefficient<0.01) coefficient=0.01;
                w*=coefficient;
            }
            weights[i]=w;
            if (isnan(w)) missing_weight=1;
            if (w>maximum) maximum=w;
        }
        double floor_weight=missing_weight ? NAN : maximum/100;
        for (Py_ssize_t i=0;i<m;i++) {
            if (isnan(floor_weight)) weights[i]=NAN;
            else if (weights[i]<floor_weight) weights[i]=floor_weight;
        }
        for (Py_ssize_t j=0;j<p;j++) {
            double value=0;
            for (Py_ssize_t i=0;i<m;i++) value+=x[i*p+j]*(weights[i]*response[i]);
            proposal[j]=value;
            for (Py_ssize_t k=0;k<p;k++) {
                value=0;
                for (Py_ssize_t i=0;i<m;i++) value+=x[i*p+j]*(weights[i]*x[i*p+k]);
                gram[j*p+k]=value; regularized[j*p+k]=value+(j==k ? ridge : 0);
            }
        }
        memcpy(factor,regularized,p*p*sizeof(double));
        if (solve(factor,proposal,(int)p)<0) {singular=1;break;}
        double sum=0, largest=0, difference=0, beta_size=0;
        for (Py_ssize_t j=0;j<p;j++) {
            sum+=proposal[j];
            if (fabs(proposal[j])>largest) largest=fabs(proposal[j]);
            if (j>=n) {double delta=proposal[j]-beta[j];difference+=delta*delta;beta_size+=proposal[j]*proposal[j];}
        }
        if (isnan(sum) || largest>100) {status=2;break;}
        memcpy(beta,proposal,p*sizeof(double)); iteration++;
        if (fabs(beta_size)<1e-9) beta_size=1;
        if (difference/beta_size<1e-9) break;
        if (iteration>1000) {status=1;break;}
    }
    Py_END_ALLOW_THREADS
    if (singular) {PyErr_SetString(PyExc_ArithmeticError,"Singular ridge system");goto cleanup;}
    result=PyTuple_New(8);
    if (!result) goto cleanup;
    PyObject *status_object=PyLong_FromLong(status);
    if (!status_object) {Py_CLEAR(result);goto cleanup;}
    PyTuple_SET_ITEM(result,0,status_object);
    double *outputs[7]={eff,beta,mu,res,weights,gram,regularized};
    Py_ssize_t sizes_out[7]={n,p,m,m,m,p*p,p*p};
    for (int i=0;i<7;i++) {
        PyObject *bytes=PyBytes_FromStringAndSize((char*)outputs[i],sizes_out[i]*sizeof(double));
        if (!bytes) {Py_CLEAR(result);goto cleanup;}
        PyTuple_SET_ITEM(result,i+1,bytes);
    }
cleanup:
    PyMem_Free(storage);
    for (int i=0;i<6;i++) if (buffers[i].obj) PyBuffer_Release(&buffers[i]);
    return result;
}

static PyMethodDef methods[]={{"loop",loop,METH_VARARGS,"Execute the validated native EM loop."},{NULL,NULL,0,NULL}};
static struct PyModuleDef module={PyModuleDef_HEAD_INIT,"_native",NULL,-1,methods};
PyMODINIT_FUNC PyInit__native(void) {
    if (initialize_special()<0) return NULL;
    return PyModule_Create(&module);
}
