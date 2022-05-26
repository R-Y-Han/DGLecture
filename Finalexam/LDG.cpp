/**
 * @file LDG.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief LDG scheme for u_t + u_x = \epsilon * u_{xxx}, P^1 and P^2
 * @version 0.1
 * @date 2022-05-26
 * 
 * @copyright Copyright (c) 2022, R.Y. Han, All Rights Reserved.
 * 
 */

#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

/**
 * @name Configration
 * @brief Some setups for the scheme
 * @{
 */
const double T = 1;   /**< final time */
double dt = 1e-4;   /**< time step length */
const double epsilon = 1;   /**< dispersion coefficients */
const int n = 20;   /**< number of cells */
const int k = 1;    /**< polynomial order of the basis */
const double pi = 3.1415926;    /**< \pi */
const double h = 2.0 * pi / n;  /**< spacial step length, uniform mesh */
double p[n+1];  /**< cell boundary's location, p_j = j * h = x_{j+1/2}, j =0,1, ..., n */

const double lobattopoint[5] = {-1.0, -0.6546536707079771, 0, 0.6546536707079771, 1.0}; /**< Gauss-Lobatto quadrature points */
const double lobattoco[5] = {0.1, 0.5444444444444444, 0.7111111111111111, 0.5444444444444444, 0.1}; /**< Gauss-Lobatto quadrature weights */
/** @} */

/**
 * @brief The initial value
 * 
 * @param x 
 * @return sin(x)
 */
double u_0(double x);

/**
 * @brief The exact solution
 * 
 * @param x 
 * @param t 
 * @return double sin(x - (1 + epsilon )t )
 */
double u_exact(double x, double t);

/**
 * @brief Legendre polynomial
 * 
 * @param l order of the l-th polynomial
 * @param x 
 * @return the l-th Legendre polynomial
 */
double phi(int l, double x);

/**
 * @brief Initialize the numerical results u_h
 * 
 * @return The L^2 projection of u_0 to the solution space
 */
double** initial();

/**
 * @brief Compute the parameters p
 * 
 * @param ut u^n
 * @return p
 */
double** getp(double **ut);

/**
 * @brief Compute the parameters q
 * 
 * @param pt p
 * @return q
 */
double** getq(double **pt);

/**
 * @brief Function to compute u^{n+1}
 * 
 * @param ut Numerical results of u^n
 * @return du/dt 
 */
double** L(double **ut);

/**
 * @brief 2-stage 2-order RK
 * 
 * @param un numerical results u^n
 * @param dtn time step length
 * @return u^{n+1}
 */
double** RK22(double** un, double dtn);

/**
 * @brief 3-stage 3-order RK
 * 
 * @param un numerical results u^n
 * @param dtn time step length
 * @return u^{n+1}
 */
double** RK33(double** un, double dtn);

/**
 * @brief RK, choose RK22 or RK33 according to k
 * 
 * @param k the order of polynomial
 * @param un u^n
 * @param dtn time step length
 * @return u^{n+1}
 */
double** RK(int k, double** un, double dtn);

int main()
{
    int i, j, l;
    double t, temp1, temp2, norm1, norm2, xi;

    double** u1 = new double* [n];
    double** u2 = new double* [n];
    for (i=0; i<n; i++)
    {
        u1[i] = new double [k+1];
        u2[i] = new double [k+1];
    }
    for (j=0; j<=n; j++)
    {
        p[j] = j * h;
    }

    u1 = initial();
    t = 0;
    while(t<T)
    {
        if ( t + dt >T)
        {
            dt = T - t;
        }
        t = t + dt;
        u2 = RK(k,u1,dt);
        u1 = u2;
        cout<<"t="<<t<<endl;
    }
    
    //¼ÆËãÎó²î
    norm1 = 0;
    norm2 = 0;
    for (j=1; j<=n; j++)
    {
        for (i=0; i<5; i++)
        {
            xi = lobattopoint[i];
            temp1 = 0;
            for (l=0; l<=k; l++)
            {
                temp1 = temp1 + u1[j-1][l] * phi(l,xi);
            }

            temp2 = h * (xi + 1) / 2. + p[j-1];
            temp1 = u_exact(temp2,T) - temp1;

            if (abs(temp1) > norm2)
            {
                norm2 = abs(temp1);
            }

            temp1 = temp1 * temp1;
            norm1 = norm1 + lobattoco[i] * temp1;
        }

    }
    norm1 = norm1 * h / 2.;
    norm1 = sqrt(norm1);
    cout<<"L2="<<norm1<<endl<<"Linf="<<norm2<<endl;//*/

//¿ªÊ¼»­Í¼
    {
        const char* fn = "DGLecture\\Finalexam\\LDG.plt";
        remove(fn);
        fstream f, f1;
        f.open(fn, ios::out | ios :: app);
        f<<"VARIABLES="<<"X"<<","<<"u_h"<<","<<"u"<<endl;
        for (j=1; j<=n; j++)
        {
            temp1 = 0;
            for (l=0; l<=k; l++)
            {
                temp1 = temp1 + u1[j-1][l] * phi(l,0);
            }
            f<<"\t"<<(p[j-1] + p[j])/2.0<<"\t"<<temp1<<"\t"<<u_exact((p[j-1] + p[j])/2.0,T)<<endl;

        }
        f.close();
    }

    for (i=0; i<n; i++)
    {
        delete[] u1[i];
        delete[] u2[i];
    }
    delete[] u1;
    delete[] u2;

    system("pause");
}

double u_0(double x)
{
    return sin(x);
}

double u_exact(double x, double t)
{
    return sin(x - (1 + epsilon) * t);
}

double phi(int l, double x)
{
    if (l==0)
    {
        return 1;
    }
    else if (l == 1)
    {
        return x;
    }
    else if (l == 2)
    {
        return (3*x*x - 1)/2;
    }
    else{
        return 0;
    }
}

double** initial()
{
    double ans, temp;
    int j, l, m;
    double** ut = new double* [n];
    double* Bt = new double [n];
    for (j=0; j<n; j++)
    {
        ut[j] = new double [k+1];
    }

    for (j=1; j<=n; j++)
    {
        for (m=0; m<=k; m++)
        {
            ans = 0;
            for (l=0; l<5; l++)
            {
                temp = h * (lobattopoint[l] + 1)/2 + p[j-1];
                ans = ans + lobattoco[l] * u_0(temp) * phi(m,lobattopoint[l]);
            }
            ans = ans / 2;
            Bt[m] = ans;
        }

        double A[3][3] = {{1,0,0},{0,3,0},{0,0,5}};
        for (m=0; m<=k ;m++)
        {
            ut[j-1][m] = 0;
            for (l=0; l<=k; l++)
            {
                ut[j-1][m] = ut[j-1][m] + A[m][l] * Bt[l];
            }
        } 
    }
    delete[] Bt;
    return ut;
}

double** getp(double **ut)
{
    int j,l,m, temp;
    double ul;
    double** ans = new double *[n];
    for (j=0; j<n; j++)
    {
        ans[j] = new double [k+1];
    }

    double A[3] = {1,3,5};
    double B[3][3] = {{0,0,0},{2,0,0},{0,2,0}};

    for (j=1; j<=n; j++)
    {
        for (m=0; m<=k; m++)
        {
            ans[j-1][m] = 0;
            
            for (l=0; l<=k; l++)
            {
                ans[j-1][m] = ans[j-1][m] - B[m][l] * ut[j-1][l];
            }

            ul = 0;
            for (l=0; l<=k; l++)
            {
                ul = ul + ut[j-1][l] * phi(l,1);
            }
            ans[j-1][m] = ans[j-1][m] + ul * phi(m,1);

            ul = 0;
            temp = j-2;
            if (temp == -1)
            {
                temp = n-1;
            }
            for (l=0; l<=k; l++)
            {
                ul = ul + ut[temp][l] * phi(l,1);
            }
            ans[j-1][m] = ans[j-1][m] - ul * phi(m,-1);
            ans[j-1][m] = ans[j-1][m] * A[m] / h;
        }
    }
    return ans;
}

double** getq(double **pt)
{
    int j,l,m, temp;
    double pl;
    double** ans = new double *[n];
    for (j=0; j<n; j++)
    {
        ans[j] = new double [k+1];
    }

    double A[3] = {1,3,5};
    double B[3][3] = {{0,0,0},{2,0,0},{0,2,0}};

    for (j=1; j<=n; j++)
    {
        for (m=0; m<=k; m++)
        {
            ans[j-1][m] = 0;
            
            for (l=0; l<=k; l++)
            {
                ans[j-1][m] = ans[j-1][m] - B[m][l] * pt[j-1][l];
            }

            pl = 0;
            for (l=0; l<=k; l++)
            {
                pl = pl + pt[j-1][l] * phi(l,1);
            }
            ans[j-1][m] = ans[j-1][m] + pl*phi(m,1);

            pl = 0;
            temp = j-2;
            if (temp == -1)
            {
                temp = n-1;
            }
            for (l=0; l<=k; l++)
            {
                pl = pl + pt[temp][l] * phi(l,1);
            }
            ans[j-1][m] = ans[j-1][m] - pl * phi(m,-1);
            ans[j-1][m] = ans[j-1][m] * A[m] / h;
        }
    }
    return ans;
}

double** L(double **ut)
{
    int j,l,m,temp;
    double ul, qr;
    double** ans = new double *[n];
    double** qt;
    double** pt;
    for (j=0; j<n; j++)
    {
        ans[j] = new double [k+1];
    }

    pt = getp(ut);
    qt = getq(pt);

    double A[3] = {1,3,5};
    double B[3][3] = {{0,0,0},{2,0,0},{0,2,0}};

    for (j=1; j<=n; j++)
    {
        for (m=0; m<=k; m++)
        {
            ans[j-1][m] = 0;
            
            for (l=0; l<=k; l++)
            {
                ans[j-1][m] = ans[j-1][m] + B[m][l] * (ut[j-1][l] - epsilon* qt[j-1][l]);
            }

            ul = 0;
            for (l=0; l<=k; l++)
            {
                ul = ul + ut[j-1][l]*phi(l,1);
            }
            ans[j-1][m] = ans[j-1][m] - ul * phi(m,1);

            ul = 0;
            temp = j-2;
            if (temp == -1)
            {
                temp = n-1;
            }
            for (l=0; l<=k; l++)
            {
                ul = ul + ut[temp][l] * phi(l,1);
            }
            ans[j-1][m] = ans[j-1][m] + ul * phi(m,-1);

            qr = 0;
            temp = j;
            if (temp == n)
            {
                temp = 0;
            }
            for (l=0; l<=k; l++)
            {
                qr = qr + qt[temp][l] * phi(l,-1);
            }
            ans[j-1][m] = ans[j-1][m] + epsilon * qr * phi(m,1);

            qr = 0;
            for (l=0; l<=k; l++)
            {
                qr = qr + qt[j-1][l] * phi(l,-1);
            }
            ans[j-1][m] = ans[j-1][m] - epsilon * qr * phi(m,-1);
            
            ans[j-1][m] = ans[j-1][m] * A[m] / h;
        }
    }

    for (j=0; j<n; j++)
    {
        delete[] pt[j];
        delete[] qt[j];
    }
    delete[] pt, qt;

    return ans;
}

double** RK22(double** un, double dtn)
{
    int j,l,m;
    double ul;
    double** ans = new double* [n];
    for (j=0; j<n; j++)
    {
        ans[j] = new double [k+1];
    }
    double** u0 = new double* [n];
    double** u1 = new double* [n];
    double** u2 = new double* [n];

    for (j=0; j<n; j++)
    {
        u0[j] = new double [k+1];
        u1[j] = new double [k+1];
        u2[j] = new double [k+1];
    }

    for (j=1; j<=n; j++)
    {
        for (l=0; l<=k; l++)
        {
            u0[j-1][l] = un[j-1][l];
        }
    }

    u1 = L(u0);
    for (j=1; j<=n; j++)
    {
        for (l=0; l<=k; l++)
        {
            u1[j-1][l] = u1[j-1][l] * dtn  + u0[j-1][l];
        }
    }
        
    u2 = L(u1);
    for (j=1; j<=n;j++)
    {
        for (l=0; l<=k; l++)
        {
            u2[j-1][l] = u0[j-1][l] * 0.5 + u1[j-1][l] * 0.5 + u2[j-1][l] * 0.5 * dtn;
            ans[j-1][l] = u2[j-1][l];
        }
    }

    delete[] u0;
    delete[] u1;
    delete[] u2;

    return ans;
}

double** RK33(double** un, double dtn)
{
    int j,l,m;
    double ul;
    double** ans = new double* [n];
    for (j=0; j<n; j++)
    {
        ans[j] = new double [k+1];
    }
    double** u0 = new double* [n];
    double** u1 = new double* [n];
    double** u2 = new double* [n];
    double** u3 = new double* [n];

    for (j=0; j<n; j++)
    {
        u0[j] = new double [k+1];
        u1[j] = new double [k+1];
        u2[j] = new double [k+1];
        u3[j] = new double [k+1];
    }

    for (j=1; j<=n; j++)
    {
        for (l=0; l<=k; l++)
        {
            u0[j-1][l] = un[j-1][l];
        }
    }

    u1 = L(u0);
    for (j=1; j<=n; j++)
    {
        for (l=0; l<=k; l++)
        {
            u1[j-1][l] = u1[j-1][l] * dtn  + u0[j-1][l];
        }
    }
    
    u2 = L(u1);
    for (j=1; j<=n;j++)
    {
        for (l=0; l<=k; l++)
        {
            u2[j-1][l] = u0[j-1][l] * 3 / 4.0 + u1[j-1][l] /4 + u2[j-1][l] * dtn / 4.0;
        }
    }

    u3 = L(u2);
    for (j=1; j<=n; j++)
    {
        for (l=0; l<=k; l++)
        {
            u3[j-1][l] = u0[j-1][l] / 3.0 + 2 * u2[j-1][l] / 3.0 + 2 * dtn * u3[j-1][l] /3.0;
            ans[j-1][l] = u3[j-1][l];
        }
    }

    delete[] u0;
    delete[] u1;
    delete[] u2;
    delete[] u3;

    return ans;
}

double** RK(int k, double** un, double dtn)
{
    double** u2;
    if ( k == 1 )
    {
        u2 = RK22(un,dtn);
    }
    else if (k == 2)
    {
        u2 = RK33(un,dtn);
    }
    else{
        ;
    }

    return u2;
}