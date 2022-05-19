#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

//先定义一些全局的变量
const int n = 160;   //划分单元个数, by CFL condition, n should larger than 10 when dt = h^2/2
const int k = 1;    //多项式最高次项次数
const double h = (double) 2/n;    //空间步长
const double dt = (double) h * h ;    //时间步长
const double pi = 3.1415926;

double p[n+1];  //节点位置，p_j = j * h = x_{j+1/2}, j =0,1, ..., n
//存储不变的系数矩阵

const double lobattopoint[5] = {-1.0, -0.6546536707079771, 0, 0.6546536707079771, 1.0};
const double lobattoco[5] = {0.1, 0.5444444444444444, 0.7111111111111111, 0.5444444444444444, 0.1};

//*************函数声明*************//
double u_0(double y);   //初值
double u_exact(double y, double t);   //真解
double f(double y); //通量函数
double phi(int l, double y);    //参考单元基函数
double** initial();   //计算初始时刻u_j
double* cellaverage(double** un); //求单元平均
double minimal(double a1, double a2, double a3);
double** modify(double** un);
double flux(double ul, double ur);  //数值通量计算
double** L(double** ut);   //用于计算RK的函数，u_t= F(u)
double** RK22(double** un);    //2步二阶RK
double** RK33(double** un);
double** RK(int k, double** un);
//************声明完毕**************//

int main()
{
    int i, j, l;
    double t, temp1, temp2, norm1, norm2, xi;
    double T = 1.5;
    double** u1 = new double* [n];
    double** u2 = new double* [n];
    for (i=0; i<n; i++)
    {
        u1[i] = new double [k+1];
        u2[i] = new double [k+1];
    }
    for (j=0; j<=n; j++)
    {
        p[j] = j * h - 1;
    }


    u1 = initial();
    t = 0;
    while(t<T - 1e-10)
    {

        t = t + dt;
        u2 = RK(k,u1);

        u1 = u2;
        cout<<t<<endl;
    }
    
    //*
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

    {
        const char* fn = "DGLecture\\homework4\\MPP.plt";
        remove(fn);
        fstream f, f1;
        f.open(fn, ios::out | ios :: app);
        f<<"VARIABLES="<<"X"<<","<<"u_h"<<","<<"u"<<endl;
        for (j=1; j<=n; j++)
        {
            temp1 = 0;
            temp2 = 0;
            for (l=0; l<=k; l++)
            {
                temp1 = temp1 + u1[j-1][l] * phi(l,0);
                temp2 = temp2 + u1[j-1][l] * phi(l,0);
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

double u_0(double y)
{
    return sin(pi*y) / 3.0 + 2.0 / 3.0 ;
}

double u_exact(double y, double t)
{
    int time=1;
    double ans, xi1=0, xi2=0.1;
    double e = 1e-6;

    //if (t==1.5)  //1.5时刻发生激波，单独算
    {
        if (y < -1+2*t/3.0)
        {
            y = y+2;
        }
    }
    
        while(abs(xi1 - xi2)>e)
        {
            xi1 = xi2;
            xi2 = xi1 + ( y - u_0(xi1) * t - xi1 ) / (t * pi * cos(pi*xi1)/3.0 + 1);
            time++;
            if (time > 10000)
            {
                cout<<"error"<<endl;
                cout<<y<<endl;
                break;
            }
        }
        ans = u_0(xi2);

    return ans;
}

double f(double y)
{
    return y * y /2;
}

double phi(int l, double y)
{
    if (l==0)
    {
        return 1;
    }
    else if (l == 1)
    {
        return y;
    }
    else if (l == 2)
    {
        return (3*y*y - 1)/2;
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

double* cellaverage(double** un)
{
    int j, l;
    double* ca = new double[n];
    double D[3] = {2, 0, 0};

    for (j=1; j<=n; j++)
    {
        /*ca[j-1] = 0;
        for (l=0; l<=k; l++)
        {
            ca[j-1] = ca[j-1] + un[j-1][l] * D[l];
        }*/
        ca[j-1] = un[j-1][0];
    }

    return ca;
}

double minimal(double a1, double a2, double a3)
{
    double ans, temp;
    
    temp = min(a1, a2);
    ans = min(temp, a3);

    return ans;
}

double** modify(double** un)
{
    int j, l, r;
    double** mod = new double* [n];
    for (j=1; j<=n; j++)
    {
        mod[j-1] = new double [k+1];
    }

    double theta, M, m, temp, a1, a2;
    double* ca = new double [n];
    
    ca = cellaverage(un);

    for (j=1; j<=n; j++)
    {
        M = ca[j-1];
        m = ca[j-1];
        
        for (l=0; l<5; l++)
        {
            temp = 0;
            for (r=0; r<=k; r++)
            {
                temp = temp + un[j-1][r] * phi(r,lobattopoint[l]);
            }
            //temp表示第l个积分节点的数值解

            if (temp > M)
            {
                M = temp;
            }
            if (temp < m)
            {
                m = temp;
            }
        }

        a1 = (1 - ca[j-1])/(M - ca[j-1]);
        a2 = (ca[j-1] - 1./3.) / (ca[j-1] - m);
        theta = minimal(a1,a2,1);

        for (l=0; l<=k; l++)
        {
            mod[j-1][l] = theta * un[j-1][l];
        }
        mod[j-1][0] = mod[j-1][0] + (1-theta) * ca[j-1];
        
    }
    delete[] ca;
    return mod;
}

double flux(double ul, double ur)
{
    double ans, alpha;
    int i, l;
    if ( ul <= ur)
    {
        //ans = f(ul);  //Godnov
        alpha = ur; //Lax-Friedrichs
    }
    else{
        //ans = f(ul);
        alpha = ul;
    }

    ans = f(ul) + f(ur) - alpha * (ur - ul);
    ans = ans * 0.5;

    return ans;
}

double** L(double** ut)
{
    int i, j, l, m, p, q;
    double ul, ur;
    double** ans = new double* [n];
    for (i=0; i<n; i++)
    {
        ans[i] = new double [k+1];
    }

    double A[3] = {1,3,5};
    double B[3][3][3] = {{{0,0,0},{0,0,0},{0,0,0}}, {{2,0,0},{0,2.0/3.0,0},{0,0,2.0/5.0}}, {{0,2,0},{2,0,4.0/5.0},{0,4.0/5.0,0}}};


    for (j=1; j<=n; j++)
    {
        for (m=0; m<=k; m++)
        {
            ans[j-1][m] = 0;

            for (p=0; p<=k; p++)
            {
                for (q=0; q<=k; q++)
                {
                    ans[j-1][m] = ans[j-1][m] + ut[j-1][p] * B[m][p][q] * ut[j-1][q];
                }
            }
            ans[j-1][m] = ans[j-1][m] / 2;

            //计算第一个数值通量
            {
                ul = 0;
                ur = 0;
                q = j;
                if (q == n)
                {
                    q = 0;
                }
                for (l=0; l<=k; l++)
                {
                    ul = ul + ut[j-1][l] * phi(l,1);
                    ur = ur + ut[q][l] * phi(l,-1);
                }
                ans[j-1][m] = ans[j-1][m] - flux(ul,ur) * phi(m,1);

            }
            
            //计算第二个数值通量
            {
                ul = 0;
                ur = 0;

                p = j-2;
                if ( p == -1)
                {
                    p = n-1;
                }
                for (l=0; l<=k; l++)
                {
                    ul = ul + ut[p][l] * phi(l,1);
                    ur = ur + ut[j-1][l] * phi(l,-1);
                }
                ans[j-1][m] = ans[j-1][m] + flux(ul,ur) * phi(m,-1);
            }

            ans[j-1][m] = ans[j-1][m] * A[m] / h;
        }
    }

    return ans;
}

double** RK22(double** un)
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
    u0 = modify(u0);

    u1 = L(u0);
    for (j=1; j<=n; j++)
    {
        for (l=0; l<=k; l++)
        {
            u1[j-1][l] = u1[j-1][l] * dt  + u0[j-1][l];
        }
    }
    u1 = modify(u1);
        

    u2 = L(u1);
    for (j=1; j<=n;j++)
    {
        for (l=0; l<=k; l++)
        {
            u2[j-1][l] = u0[j-1][l] * 0.5 + u1[j-1][l] * 0.5 + u2[j-1][l] * 0.5 * dt;
            ans[j-1][l] = u2[j-1][l];
        }
    }
    ans = modify(ans);
    

    delete[] u0;
    delete[] u1;
    delete[] u2;

    return ans;
}

double** RK33(double** un)
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
    u0 = modify(u0);

    u1 = L(u0);
    for (j=1; j<=n; j++)
    {
        for (l=0; l<=k; l++)
        {
            u1[j-1][l] = u1[j-1][l] * dt  + u0[j-1][l];
        }
    }
    u1 = modify(u1);
        

    u2 = L(u1);
    for (j=1; j<=n;j++)
    {
        for (l=0; l<=k; l++)
        {
            u2[j-1][l] = u0[j-1][l] * 3 / 4 + u1[j-1][l] /4 + u2[j-1][l] * dt / 4;
        }
    }
    u2 = modify(u2);

    u3 = L(u2);
    for (j=1; j<=n; j++)
    {
        for (l=0; l<=k; l++)
        {
            u3[j-1][l] = u0[j-1][l] / 3 + 2 * u2[j-1][l] / 3 + 2 * dt * u3[j-1][l] /3;
            ans[j-1][l] = u3[j-1][l];
        }
    }
    ans = modify(ans);
    

    delete[] u0;
    delete[] u1;
    delete[] u2;
    delete[] u3;

    return ans;
}

double** RK(int k, double** un)
{
    double** u2;
    if ( k == 1 )
    {
        u2 = RK22(un);
    }
    else if (k == 2)
    {
        u2 = RK33(un);
    }
    else{
        ;
    }

    return u2;
}