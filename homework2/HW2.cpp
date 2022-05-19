#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

//先定义一些全局的变量
const int n = 100;   //划分单元个数
const int k = 2;    //多项式最高次项次数
const double h = (double) 1/n;    //空间步长
double a = 0.2;
const double dt = (double) h * a;    //时间步长
const double pi = 3.1415926;

double p[n+1];  //节点位置，p_j = j * h = x_{j+1/2}, j =0,1, ..., n
//存储不变的系数矩阵

const double Legendre[5] = {-0.9061798459,-0.5384693101,0.,0.5384693101,0.9061798459};
const double LegendreCo[5] = {0.2369268851,0.4786286705,0.5688888889,0.4786286705,0.2369268851};

//*************函数声明*************//
double u_exact(double y);   //真解
double phi(int l, double y);    //参考单元基函数
double** initial();   //计算初始时刻u_j
double** F(double** ui);   //用于计算RK的函数，u_t= F(u)
double** RK22(double** un);    //2步二阶RK
double** RK33(double** un);
double** RK(int k, double** un);
//************声明完毕**************//

int main()
{
    int i, j, l;
    double t, temp1, temp2, norm1, norm2;
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
    while(t<2 - 1e-10)
    {
/*
for (i=0; i<n;i++)
{
    for (j=0; j<=k;j++)
    {
        cout<<u1[i][j]<<"\t";
    }
}
cout<<endl;//*/

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
        temp1 = 0;
        temp2 = 0;
        for (l=0; l<=k; l++)
        {
            temp1 = temp1 + u1[j-1][l] * phi(l,-1);
            temp2 = temp2 + u1[j-1][l] * phi(l,1);
        }
        temp1 = u_exact(p[j-1]+ 1e-10) - temp1;
        temp2 = u_exact(p[j]-1e-10) - temp2;

        if (abs(temp1) > norm2)
        {
            norm2 = abs(temp1);
        }
        if (abs(temp2) > norm2)
        {
            norm2 = abs(temp2);
        }

        temp1 = temp1 * temp1;
        temp2 = temp2 * temp2;

        norm1 = norm1 + (temp1 + temp2)*h/2;

    }
    norm1 = sqrt(norm1);
cout<<norm1<<endl<<norm2<<endl;//*/

    {
        const char* fn = "DGLecture\\homework2\\output2.plt";
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
                temp1 = temp1 + u1[j-1][l] * phi(l,-1);
                temp2 = temp2 + u1[j-1][l] * phi(l,1);
            }
            f<<"\t"<<p[j-1]<<"\t"<<temp1<<"\t"<<u_exact(p[j-1])<<endl;
            f<<"\t"<<p[j]<<"\t"<<temp2<<"\t"<<u_exact(p[j])<<endl;
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

double u_exact(double y)
{
    //return sin(2*pi*y);

    //*
    if (y>=0.25 && y<= 0.75)
    {
        return 1;
    }
    else{
        return 0;
    }//*/
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
                temp = h * (Legendre[l] + 1)/2 + p[j-1];
                ans = ans + LegendreCo[l] * u_exact(temp) * phi(m,Legendre[l]);
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

double** F(double** ui)
{
    double** ans = new double* [n];
    double ul;
    int i,j,l,m;
    double C[3] = {1,-3,5};

    double B[3][3] = {{-1,-1,-1},{3,-3,-3},{-5,5,-5}};

    for (i=0; i<n; i++)
    {
        ans[i] = new double [k+1];
    }


    for (j=1; j<=n; j++)
    {
        //先计算ul
        ul = 0;
        m = j-2;
        if (m == -1)
        {
            m = n-1;
        }
        for (l=0; l<=k; l++)
        {
            ul = ul + ui[m][l] * phi(l,1);
        }

        for (l=0; l<=k;l++)
        {
            ans[j-1][l] = 0;
            for (m=0; m<=k;m++)
            {
                ans[j-1][l] = ans[j-1][l] + B[l][m] * ui[j-1][m];
            }
            ans[j-1][l] = ans[j-1][l] + ul * C[l];
            ans[j-1][l] = ans[j-1][l] / h;
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

    u1 = F(u0);
    for (j=1; j<=n; j++)
    {
        for (l=0; l<=k; l++)
        {
            u1[j-1][l] = u1[j-1][l] * dt  + u0[j-1][l];
        }
    }
        

    u2 = F(u1);
    for (j=1; j<=n;j++)
    {
        for (l=0; l<=k; l++)
        {
            u2[j-1][l] = u0[j-1][l] * 0.5 + u1[j-1][l] * 0.5 + u2[j-1][l] * 0.5 * dt;
            ans[j-1][l] = u2[j-1][l];
        }
    }
    

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

    u1 = F(u0);
    for (j=1; j<=n; j++)
    {
        for (l=0; l<=k; l++)
        {
            u1[j-1][l] = u1[j-1][l] * dt  + u0[j-1][l];
        }
    }
        

    u2 = F(u1);
    for (j=1; j<=n;j++)
    {
        for (l=0; l<=k; l++)
        {
            u2[j-1][l] = u0[j-1][l] * 3 / 4 + u1[j-1][l] /4 + u2[j-1][l] * dt / 4;
        }
    }

    u3 = F(u2);
    for (j=1; j<=n; j++)
    {
        for (l=0; l<=k; l++)
        {
            u3[j-1][l] = u0[j-1][l] / 3 + 2 * u2[j-1][l] / 3 + 2 * dt * u3[j-1][l] /3;
            ans[j-1][l] = u3[j-1][l];
        }
    }
    

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