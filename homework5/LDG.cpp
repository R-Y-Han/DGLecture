#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;


//************部分参数**************//
const int n = 320;   //划分单元个数
const int k = 2;    //多项式最高次项次数

const double pi = 3.1415926;
const double h = 2 * pi / n;    //空间步长
double dt = 1e-5 ;    //时间步长
double p[n+1];  //节点位置，p_j = j * h = x_{j+1/2}, j =0,1, ..., n

const double lobattopoint[5] = {-1.0, -0.6546536707079771, 0, 0.6546536707079771, 1.0};
const double lobattoco[5] = {0.1, 0.5444444444444444, 0.7111111111111111, 0.5444444444444444, 0.1};
//*************参数定义完毕*********//

//*************函数声明*************//
double u_0(double x);   //初值
double u_exact(double x, double t);   //真解
double phi(int l, double x);    //参考单元基函数
double phix(int l, double x);   //基函数一阶导
double** initial();   //计算初始时刻u_j
double fluxu(double ul, double ur);  //数值通量计算
double fluxq(double ql, double qr);
void getq(double ut);
double** L(double** ut);   //用于计算RK的函数，u_t= F(u)
double** RK22(double** un, double dtn);    //2步二阶RK
double** RK33(double** un, double dtn);
double** RK(int k, double** un, double dtn);
//************声明完毕**************//


int main()
{
    int i, j, l;
    double t, temp1, temp2, norm1, norm2, xi;
    double T = 1;
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
    
    //计算误差
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

    //开始画图
    {
        const char* fn = "DGLecture\\homework5\\LDG.plt";
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

//**********函数定义***********//

double u_0(double x)
{
    return sin(x);
}

double u_exact(double x, double t)
{
    double ans;
    ans = exp(-t) * sin(x);
    return ans;
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

double fluxu(double ul, double ur)
{
    double ans;
    //average flux
    //ans = 0.5 * (ul + ur);

    //alternating flux
    ans = ul;

    return ans;
}

double fluxq(double ql, double qr)
{
    double ans;
    //average flux
    //ans = 0.5 * (ql + qr);

    //alternating flux
    ans = qr;

    return ans;
}

double** getq(double** ut)
{
    int j, l, m, r;
    double ul, ur;
    double** q = new double* [n];
    for (j=0; j<n; j++)
    {
        q[j] = new double [k+1];
    }
    
    double A[3] = {1,3,5};
    double B[3][3] = {{0,0,0},{2,0,0},{0,2,0}};

    for ( j=1; j<=n; j++)
    {
        for (m=0; m<=k; m++)
        {
            q[j-1][m] = 0;

            for (l=0; l<=k; l++)
            {
                q[j-1][m] = q[j-1][m] - B[m][l] * ut[j-1][l];
            }
            //第一个数值通量
            {
                ul = 0;
                ur = 0;

                r = j;
                if (r == n)
                {
                    r = 0;
                }
                for (l=0; l<=k; l++)
                {
                    ul = ul + ut[j-1][l] * phi(l,1);
                    ur = ur + ut[r][l] * phi(l,-1);
                }

                q[j-1][m] = q[j-1][m] + fluxu(ul, ur) * phi(m,1);
            }
            //第二个数值通量
            {
                ul = 0;
                ur = 0;

                r = j-2;
                if (r == -1)
                {
                    r = n-1;
                }
                for (l=0; l<=k; l++)
                {
                    ul = ul + ut[r][l] * phi(l,1);
                    ur = ur + ut[j-1][l] * phi(l,-1);
                }

                q[j-1][m] = q[j-1][m] - fluxu(ul, ur) * phi(m,-1);
            }
            q[j-1][m] = q[j-1][m] * A[m] / h;
        }
    }

    return q;
}

double** L(double** ut)
{
    int i, j, l, m, p;
    double ql, qr;
    double** ans = new double* [n];
    double** q = new double* [n];
    for (i=0; i<n; i++)
    {
        ans[i] = new double [k+1];
        q[i] = new double [k+1];
    }
    q = getq(ut);

    double A[3] = {1,3,5};
    double B[3][3] = {{0,0,0},{2,0,0},{0,2,0}};

    for (j=1; j<=n; j++)
    {
        for (m=0; m<=k; m++)
        {
            ans[j-1][m] = 0;

            for (l=0; l<=k; l++)
            {
                ans[j-1][m] = ans[j-1][m] - B[m][l] * q[j-1][l];
            }

            //第一个数值通量
            {
                p = j;
                if (p == n)
                {
                    p = 0;
                }

                ql = 0;
                qr = 0;

                for (l=0; l<=k; l++)
                {
                    ql = ql + q[j-1][l] * phi(l,1);
                    qr = qr + q[p][l] * phi(l,-1);
                }

                ans[j-1][m] = ans[j-1][m] + fluxq(ql,qr) * phi(m,1);
            }
            //第二个数值通量
            {
                p = j-2;
                if (p == -1)
                {
                    p = n-1;
                }

                ql = 0;
                qr = 0;

                for (l=0; l<=k; l++)
                {
                    ql = ql + q[p][l] * phi(l,1);
                    qr = qr + q[j-1][l] * phi(l,-1);
                }

                ans[j-1][m] = ans[j-1][m] - fluxq(ql,qr) * phi(m,-1);
            }
            ans[j-1][m] = ans[j-1][m] * A[m] / h;
        }
        
    }
    for (i=0; i<n; i++)
    {
        delete[] q[i];
    }
    delete[] q;
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