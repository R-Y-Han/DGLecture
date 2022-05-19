#include <iostream>
#include <cmath>

using namespace std;

int d = 5;  //网格划分从10份到320份
double pi = 3.14159;

double u_exact(double x)    //真解
{
    return sin(2 * pi * x);
}

double f (double x) //右端项
{
    return 2 * pi * cos(2 * pi * x);
}

double phi(int k, double x) //参考单元I上的基函数
{
    if (k == 0)
    {
        return 1;
    }
    else if (k == 1)
    {
        return x;
    }
    else if (k == 2)
    {
        return x * x;
    }
    else{
        return 0;
    }
}

int main()
{
    int t=0;
    int n = 10; //从剖分10份开始
    int k = 2;

    cout<<"k="<<k<<endl;

    //double A_inv[1][1] = {{1}};
    //double A_inv[2][2] = {{1, -2}, {0, 2}}; //k=1
    double A_inv[3][3] = {{1,-6,6}, {0,18, -24}, {0, -12, 18}};    //线性方程组的系数逆矩阵，由于简单直接计算出

    int i,j,l;

    do{
        cout<<"剖分个数:"<<n<<endl; 
        double h = 1 /(double) n;   //步长，取均匀剖分
        cout<<"步长h:"<<h<<endl;
        
        double* x;  //记录剖分节点
        x = new double [ n + 1 ];
        for (j=0; j< n+1; j++)
        {
            x[j] = j * h;   //x_j 是第(j+1)个网格的左端点
        }

        double u_h[n][k+1];   //u_h(j,l)表示第j+1个网格上多项式的系数，注意此时j是从0开始的
        double B[k+1];  //存储右端项

        //从左到右计算
        for (j = 0; j< n; j++)
        {   
            
            double ul;  //ul是第j次运算时用到的左端量
            if (j==0)
            {
                ul = u_exact(0);
            }
            else{
                ul = 0;
                for (l=0; l<=k; l++)
                {
                    ul = ul + u_h[j-1][l] * phi(l, 1);
                }
            }

            //给B赋值
            for (l=0; l<=k; l++)
            {
                B[l] = 0;
                int ntemp = 100;  //计算积分用的参数
////////////////////////ntemp取不同值时在k=3的情况下得到的结果不同，但是B又是对的/////////////////////
                double htemp = 1/ (double) ntemp;
            
                for (i=0; i<ntemp; i++)
                {
                    double xi1 = i * htemp;
                    double xi2 = ( i + 1 ) * htemp;

                    B[l] = B[l] + ( f(xi1*h+x[j]) * phi(l, xi1) + f(xi2*h+x[j]) * phi(l, xi2) ) * htemp * 0.5;
                }   //积分
                B[l] = B[l] * h + ul * phi(l,0);
                //cout<<B[l]<<endl;
            }

            for (i = 0; i<=k ; i++)
            {
                u_h[j][i] = 0;
                for (l=0; l<=k; l++)
                {
                    u_h[j][i] = u_h[j][i] + A_inv[i][l] * B[l];
                }

                
                //cout<<u_h[j][i]<<"\t";
            }
            //cout<<endl;
        }
        
        //计算L2误差
        double err=0, u1, u2;   //err为L^2范数误差，u1, u2 为数值积分的两个节点上的量

        for (j=1; j <= n; j++)  //在每个I_j上积分，注意I_j从I_1开始，而程序j从0开始数
        {
            u1 = 0;
            u2 = 0;

            for (l = 0; l <= k; l++)
            {
                u1 = u1 + u_h[j-1][l] * phi(l, 0);
                u2 = u2 + u_h[j-1][l] * phi(l, 1);    //I_j上的基函数换到了参考单元I上的基函数
            }

            cout<<x[j-1]<<"\t"<<u1<<"\t"<<x[j]<<"\t"<<u2<<endl;

            u1 = u_exact(x[j-1]) - u1;
            u2 = u_exact(x[j]) - u2;


            err = err + (u2 * u2 + u1 * u1 ) * h * 0.5;
        }
        err = sqrt(err);

        n = n * 2;
        t++;
        cout<<"L^2 norm = "<<err<<endl;
        cout<<endl<<endl;

    }while(t<=d);
    system("pause");
}