#include <iostream>
#include <cmath>

using namespace std;

int d = 5;  //���񻮷ִ�10�ݵ�320��
double pi = 3.14159;

double u_exact(double x)    //���
{
    return sin(2 * pi * x);
}

double f (double x) //�Ҷ���
{
    return 2 * pi * cos(2 * pi * x);
}

double phi(int k, double x) //�ο���ԪI�ϵĻ�����
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
    int n = 10; //���ʷ�10�ݿ�ʼ
    int k = 2;

    cout<<"k="<<k<<endl;

    //double A_inv[1][1] = {{1}};
    //double A_inv[2][2] = {{1, -2}, {0, 2}}; //k=1
    double A_inv[3][3] = {{1,-6,6}, {0,18, -24}, {0, -12, 18}};    //���Է������ϵ����������ڼ�ֱ�Ӽ����

    int i,j,l;

    do{
        cout<<"�ʷָ���:"<<n<<endl; 
        double h = 1 /(double) n;   //������ȡ�����ʷ�
        cout<<"����h:"<<h<<endl;
        
        double* x;  //��¼�ʷֽڵ�
        x = new double [ n + 1 ];
        for (j=0; j< n+1; j++)
        {
            x[j] = j * h;   //x_j �ǵ�(j+1)���������˵�
        }

        double u_h[n][k+1];   //u_h(j,l)��ʾ��j+1�������϶���ʽ��ϵ����ע���ʱj�Ǵ�0��ʼ��
        double B[k+1];  //�洢�Ҷ���

        //�����Ҽ���
        for (j = 0; j< n; j++)
        {   
            
            double ul;  //ul�ǵ�j������ʱ�õ��������
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

            //��B��ֵ
            for (l=0; l<=k; l++)
            {
                B[l] = 0;
                int ntemp = 100;  //��������õĲ���
////////////////////////ntempȡ��ֵͬʱ��k=3������µõ��Ľ����ͬ������B���ǶԵ�/////////////////////
                double htemp = 1/ (double) ntemp;
            
                for (i=0; i<ntemp; i++)
                {
                    double xi1 = i * htemp;
                    double xi2 = ( i + 1 ) * htemp;

                    B[l] = B[l] + ( f(xi1*h+x[j]) * phi(l, xi1) + f(xi2*h+x[j]) * phi(l, xi2) ) * htemp * 0.5;
                }   //����
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
        
        //����L2���
        double err=0, u1, u2;   //errΪL^2������u1, u2 Ϊ��ֵ���ֵ������ڵ��ϵ���

        for (j=1; j <= n; j++)  //��ÿ��I_j�ϻ��֣�ע��I_j��I_1��ʼ��������j��0��ʼ��
        {
            u1 = 0;
            u2 = 0;

            for (l = 0; l <= k; l++)
            {
                u1 = u1 + u_h[j-1][l] * phi(l, 0);
                u2 = u2 + u_h[j-1][l] * phi(l, 1);    //I_j�ϵĻ����������˲ο���ԪI�ϵĻ�����
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