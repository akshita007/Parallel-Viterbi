#include<bits/stdc++.h>
#include<omp.h>
#include<sys/stat.h>
#include<sys/types.h>
#include<getopt.h>
#include<time.h>
#include<fcntl.h>
#include<unistd.h>
using namespace std;

#include "viterbi.h"

vector <int> seq_LTDP()
{
    for(int i=0;i<m;i++)
    V[0][i]=init[i]*B[i][observe[0]];
    vector <int> result;
    // Forward Phase
    for(int j=1;j<t;j++)
    {
        for(int i=0;i<m;i++)
        {
            double mx = INT_MIN;
            int ind_max=0;
            for(int k=0;k<m;k++)
            {
	        double temp = V[j-1][k]*A[k][i]*B[i][observe[j]];
                if(temp > mx)
                {
                    mx = temp;
                    ind_max = k;
                }
            }
            V[j][i] = mx;
            argmnt[j][i] = ind_max;
        }
    }
    // Backward Phase
    double mx = INT_MIN;
    int indx = 0;
    result.resize(t,0);
    for(int i=0;i<m;i++)
    {
		
        if(V[t-1][i]>mx)
        {
            mx = V[t-1][i];
            indx=i;
        }
    }
    int q;
    q=result[t-1]=indx;
    //cout<<q<<endl;
    for(int i=t-2;i>=0;i--)
    {
        q=result[i] = argmnt[i+1][q];
        
    }
    return result;
}

