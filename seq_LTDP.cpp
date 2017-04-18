#include<bits/stdc++.h>
#include<omp.h>
using namespace std;

/* contants defined
* n : size of sequence of observations
* m : number of hidden states
* t : size of the current observation
* A : Transmission matrix assumed to be time invariant (m*m)
* V : The Viterbi output matrix of size (m*t)
* B : Emission probablilty (m*n)
* init : Initial probability matrix(m)
*/


const int n = 40, m = 30, t = 10;

double B[m][n], A[m][m], V[m][t];
double init[m];

int argmnt[m][t]; // to store precedence product

vector <int> S(n); // to store solution vector
vector <int> pred[n];
vector <int> observe(t);

vector <int> seq_LTDP()
{
    for(int i=0;i<m;i++)
    V[i][0]=init[i]*B[i][observe[0]];
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
	        double temp = V[k][j-1]*A[k][i]*B[i][observe[j]];
                if(temp > mx)
                {
                    mx = temp;
                    ind_max = k;
                }
            }
            V[i][j] = mx;
            argmnt[i][j] = ind_max;
        }
    }
    cout<<"The Viterbi Matrix "<<endl;
    for(int i=0;i<m;i++)
    {
		for(int j=0;j<t;j++)
		{
			cout<<V[i][j]<<" ";
		}
		cout<<endl;
	}

    // Backward Phase
    double mx = INT_MIN;
    int indx = 0;
    result.resize(t,0);
    for(int i=0;i<m;i++)
    {
		
        if(V[i][t-1]>mx)
        {
            mx = V[i][t-1];
            indx=i;
        }
        //cout<<V[i][t-1]<<" "<<mx<<" "<<indx<<endl;
    }
    int q;
    q=result[t-1]=indx;
     //cout<<q<<endl;
    for(int i=t-2;i>=0;i--)
    {
        q=result[i] =  argmnt[q][i+1];
    }
    return result;
}

vector <int> parallel_LTDP();// prototype here

int main()
{
    // Initialization
    //cin>>n>>m>>t;
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<m;j++)
        {
            /*int r = rand()%20;
            double val = r/100;
            A[i][j] = val;*/
            cin>>A[i][j];
        }
    }
    for(int i=0;i<m;i++){
		cin>>init[i];
	}
		//init[i]=((double)(rand()%20))/100;
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
        {
            /* int r = rand()%20;
            double val = r/100;
            B[i][j] = val; */
            cin>>B[i][j];
        }
    } 
    //Input the observation vector(must have numbers between 0...n-1)   
    /*for(int i=0;i<n;i++)
        S[i] = i;*/
	  for(int i=0;i<t;i++)
		cin>>observe[i];

    vector <int> LTDP_result = seq_LTDP();
    //vector <int> TLDP_parallel_result = parallel_LTDP();
    int s=LTDP_result.size();
    cout<<"The resultant hidden sequence"<<endl;
    for(int i=0;i<s;i++)
        cout << LTDP_result[i] << " ";
    cout<<endl;
    return 0;
}


// Parallel LTDP

/*vector <int> parallel_LTDP()
{
    vector <int> S[n];
    vector <bool> conv;
    bool converge = 0;

    int P = omp_get_num_threads();
    #pragma omp parallel
    {
        vector <int> s;
        int p = omp_get_thread_num();
        int l = n/P*(p-1);
        int r = n/P*(p);

        for(int i=l;i<=r;i++)
        {
            s = S[i] = dot_product(A,s);
            pred[i] = get_max(A,s);
        }
    }
    #pragma omp barrier
    //fix up loop

    do
    {
        #pragma omp parallel for
        for(int p=2;p<=P;p++)
        {
            conv[p] = false;
            int l = n/P*(p-1);
            int r = n/P*(p);
            s = s[l];
            for(int i=l;i<=r;i++)
            {
                s = dot_product(A,s);
                pred[i] = get_max(A,s);
                if(s == S[i])
                {
                    conv[p] = true;
                    break;
                }
                S[i] = s;
            }
            #pragma omp barrier
            converge = converge & conv[p];
        }
    }while(!converge);
    return Backward_Phase();
}

vector <int> Backward_Phase()
{
    vector <int> res;
    vector <bool> conv(n);
    int P = omp_get_num_threads();
    #pragma omp parallel
    {
        vector <int> s;
        int p = omp_get_thread_num();
        int l = n/P*(p-1);
        int r = n/P*(p);

        int x = 0; // local x
        for(int i=r;i>l;i--)
        {
            x = res[i] = pred[i][x];
        }
    }
    #pragma omp barrier
    do // fix up loop
    {
        #pragma omp parallel for
        for(int p=P-1;p>=1;p--)
        {
            conv[p] = false;
            int l = n/P*(p-1);
            int r = n/P*(p);
            itn x = res[r+1];
            for(int i=l;i<=r;i++)
            {
                x = pred[i][x];
                if(x == S[i])
                {
                    conv[p] = true;
                    break;
                }
                res[i] = x;
            }
            #pragma omp barrier
            converge = converge & conv[p];
        }
    }while(!converge);
}*/
