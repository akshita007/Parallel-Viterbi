#include<bits/stdc++.h>
#include<omp.h>
#include<sys/stat.h>
#include<sys/types.h>
#include<getopt.h>
#include<time.h>
#include<fcntl.h>
#include<unistd.h>
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

int n ,m ,t;
int P;

vector < vector <double > > B;
vector < vector <double > > A;
vector < vector <double > > V;
vector < double > init ;
vector < vector <int  > > argmnt;
vector <int> observe;

void Allocate_Memory()
{
	cin>>n>>m>>t;
	A.resize(m);
	for (int i = 0; i < m; ++i)
    A[i].resize(m);
    
    B.resize(m);
	for (int i = 0; i < m; ++i)
    B[i].resize(n);
    
    V.resize(t);
	for (int i = 0; i < t; ++i)
    V[i].resize(m);
    
    argmnt.resize(t);
	for (int i = 0; i < t; ++i)
    argmnt[i].resize(m);
    
    init.resize(m);
    
    observe.resize(t);
    
    for(int i=0;i<m;i++)
        for(int j=0;j<m;j++)
			cin>>A[i][j];

    for(int i=0;i<m;i++)
		cin>>init[i];
	
    for(int i=0;i<m;i++) 
       for(int j=0;j<n;j++)
            cin>>B[i][j];

    for(int i=0;i<t;i++)
		cin>>observe[i];
}

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




vector <int> parallel_LTDP();// prototype here


static void usage(char *argv0) {
    const char *help =
        "Usage: %s [switches]  -p num_processors\n"
        "       -n Number of processors   : No of input data points\n" ;
    fprintf(stderr, help, argv0);
    exit(-1);
}
int main( int argc, char *argv[])
{
    // Initialization
    Allocate_Memory();
    int option;
    while ((option = getopt(argc, argv,"p:")) != -1) {
        		switch (option) {
             			case 'p' : P=atoi(optarg);
						break;
             			default: usage(argv[0]);
                 		exit(EXIT_FAILURE);
        		}
    		}
   	clock_t startTime = clock();
	vector <int> LTDP_result = seq_LTDP();
	
	cout<<n<<"\t"<<m<<"\t"<<t<<"\t"<<P<<"\t";
	cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<<"\t";
	startTime=clock();
    vector <int> Parallel_result = parallel_LTDP();
    cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< endl;
    int s=Parallel_result.size();
    /* cout<<"The resultant sequence(from sequential)"<<endl;
    for(int i=0;i<s;i++)
        cout << LTDP_result[i] << " ";
   cout<<endl;
    cout<<"The resultant sequence(from rank convergence)"<<endl;
    for(int i=0;i<s;i++)
        cout << Parallel_result[i] << " ";
    cout<<endl;
    */
    return 0;
}	


// Parallel LTDP
vector<double> dot_product(vector <double> s,int j,vector<int> &pre)
{
	vector <double> temp(m);
	//cout<<j<<" "<<s.size()<<endl;;
	pre.clear();
	for(int i=0;i<m;i++)
	{
		double mx=INT_MIN;
		int indx=0;
		for(int k=0;k<m;k++)
		{
			double val =s[k]*A[k][i]* B[i][observe[j]]; ;
			if(val >mx)
			{
				mx=val;
				indx=k;
			}
		}
		temp[i] = mx;
		pre.push_back(indx);
	}
	return temp;		
}

int is_parallel(vector<double> s1,vector < double> s2)
{
	int i;
	double diff=s1[0]-s2[0];
	for(i=1;i<m;i++)
	{
		if(s1[i]-s2[i]!=diff)
		return 0;
	}
	return 1;
}

vector <int> parallel_LTDP()
{
    vector <double> S[t];
    vector <int> pred[t];
    
    bool converge = 0;
    vector <bool> conv(P);

    vector<double> s(m);
    
    for(int i=0;i<m;i++)
    s[i]=(init[i]*B[i][observe[0]]);
    
    int lower[P],upper[P];
    int count=0;
    
    for (int i = 0; i < P; ++i)
    {
         count =t/P;
         if (i < t%P)
            count++;
		 if(i==0)
		 lower[i]=0;
		 else 
		 lower[i]=upper[i-1]+1;
		 upper[i]=lower[i]+count-1;
	 }
	 lower[0]=1;
	#pragma omp parallel for firstprivate(s)
	for(int p=0;p<P;p++)
	{
		int l=lower[p];
		int r=upper[p];
		for(int j=l;j<=r;j++)
		{
			s = S[j]=dot_product(s,j,pred[j]);
		}
	}
	// fix up loop
	int iter=1;
    do
    {
        	#pragma omp parallel for private(s)
			for(int p=iter;p<P;p++)
			{
				conv[p] = false;
				int l=lower[p];
				int r=upper[p];
				s = S[l-1];
				for(int j=l;j<=r;j++)
				{
					s=dot_product(s,j,pred[j]);
					if(is_parallel( s,S[j]))
					{
						conv[p] = true;
						break;
					}
					S[j] = s;
				}
			}
		converge=1;
		for(int p=iter;p<P;p++)
           converge = converge & conv[p];
           iter++;
    }while(!converge );

    double mx = INT_MIN;
    int indx = 0;
    vector<int> result;
    result.resize(t,0);
    for(int i=0;i<m;i++)
    {
		if(pred[t-1][i]>mx)
        {
            mx = pred[t-1][i];
            indx=i;
        }
    }
    int q;
    q=result[t-1]=indx;
    for(int i=t-2;i>=0;i--)
    {
        q=result[i] =  pred[i+1][q];
    }
    return result;
}

vector<int> block_cyclic()
{
	 int P = 5;
    vector <int> stages[P];
    vector <double> S[t];
    vector <int> pred[t];
    vector<double> s(m);
    bool converge=0;
     int ptr=0;
     //0th processor will exxecute stage 0 so 0 will not be assigned specifically
     //Assigned the stages to each
     for(int i=0;i<P;i++)
     {
		 stages[ptr].push_back(i);
		 ptr=(ptr+1)%P;
	 }
	 int flag[t];
	 omp_lock_t flag[t];
	 //Variable whose value determines whether the current stage can be executed or not
	 for(int i=0;i<t;i++)
	 flag[i]=0;
	 //Set the value to 1 for all stages in the first processor
	 for(int i:stages[0])
		flag[i]=1;
	 for(int i=0;i<m;i++)
		s[i]=(init[i]*B[i][observe[0]]);
	for(int i=0;i<t;i++)
	S[i]=s;
	flag[1]=1;
	flag[0]=1;
	 #pragma omp parallel for
	 for(int p=0; p<P; p++)
		{
			for(int i: stages[p])
			{
				if(flag[i]==1)
				{
					S[i]=dot_product(S[i-1],i,pred[i]);
					flag[i]=0;
					flag[(i+1)%t]=1;					
				}
			}
		}
	vector <bool> conv(n/P);
		 //instead of p we have n/p 
	//fix up loop
	int iter=1;
	do
	{
		#pragma omp parallel for private(s)
		for(int p=0; p<P; p++)
		{
			
			for(int j=iter;j<stages[p].size();j++)
			{
				int i=stages[p][j];
				if(p==0)
				conv[i/P]=false;
				#pragma omp critical
				{
					if(flag[i]==1)
					{
						s=dot_product(S[i-1],i,pred[i]);
						flag[i]=0;
						if(is_parallel( s,S[i]))
						{
							//if convergence occurs for this then we need to set for i/P
							conv[i/P] = true;
							flag[(i/P +1 )*P] =1;
							// Set the flag for the next set of stages
						}
						else
						flag[(i+1)%t]=1;
						
						S[i] = s;
						
					}
				}
			}
		}
		converge=1;
		for(int p=iter;p<n/P;p++)
           converge = converge & conv[p];
           iter++;
    }while(!converge );
    
    //backward phase
    double mx = INT_MIN;
    int indx = 0;
    vector<int> result;
    result.resize(t,0);
    for(int i=0;i<m;i++)
    {
		if(pred[t-1][i]>mx)
        {
            mx = pred[t-1][i];
            indx=i;
        }
    }
    int q;
    q=result[t-1]=indx;
    for(int i=t-2;i>=0;i--)
    {
        q=result[i] =  pred[i+1][q];
    }
    return result;
    
}
	

vector<int> block_cyclic_x()
{
	 int P = 5;
    vector <int> stages[P];
    vector <double> S[t];
    vector <int> pred[t];
    vector<double> s(m);
    bool converge=0;
     int ptr=0;
     //0th processor will exxecute stage 0 so 0 will not be assigned specifically
     //Assigned the stages to each
     for(int i=0;i<P;i++)
     {
		 stages[ptr].push_back(i);
		 ptr=(ptr+1)%P;
	 }
	 int flag[t];
	 omp_lock_t flag[t];
	 //Variable whose value determines whether the current stage can be executed or not
	 for(int i=0;i<t;i++)
	 flag[i]=0;
	 //Set the value to 1 for all stages in the first processor
	 for(int i:stages[0])
		flag[i]=1;
	 for(int i=0;i<m;i++)
		s[i]=(init[i]*B[i][observe[0]]);
	for(int i=0;i<t;i++)
	S[i]=s;
	flag[1]=1;
	flag[0]=1;
	 #pragma omp parallel for
	 for(int p=0; p<P; p++)
		{
			for(int i: stages[p])
			{
				if(flag[i]==1)
				{
					S[i]=dot_product(S[i-1],i,pred[i]);
					flag[i]=0;
					flag[(i+1)%t]=1;					
				}
			}
		}
	vector <bool> conv(n/P);
		 //instead of p we have n/p 
	//fix up loop
	int iter=1;
	do
	{
		#pragma omp parallel for private(s)
		for(int p=0; p<P; p++)
		{
			
			for(int j=iter;j<stages[p].size();j++)
			{
				int i=stages[p][j];
				if(p==0)
				conv[i/P]=false;
				#pragma omp critical
				{
					if(flag[i]==1)
					{
						s=dot_product(S[i-1],i,pred[i]);
						flag[i]=0;
						if(is_parallel( s,S[i]))
						{
							//if convergence occurs for this then we need to set for i/P
							conv[i/P] = true;
							flag[(i/P +1 )*P] =1;
							// Set the flag for the next set of stages
						}
						else
						flag[(i+1)%t]=1;
						
						S[i] = s;
						
					}
				}
			}
		}
		converge=1;
		for(int p=iter;p<n/P;p++)
           converge = converge & conv[p];
           iter++;
    }while(!converge );
    
    //backward phase
    double mx = INT_MIN;
    int indx = 0;
    vector<int> result;
    result.resize(t,0);
    for(int i=0;i<m;i++)
    {
		if(pred[t-1][i]>mx)
        {
            mx = pred[t-1][i];
            indx=i;
        }
    }
    int q;
    q=result[t-1]=indx;
    for(int i=t-2;i>=0;i--)
    {
        q=result[i] =  pred[i+1][q];
    }
    return result;
    
}
	
