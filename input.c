#include "viterbi.h"

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
