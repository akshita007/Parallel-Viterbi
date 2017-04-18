#include<bits/stdc++.h>
using namespace std;
int main(){
    int i,j,m,n,t;
    double r;
    ofstream myfile;
    myfile.open("Large_Input");
    cin>>n>>m>>t;
    myfile<<n<<" "<<m<<" "<<t<<endl;
    for(i=0;i<m;i++){
        for(j=0;j<m;j++){
            r = rand()%99+1.0;
            r/=100;
            cout<<r<<" ";
            myfile<<r<<" ";
        }
        cout<<endl;
        myfile<<endl;
    }
    cout<<endl;
    myfile<<endl;
    for(j=0;j<m;j++){
        r = rand()%99+1.0;
        r/=100;
        cout<<r<<" ";
        myfile<<r<<" ";
    }
    cout<<endl<<endl;
    myfile<<endl<<endl;
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            r = rand()%99+1.0;
            r/=100;
            cout<<r<<" ";
            myfile<<r<<" ";
        }
        cout<<endl;
        myfile<<endl;
    }
    cout<<endl;
    myfile<<endl;
    int rr;
    for(j=0;j<t;j++){
        rr=rand()%n;
        cout<<rr<<" ";
        myfile<<rr<<" ";
    }
    cout<<endl<<endl;
    myfile<<endl;
    myfile.close();
    return 0;
}
