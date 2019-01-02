#include <bits/stdc++.h>
using namespace std;
const int M = 1e4 + 9;


int main (){

    #ifndef ONLINE_JUDGE
        freopen("i.txt", "r", stdin);
    #endif
    int a, b, i=1;
    while(scanf("%d %d", &a, &b)==2 & (a!=0 && b!=0)){
        printf("Case %d: ", i++);
        for(int j=1; j<=27; j++){
            if(a <= j*b){
                printf("%d\n", j-1);
                goto lab;
            }
        }
        puts("impossible");
        lab:
            int fraya;
    }
}
