#include <bits/stdc++.h>
using namespace std;

int a[109][5], dp[109][5];

int rec (int i, int p, int n){
    if(i==n)return 0;
    if(dp[i][p]!=-1)return dp[i][p];
    int mn = INT_MAX;
    for(int j=0; j<3; j++)if(j!=p)
        mn = min(mn, a[i][j] + rec(i+1, j, n));
    return dp[i][p] = mn;
}

void solve (int t){
    memset(dp, -1, sizeof dp);
    int n;
    scanf("%d", &n);
    for(int i=0; i<n; i++)scanf("%d %d %d", &a[i][0], &a[i][1], &a[i][2]);
    printf("Case %d: %d\n", t, rec(0, 3, n));
}

int main (){
    #ifndef ONLINE_JUDGE
        freopen("i.txt", "r", stdin);
    #endif
    int t;
    scanf("%d", &t);
    for(int i=1; i<=t; i++)solve(i);
}
