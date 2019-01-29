/**
 * Author: 
 * Description: polynomial exp/log
 */
const int mo=998244353,g=3;
 int qp(int a,int b){
   int ans=1;
   do{if(b&1)ans=1ll*ans*a%mo;a=1ll*a*a%mo;}while(b>>=1);
   return ans;
 }
 int w[2][N+1];
 void dft(int *a,int n,bool v){//0 <= a[i] < mo
   int j=0;
   rep(i,0,n){
     if(i>j)swap(a[i],a[j]);
     for (int v=n>>1;(j^=v)<v;v>>=1);
   }
   for (int i=2;i<=n;i<<=1)
     for (int j=0,s=n/i;j<n;j+=i)
       rep(v,0,i>>1){
         int t=1ll*a[j+v+(i>>1)]*w[v][s*v]%mo;
         a[j+v+(i>>1)]=(a[j+v]+mo-t)%mo;
         a[j+v]=(a[j+v]+t)%mo;
} if(v){
     int y=qp(n,mo-2);
     rep(i,0,n)a[i]=1ll*a[i]*y%mo;
  }
}

void init(int n){
   int ww=qp(g,(mo-1)/n);
   w[0][0]=1;
   rep(i,1,n+1)w[0][i]=1ll*w[0][i-1]*ww%mo;
   rep(i,0,n+1)w[1][i]=w[0][n-i];
 }
 void mul(int *a,int *b,int n){
   static int x[N];
   rep(i,0,2*n)x[i]=b[i];
   init(2*n);
   dft(x,2*n,0);
   dft(a,2*n,0);
   rep(i,0,2*n)a[i]=1ll*x[i]*a[i]%mo;
   dft(a,2*n,1);
   rep(i,n,2*n)a[i]=0;
 }

 void inv(int *a,int n,int *b){
   static int x[N];
   b[0]=qp(a[0],mo-2);b[1]=0;
   for (int m=2;m<=n;m<<=1){
     rep(i,0,m) x[i]=a[i],x[i+m]=b[i+m]=0;
     init(2*m);
     dft(x,2*m,0);
     dft(b,2*m,0);
     rep(i,0,2*m) b[i]=1ll*b[i]*(2-1ll*x[i]*b[i]%mo+mo)%mo;
     dft(b,2*m,1);
     rep(i,m,2*m)b[i]=0;
  }
}

void Ln(int *a,int n){
  static int x[N];
  a[0]=1;inv(a,n,x);
  rep(i,0,n-1)a[i]=1ll*(i+1)*a[i+1]%mo;
  mul(a,x,n);
  per(i,1,n)a[i]=1ll*a[i-1]*qp(i,mo-2)%mo;
  a[0]=0;
}

void Exp(int *a,int n,int *r){
  static int x[N];
  r[0]=1;r[1]=0;
  for (int m=2;m<=n;m<<=1){
    rep(i,0,m)x[i]=r[i];
    Ln(x,m);
    rep(i,0,m)x[i]=(a[i]-x[i]+mo)%mo;
    x[0]=(x[0]+1)%mo;
    rep(i,m,2*m) r[i]=x[i]=0;
    mul(r,x,m);
    rep(i,m,2*m)r[i]=0;
  }
}
