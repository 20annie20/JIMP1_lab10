#include <stdio.h>
#include <conio.h>

int n;
float a[100],b[100],c[100];
float x;
  
float wart(int k)            //wartosc wielomianu
{
    if (k==n){
        return b[n];
    }
    else {
        return wart(k+1)*x+b[k];
    }
}
 
float wspol(int k)           //wspolczynniki wielomianu po podzieleniu
{
    if (k==n) {
        c[k-1]=b[k]; 
        return b[k];
    }
    else {
        if (k>0) {
            c[k-1]=wspol(k+1)*x+b[k];
        } 
        return wspol(k+1)*x+b[k]; 
    }
}
 
double horner(int stopien, double val, int pochodna){
    int i,j,s;
    float tmp;
    n = stopien;
    if (n>100)
        return(1); 
    //kolejne wspolczynniki 
    for(i=n; i>=0; i--){
        printf("a%d ", i);
        scanf("%e", &a[i]);
    }
    x = val;
    s = pochodna;
    if ((s<0)||(s>n)) {
        printf("zle wpisany stopien pochodnej"); 
        getche(); 
        return(1); 
    }
 
    for (i=0; i<=n; i++) {
        b[i]=a[i];
    }
    for (j=1; j<=s; j++) {
        tmp=wspol(0);
        for (i=0; i<=n; i++) {
            b[i]=c[i];
        }
    }
 
    printf("Wartosc pochodnej znormalizowanej %d stopnia wynosi: %f",s,wart(0));
 
    getche();
    return wart(0);
}