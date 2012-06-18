#include "matriz.h"

/************************************Matriz**********************************/
//constructor por defecto
matriz :: matriz()
{
    e=NULL;
	n=0; m=0;
}

// constructor general
matriz :: matriz(int l, int k)
{
	n=l;
	m=k;
    e = new double*[n];
	for(int i=0;i<n;i++){
		e[i]=new double[m];
	}
}

//Destructor
matriz :: ~matriz (void)
{
	delete e;
	n=0; m=0;
}

//Crear Matriz con sus dimenciones
void matriz :: crear(int l, int k)
{
	n=l;
	m=k;
    e = new double*[n];
	for(int i=0;i<n;i++){
		e[i]=new double[m];
	}
}

//Obtiene el tamaño de n
int matriz :: getN()
{
    return(n);
}

//Obtiene el tamaño de m
int matriz :: getM()
{
    return(m);
}

//Da valor al elemento en la pocicion i,j
void matriz :: setElemento(int i, int j, double elem)
{
    e[i][j]=elem;
    return;
}

//Obtiene el elemento en la pocicion i,j 
double matriz :: getElemento(int i, int j)
{
    return(e[i][j]);
}

//cambio de filas entre p e i
void matriz :: swap(int p, int i)
{
	double aux;
	for(int j=0; j<=n; j++){
		aux=e[p][j];
		e[p][j]=e[i][j];
		e[i][j]=aux;
	}
}

//Operacion Empleada en las eliminaciones gaucianas.
void matriz :: restaFilas(int j, int i, double m)
{
	double aux;
	
	for(int k=0; k<=n; k++){
		aux=e[j][k]-(m*e[i][k]);
		e[j][k]=aux;		
	}
	
}
//Sumatoria de los elementos de una fila multiplicados por un Xi
double matriz :: sumatoria(int i, double *x)
{
	double aux=0;
	for(int j=i; j<n; j++){
		aux+=e[i][j]*x[j];
	}
	return aux;
	
}

//Sumatoria de los elementos de una fila multiplicados por un Xi
//para su uso con Nrow
double matriz :: sumatoriaNROW(int i, int *nrow, double *x)
{
	double aux=0;
	for(int j=i; j<n; j++){
		aux+=e[nrow[i]][j]*x[j];
	}
	return aux;
	
}

//Eliminacion Gaussiana
bool matriz :: eliGauss(double *x)
{
  
	//Paso 1: eliminacion{paso 2-4}
	for (int i=0; i<n-1;i++){
		int p=0;
		for(p=i; p<n; p++){ //Paso 2: entero p mas pequeño !=0
			if(e[p][i]!=0)  break;
		}
		if(p!=i) swap(p,i); //Paso 3: intercambio filas
		for(int j=i+1; j<n;j++){ //Paso 4: Eliminacion{5-6} 
			double l=0;
			l=e[j][i]/e[i][i]; //5
			restaFilas(j,i,l); //6
		}
	}
	if(e[n-1][n-1]==0) return false;//Paso 7
	
	/*************************
		for(int i=0;i<n;i++){
		for(int j=0; j<m;j++){
			cout<< e[i][j]<<"    ";
		}
		cout<<endl;
		}
		cout<<endl;
		cout<<endl;
	/*************************/
	
	for(int i=0;i<n;i++) x[i]=0;
	
	//Paso 8 comienza sustitucion hacia atras		
	x[n-1]=e[n-1][n]/e[n-1][n-1];
	
	for(int i=n-2; i>=0; i--){
		x[i]=(e[i][n]-sumatoria(i,x))/e[i][i];//Paso 9
		
	}
			
	//for(int i=0;i<n;i++) cout<<x[i]<<" ";//Paso 10: Salida
    return true;

}

//Eliminacion Gaussiana con Pivoteo Parcial
bool matriz :: eliGaussPivotPar(double *x)
{
    int *nrow;
	nrow=new int [n];
	
	for(int i=0; i<n; i++){//Paso 1: iniciar apuntadores
		nrow[i]=i;
	}
	
	//Paso 2 {3-6}: eliminacion
	for(int i=0;i<n-1;i++){
		//paso 3
		double mayor=0;
		int p=-1;
		
		for(int k=i; k<n; k++){ //Paso 2: entero p mas pequeño !=0
			if(mayor<e[nrow[k]][i]){
				mayor=e[nrow[k]][i];
				p=k;
			}
		}
		//cout<<e[nrow[p]][i]<<endl;
		
		if(e[nrow[p]][i]==0) return false;//Paso 4
			
		if(nrow[p]!=nrow[i]){//Paso 5
			int ncopy;
			ncopy=nrow[p];
			nrow[p]=nrow[i];
			nrow[i]=ncopy;
		}
		
		for(int j=i+1; j<n;j++){//Paso 6{7-8}: eliminacion
			double l=0;
			l=e[nrow[j]][i]/e[nrow[i]][i];//Paso 7
			restaFilas(nrow[j],nrow[i],l);//Paso 8
		}
		
	}
	if(e[nrow[n-1]][n-1]==0){ return 0;}//Paso 9
	
	/*************************
		for(int i=0;i<n;i++){
		for(int j=0; j<m;j++){
			cout<< e[i][j]<<"    ";
		}
		cout<<endl;
		}
		cout<<endl;
		cout<<endl;
	/*************************/
	
	for(int i=0;i<n;i++)x[i]=0;
	
	//Paso 10 sustitucion hacia atras		
	x[n-1]=e[nrow[n-1]][n]/e[nrow[n-1]][n-1];
	
	for(int i=n-2; i>=0; i--){
		x[i]=(e[nrow[i]][n]-sumatoriaNROW(i, nrow, x))/e[nrow[i]][i];//Paso 11
		
	}
			
	//for(int i=0;i<n;i++) cout<<x[i]<<" ";//Paso 12: Salida
    return true;
}

//Eliminacion Gaussiana con Pivoteo Parcial Escalado
bool matriz :: eliGaussPivotParEsc(double *x)
{
    int *nrow;
	double mayor;
	double *S=new double[n];
    double *X=new double[n];
	nrow=new int [n];
	
	for(int i=0;i<n;i++){//paso 1
		nrow[i]=0;
	  	mayor=0;
	  	for(int j=0;j<n;j++){
       		if(mayor<fabs(e[i][j]))
            	mayor=e[i][j];
        }
        S[i]=fabs(mayor);
		if(S[i]==0) return 0;
        nrow[i]=i;
    }
	
	//Paso 2 {3-6}: eliminacion
	for(int i=0;i<n-1;i++){
		//paso 3
		mayor=0;
		int p=-1;
		
		for(int j=i;j<n;j++){
        	if(mayor<(fabs(e[nrow[j]][i])/S[nrow[j]])){
            	mayor=fabs(e[nrow[j]][i])/S[nrow[j]];
                p=j;
            }
        }
		
		if(e[nrow[p]][i]==0) return false;//Paso 4
		if(nrow[p]!=nrow[i]){//Paso 5
			int ncopy;
			ncopy=nrow[p];
			nrow[p]=nrow[i];
			nrow[i]=ncopy;
		}
		
		for(int j=i+1; j<n;j++){//Paso 6{7-8}: eliminacion
			double l=0;
			l=e[nrow[j]][i]/e[nrow[i]][i];//Paso 7
			restaFilas(nrow[j],nrow[i],l);//Paso 8
		}
		
	}
	if(e[nrow[n-1]][n-1]==0) return false;//Paso 9
	
	/*************************
		for(int i=0;i<n;i++){
		for(int j=0; j<m;j++){
			cout<< e[i][j]<<"    ";
		}
		cout<<endl;
		}
		cout<<endl;
		cout<<endl;
	/*************************/
	
	for(int i=0;i<n;i++) x[i]=0;
	
	//Paso 10 sustitucion hacia atras		
	x[n-1]=e[nrow[n-1]][n]/e[nrow[n-1]][n-1];
	
	for(int i=n-2; i>=0; i--){
		x[i]=(e[nrow[i]][n]-sumatoriaNROW(i, nrow, x))/e[nrow[i]][i];//Paso 11
		
	}
			
	//for(int i=0;i<n;i++) cout<<x[i]<<" ";//Paso 12: Salida
    return true;
	
}

//
double sumatoriaLU(int i, matriz *l, matriz *u)
{
	double res=0;
	for(int k=0;k<=i-1; k++){
		res+=l->getElemento(i,k)*u->getElemento(k,i);
	}
	return res;
}


//Factorizacion LU y resolucion Ax=b
bool matriz :: factLU(matriz *b, double *X)
{
    int i,j,k;
    matriz l(n,n),u(n,n);
    double *Y,res,res2;
    Y=new double[n*n];

    for(j=0;j<n;j++){
    	for(i=0;i<n;i++){
        	if(i==j)
            l.e[i][j]=1;
        }
    }

    u.e[0][0]=e[0][0];

    if(u.e[0][0]==0) return false;

    for(j=1;j<n;j++){
    	u.e[0][j]= e[0][j];
        l.e[j][0]= e[j][0]/u.e[0][0];
    }

    for(i=1;i<n-1;i++){
    	u.e[i][i]= e[i][i]-sumatoriaLU(i,&l,&u);
        if(u.e[i][i]==0) return false;

        for(j=i+1;j<n;j++){
        	res=0;
        	for(k=0;k<=i-1;k++){
            	res= res + (l.e[i][k]*u.e[k][j]);
            }
            u.e[i][j]=(e[i][j]-res)/l.e[i][i];
            res2=0;
            for(k=0;k<=i-1;k++){
            	res2= res2 + (l.e[j][k]*u.e[k][i]);
            }
            l.e[j][i]=(e[j][i]-res2)/u.e[i][i];
        }
	}
    res2=0;
    for(k=0;k<n-1;k++)
    	res2= res2 + (l.e[n-1][k]*u.e[k][n-1]);

    u.e[n-1][n-1]= e[n-1][n-1]-res2;

    Y[0]=b->e[0][0];
    for(i=1;i<n;i++){
    	res=0;
        for(k=0;k<=i-1;k++)
        	res = res +(l.e[i][k]*Y[k]);
        Y[i]=(b->e[i][0]-res);
    }
    X[n-1]= Y[n-1]/u.e[n-1][n-1];

    for(i=n-1;i>=0;i--){
    	res=0;
        for(k=i+1;k<n;k++)
        	res = res +(u.e[i][k]*X[k]);
        
		X[i]=(Y[i]-res)/u.e[i][i];
    }
     
	return true;
}

//Jacobi
bool matriz :: jacobi(matriz *b, double *X, int N, float tol)
{
    int i,k=0,j,a;
	bool band=false;
    double res,mayor,aux5,aux0;
    double *Xo=new double[n];
    double *aux=new double[n];
    
    for(i=0;i<n;i++)//inicializacion
    	Xo[i]=0;
	
    while (k<N && band==0){//Paso 2 {3-6}
    	for(i=0;i<n;i++){//paso 3
        	res=0;
            for(j=0;j<n;j++){//Sumatoria
				if(i!=j){
					res=res + (e[i][j]*Xo[j]);
                }
            }
			X[i]=((res*-1)+b->e[i][0])/e[i][i];
		}
        for(a=0;a<n;a++){//paso 4
        	aux[a]=X[a]-Xo[a];
        }

        mayor=0;
        for(a=0;a<n;a++){//Norma infinito
        	if(fabs(aux[a])>mayor)
            	mayor=fabs(aux[a]);
        }
        if(mayor<tol){
			band=true;
        }
        k++;//Paso 5

        for(a=0;a<n;a++)//Paso 6
			Xo[a]=X[a];
	}
    return band;
}

//Gauss-Seidel
bool matriz :: gaussSeidel(matriz *b, double *X, int N, float tol) 
{
	int i,k=0,j,a;
	bool band=false;
	double res,res2,mayor;
	double *Xo=new double[n];
	double *aux=new double[n];

    for(i=0;i<n;i++)
		Xo[i]=0;
	while (k<N && !band){//Paso 2
		for(i=0;i<n;i++){//Paso 3
			res=0;
			res2=0;
			for(j=0;j<=i-1;j++){//Sumatoria 1
				res=res + (e[i][j]*X[j]);
			}
			res=(-1)*res;
			for(j=i+1;j<n;j++){//Sumatoria 2
				res2=res2 + (e[i][j]*Xo[j]);
			}
			X[i]=(res-res2+b->e[i][0])/e[i][i];
		}
		for(a=0;a<n;a++){//paso 4
			aux[a]=X[a]-Xo[a];
		}
		mayor=0;
		for(a=0;a<n;a++){//norma infinito
			if(fabs(aux[a])>mayor)
				mayor=fabs(aux[a]);
		}
		if(mayor<tol){
			band=true;
		}
		k++;//Paso 5			
		
		for(a=0;a<n;a++)//Paso 6
			Xo[a]=X[a];
		
	}
	return band;
}

//Successive over-relaxation
bool matriz :: sor(double *X, int N, float tol, float w)
{
	int i, j, k=0;
   	double* Xo=new double[n];
   	double sum;
   	double sep=1.0;
	
	for(j=0;j<n;++j)
		X[j]=Xo[j]=0;
	
	while(k++<N && sep>tol){//paso 2
		for(i=0;i<n;i++){//Paso 3
			sum=e[i][n];
		 	for(j=0;j<n;j++)//Sumatoria 1
		  		if(j!=i)
		   			sum-=e[i][j]*X[j];
		   	X[i]=sum/e[i][i];
	   	}
	   	for(j=0;j<n;j++)//Sumatoria 2
			X[j]=w*X[j]+(1-w)*Xo[j];
		
		sep=0;
	   	for(j=0;j<n;j++)//Paso 4
			sep+=fabs(X[j]-Xo[j]);
		for(j=0;j<n;j++)//Paso 6
			Xo[j]=X[j];
	}
	if(sep<tol){
		return true;
	}
	else{
		return false;
	}
}

//Diferencias Divididas interpolantes de Newton
void matriz :: diferenciasDivididas(matriz *Xi)
{
	/* Paso 1 */
      for (int i=1; i<n; i++)
         for (int j=1; j<=i; j++)
            e[i][j] = (e[i][j-1] - e[i-1][j-1]) / (Xi->e[i][0] - Xi->e[i-j][0]);
}

//Interpolacion iterada de Neville
void matriz :: neville(matriz *Xi, double x)
{
    /* Paso 1 */
	for(int i=1; i<n; i++){		
        for(int j=1; j<=i; j++){			
            double aux=x-Xi->e[i-j][0], aux2=x-Xi->e[i][0];			
            if(j>0&&i>0){
				double W;
                W=(aux * e[i][j-1])-(aux2 * e[i-1][j-1]);				
				e[i][j]=W/(aux2 - aux);
				if(e[i][j]<0)e[i][j]*=-1;								
            }					
        }
    }    
}

//Interpolacion de Hermite
void matriz :: hermite(matriz *Mat)
{
	
	double *z=new double[n];	
	int i, j, f;
	
	for(i=0; i<n; i++){
		for(j=0; j<n; j++)	
			e[i][j]=0;
	}
	//Paso 1
	for(i=0; i<Mat->getN(); i++){
	
		f=2*i;
		//Paso2
		z[f]=z[f+1]=Mat->e[i][0];	
		e[f][0]=e[f+1][0]=Mat->e[i][1];	
		e[f+1][1]=Mat->e[i][2];
		//Paso 3
		if(i!=0)	
			e[f][1]=(e[f][0]-e[f-1][0])/(z[f]-z[f-1]);	
	}
	//Paso 4
	for(i=2; i<n; i++){
		for(j=2; j<=i; j++)	
			e[i][j]= (e[i][j-1]-e[i-1][j-1])/(z[i]-z[i-j]);
	}	
}
      
//Interpolacion de Lagrange
double matriz :: lagrange(matriz *f, double x)
{
	double sum=0, prod=0;
	for(int k=0; k<n; k++){//Sumatoria
		prod=1;
		for(int i=0; i<n; i++){//Productoria
			if(k!=i){
				prod*=(x-e[i][0])/(e[k][0]-e[i][0]);				
			}			
		}
		sum+=f->e[k][0]*prod;
	}
	return(sum);
}

// operador friend << sobrecargado
ostream& operator << (ostream& co, const matriz &mat){
    co << mat.n <<" "<< mat.m<<endl;
	for(int i=0;i<mat.n;i++){
		for(int j=0; j<mat.m;j++){
			co<< mat.e[i][j]<<" ";
		}
		co<<endl;
	}
    return co;
}
