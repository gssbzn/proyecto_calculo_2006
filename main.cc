
#include "matriz.h"
#include "matriz.cc"

using namespace std;
//Menu de la Aplicacion
int menu();
//Para cargar "matriz.in" en 2 matrices, una con los valores de x y otra con los
//f(x)
bool cargaMatriz(matriz *, matriz *);
//Para cargar "matriz.in" en 1 matriz, con los valores de x, f(x) y f'(x)
bool cargaMatriz2(matriz *);

int main()
{
	int o;
	float tol, w;
	double x;
	double *xi;
	matriz Q, Xi;
	matriz A, b;
	system("clear");
	do{
		o=menu();
		
		switch(o){
			case 1:
				cout<<"Eliminacion Gaussiana\n\n";
				cargaMatriz2(&A);
				x=new double[A.getN()];
				if(A.eliGauss(x)){
					for(int i=0; i<A.getN(); i++)
						cout<<"X["<<i<<"]= "<<xi[i]<<" ";
				}
				else
					cout<<"Error\n\n";
				free(x);
			break;
			case 2:
				cargaMatriz2(&A);
				cout<<"Eliminacion Gaussiana con Pivoteo Parcial\n\n";
				
				x=new double[A.getN()];
				if(A.eliGaussPivotPar(x)){
					for(int i=0; i<A.getN(); i++)
						cout<<"X["<<i<<"]= "<<xi[i]<<" ";
				}
				else
					cout<<"Error\n\n";
				free(x);
			break;
			case 3:
				cout<<"Eliminacion Gaussiana con Pivoteo Parcial Escalado\n\n";
				cargaMatriz2(&A);
			
				x=new double[A.getN()];
				if(A.eliGaussPivotParEsc(x)){
					for(int i=0; i<A.getN(); i++)
						cout<<"X["<<i<<"]= "<<xi[i]<<" ";
				}
				else
					cout<<"Error\n\n";
				free(x);
			break;
			case 4:
				cargaMatriz(&A, &b);
				cout<<"Factorizacion LU\n\n";
				A.factLU(&b);
			break;
			case 5:
				cargaMatriz(&A, &b);
				cout<<"Jacobi\n\n";
				cout<<"Numero de iteraciones:  ";
				cin>>N;
				cout<<"Tolerancia: ";
				cin>>tol;
			
				x=new double[A.getN()];
				A.jacobi(&b, x, N, tol);
				for(int i=0; i<A.getN(); i++)
					cout<<"X["<<i<<"]= "<<xi[i]<<" ";
				free(x);
			break;
			case 6:
				cargaMatriz(&A, &b);
				cout<<"Gauss-Seidel\n\n";
				cout<<"Numero de iteraciones: ";
				cin>>N;
				cout<<"Tolerancia: ";
				cin>>tol;
			
				x=new double[A.getN()];
				A.gaussSeidel(&b, x, N, tol);
				for(int i=0; i<A.getN(); i++)
					cout<<"X["<<i<<"]= "<<xi[i]<<" ";
				free(x);
			break;
			case 7:
				cargaMatriz2(&A);
				cout<<"SOR\n\n";
				cout<<"Numero de iteraciones: ";
				cin>>N;
				cout<<"Tolerancia: ";
				cin>>tol;
				cout<<"W: ";
				cin>>w;
				A.sor(N,tol,w);
			break;
			case 8:
				cout<<"Metodo de Diferencias Divididas\n\n";
				if(cargaMatriz(&Xi, &Q)){
					Q.diferenciasDivididas(&Xi);
					cout<<Q<<endl;
				}
				else
					cerr << "\n\n\n***ERROR*** No se ha podido encontrar el archivo\n\b";
			break;
			case 9:
				cout<<"Metodo de Neville\n\n";
				if(cargaMatriz(&Xi, &Q)){
					cout<<"Introduzca el valor de x para el cual\naproximar f(x): ";
					cin>>x;
					Q.neville(&Xi,x);
					cout<<Q<<endl;
				}
				else
					cerr << "\n\n\n***ERROR*** No se ha podido encontrar el archivo\n\b";
			break;
			case 10:
				cout<<"Metodo de Hermite\n\n";
				if(cargaMatriz2(&Xi)){
					Q.crear(2*Xi.getN(), 2*Xi.getN());
					Q.hermite(&Xi);
					cout<<Q<<endl;
				}
				else
					cerr << "\n\n\n***ERROR*** No se ha podido encontrar el archivo\n\b";
			break;
			case 11:
				cout<<"Metodo de Lagrange\n\n";
				if(cargaMatriz(&Q, &Xi)){
					cout<<"Introduzca el valor de x para el cual\naproximar f(x): ";
					cin>>x;
					cout<<"f("<<x<<")= "<<Q.lagrange(&Xi,x)<<endl;
				}
				else
					cerr << "\n\n\n***ERROR*** No se ha podido encontrar el archivo\n\b";
			break;
		
			case 0:
				cout<<endl<<"Adios"<<endl;
			break;
			default: cerr<<"opcion erronea\n\n"; break;
		}
	}while(o!=5);
	
	return EXIT_SUCCESS;
}

/******************************************************************************/

int menu()
{
	char res[20];
	cout<<"\n\nEscoja el metodo para la solucion del sistema\n\n\n";
	cout<<"1.-Eliminacion Gaussiana\n";
	cout<<"2.-Eliminacion Gaussiana con Pivoteo Parcial\n";
	cout<<"3.-Eliminacion Gaussiana con Pivoteo Parcial Escalado\n";
	cout<<"4.-Descomposicion L.U.\n";
	cout<<"5.-Jacobi\n";
	cout<<"6.-Gauss-Seidel\n";
	cout<<"7.-SOR\n";
	cout<<"8.-Metodo de Diferencias Divididas\n";
	cout<<"9.-Metodo de Neville\n";
	cout<<"10.-Metodo de Hermite\n";
	cout<<"11.-Metodo de Lagrange\n\n";

	cout<<"0.-Salir\n";
	cin>>res;
	system("clear");

	return atoi(res);
}

/******************************************************************************/

bool cargaMatriz(matriz *A, matriz *b){
	ifstream entrada;
	int n,m;
	
	//Carga del archivo de entrada matriz.in
	entrada.open("matriz.in", ios:: in);
	if (entrada==NULL)/*SI NO SE PUEDE ABRIR EL FICHERO*/
		return false;
	else {
		entrada >> n >> m;
		//Lectura A
			A->crear(n,1);
			    for(int i=0;i<A->getN();i++)
                {
					double num;
                	entrada >> num;
					A->setElemento(i,0,num);
                }
        	//Lectura b
			b->crear(n,n);
			for(int i=0;i<b->getN();i++)
            {
				double num;
                entrada >> num;
				b->setElemento(i,0,num);
          	}
		
	}
	entrada.close();
	return true;
}

/******************************************************************************/

bool cargaMatriz2(matriz *A){
	ifstream entrada;
	int n,m;
	
	//Carga del archivo de entrada matriz.in
	entrada.open("matriz.in", ios:: in);
	if (entrada==NULL)/*SI NO SE PUEDE ABRIR EL FICHERO*/
		return false;
	else {
		entrada >> n >> m;
		//Lectura A
			A->crear(n,m);
			for(int j=0;j<A->getM();j++)
            {
            	for(int i=0;i<A->getN();i++)
                {
					double num;
                	entrada >> num;
					A->setElemento(i,j,num);
                }
        	}		
	}
	entrada.close();
	return true;
}

/******************************************************************************/
