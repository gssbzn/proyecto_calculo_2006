#ifndef __MATRIZ_H__
#define __MATRIZ_H__

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

class matriz{
	private:
		int n, m;
        double **e;
		
       
    public:
        //constructor por defecto
        matriz();
        //constructor general
        matriz(int, int); 
        //Destructor
		~matriz (void);
	
		//Crear Matriz con sus dimenciones
		void crear(int, int);
		//Obtiene el tamaño de n
		int getN();
		//Obtiene el tamaño de m
	 	int getM();
		//Da valor al elemento en la pocicion i,j
        void setElemento(int, int, double);
		//Obtiene el elemento en la pocicion i,j 
		double getElemento(int, int);
		
		/***************Operaciones Sobre matrices****************/
		//cambio de filas entre las filas dadas por el indice
		void swap(int, int);
		//Operacion Empleada en las eliminaciones gaucianas.
		void restaFilas(int, int, double);
		//Sumatoria de los elementos de una fila multiplicados por un Xi
		double sumatoria(int, double*);
		//Eliminacion Gaussiana
		bool eliGauss(double *);
		//Sumatoria de los elementos de una fila multiplicados por un Xi
		//para su uso con Nrow
		double sumatoriaNROW(int , int *, double *);
		//Eliminacion Gaussiana con Pivoteo Parcial
		bool eliGaussPivotPar(double *);
		//Eliminacion Gaussiana con Pivoteo Parcial Escalado
		bool eliGaussPivotParEsc(double *);
		//Factorizacion LU
		bool factLU(matriz *, double *);
		//Jacobi
		bool jacobi(matriz *, double *, int , float );
		//Gauss-Seidel
		bool gaussSeidel(matriz *, double*, int , float );
		//Successive over-relaxation
		bool sor(double *, int, float, float);
		//Diferencias Divididas
		void diferenciasDivididas(matriz * );
		//neville
		void neville(matriz * , double );
		//hermite
		void hermite(matriz *);
		//lagrange
		double lagrange(matriz *, double );
		
		//Operador << Sobrecarga
		friend ostream& operator<< (ostream&, const matriz&);
		
};

#endif
