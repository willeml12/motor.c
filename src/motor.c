#include "motor.h"
#include <stdio.h>
#include <math.h>
#include "glfem.h"
#include"fem.h"

/*******************************************************
Code : motor.c contient les fonctions à réaliser dans le
	cadre du projet d'élément fini 2020-2021 ainsi que
	les sous-fonctions liée.
Auteur : Laureline Willem, 21381700, basé sur le canevas
	fourni par l'équipe LEPL1110 et les devoirs réalisé
	précédements.
Date : 08/05/2021
********************************************************/

void motorAdaptMesh(motor* theMotor, double delta)
{
	motorMesh* theMesh = theMotor->mesh;

	double x, y;
	for (int i = 0; i < theMesh->nNode; ++i) {
		if (theMotor->movingNodes[i] == 1) {
			x = theMesh->X[i] * cos(delta) - theMesh->Y[i] * sin(delta);
			y = theMesh->X[i] * sin(delta) + theMesh->Y[i] * cos(delta);
			theMesh->X[i] = x;
			theMesh->Y[i] = y;
		}
	}
	theMotor->theta += delta;
}

double motorComputeCouple(motor* theMotor)
{
	return 0.0;

}

void motorComputeCurrent(motor* theMotor)
{
	return;
}

/*******************************************************
IMPORT ET EDITION DEVOIR 4
Import des fonctions du devoir 4 qui résoud l'équation
de Poisson. Des portions du code sont donc rédigée par
l'équipe didactique.
De petit édition on été réalisée pour pouvoir appliquer
celui-ci à la nouvelle structure.
Ce fragment de code sert de base au calcul du potentiel
magnetique.
********************************************************/

typedef struct {
	femMesh* mesh;
	femEdges* edges;
	femDiscrete* space;
	femIntegration* rule;
	femFullSystem* system;
} femPoissonProblem;

# ifndef NOPOISSONCREATE

femPoissonProblem* femPoissonCreate(motor* theMotor)
{
	femPoissonProblem* theProblem = malloc(sizeof(femPoissonProblem));
	theProblem->mesh = (femMesh*)theMotor->mesh;
	theProblem->edges = femEdgesCreate(theProblem->mesh);
	theProblem->space = femDiscreteCreate(3, FEM_TRIANGLE);
	theProblem->rule = femIntegrationCreate(3, FEM_TRIANGLE);
	theProblem->system = femFullSystemCreate(theProblem->mesh->nNode);
	return theProblem;
}

# endif

# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem* theProblem)
{
	femFullSystemFree(theProblem->system);
	femIntegrationFree(theProblem->rule);
	femDiscreteFree(theProblem->space);
	femEdgesFree(theProblem->edges);
	free(theProblem);
}


# endif

# ifndef NOPOISSONSOLVE

void femPoissonSolve(femPoissonProblem* theProblem, motor* theMotor)
{
	int n = theProblem->space->n;
	int nLocal = theProblem->mesh->nLocalNode;
	double** A = theProblem->system->A;
	double* B = theProblem->system->B;
	// Si on ajoute des règles pour des nombres > 4, alors il suffit de remplacer les 4 par la
	// nouvelle valeur maximale.
	if (n > 4) Error("Unexpected discrete space size !");
	double x[4];
	double y[4];
	double phi[4];
	double dphidxsi[4];
	double dphideta[4];
	double dphidx[4];
	double dphidy[4];
	int map[4];
	for (int elem = 0; elem < theProblem->mesh->nElem; elem++) {
		femMeshLocal(theProblem->mesh, elem, map, x, y);
		for (int currentPoint = 0; currentPoint < n; currentPoint++) {
			double xsi = theProblem->rule->xsi[currentPoint];
			double eta = theProblem->rule->eta[currentPoint];
			double weight = theProblem->rule->weight[currentPoint];
			femDiscretePhi2(theProblem->space, xsi, eta, phi);
			femDiscreteDphi2(theProblem->space, xsi, eta, dphidxsi, dphideta);
			double xloc = 0;
			double yloc = 0;
			double dxdxsi = 0;
			double dxdeta = 0;
			double dydxsi = 0;
			double dydeta = 0;
			for (int i = 0; i < n; i++) {
				xloc += x[i] * phi[i];
				yloc += y[i] * phi[i];
				dxdxsi += x[i] * dphidxsi[i];
				dxdeta += x[i] * dphideta[i];
				dydxsi += y[i] * dphidxsi[i];
				dydeta += y[i] * dphideta[i];
			}
			double J = dxdxsi * dydeta - dxdeta * dydxsi;
			J = fabs(J);
			for (int i = 0; i < n; i++) {
				dphidx[i] = dphidxsi[i] * dydeta - dphideta[i] * dydxsi;
				dphidx[i] = dphidx[i] / J;
				dphidy[i] = dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta;
				dphidy[i] = dphidy[i] / J;
			}
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					A[map[i]][map[j]] += (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * J * weight / theMotor->mu[theMotor->mesh->domain[elem]];
				}
			}
			for (int i = 0; i < n; i++) {
				B[map[i]] += weight * J * phi[i] * theMotor->js[theMotor->mesh->domain[elem]];
			}
		}
	}
	for (int ed = 0; ed < theProblem->edges->nEdge; ed++) {
		if (theProblem->edges->edges[ed].elem[1] < 0) {
			femFullSystemConstrain(theProblem->system, theProblem->edges->edges[ed].node[0], 0);
			femFullSystemConstrain(theProblem->system, theProblem->edges->edges[ed].node[1], 0);
		}
	}
	femFullSystemEliminate(theProblem->system);
}

# endif

/*******************************************************
FIN IMPORT ET EDITION DEVOIR 4
********************************************************/

/*******************************************************
Résoud l'équation de Poisson pour le moteur theMotor, et
place le résultat dans theMotor->a.
********************************************************/
void motorComputeMagneticPotential(motor* theMotor)
{
	femPoissonProblem* myProblem = femPoissonCreate(theMotor);
	femPoissonSolve(myProblem, theMotor);
	memcpy(theMotor->a, myProblem->system->B, sizeof(double) * theMotor->mesh->nNode);
	femPoissonFree(myProblem);
	return;
}

void motorFree(motor* theMotor) {

}