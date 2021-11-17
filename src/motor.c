#include "motor.h"
#include <stdio.h>
#include <math.h>
#include "glfem.h"
#include"fem.h"

/*******************************************************
Code : motor.c contient les fonctions a realiser dans le
	cadre du projet d'element fini 2020-2021 ainsi que
	les sous-fonctions liee.
Auteur : Laureline Willem, 21381700, base sur le canevas
	fourni par l'equipe LEPL1110 et les devoirs realise
	precedements.
Date : 16/08/2021
********************************************************/

// VARIABLE GLOBALE.
femDiscrete* theGlobalSpace = NULL;

/*******************************************************
Calcule la mediane du triangle
********************************************************/
double computeMedian(double x[2], double y[2], double x3, double y3) {
	double xmedian = x[0] + x[1];
	xmedian = xmedian / 2.0;
	double ymedian = y[0] + y[1];
	ymedian = ymedian / 2.0;
	double distx = xmedian - x3;
	double disty = ymedian - y3;
	return sqrt((distx * distx) + (disty * disty));
}

/*******************************************************
Effectue la rotation du moteur et remaille la zone de
l'entrefer pour garder des triangles optimaux.
********************************************************/
void motorAdaptMesh(motor* theMotor, double delta)
{
	femMesh* theMesh = malloc(sizeof(femMesh));
	theMesh->elem = theMotor->mesh->elem;
	theMesh->X = theMotor->mesh->X;
	theMesh->Y = theMotor->mesh->Y;
	theMesh->nElem = theMotor->mesh->nElem;
	theMesh->nNode = theMotor->mesh->nNode;
	theMesh->nLocalNode = theMotor->mesh->nLocalNode;

	// Rotation
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

	// Remaillage

	// Creation d'une structure edges pour pouvoir cibler les arretes qui nous interesse.
	femEdges* theEdges = femEdgesCreate((femMesh*)theMesh);
	int starting = 0;
	for (int i = 0; i < 11; i++) {
		starting += theMotor->mesh->nElemDomain[i];
	}
	for (int currEdge = 0; currEdge < theEdges->nEdge; currEdge++) {
		int elem1 = theEdges->edges[currEdge].elem[0];
		int elem2 = theEdges->edges[currEdge].elem[1];
		int elemDomain1 = theMotor->mesh->domain[elem1];
		int elemDomain2 = 0;
		if (elem2 != -1) {
			elemDomain2 = theMotor->mesh->domain[elem2];
		}
		// remaille la frontiere avec le domaine stator_gap avec le meilleur noeud sur la frontiere avec le domaine rotor_gap.
		if ((elemDomain1 == 7 && elemDomain2 == 11) || (elemDomain1 == 11 && elemDomain2 == 7)) {
			double bestMedian = 99999;
			int bestNode;
			for (int currEdgebis = 0; currEdgebis < theEdges->nEdge; currEdgebis++) {
				int elem3 = theEdges->edges[currEdgebis].elem[0];
				int elem4 = theEdges->edges[currEdgebis].elem[1];
				int elemDomain3 = theMotor->mesh->domain[elem3];
				int elemDomain4 = 0;
				if (elem4 != -1) {
					elemDomain4 = theMotor->mesh->domain[elem4];
				}
				if ((elemDomain3 == 10 && elemDomain4 == 11) || (elemDomain3 == 11 && elemDomain4 == 10)) {
					double x[2];
					x[0] = theMesh->X[theEdges->edges[currEdge].node[0]];
					x[1] = theMesh->X[theEdges->edges[currEdge].node[1]];
					double y[2];
					y[0] = theMesh->Y[theEdges->edges[currEdge].node[0]];
					y[1] = theMesh->Y[theEdges->edges[currEdge].node[1]];
					for (int i = 0; i < 2; i++) {
						double x3 = theMesh->X[theEdges->edges[currEdgebis].node[i]];
						double y3 = theMesh->Y[theEdges->edges[currEdgebis].node[i]];
						double median = computeMedian(x, y, x3, y3);
						if (median < bestMedian) {
							bestMedian = median;
							bestNode = theEdges->edges[currEdgebis].node[i];
						}
					}
				}
			}if (elemDomain1 == 7 && elemDomain2 == 11) {
				theMesh->elem[3 * starting + 0] = theEdges->edges[currEdge].node[0];
				theMesh->elem[3 * starting + 2] = theEdges->edges[currEdge].node[1];
				theMesh->elem[3 * starting + 1] = bestNode;
			}
			else {
				theMesh->elem[3 * starting + 2] = theEdges->edges[currEdge].node[0];
				theMesh->elem[3 * starting + 0] = theEdges->edges[currEdge].node[1];
				theMesh->elem[3 * starting + 1] = bestNode;
			}
			starting++;
		}
		// remaille la frontiere avec le domaine rotor_gap avec le meilleur noeud sur la frontiere avec le domaine stator_gap.
		if ((elemDomain1 == 10 && elemDomain2 == 11) || (elemDomain1 == 11 && elemDomain2 == 10)) {
			double bestMedian = 99999;
			int bestNode;
			for (int currEdgebis = 0; currEdgebis < theEdges->nEdge; currEdgebis++) {
				int elem3 = theEdges->edges[currEdgebis].elem[0];
				int elem4 = theEdges->edges[currEdgebis].elem[1];
				int elemDomain3 = theMotor->mesh->domain[elem3];
				int elemDomain4 = 0;
				if (elem4 != -1) {
					elemDomain4 = theMotor->mesh->domain[elem4];
				}
				if ((elemDomain3 == 7 && elemDomain4 == 11) || (elemDomain3 == 11 && elemDomain4 == 7)) {
					double x[2];
					x[0] = theMesh->X[theEdges->edges[currEdge].node[0]];
					x[1] = theMesh->X[theEdges->edges[currEdge].node[1]];
					double y[2];
					y[0] = theMesh->Y[theEdges->edges[currEdge].node[0]];
					y[1] = theMesh->Y[theEdges->edges[currEdge].node[1]];
					for (int i = 0; i < 2; i++) {
						double x3 = theMesh->X[theEdges->edges[currEdgebis].node[i]];
						double y3 = theMesh->Y[theEdges->edges[currEdgebis].node[i]];
						double median = computeMedian(x, y, x3, y3);
						if (median < bestMedian) {
							bestMedian = median;
							bestNode = theEdges->edges[currEdgebis].node[i];
						}
					}
				}
			}
			if (elemDomain1 == 10 && elemDomain2 == 11) {
				theMesh->elem[3 * starting + 0] = theEdges->edges[currEdge].node[0];
				theMesh->elem[3 * starting + 2] = theEdges->edges[currEdge].node[1];
				theMesh->elem[3 * starting + 1] = bestNode;
			}
			else {
				theMesh->elem[3 * starting + 2] = theEdges->edges[currEdge].node[0];
				theMesh->elem[3 * starting + 0] = theEdges->edges[currEdge].node[1];
				theMesh->elem[3 * starting + 1] = bestNode;
			}
			starting++;
		}
	}

	// Liberation de la memoire
	femEdgesFree(theEdges);
	free(theMesh);
}

/*******************************************************
Calcule le couple du moteur et en retourne la valeur.
********************************************************/
double motorComputeCouple(motor* theMotor)
{
	double x[3];
	double y[3];
	double phi[3];
	double dphidxsi[3];
	double dphideta[3];
	double dphidx[3];
	double dphidy[3];
	int map[3];
	double r;
	double theta;
	double xsitab[3] = { 1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0 };
	double etatab[3] = { 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 };
	double dphidtheta[3];
	double dphidr[3];

	double dadr;
	double dadtheta;
	double C = 0;

	motorMesh* mesh = theMotor->mesh;

	// pour commencer le for sur le premier element de Stator_gap.
	int starting = 0;
	for (int i = 0; i < 10; i++) {
		starting += mesh->nElemDomain[i];
	}

	double rmax = 0;
	double rmin = 99999999;
	double d;

	if (theGlobalSpace == NULL) {
		theGlobalSpace = femDiscreteCreate(3, FEM_TRIANGLE);
	}
	femDiscrete* theSpace = theGlobalSpace;

	for (int i = 0; i < 3; i++) {
		double xr = mesh->X[mesh->elem[3 * starting + i]];
		double yr = mesh->Y[mesh->elem[3 * starting + i]];
		double currentr = sqrt(xr * xr + yr * yr);
		if (rmax < currentr) {
			rmax = currentr;
		}
		if (rmin > currentr) {
			rmin = currentr;
		}
	}
	d = rmax - rmin;

	// Edite depuis le devoir 4.

	for (int elem = starting; elem < starting + mesh->nElemDomain[10]; elem++) {
		femMeshLocal(theMotor->mesh, elem, map, x, y);
		for (int currentPoint = 0; currentPoint < 3; currentPoint++) {
			double xsi = xsitab[currentPoint];
			double eta = etatab[currentPoint];
			double weight = 1.0 / 6.0;
			femDiscretePhi2(theSpace, xsi, eta, phi);
			femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);
			double xloc = 0;
			double yloc = 0;
			double dxdxsi = 0;
			double dxdeta = 0;
			double dydxsi = 0;
			double dydeta = 0;
			for (int i = 0; i < 3; i++) {
				xloc += x[i] * phi[i];
				yloc += y[i] * phi[i];
				dxdxsi += x[i] * dphidxsi[i];
				dxdeta += x[i] * dphideta[i];
				dydxsi += y[i] * dphidxsi[i];
				dydeta += y[i] * dphideta[i];
			}
			r = sqrt(xloc * xloc + yloc * yloc);
			theta = atan2(yloc, xloc);
			double dxdr = cos(theta);
			double dydr = sin(theta);
			double dxdtheta = -r * sin(theta);
			double dydtheta = r * cos(theta);
			double J = dxdxsi * dydeta - dxdeta * dydxsi;
			J = fabs(J);
			for (int i = 0; i < 3; i++) {
				dphidx[i] = dphidxsi[i] * dydeta - dphideta[i] * dydxsi;
				dphidx[i] = dphidx[i] / J;
				dphidy[i] = dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta;
				dphidy[i] = dphidy[i] / J;
				dphidtheta[i] = dphidx[i] * dxdtheta + dphidy[i] * dydtheta;
				dphidr[i] = dphidx[i] * dxdr + dphidy[i] * dydr;
			}
			dadr = 0;
			dadtheta = 0;
			for (int i = 0; i < 3; i++) {
				dadr += theMotor->a[map[i]] * (dphidr[i]);
				dadtheta += theMotor->a[map[i]] * (dphidtheta[i]);
			}
			C += weight * dadr * dadtheta * J;
		}
	}
	double coef = -theMotor->L;
	coef = coef / (4 * 3.141592653589793 * 1e-7);
	coef = coef / d;
	return coef * C;
}

/*******************************************************
Dztermine quelle bobine allumer pour optimiser le couple
du moteur.
********************************************************/
void motorComputeCurrent(motor* theMotor)
{
	for (int i = 1; i < 7; i++) {
		theMotor->js[i] = 0;
	}
	double js = 8.8464 * 1e5;
	double pi = 3.14159265358979323846;
	double theta = fmod(fabs(theMotor->theta), pi / 2.0);

	if (theta < 0.32) {
		theMotor->js[3] = js;
		theMotor->js[4] = -js;
	}
	else if (theta < 0.82) {
		theMotor->js[1] = js;
		theMotor->js[2] = -js;
	}
	else if (theta < 1.35) {
		theMotor->js[5] = js;
		theMotor->js[6] = -js;
	}
	else {
		theMotor->js[3] = js;
		theMotor->js[4] = -js;
	}
	return;
}

/*******************************************************
IMPORT ET EDITION DEVOIR 4
Import des fonctions du devoir 4 qui resoud l'equation
de Poisson. Des portions du code sont donc redigee par
l'equipe didactique.
De petit edition on ete realisee pour pouvoir appliquer
celui-ci a la nouvelle structure.
Ce fragment de code sert de base au calcul du potentiel
magnetique motorComputeMagneticPotentialBasic.
********************************************************/

typedef struct {
	femMesh* mesh;
	femEdges* edges;
	femDiscrete* space;
	femIntegration* rule;
	femFullSystem* system;
} femPoissonProblem;

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

void femPoissonFree(femPoissonProblem* theProblem)
{
	femFullSystemFree(theProblem->system);
	femIntegrationFree(theProblem->rule);
	femDiscreteFree(theProblem->space);
	femEdgesFree(theProblem->edges);
	free(theProblem);
}

void femPoissonSolve(femPoissonProblem* theProblem, motor* theMotor)
{
	int n = theProblem->space->n;
	int nLocal = theProblem->mesh->nLocalNode;
	double** A = theProblem->system->A;
	double* B = theProblem->system->B;
	double x[3];
	double y[3];
	double phi[3];
	double dphidxsi[3];
	double dphideta[3];
	double dphidx[3];
	double dphidy[3];
	int map[3];
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

/*******************************************************
FIN IMPORT ET EDITION DEVOIR 4
********************************************************/

/*******************************************************
Resoud l'zquation de Poisson pour le moteur theMotor, et
place le resultat dans theMotor->a. Il s'agit de la
version basique par solveur lineaire.
********************************************************/
void motorComputeMagneticPotentialBasic(motor* theMotor)
{
	femPoissonProblem* myProblem = femPoissonCreate(theMotor);
	femPoissonSolve(myProblem, theMotor);
	memcpy(theMotor->a, myProblem->system->B, sizeof(double) * theMotor->mesh->nNode);
	// Liberation de la memoire
	femPoissonFree(myProblem);
	return;
}

/*******************************************************
IMPORT ET EDITION FONCTION BAND
Import des fonctions de fem.c qui permettent de creer un
solveur Bande.  Des portions du code sont donc redigee par
l'equipe didactique.
De petit edition on ete realisee pour pouvoir appliquer
celui-ci a la nouvelle structure.
Ce fragment de code sert de base au calcul du potentiel
magnetique motorComputeMagneticPotentialFem.
********************************************************/

femDiffusionProblem* femDiffusionCreateMotor(motor* theMotor, femSolverType solverType, femRenumType renumType)
{
	int i, band;

	femDiffusionProblem* theProblem = malloc(sizeof(femDiffusionProblem));
	femMesh* theMesh = malloc(sizeof(femMesh));
	theProblem->mesh = theMesh;
	theMesh->elem = theMotor->mesh->elem;
	theMesh->X = theMotor->mesh->X;
	theMesh->Y = theMotor->mesh->Y;
	theMesh->nElem = theMotor->mesh->nElem;
	theMesh->nNode = theMotor->mesh->nNode;
	theMesh->nLocalNode = theMotor->mesh->nLocalNode;
	theProblem->mesh->number = malloc(sizeof(int) * theProblem->mesh->nNode);
	for (int i = 0; i < theProblem->mesh->nNode; i++) {
		theProblem->mesh->number[i] = i;
	}
	if (theProblem->mesh->nLocalNode == 4) {
		theProblem->space = femDiscreteCreate(4, FEM_QUAD);
		theProblem->rule = femIntegrationCreate(4, FEM_QUAD);
	}
	else if (theProblem->mesh->nLocalNode == 3) {
		theProblem->space = femDiscreteCreate(3, FEM_TRIANGLE);
		theProblem->rule = femIntegrationCreate(3, FEM_TRIANGLE);
	}
	theProblem->size = theProblem->mesh->nNode;
	theProblem->sizeLoc = theProblem->mesh->nLocalNode;
	femMeshRenumber(theProblem->mesh, renumType);
	theProblem->sourceValue = 1.0;
	theProblem->dirichletValue = 1.0;


	theProblem->dirichlet = malloc(sizeof(int) * theProblem->size);
	for (i = 0; i < theProblem->size; i++)
		theProblem->dirichlet[i] = 0;
	femEdges* theEdges = femEdgesCreate(theProblem->mesh);
	for (i = 0; i < theEdges->nEdge; i++) {
		if (theEdges->edges[i].elem[1] < 0) {
			theProblem->dirichlet[theEdges->edges[i].node[0]] = 1;
			theProblem->dirichlet[theEdges->edges[i].node[1]] = 1;
		}
	}
	femEdgesFree(theEdges);

	switch (solverType) {
	case FEM_FULL:
		theProblem->solver = femSolverFullCreate(theProblem->size,
			theProblem->sizeLoc); break;
	case FEM_BAND:
		band = femMeshComputeBand(theProblem->mesh);
		theProblem->solver = femSolverBandCreate(theProblem->size,
			theProblem->sizeLoc, band); break;
	case FEM_ITER:
		theProblem->solver = femSolverIterativeCreate(theProblem->size,
			theProblem->sizeLoc); break;
	default: Error("Unexpected solver option");
	}

	theProblem->soluce = malloc(sizeof(double) * theProblem->size);
	for (i = 0; i < theProblem->size; i++)
		theProblem->soluce[i] = 0;

	return theProblem;
}

void femDiffusionFreeMotor(femDiffusionProblem* theProblem)
{
	femIntegrationFree(theProblem->rule);
	femDiscreteFree(theProblem->space);
	free(theProblem->mesh->number);
	free(theProblem->mesh);
	femSolverFree(theProblem->solver);
	free(theProblem->dirichlet);
	free(theProblem->soluce);
	free(theProblem);
}

void femDiffusionComputeMotor(femDiffusionProblem* theProblem, motor* theMotor)
{
	femMesh* theMesh = theProblem->mesh;
	femIntegration* theRule = theProblem->rule;
	femDiscrete* theSpace = theProblem->space;
	femSolver* theSolver = theProblem->solver;
	int* number = theProblem->mesh->number;
	double source = theProblem->sourceValue;
	double dirichlet = theProblem->dirichletValue;

	if (theSpace->n > 4) Error("Unexpected discrete space size !");
	double Xloc[4], Yloc[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
	double Uloc[4];
	int iEdge, iElem, iInteg, i, j, map[4], ctr[4];
	double** A = theSolver->local->A;
	double* Aloc = theSolver->local->A[0];
	double* Bloc = theSolver->local->B;

	for (iElem = 0; iElem < theMesh->nElem; iElem++) {
		for (i = 0; i < theSpace->n; i++)      Bloc[i] = 0;
		for (i = 0; i < (theSpace->n) * (theSpace->n); i++) Aloc[i] = 0;
		femDiffusionMeshLocal(theProblem, iElem, map, ctr, Xloc, Yloc, Uloc);
		for (iInteg = 0; iInteg < theRule->n; iInteg++) {
			double xsi = theRule->xsi[iInteg];
			double eta = theRule->eta[iInteg];
			double weight = theRule->weight[iInteg];
			femDiscretePhi2(theSpace, xsi, eta, phi);
			femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);
			double dxdxsi = 0;
			double dxdeta = 0;
			double dydxsi = 0;
			double dydeta = 0;
			for (i = 0; i < theSpace->n; i++) {
				dxdxsi += Xloc[i] * dphidxsi[i];
				dxdeta += Xloc[i] * dphideta[i];
				dydxsi += Yloc[i] * dphidxsi[i];
				dydeta += Yloc[i] * dphideta[i];
			}
			double jac = dxdxsi * dydeta - dxdeta * dydxsi;
			for (i = 0; i < theSpace->n; i++) {
				dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
				dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
			}
			for (i = 0; i < theSpace->n; i++) {
				for (j = 0; j < theSpace->n; j++) {
					A[i][j] += (dphidx[i] * dphidx[j]
						+ dphidy[i] * dphidy[j]) * jac * weight / theMotor->mu[theMotor->mesh->domain[iElem]];
				}
			}
			for (i = 0; i < theSpace->n; i++) {
				Bloc[i] += phi[i] * jac * source * weight * theMotor->js[theMotor->mesh->domain[iElem]];
			}
		}
		for (i = 0; i < theSpace->n; i++)
			if (ctr[i] == 1) femFullSystemConstrain(theSolver->local, i, 0);
		femSolverAssemble(theSolver, Aloc, Bloc, Uloc, map, theSpace->n);
	}

	double* soluce = femSolverEliminate(theSolver);
	for (i = 0; i < theProblem->size; i++)
		theProblem->soluce[i] += soluce[number[i]];
}

/*******************************************************
FIN IMPORT ET EDITION FONCTION BAND
********************************************************/

/*******************************************************
Resoud l'equation de Poisson pour le moteur theMotor, et
place le resultat dans theMotor->a. Il s'agit de la
version en solver Bande avec renumerotation selon la
coordonnee Y.
********************************************************/
void motorComputeMagneticPotentialFem(motor* theMotor)
{
	femSolverType solverType = FEM_BAND;
	femRenumType  renumType = FEM_YNUM;
	femDiffusionProblem* myProblem = femDiffusionCreateMotor(theMotor, solverType, renumType);
	femDiffusionComputeMotor(myProblem, theMotor);
	memcpy(theMotor->a, myProblem->soluce, sizeof(double) * theMotor->mesh->nNode);
	// Liberation de la memoire
	femDiffusionFreeMotor(myProblem);
	return;
}

/*******************************************************
Resoud l'equation de Poisson pour le moteur theMotor, et
place le resultat dans theMotor->a. On peut ici appeler
la version de son choix Basic ou Fem (solveur Bande).
********************************************************/
void motorComputeMagneticPotential(motor* theMotor)
{
	motorComputeMagneticPotentialFem(theMotor);
	return;
}

/*******************************************************
Libere la memoire partagee par toute les iterations &
fonctions a la fin de la main.
********************************************************/
void motorFree(motor* theMotor) {
	if (theGlobalSpace != NULL) {
		femDiscreteFree(theGlobalSpace);
	}
	free(theMotor->a);
	free(theMotor->movingNodes);
	free(theMotor->js);
	free(theMotor->mu);
	free(theMotor);
}