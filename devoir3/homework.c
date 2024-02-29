#include "fem.h"

# ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename)
{
    femGeo* theGeometry = geoMeshCreate(filename);
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->geo  = theGeometry;
    femMesh *theMesh = theGeometry->theElements;
    if (theMesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4,FEM_QUAD); }
    else if (theMesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE); }
    theProblem->system = femFullSystemCreate(theMesh->nodes->nNodes);
    return theProblem;
}

# endif
# ifndef NOPOISSONBOUNDARY

void femPoissonFindBoundaryNodes(femPoissonProblem *theProblem)
{
    
    femGeo* theGeometry = theProblem->geo;  
    femMesh* theEdges = theGeometry->theEdges; 
    int nBoundary = 0;
    
    //  A completer :-)
    for(int i=0;i<theEdges->nElem;i++){
        if(theEdges->elem[i*2+1]==0){
            nBoundary++;
        }
    }

    femDomain *theBoundary = malloc(sizeof(femDomain));
    theGeometry->nDomains++;
    theGeometry->theDomains = realloc(theGeometry->theDomains,theGeometry->nDomains*sizeof(femDomain*));
    theGeometry->theDomains[theGeometry->nDomains-1] = theBoundary;
    theBoundary->nElem = nBoundary;
    theBoundary->elem = malloc(nBoundary*sizeof(int));
    theBoundary->mesh = NULL;
    sprintf(theBoundary->name,"Boundary");

    
 
    // A completer :-)
    int j=0;
    for(int i=0;i<theEdges->nElem;i++){
        if(theEdges->elem[i*2+1]==0){
            theBoundary->elem[j]=theEdges->elem[i*2];
            j++;
        }
    }
}
    
# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{

    femFullSystemFree(theProblem->system);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    geoMeshFree(theProblem->geo);
    free(theProblem);
}
    
# endif
# ifndef NOPOISSONLOCAL

void femPoissonLocal(femPoissonProblem *theProblem, const int iElem, int *map, double *x, double *y)
{
    femMesh *theMesh = theProblem->geo->theElements;
    int nLocal = theMesh->nLocalNode;
    for(int j=0;j<nLocal;j++){
        map[j] = theMesh->elem[iElem*nLocal+j];
        x[j] = theMesh->nodes->X[map[j]];
        y[j] = theMesh->nodes->Y[map[j]];
    }


}

# endif
# ifndef NOPOISSONSOLVE

void femPoissonSolve(femPoissonProblem *theProblem)
{
    femMesh *mesh = theProblem->geo->theElements;
    femIntegration *rule = theProblem->rule;
    femDiscrete *space = theProblem->space;

    int nlocalNode = mesh->nLocalNode;

    double **A=theProblem->system->A;
    double *B=theProblem->system->B;

    //linear system:
    for(int iElement=0;iElement<mesh->nElem;iElement++){
        double *x=malloc(nlocalNode*sizeof(double));
        double *y=malloc(nlocalNode*sizeof(double));
        int *map=malloc(nlocalNode*sizeof(int));
        femPoissonLocal(theProblem,iElement,map,x,y);
        for(int iIntegral=0;iIntegral<rule->n;iIntegral++){
            double weight=rule->weight[iIntegral];
            double xsi = rule->xsi[iIntegral];
            double eta = rule->eta[iIntegral];

            double *phi=malloc(nlocalNode*sizeof(double));
            femDiscretePhi2(space,xsi,eta,phi);

            double *dphidxsi=malloc(nlocalNode*sizeof(double));
            double *dphideta=malloc(nlocalNode*sizeof(double));
            femDiscreteDphi2(space,xsi,eta,dphidxsi,dphideta);

            double dxdxsi=0;                double dydxsi =0;
            double dxdeta=0;                double dydeta =0;
            for(int i=0;i<space->n;i++){
                dxdxsi+=dphidxsi[i]*x[i]; 
                dydxsi+=dphidxsi[i]*y[i];
                dxdeta+=dphideta[i]*x[i];
                dydeta+=dphideta[i]*y[i];
            }

            double jacobian = dxdxsi*dydeta-dxdeta*dydxsi;
            if(jacobian<0){
                int node = mesh->elem[iElement*nlocalNode];
                mesh->elem[nlocalNode*iElement]=mesh->elem[nlocalNode*iElement+2];
                mesh->elem[nlocalNode*iElement+2]=node;
            }
            jacobian=fabs(jacobian);
            double *dphidx=malloc(space->n*sizeof(double));
            double *dphidy=malloc(space->n*sizeof(double));
            for(int i=0;i<space->n;i++){
                dphidx[i]=(dphidxsi[i]*dydeta-dphideta[i]*dydxsi)/jacobian;
                dphidy[i]=(dphideta[i]*dxdxsi+dphidxsi[i]*dxdeta)/jacobian;
            }

            for(int i=0;i<space->n;i++){
                B[map[i]]+=jacobian*weight*phi[i];
                for(int j=0;j<space->n;j++){
                    A[map[i]][map[j]]+=jacobian*weight*(dphidx[i]*dphidx[j]+dphidy[i]*dphidy[j]);
                }
            }
            free(dphidx);
            free(dphidy);
            free(dphidxsi);
            free(dphideta);
            free(phi);
        }
        free(x);
        free(y);
        free(map);
    }
    femMesh *Edges = theProblem->geo->theEdges;
    for(int iEdge=0;iEdge<Edges->nElem;iEdge++){
        double value=0.0;
        for (int i = 0; i < 2; i++){
            femFullSystemConstrain(theProblem->system,Edges->elem[i],value);
        }
        
    }

    
    femFullSystemEliminate(theProblem->system);

}

# endif



