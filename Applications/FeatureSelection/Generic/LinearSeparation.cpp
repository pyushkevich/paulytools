#include "LinearSeparation.h"
#include "globals.h"

LinearSeparation::LinearSeparation() 
    : AbstractLPClassification()
{
	m_provider = NULL;
}

LinearSeparation::~LinearSeparation() {
}

double LinearSeparation::getSeparatingPlane(Vec &w) {
	w = wResult;
	return gResult;
}

void LinearSeparation::initialize(const Mat &A, const Mat &B) 
{
    // Is there a provider?
    assert(m_provider);

    // Get the basic parameters
    m = A.rows();
    k = B.rows();
    n = A.columns();

    // Construct the linear programming matrix
    Mat M(m+k,m+k+n+1);

    // Insert the data matrices and window matrix
    M.insertMatrix(0,0,A);
    M.insertMatrix(m,0,-B);

    // Insert all the identities
    insertDiagonal(M,0     ,n     ,m+k ,1);

    // Insert the vertical unit vectors
    insertVertical(M,0,n+m+k,m,-1);
    insertVertical(M,m,n+m+k,k,1);

    // Now, specify the vector B
    Vec b(m+k);
    b.setAll(1);

    // Construct the basic vector C
    c.setSize(n+m+k+1);
    insertVertical(c,n,0,m,1.0/m);
    insertVertical(c,n+m,0,k,1.0/k);

    // Construct the upper and lower limits
    Vec upper(n+m+k+1), lower(n+m+k+1);
    upper.setAll(1e100);
    insertVertical(lower,0,0,n,-1e100);
    lower(n+m+k) = -1e100;

    // Create a soplex problem
    m_provider->setProblem(c,M,b,lower,upper);
    state = INIT;
}

bool LinearSeparation::run() 
{
    assert(state == INIT && m_provider);

    // Solve the problem
    Vec result;
	result.setSize(c.size());
    bool status = m_provider->solve(result);

    // Now, send back some results
    if(!status) {
    	result.setAll(0);
    }
    
    wResult.setSize(n);
    result.extractMatrix(0,0,wResult);
    gResult = result(m+k+n);
    objResult = c.dotProduct(result);

    state = status ? SUCCESS : FAILURE;
    return status;
}

