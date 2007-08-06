#include "BruteForceSubdivisionMedialModel.h"

BruteForceSubdivisionMedialModel::BruteForceSubdivisionMedialModel() :
  SubdivisionMedialModel()
{
  xAtoms = NULL;
}

BruteForceSubdivisionMedialModel::~BruteForceSubdivisionMedialModel()
{
}

void
BruteForceSubdivisionMedialModel
::SetMesh(const MeshLevel &mesh, 
  const Vec &C, const Vec &u, const Vec &v,
  size_t nAtomSubs, size_t nCoeffSubs)
{
  // Call the parent method
  SubdivisionMedialModel::SetMesh(mesh, C, u, v, nAtomSubs, nCoeffSubs);

  // Pass the mesh to the loop scheme
  xLoopScheme.SetMeshLevel(&mlAtom);

  // Compute the sparse matrix that gives u and v derivatives of 
  // functions/vectors over the mesh. These are going to be derivatives
  // with respect to the U and V arrays passed in, not to some arbitrary
  // local parameterization (which is referred to as p, q)
  SparseMat::STLSourceType srcWuv;
  for(size_t i = 0; i < mlAtom.nVertices; i++)
    {
    // Generate the row
    SparseMat::STLRowType row;

    // Compute the weigths for this row. The weights arise in the following:
    //   dF_i/dp = Sum_{j \in N(i)} W^p_j F_j
    //   dF_i/du = Sum_{j \in N(i)} (dp/du W^p_j + dq/du W^p_j) F_j
    // So we are here computing (dp/du W^p_j + dq/du W^p_j) and storing them
    // as elements of the sparse matrix Wu. We also do the same for Wv
    
    // First, we are going to compute du/dp and dv/dp
    double dudp = xLoopScheme.GetPartialDerivative(0, i, uAtom.data_block());
    double dudq = xLoopScheme.GetPartialDerivative(1, i, uAtom.data_block());
    double dvdp = xLoopScheme.GetPartialDerivative(0, i, vAtom.data_block());
    double dvdq = xLoopScheme.GetPartialDerivative(1, i, vAtom.data_block());

    // Now compute the determinant of the Jacobian
    double detj = 1.0 / (dudp * dvdq - dudq * dvdp);

    // Now compute the inverse of the Jacobian
    double dpdu =   detj * dvdq;
    double dpdv = - detj * dudq;
    double dqdu = - detj * dvdp;
    double dqdv =   detj * dudp; 

    // Compute how much the vertex itself contributes to the weights
    row.push_back(
      make_pair(i, 
        make_pair(
          xLoopScheme.GetOwnWeight(0, i) * dpdu + 
          xLoopScheme.GetOwnWeight(1, i) * dqdu,
          xLoopScheme.GetOwnWeight(0, i) * dpdv +
          xLoopScheme.GetOwnWeight(1, i) * dqdv)));

    // Now, compute the contributing weights for each vertex
    for(EdgeWalkAroundVertex w(&mlAtom, i); !w.IsAtEnd(); ++w)
      {
      row.push_back(
        make_pair(w.MovingVertexId(), 
          make_pair(
            xLoopScheme.GetNeighborWeight(0, w) * dpdu + 
            xLoopScheme.GetNeighborWeight(1, w) * dqdu,
            xLoopScheme.GetNeighborWeight(0, w) * dpdv + 
            xLoopScheme.GetNeighborWeight(1, w) * dqdv)));
      }

    // Store the rows
    srcWuv.push_back(row);
    }

  // Set the derivative matrices
  Wuv.SetFromSTL(srcWuv, mlAtom.nVertices);
}

void
BruteForceSubdivisionMedialModel
::ComputeAtoms(const double *xHint)
{
  size_t i;

  // The first step is to compute the X and rho of the atoms based on the
  // coefficients. This step is performed in this class, not in the solver
  for(i = 0; i < mlAtom.nVertices; i++)
    {
    // Set up i-th atom
    MedialAtom &a = xAtoms[i]; a.X.fill(0.0); a.R = 0.0;

    // Compute the weighted sum of the coefficients
    ImmutableSparseMatrix<double>::RowIterator it = mlAtom.weights.Row(i);
    for( ; !it.IsAtEnd(); ++it)
      {
      size_t j = it.Column() << 2; double w = it.Value();
      a.X += w * xCoefficients.extract(3, j);
      a.R += w * xCoefficients[j+3];
      }

    // Set F = R^2
    a.F = a.R * a.R;
    }

  // We must now compute the derivatives of X and R to make real atoms
  for(i = 0; i < mlAtom.nVertices; i++)
    {
    // Set up i-th atom
    MedialAtom &a = xAtoms[i];

    // Compute the partial derivatives of F and X
    a.Xu.fill(0.0); a.Xv.fill(0.0); a.Fu = 0; a.Fv = 0;
    for(SparseMat::RowIterator it = Wuv.Row(i); !it.IsAtEnd(); ++it)
      {
      MedialAtom &anbr = xAtoms[it.Column()];
      double wu = it.Value().first;
      double wv = it.Value().second;

      a.Xu += wu * anbr.X;
      a.Xv += wv * anbr.X;
      a.Fu += wu * anbr.F;
      a.Fv += wv * anbr.F;
      }

    /*
    // Compute partial derivatives
    /// a.Fu = xLoopScheme.Fu(i, xAtoms);
    /// a.Fv = xLoopScheme.Fv(i, xAtoms);
    /// a.Xu = xLoopScheme.Xu(i, xAtoms);
    /// a.Xv = xLoopScheme.Xv(i, xAtoms);

    // Compute the differential geometry of the atom
    a.G.SetOneJet(a.X.data_block(), a.Xu.data_block(), a.Xv.data_block());

    // Compute the normal vector
    a.ComputeNormalVector();

    // Compute all the components of the atoms
    a.ComputeBoundaryAtoms(!mlAtom.IsVertexInternal(i));
    */

    // Compute some form of strain tensor. Loop over all neighbors
    /*
    SMLVec3d XU, XV; XU.fill(0); XV.fill(0);
    EdgeWalkAroundVertex it(&mlAtom, i);

    // Flag all neighbors
    for(; !it.IsAtEnd(); ++it)
      {
      MedialAtom &aj = xAtoms[it.MovingVertexId()];
      XU += (aj.X - a.X) / (aj.u - a.u);
      XV += (aj.X - a.X) / (aj.v - a.v);
      }

    double E = dot_product(XU, XU);
    double F = dot_product(XU, XV);
    double G = dot_product(XV, XV);
    a.G.xChristoffelSecond[0][0][0] = E + G;
    */
    

    // Compute the Riemannian gradient of the parameterization (for stretch tensors)
    /*
    double udu = xLoopScheme.Uu(i, xAtoms);
    double udv = xLoopScheme.Uv(i, xAtoms);
    double vdu = xLoopScheme.Vu(i, xAtoms);
    double vdv = xLoopScheme.Vv(i, xAtoms);
    double xMagGradU = 
      udu * udu * a.G.xContravariantTensor[0][0] + 
      2 * udu * udv * a.G.xContravariantTensor[0][1] + 
      udv * udv * a.G.xContravariantTensor[1][1];
    double xMagGradV = 
      vdu * vdu * a.G.xContravariantTensor[0][0] + 
      2 * vdu * vdv * a.G.xContravariantTensor[0][1] + 
      vdv * vdv * a.G.xContravariantTensor[1][1];
    double xStretch = xMagGradV + xMagGradU;
    a.G.xChristoffelSecond[0][0][0] = xStretch;
    */
    }

  // On the third pass, we compute the second order partial derivatives
  for(i = 0; i < mlAtom.nVertices; i++)
    {
    // Set up i-th atom
    MedialAtom &a = xAtoms[i];

    // Compute the partial derivatives of F and X
    a.Xuu.fill(0.0); a.Xuv.fill(0.0); a.Xvv.fill(0.0);
    a.Fuu = 0.0; a.Fuv = 0.0; a.Fvv = 0.0;
    SMLVec3d Xvu; Xvu.fill(0.0);
    for(SparseMat::RowIterator it = Wuv.Row(i); !it.IsAtEnd(); ++it)
      {
      MedialAtom &anbr = xAtoms[it.Column()];
      double wu = it.Value().first;
      double wv = it.Value().second;

      a.Xuu += wu * anbr.Xu;
      a.Xuv += 0.5 * (wu * anbr.Xv + wv * anbr.Xu);
      a.Xvv += wv * anbr.Xv;
      a.Fuu += wu * anbr.Fu;
      a.Fuv += 0.5 * (wu * anbr.Fv + wv * anbr.Fu);
      a.Fvv += wv * anbr.Fv;
      }

    // Compute the atom's differential geometry, normal and boundary
    a.ComputeDifferentialGeometry();
    a.ComputeNormalVector();
    a.ComputeBoundaryAtoms(!mlAtom.IsVertexInternal(i));

    // Compute second order differential geometry
    // a.ComputeDifferentialGeometry();
    
    // The second order derivatives are set to zero
    /// a.Xuu.fill(0.0); a.Xuv.fill(0.0); a.Xvv.fill(0.0);
    
    // Compute the derivatives of the normal vector
    /// SMLVec3d Nu = xLoopScheme.Nu(i, xAtoms);
    /// SMLVec3d Nv = xLoopScheme.Nv(i, xAtoms);

    // Compute the derivative of the area element
    /// double Au = xLoopScheme.AreaEltDu(i, xAtoms);
    /// double Av = xLoopScheme.AreaEltDv(i, xAtoms);
    /// double AGradMag2 = Au * Au + Av * Av;
    /// a.Ruv = AGradMag2;

    // Compute the principal curvatures
    /// double ee = - dot_product(Nv, a.Xv);
    /// double ff = -0.5 * (dot_product(Nu, a.Xv) + dot_product(Nv, a.Xu));
    /// double gg = - dot_product(Nu, a.Xu);
    /// a.xMeanCurv = 0.5 * (
    ///   ee * a.G.xContravariantTensor[0][0] +
    ///   2 * ff * a.G.xContravariantTensor[1][0] + 
    ///   gg * a.G.xContravariantTensor[1][1]);
    /// a.xGaussCurv = (ee * gg - ff * ff) * a.G.gInv;
    /// a.Xuu = Nu;
    /// a.Xvv = Nv;
    }
}

void 
printvec(SMLVec3d &x)
{
  cout << "{" << x[0] << "," << x[1] << "," << x[2] << "}";
}

void 
BruteForceSubdivisionMedialModel::
ComputeBoundaryCurvature(Vec &xMeanCurv, Vec &xGaussCurv)
{
  // The curvatures are computed at non-crest locations. At the crest, the rules are
  // different and I don't feel like screwing around with them
  for(MedialBoundaryPointIterator bip(this->GetIterationContext());
    !bip.IsAtEnd(); ++bip)
    {
    // Get the medial atom and side
    size_t ib = bip.GetIndex();
    size_t ia = bip.GetAtomIndex();
    size_t side = bip.GetBoundarySide();
    MedialAtom &a = xAtoms[ia];

    // Just set up an edge walk
    xMeanCurv[ib] = 0.0;
    xGaussCurv[ib] = 0.0;

    // Just set up an edge walk
    if(a.flagValid && !a.flagCrest)
      {
      // For debugging purposes
      /*
      if(ib!=17) continue;
      cout << "*** IB = " << ib << endl;
      cout << "{X,Xu,Xv,Xuu,Xuv,Xvv} = {"; 
                          printvec(a.X);
      cout << " , ";      printvec(a.Xu);
      cout << " , ";      printvec(a.Xv);
      cout << " , ";      printvec(a.Xuu);
      cout << " , ";      printvec(a.Xuv);
      cout << " , ";      printvec(a.Xvv);
      cout << "};" << endl;
      cout << "{F,Fu,Fv,Fuu,Fuv,Fvv} = {"
        << a.F << ", " << a.Fu << ", " << a.Fv << ", "
        << a.Fuu << ", " << a.Fuv << ", " << a.Fvv << "};" << endl;
      */

      // Create arrays to allow easier looping
      SMLVec3d Xi[2] = {a.Xu, a.Xv};
      SMLVec3d Xij[2][2] = {{a.Xuu, a.Xuv},{a.Xuv, a.Xvv}};
      double Fi[2] = {a.Fu, a.Fv};
      double Fij[2][2] = {{a.Fuu, a.Fuv},{a.Fuv, a.Fvv}};

      // Compute the partial derivatives of GradF
      SMLVec3d gradF_i[2], gradR_i[2], Ni[2], NB_i[2], XB_i[2];

      size_t i,j,k,l;
      for(i=0;i<2;i++) 
        {
        // Derivative of gradF wrt U^i
        gradF_i[i].fill(0.0);

        // Sum over j, k 
        for(j=0;j<2;j++) for(k=0;k<2;k++)
          {
          // Compute dg^jk/du^i
          double dg_jk_d_ui = 0.0;
          for(l=0;l<2;l++)
            {
            dg_jk_d_ui += -(
              a.G.xContravariantTensor[l][k] * a.G.xChristoffelSecond[l][i][j] +
              a.G.xContravariantTensor[j][l] * a.G.xChristoffelSecond[l][i][k]);
            }

          // Compute contribution to gradFi
          gradF_i[i] += 
            (dg_jk_d_ui * Fi[k]) * Xi[j]
            + (a.G.xContravariantTensor[j][k] *  Fi[k]) * Xij[j][i] 
            + (a.G.xContravariantTensor[j][k] * Fij[k][i]) * Xi[j];
          }

        // Compute the partial derivative of gradR
        gradR_i[i] = (gradF_i[i] * a.R - a.xGradR * Fi[i]) / (2.0 * a.F);

        // Compute the partial derivative of the unit normal to the medial axis
        double dg_di = 2 * a.G.g * (
          a.G.xChristoffelSecond[0][i][0] + a.G.xChristoffelSecond[1][i][1]);
        SMLVec3d dN0_di = 
          vnl_cross_3d(Xij[0][i], Xi[1]) + vnl_cross_3d(Xi[0], Xij[1][i]);
        Ni[i] = (dN0_di * a.aelt - 0.5 * a.N * dg_di) / (a.G.g);

        // Compute the boundary node's normal
        NB_i[i] = 
          - gradR_i[i] + (side ? 1.0 : -1.0) * (
            - dot_product(a.xGradR, gradR_i[i]) * a.N / a.xNormalFactor +
            a.xNormalFactor * Ni[i]);

        // Compute the boundary node's position
        XB_i[i] = Xi[i] + (Fi[i] / (2 * a.R)) * a.xBnd[side].N + a.R * NB_i[i];
        }

      // Now that we have the first derivatives of the normal and X on the boundary, curvature
      // computation is easy
      SMLVec3d &Xu = XB_i[0], &Xv = XB_i[1], &Nu = NB_i[0], &Nv = NB_i[1];

      /* 
      cout << "BND[" << side << "] : " << endl;
      cout << "X = " << a.xBnd[side].X << endl;
      cout << "N = " << a.xBnd[side].N << endl;
      cout << "Xu = " << Xu << endl;
      cout << "Xv = " << Xv << endl;
      cout << "Nu = " << Nu << endl;
      cout << "Nv = " << Nv << endl;
      */

      double ee = - dot_product(Nv,Xv);
      double ff = -0.5 * (dot_product(Nu, Xv) + dot_product(Nv, Xu));
      double gg = - dot_product(Nu, Xu);
      double E = dot_product(Xu, Xu);
      double F = dot_product(Xu, Xv);
      double G = dot_product(Xv, Xv);
      double det = (E*G - F*F);

      xMeanCurv[ib] = (ee * G - 2 * ff * F + gg * E) / (2 * det);
      xGaussCurv[ib] = (ee * gg - ff * ff) / det;
      }
    }

}
/*
void 
BruteForceSubdivisionMedialModel::
ComputeBoundaryCurvature(Vec &xMeanCurv, Vec &xGaussCurv)
{
  // The curvatures are computed at non-crest locations. At the crest, the rules are
  // different and I don't feel like screwing around with them
  for(MedialBoundaryPointIterator bip(this->GetIterationContext());
    !bip.IsAtEnd(); ++bip)
    {
    // Get the medial atom and side
    size_t ib = bip.GetIndex();
    size_t ia = bip.GetAtomIndex();
    size_t side = bip.GetBoundarySide();

    // Just set up an edge walk
    xMeanCurv[ib] = 0.0;
    xGaussCurv[ib] = 0.0;

    // Just set up an edge walk
    EdgeWalkAroundVertex it(&mlAtom, ia);
    if(!it.IsOpen()) 
      {
      // The vector whose magnitude is the mean curvature
      double area = 0.0;
      SMLVec3d kvec(0.0);

      // The sum of the angles around I
      double sum_theta = 0.0;

      double scope = 0;

      // Loop around vertex
      for( ; !it.IsAtEnd(); ++it)
        {
        // Naming convention:
        // I - central vertex; J - opposite vertex; A - ahead vertex; B - behind vertex
        // AI - vector from A to I, etc.
        SMLVec3d &I = xAtoms[ia].xBnd[side].X;
        SMLVec3d &J = xAtoms[it.MovingVertexId()].xBnd[side].X;
        SMLVec3d &A = xAtoms[it.VertexIdAhead()].xBnd[side].X;
        SMLVec3d &B = xAtoms[it.VertexIdBehind()].xBnd[side].X;
        SMLVec3d AI = I - A;
        SMLVec3d AJ = J - A; 
        SMLVec3d BI = I - B;
        SMLVec3d BJ = J - B;
        SMLVec3d IA = A - I;
        SMLVec3d IJ = J - I;
    
        // Compute cot(a)
        double cos_a = dot_product(AI, AJ);
        double sin_a = 
          sqrt(dot_product(AI, AI) * dot_product(AJ, AJ) - cos_a * cos_a);
        double cot_a = cos_a / sin_a;

        // Compute cot(b)
        double cos_b = dot_product(BI, BJ);
        double sin_b = 
          sqrt(dot_product(BI, BI) * dot_product(BJ, BJ) - cos_b * cos_b);
        double cot_b = cos_b / sin_b;

        // Compute the contribution to mean curvature
        kvec += - 0.5 * (cot_a + cot_b) * IJ;
        scope += kvec[0];

        // Compute the area fraction
        area += TriangleArea(I, J, A);

        // Compute the contribution to the Gaussian curvature
        double num = dot_product(IA, IJ);
        double den = sqrt(dot_product(IA,IA) * dot_product(IJ,IJ));
        double cos_t = num / den;
        double theta = acos(cos_t);
        sum_theta += theta;
        }

      // Compute the true values
      SMLVec3d kvec_scale = 3.0 * kvec / area;
      double angle_def = sum_theta - 2.0 * M_PI;

      double kvec_scale_mag2 = dot_product(kvec_scale, kvec_scale);
      xMeanCurv[ib] = 0.5 * sqrt(kvec_scale_mag2);
      xGaussCurv[ib] = 3.0 * angle_def / area;

      // TEMP
      // xMeanCurv[ib] = scope;
      // xGaussCurv[ib] = sum_theta;
      }
    }
*/

    // Compute differential geometry for each atom
    /* 
    // Compute partials of GradR
    SMLVec3d gradPhiDk[2];
    SMLVec3d gradR_D_Uk[2];
    SMLVec3d D_N_side_Dk[2][2];
    SMLVec3d D_X_side_Dk[2][2];
    for(size_t k = 0; k < 2; k++)
      {
      gradPhiDk[k].fill(0.0);
      for(size_t i = 0; i < 2; i++)
        {
        for(size_t j = 0; j < 2; j++)
          {
          double D_Gij_D_Uk = 0.0;
          for(size_t z = 0; z < 2; z++)
            {
            D_Gij_D_Uk += 
              - a.xContravariantTensor[z][j] * a.xChristoffelSecond[z][k][i]
              - a.xContravariantTensor[i][z] * a.xChristoffelSecond[z][k][j];
            }
          
          gradPhiDk[k] += 
              D_Gij_D_Uk * D_X_D_Ui[i] * D_F_D_Uj[j] 
            + a.xContravariantTensor[i][j] * D_X_D_UiUj[i][k] * D_F_D_Uj[j]
            + a.xContravariantTensor[i][j] * D_X_D_Ui[i] * D_F_D_UiUj[i][k];
          }
        }

      // Compute the partials of GradR
      gradR_D_Uk[k] = 0.5 * (a.gradPhiDk[k] * a.R - a.gradR * D_F_D_Uj[k]) / a.F;

      // Compute the partials of boundary bormal
      SMLVec3d D_N_D_Uk = 
      
      }



    // Compute the partials of GradR




    // Compute the partial derivatives of the normal vector 
    SMLVec3d Xu(0.0), Xv(0.0), Nu(0.0), Nv(0.0);
    for(SparseMat::RowIterator it = Wuv.Row(ia); !it.IsAtEnd(); ++it)
      {
      MedialAtom &anbr = xAtoms[it.Column()];
      size_t s = (anbr.flagCrest) ? 0 : side;
      double wu = it.Value().first;
      double wv = it.Value().second;
      Nu += wu * anbr.xBnd[s].N;
      Nv += wv * anbr.xBnd[s].N;
      Xu += wu * anbr.xBnd[s].X;
      Xv += wv * anbr.xBnd[s].X;
      }

    // Compute the first and second fundamental forms
    double E = dot_product(Xu, Xu);
    double F = dot_product(Xu, Xv);
    double G = dot_product(Xv, Xv);
    double det = E*G-F*F;
    
    double ee = - dot_product(Nv, Xv);
    double ff = -0.5 * (dot_product(Nu, Xv) + dot_product(Nv, Xu));
    double gg = - dot_product(Nu, Xu);

    cout << "ia = " << ia << "; side = " << side << endl;
    cout << "E = " << E << endl;
    cout << "F = " << F << endl;
    cout << "G = " << G << endl;
    cout << "ee = " << ee << endl;
    cout << "ff = " << ff << endl;
    cout << "gg = " << gg << endl;

    xMeanCurv[bip.GetIndex()] = 
      0.5 * (ee * G - 2 * ff * F + gg * E) / det;
    xGaussCurv[bip.GetIndex()] = 
      (ee * gg - ff * ff) / det;
    }
  */
//}

void 
BruteForceSubdivisionMedialModel::
ComputeBoundaryCurvaturePartial(
  Vec &dMeanCurv, Vec &dGaussCurv, MedialAtom *dAtoms)
{
  // The curvatures are computed at non-crest locations. At the crest, the rules are
  // different and I don't feel like screwing around with them
  for(MedialBoundaryPointIterator bip(this->GetIterationContext());
    !bip.IsAtEnd(); ++bip)
    {
    // Get the medial atom and side
    size_t ib = bip.GetIndex();
    size_t ia = bip.GetAtomIndex();
    size_t side = bip.GetBoundarySide();
    MedialAtom &a = xAtoms[ia];
    MedialAtom &da = dAtoms[ia];

    // Just set up an edge walk
    dMeanCurv[ib] = 0.0;
    dGaussCurv[ib] = 0.0;

    // Just set up an edge walk
    if(a.flagValid && !a.flagCrest)
      {
      /* 
      cout << "*** IB = " << ib << endl;
      cout << "{X,Xu,Xv,Xuu,Xuv,Xvv} -> \n{"; 
                          printvec(a.X);
      cout << " , \n";      printvec(a.Xu);
      cout << " , \n";      printvec(a.Xv);
      cout << " , \n";      printvec(a.Xuu);
      cout << " , \n";      printvec(a.Xuv);
      cout << " , \n";      printvec(a.Xvv);
      cout << "};" << endl;
      cout << "{F,Fu,Fv,Fuu,Fuv,Fvv} -> {"
        << a.F << ", " << a.Fu << ", " << a.Fv << ", \n"
        << a.Fuu << ", " << a.Fuv << ", " << a.Fvv << "};" << endl;

      cout << "{dX,dXu,dXv,dXuu,dXuv,dXvv} -> \n{"; 
                          printvec(da.X);
      cout << " , \n";      printvec(da.Xu);
      cout << " , \n";      printvec(da.Xv);
      cout << " , \n";      printvec(da.Xuu);
      cout << " , \n";      printvec(da.Xuv);
      cout << " , \n";      printvec(da.Xvv);
      cout << "};" << endl;
      cout << "{dF,dFu,dFv,dFuu,dFuv,dFvv} -> {"
        << da.F << ", " << da.Fu << ", " << da.Fv << ", \n"
        << da.Fuu << ", " << da.Fuv << ", " << da.Fvv << "};" << endl;
      */

      // Create arrays to allow easier looping
      SMLVec3d Xi[2] = {a.Xu, a.Xv};
      SMLVec3d Xij[2][2] = {{a.Xuu, a.Xuv},{a.Xuv, a.Xvv}};
      SMLVec3d DXi[2] = {da.Xu, da.Xv};
      SMLVec3d DXij[2][2] = {{da.Xuu, da.Xuv},{da.Xuv, da.Xvv}};

      double Fi[2] = {a.Fu, a.Fv};
      double Fij[2][2] = {{a.Fuu, a.Fuv},{a.Fuv, a.Fvv}};
      double DFi[2] = {da.Fu, da.Fv};
      double DFij[2][2] = {{da.Fuu, da.Fuv},{da.Fuv, da.Fvv}};

      // Compute the partial derivatives of GradF
      SMLVec3d gradF_i[2], gradR_i[2], Ni[2], NB_i[2], XB_i[2];
      SMLVec3d d_gradF_i[2], d_gradR_i[2], d_Ni[2], d_NB_i[2], d_XB_i[2];

      size_t i,j,k,l;
      for(i=0;i<2;i++) 
        {
        // Derivative of gradF wrt U^i
        gradF_i[i].fill(0.0);
        d_gradF_i[i].fill(0.0);

        // Sum over j, k 
        for(j=0;j<2;j++) for(k=0;k<2;k++)
          {
          // Compute dg^jk/du^i
          double dg_jk_d_ui = 0.0;
          double d_dg_jk_d_ui = 0.0;
          for(l=0;l<2;l++)
            {
            dg_jk_d_ui += -(
              a.G.xContravariantTensor[l][k] * a.G.xChristoffelSecond[l][i][j] +
              a.G.xContravariantTensor[j][l] * a.G.xChristoffelSecond[l][i][k]);

            d_dg_jk_d_ui += -(
              da.G.xContravariantTensor[l][k] * a.G.xChristoffelSecond[l][i][j] +
              a.G.xContravariantTensor[l][k] * da.G.xChristoffelSecond[l][i][j] +
              da.G.xContravariantTensor[j][l] * a.G.xChristoffelSecond[l][i][k] +
              a.G.xContravariantTensor[j][l] * da.G.xChristoffelSecond[l][i][k]);
            }

          // Compute contribution to gradFi
          gradF_i[i] += 
            (dg_jk_d_ui * Fi[k]) * Xi[j]
            + (a.G.xContravariantTensor[j][k] *  Fi[k]) * Xij[j][i] 
            + (a.G.xContravariantTensor[j][k] * Fij[k][i]) * Xi[j];

          // Holy cow, what a mess!
          d_gradF_i[i] += 
            (d_dg_jk_d_ui * Fi[k]) * Xi[j]
            + (dg_jk_d_ui * DFi[k]) * Xi[j]
            + (dg_jk_d_ui * Fi[k]) * DXi[j]
            + (da.G.xContravariantTensor[j][k] *  Fi[k]) * Xij[j][i] 
            + (a.G.xContravariantTensor[j][k] *  DFi[k]) * Xij[j][i] 
            + (a.G.xContravariantTensor[j][k] *  Fi[k]) * DXij[j][i] 
            + (da.G.xContravariantTensor[j][k] * Fij[k][i]) * Xi[j]
            + (a.G.xContravariantTensor[j][k] * DFij[k][i]) * Xi[j]
            + (a.G.xContravariantTensor[j][k] * Fij[k][i]) * DXi[j];
          }

        // cout << "d_gradF_i[" << i << " = " << d_gradF_i[i] << endl;

        // Compute the partial derivative of gradR
        gradR_i[i] = (gradF_i[i] * a.R - a.xGradR * Fi[i]) / (2.0 * a.F);
        d_gradR_i[i] = 
          ((d_gradF_i[i] * a.R + gradF_i[i] * da.R 
           - (da.xGradR * Fi[i] + a.xGradR * DFi[i])) * (a.F)
           - (gradF_i[i] * a.R - a.xGradR * Fi[i]) * (da.F))
          / (2.0 * a.F * a.F);
             
        // cout << "d_gradR_i[" << i << " = " << d_gradR_i[i] << endl;

        // Compute the partial derivative of the unit normal to the medial axis
        double dg_di = 2 * a.G.g * (
          a.G.xChristoffelSecond[0][i][0] + a.G.xChristoffelSecond[1][i][1]);
        double d_dg_di = 
          2 * da.G.g * (
            a.G.xChristoffelSecond[0][i][0] + a.G.xChristoffelSecond[1][i][1])
          + 2 * a.G.g * (
            da.G.xChristoffelSecond[0][i][0] + da.G.xChristoffelSecond[1][i][1]);
            
        SMLVec3d dN0_di = 
          vnl_cross_3d(Xij[0][i], Xi[1]) + vnl_cross_3d(Xi[0], Xij[1][i]);
        SMLVec3d d_dN0_di = 
          vnl_cross_3d(DXij[0][i], Xi[1]) 
          + vnl_cross_3d(Xij[0][i], DXi[1]) 
          + vnl_cross_3d(DXi[0], Xij[1][i])
          + vnl_cross_3d(Xi[0], DXij[1][i]);

        Ni[i] = (dN0_di * a.aelt - 0.5 * a.N * dg_di) / (a.G.g);
        d_Ni[i] = (
          (d_dN0_di * a.aelt + dN0_di * da.aelt - 0.5 * 
            (da.N * dg_di + a.N * d_dg_di)) * a.G.g
          - (dN0_di * a.aelt - 0.5 * a.N * dg_di) * da.G.g) / (a.G.g * a.G.g);

        // cout << "d_Ni[" << i << " = " << d_Ni[i] << endl;

        // Compute some intermediate terms first
        double tmp1 = - dot_product(a.xGradR, gradR_i[i]);
        double d_tmp1 = - dot_product(da.xGradR, gradR_i[i]) 
          - dot_product(a.xGradR, d_gradR_i[i]);

        SMLVec3d tmp2 = a.N / a.xNormalFactor;
        SMLVec3d d_tmp2 = (da.N * a.xNormalFactor - a.N * da.xNormalFactor) /
          (a.xNormalFactor * a.xNormalFactor);

        // Compute the boundary node's normal
        NB_i[i] = 
          - gradR_i[i] + (side ? 1.0 : -1.0) * (
            tmp1 * tmp2 + a.xNormalFactor * Ni[i]);

        d_NB_i[i] = 
          - d_gradR_i[i] + (side ? 1.0 : -1.0) * (
            d_tmp1 * tmp2 + tmp1 * d_tmp2 
            + da.xNormalFactor * Ni[i] + a.xNormalFactor * d_Ni[i]);

        // cout << "a.xNF " << a.xNormalFactor << endl;
        // cout << "da.xNF " << da.xNormalFactor << endl;
        // cout << "d_NB[" << i << " = " << d_NB_i[i] << endl;

        // Compute the boundary node's position
        double Ri = 0.5 * Fi[i] / a.R;
        double d_Ri = 0.5 * (DFi[i] * a.R - Fi[i] * da.R) / a.F;

        XB_i[i] = Xi[i] + Ri * a.xBnd[side].N + a.R * NB_i[i];
        d_XB_i[i] = 
          DXi[i] 
          + d_Ri * a.xBnd[side].N 
          + Ri * da.xBnd[side].N
          + da.R * NB_i[i]
          + a.R * d_NB_i[i];

        // cout << "d_XB[" << i << " = " << d_XB_i[i] << endl;
        }

      // Now that we have the first derivatives of the normal and X on the boundary, curvature
      // computation is easy
      SMLVec3d &Xu = XB_i[0], &Xv = XB_i[1], &Nu = NB_i[0], &Nv = NB_i[1];
      SMLVec3d &dXu = d_XB_i[0], &dXv = d_XB_i[1], &dNu = d_NB_i[0], &dNv = d_NB_i[1];

      /* 
      cout << "BND[" << side << "] : " << endl;
      cout << "X = " << a.xBnd[side].X << endl;
      cout << "N = " << a.xBnd[side].N << endl;
      cout << "Xu = " << Xu << endl;
      cout << "Xv = " << Xv << endl;
      cout << "Nu = " << Nu << endl;
      cout << "Nv = " << Nv << endl;
      */

      double ee = - dot_product(Nv,Xv);
      double d_ee = - dot_product(dNv,Xv) - dot_product(Nv, dXv);
      double ff = -0.5 * (dot_product(Nu, Xv) + dot_product(Nv, Xu));
      double d_ff = -0.5 * (
          dot_product(dNu, Xv) 
        + dot_product(Nu, dXv) 
        + dot_product(dNv, Xu)
        + dot_product(Nv, dXu));
      double gg = - dot_product(Nu, Xu);
      double d_gg = - dot_product(dNu, Xu) - dot_product(Nu, dXu);
      
      double E = dot_product(Xu, Xu);
      double dE = 2.0 * dot_product(Xu, dXu);
      double F = dot_product(Xu, Xv);
      double dF = dot_product(dXu, Xv) + dot_product(Xu, dXv);
      double G = dot_product(Xv, Xv);
      double dG = 2.0 * dot_product(Xv, dXv);
      double det = (E*G - F*F);
      double d_det = (dE*G + E*dG - 2*F*dF);

      // cout << "Mean Curv " << (ee * G - 2 * ff * F + gg * E) / (2 * det);
      // cout << "Gauss Crv " << (ee * gg - ff * ff) / det;
      dMeanCurv[ib] = 0.5 *
        ((d_ee * G + ee * dG - 2 * (d_ff * F + ff * dF) + d_gg * E + gg * dE) * det 
         - (ee * G - 2 * ff * F + gg * E) * d_det) / (det * det);

      dGaussCurv[ib] = (
        (d_ee * gg + ee * d_gg - 2 * ff * d_ff) * det
        - (ee * gg - ff * ff) * d_det) / (det * det);

      // cout << "D Mean Curv " << dMeanCurv[ib] << endl; 
      // cout << "D Gauss Crv " << dGaussCurv[ib] << endl; 
      }
    }

/*
void 
BruteForceSubdivisionMedialModel::
ComputeBoundaryCurvaturePartial(
  Vec &dMeanCurv, Vec &dGaussCurv, MedialAtom *dAtoms)
{
  // The curvatures are computed at non-crest locations. At the crest, the rules are
  // different and I don't feel like screwing around with them
  for(MedialBoundaryPointIterator bip(this->GetIterationContext());
    !bip.IsAtEnd(); ++bip)
    {
    // Get the medial atom and side
    size_t ib = bip.GetIndex();
    size_t ia = bip.GetAtomIndex();
    size_t side = bip.GetBoundarySide();

    // Just set up an edge walk
    dMeanCurv[ib] = 0.0;
    dGaussCurv[ib] = 0.0;

    double dscope = 0.0;

    // Just set up an edge walk
    EdgeWalkAroundVertex it(&mlAtom, ia);
    if(!it.IsOpen()) 
      {
      // The vector whose magnitude is the mean curvature
      double area = 0.0, d_area = 0.0;
      SMLVec3d kvec(0.0), d_kvec(0.0);

      // The sum of the angles around I
      double sum_theta = 0.0, d_sum_theta = 0.0;

      // Loop around vertex
      for( ; !it.IsAtEnd(); ++it)
        {
        // Naming convention:
        // I - central vertex; J - opposite vertex; A - ahead vertex; B - behind vertex
        // AI - vector from A to I, etc.
        SMLVec3d &I = xAtoms[ia].xBnd[side].X;
        SMLVec3d &J = xAtoms[it.MovingVertexId()].xBnd[side].X;
        SMLVec3d &A = xAtoms[it.VertexIdAhead()].xBnd[side].X;
        SMLVec3d &B = xAtoms[it.VertexIdBehind()].xBnd[side].X;
        SMLVec3d &dI = dAtoms[ia].xBnd[side].X;
        SMLVec3d &dJ = dAtoms[it.MovingVertexId()].xBnd[side].X;
        SMLVec3d &dA = dAtoms[it.VertexIdAhead()].xBnd[side].X;
        SMLVec3d &dB = dAtoms[it.VertexIdBehind()].xBnd[side].X;
        SMLVec3d AI = I - A; SMLVec3d d_AI = dI - dA;
        SMLVec3d AJ = J - A; SMLVec3d d_AJ = dJ - dA; 
        SMLVec3d BI = I - B; SMLVec3d d_BI = dI - dB;
        SMLVec3d BJ = J - B; SMLVec3d d_BJ = dJ - dB;
        SMLVec3d IA = A - I; SMLVec3d d_IA = dA - dI;
        SMLVec3d IJ = J - I; SMLVec3d d_IJ = dJ - dI;
    
        // Compute cot(a)
        double cos_a = dot_product(AI, AJ);
        double sin_a = 
          sqrt(dot_product(AI, AI) * dot_product(AJ, AJ) - cos_a * cos_a);
        double cot_a = cos_a / sin_a;
        double d_cos_a = dot_product(d_AI, AJ) + dot_product(AI, d_AJ);
        double d_sin_a = 
          (dot_product(AI, d_AI) * dot_product(AJ, AJ) + 
           dot_product(AI, AI) * dot_product(AJ, d_AJ) - 
           cos_a * d_cos_a) / sin_a;
        double d_cot_a = (d_cos_a * sin_a - cos_a * d_sin_a) / (sin_a * sin_a);

        // Compute cot(b)
        double cos_b = dot_product(BI, BJ);
        double sin_b = 
          sqrt(dot_product(BI, BI) * dot_product(BJ, BJ) - cos_b * cos_b);
        double cot_b = cos_b / sin_b;
        double d_cos_b = dot_product(d_BI, BJ) + dot_product(BI, d_BJ);
        double d_sin_b = 
          (dot_product(BI, d_BI) * dot_product(BJ, BJ) + 
           dot_product(BI, BI) * dot_product(BJ, d_BJ) - 
           cos_b * d_cos_b) / sin_b;
        double d_cot_b = (d_cos_b * sin_b - cos_b * d_sin_b) / (sin_b * sin_b);

        // Compute the contribution to mean curvature
        kvec += - 0.5 * (cot_a + cot_b) * IJ;
        d_kvec += - 0.5 * ((d_cot_a + d_cot_b) * IJ + (cot_a + cot_b) * d_IJ);
        dscope += d_kvec[0];

        // Compute the area fraction
        double area_chunk = TriangleArea(I, J, A);
        double d_area_chunk = TriangleAreaPartialDerivative(I, J, A, dI, dJ, dA, area_chunk);
        area += area_chunk;
        d_area += d_area_chunk;

        // Compute the contribution to the Gaussian curvature
        double num = dot_product(IA, IJ);
        double d_num = dot_product(IA, d_IJ) + dot_product(d_IA, IJ);
        double den = sqrt(dot_product(IA,IA) * dot_product(IJ,IJ));
        double d_den = (
          dot_product(IA,d_IA) * dot_product(IJ,IJ) +
          dot_product(IA,IA) * dot_product(IJ,d_IJ)) / den;
        double cos_t = num / den;
        double d_cos_t = (d_num * den - num * d_den) / (den * den);
        double theta = acos(cos_t);
        double d_theta = - d_cos_t / sqrt(1.0 - cos_t * cos_t);
        sum_theta += theta;
        d_sum_theta += d_theta;
        }

      // Compute the true values
      SMLVec3d kvec_scale = 3.0 * kvec / area;
      SMLVec3d d_kvec_scale = 3.0 * (d_kvec * area - kvec * d_area) / (area * area);
      double angle_def = sum_theta - 2.0 * M_PI;
      double d_angle_def = d_sum_theta;

      double kvec_scale_mag2 = dot_product(kvec_scale, kvec_scale);
      double d_kvec_scale_mag2 = 2.0 * dot_product(d_kvec_scale, kvec_scale);

      // xMeanCurv[ib] = 0.5 * sqrt(kvec_scale_mag2);
      dMeanCurv[ib] = 0.25 * d_kvec_scale_mag2 / sqrt(kvec_scale_mag2);
      // xGaussCurv[ib] = 3.0 * angle_def / xAreaFrac;
      dGaussCurv[ib] = 3.0 * (d_angle_def * area - angle_def * d_area) / (area * area);


      // dMeanCurv[ib] = dscope;
      // dGaussCurv[ib] = d_sum_theta;
      }
    }
*/
  /*
        // Compute the contribution to mean curvature
        double IAJA = dot_product(J-A, I-A);
        double IBJB = dot_product(J-B, I-B);
        double d_IAJA = dot_product(dJ-dA, I-A) + dot_product(J-A, dI-dA);
        double d_IBJB = dot_product(dJ-dB, I-B) + dot_product(J-B, dI-dB);

        double cos_a = IAJA;
        double d_cos_a = d_IAJA;
        double sin_a = sqrt(dot_product(J-A,J-A)*dot_product(I-A,I-A) - cos_a * cos_a);
        double d_sin_a = 
          ( dot_product(J-A,dJ-dA) * dot_product(I-A,I-A) + 
            dot_product(J-A,J-A) * dot_product(I-A,dI-dA) - 
            cos_a * d_cos_a ) / sin_a;
        double cot_a = cos_a / sin_a;
        double d_cot_a = (d_cos_a * sin_a - cos_a * d_sin_a) / (sin_a * sin_a);

        double cos_b = IBJB;
        double d_cos_b = d_IBJB;
        double sin_b = sqrt(dot_product(J-B,J-B)*dot_product(I-B,I-B) - cos_b * cos_b);
        double d_sin_b = 
          ( dot_product(J-B,dJ-dB) * dot_product(I-B,I-B) + 
            dot_product(J-B,J-B) * dot_product(I-B,dI-dB) - 
            cos_b * d_cos_b ) / sin_b;
        double cot_b = cos_b / sin_b;
        double d_cot_b = (d_cos_b * sin_b - cos_b * d_sin_b) / (sin_b * sin_b);

        xMeanCurvNorm += 0.5 * (cot_a + cot_b) * (I - J);
        dMeanCurvNorm += 0.5 * ((d_cot_a + d_cos_b) * (I - J) + (cot_a + cot_b) * (dI - dJ));

        // Compute the area fraction
        xAreaFrac += TriangleArea(I, J, A);
        dAreaFrac += TriangleAreaPartialDerivative(I, J, A, dI, dJ, dA, xAreaFrac);

        // Compute the contribution to the Gaussian curvature
        double ab = dot_product(A-I,J-I);
        double aa = dot_product(A-I,A-I);
        double bb = dot_product(J-I,J-I);
        double d_ab = dot_product(dA-dI,J-I) + dot_product(A-I,dJ-dI); 
        double d_aa = 2 * dot_product(dA-dI,A-I);
        double d_bb = 2 * dot_product(dJ-dI,J-I);
  
        double nab = sqrt(aa*bb);
        double d_nab = 0.5 * (d_aa * bb + aa * d_bb) / nab;
        double cos_t = ab / nab;
        double d_cos_t = (d_ab * nab - ab * d_nab) / (nab * nab);
        double theta = acos(cos_t);
        double d_theta = -d_cos_t / sqrt(1.0 - cos_t * cos_t);
        xGaussCurv += theta;
        dGaussCurv[ib] += d_theta;
        }

      // Compute the true values
      SMLVec3d xMeanCurvNormFull = 3.0 * xMeanCurvNorm / xAreaFrac;
      SMLVec3d dMeanCurvNormFull =  
        3.0 * (dMeanCurvNorm * xAreaFrac - xMeanCurvNorm * dAreaFrac) 
        / (xAreaFrac * xAreaFrac);
      
      double zz = xMeanCurvNormFull.magnitude();
      xMeanCurv = 0.5 * zz;
      dMeanCurv[ib] = 0.5 * dot_product(xMeanCurvNormFull, dMeanCurvNormFull) / zz;

      double xGaussCurvFull = -(xGaussCurv - 2 * M_PI) * 3.0 / xAreaFrac;
      double dGaussCurvFull = - 3.0 * 
        (dGaussCurv[ib] * xAreaFrac - (xGaussCurv - 2 * M_PI) * dAreaFrac) 
        / (xAreaFrac * xAreaFrac);

      dGaussCurv[ib] = dGaussCurvFull;
      }
    }
    */
}

BruteForceSubdivisionMedialModel::Vec
BruteForceSubdivisionMedialModel::GetHintArray() const
{
  Vec xHint(1, 0.0);
  return xHint;
}

void
BruteForceSubdivisionMedialModel
::PrepareAtomsForVariationalDerivative(
  const Vec &xVariation, MedialAtom *dAtoms) const
{
  size_t i;

  // First, we must label the dependency structure. To begin with, for each atom
  // we must mark whether is is affected by the variation or not. At the same time
  // we can compute the derivative of X and R with respect to the variation 
  for(i = 0; i < mlAtom.nVertices; i++)
    {
    // Set up i-th atom
    MedialAtom &da = dAtoms[i]; da.X.fill(0.0); da.R = 0.0;

    // Initialize the user_data flag to indicate no dependence on the variation
    da.user_data = 0;

    // Compute the weighted sum of the coefficients
    ImmutableSparseMatrix<double>::ConstRowIterator it = mlAtom.weights.Row(i);
    for( ; !it.IsAtEnd(); ++it)
      {
      size_t j = it.Column() << 2; 
      double w = it.Value();

      SMLVec3d dXdVar = w * xVariation.extract(3, j);
      double dRdVar = w * xVariation[j+3];

      da.X += dXdVar;
      da.R += dRdVar;

      // Check if there is a dependency
      if(dXdVar.squared_magnitude() > 0.0 || dRdVar != 0.0)
        da.user_data = 1;
      }
    }

  // Propagate the dependency flag (user_data) to the 2-nd ring of neighbors. This
  // is because the neighbors are used in the computation of Xu and Xv
  for(size_t level = 0; level < 2; level++)
    {
    for(i = 0; i < mlAtom.nVertices; i++)
      {
      if(dAtoms[i].user_data == 1)
        {
        // Create a walk around the vertex
        EdgeWalkAroundVertex it(&mlAtom, i);

        // Flag all neighbors
        for(; !it.IsAtEnd(); ++it)
          dAtoms[it.MovingVertexId()].user_data = 1;
        }
      }
    }

  // Next, precompute the first partial derivatives of X and R
  for(i = 0; i < mlAtom.nVertices; i++) 
    {
    // Initialize the derivatives to zero
    MedialAtom &da = dAtoms[i];
    da.Xu.fill(0.0); da.Xv.fill(0.0); da.Fu = 0; da.Fv = 0;

    // Compute the partial derivatives of F and X
    for(SparseMat::ConstRowIterator it = Wuv.Row(i); !it.IsAtEnd(); ++it)
      {
      MedialAtom &danbr = dAtoms[it.Column()];
      double wu = it.Value().first;
      double wv = it.Value().second;

      da.Xu += wu * danbr.X;
      da.Xv += wv * danbr.X;
      }
    }

  // Now, go through and compute the second derivative terms
  for(i = 0; i < mlAtom.nVertices; i++) 
    {
    // Initialize the derivatives to zero
    MedialAtom &da = dAtoms[i];
    da.Xuu.fill(0.0); da.Xuv.fill(0.0); da.Xvv.fill(0.0);

    // Compute the partial derivatives of F and X
    for(SparseMat::ConstRowIterator it = Wuv.Row(i); !it.IsAtEnd(); ++it)
      {
      MedialAtom &danbr = dAtoms[it.Column()];
      double wu = it.Value().first;
      double wv = it.Value().second;

      da.Xuu += wu * danbr.Xu;
      da.Xuv += 0.5 * (wu * danbr.Xv + wv * danbr.Xu);
      da.Xvv += wv * danbr.Xv;
      }
    }

  // Set the derivative atoms unaffected by the computation to zeros
  for(i = 0; i < mlAtom.nVertices; i++) 
    {
    if(dAtoms[i].user_data == 0)
      dAtoms[i].SetAllDerivativeTermsToZero();
    }
}

void
BruteForceSubdivisionMedialModel
::ComputeAtomGradient(std::vector<MedialAtom *> &xVariations)
{
  size_t i, q;

  // Precompute common terms for atom derivatives
  MedialAtom::DerivativeTerms *dt = 
    new MedialAtom::DerivativeTerms[mlAtom.nVertices];
  for(i = 0; i < mlAtom.nVertices; i++)
    xAtoms[i].ComputeCommonDerivativeTerms(dt[i]);

  // Compute the derivative for each variation
  for(q = 0; q < xVariations.size(); q++)
    {
    // Get the derivative atoms for this variation
    MedialAtom *dAtoms = xVariations[q];

    // Compute F at each atom
    for(i = 0; i < mlAtom.nVertices; i++) 
      if(dAtoms[i].user_data)
        dAtoms[i].F = 2.0 * dAtoms[i].R * xAtoms[i].R;

    // Compute dF/dv for each atom
    for(i = 0; i < mlAtom.nVertices; i++) 
      {
      // Only consider atoms affected by the variation
      if(dAtoms[i].user_data)
        {
        // Get the current atom and the derivative (which we are computing)
        MedialAtom &a = xAtoms[i];
        MedialAtom &da = dAtoms[i];

        // Compute the partial derivatives of Fu
        da.Fu = da.Fv = 0.0;
        for(SparseMat::RowIterator it = Wuv.Row(i); !it.IsAtEnd(); ++it)
          {
          MedialAtom &danbr = dAtoms[it.Column()];
          double wu = it.Value().first;
          double wv = it.Value().second;

          da.Fu += wu * danbr.F;
          da.Fv += wv * danbr.F;
          }
      
        // Compute the metric tensor derivatives of the atom
        xAtoms[i].ComputeMetricTensorDerivatives(da);
        xAtoms[i].ComputeChristoffelDerivatives(da);

        // Compute the derivatives of the boundary nodes
        xAtoms[i].ComputeBoundaryAtomDerivatives(da, dt[i]);
        }
      }

    for(i = 0; i < mlAtom.nVertices; i++) 
      {
      // Only consider atoms affected by the variation
      if(dAtoms[i].user_data)
        {
        // Get the current atom and the derivative (which we are computing)
        MedialAtom &a = xAtoms[i];
        MedialAtom &da = dAtoms[i];

        // Compute the derivative of mean and gauss curvatures
        da.Fuu = da.Fuv = da.Fvv = 0.0;
        for(SparseMat::RowIterator it = Wuv.Row(i); !it.IsAtEnd(); ++it)
          {
          MedialAtom &danbr = dAtoms[it.Column()];
          double wu = it.Value().first;
          double wv = it.Value().second;

          da.Fuu += wu * danbr.Fu;
          da.Fvv += wv * danbr.Fv;
          da.Fuv += 0.5 * (wu * danbr.Fv + wv * danbr.Fu);
          }

        // Compute things in the atom that depend on second derivatives
        }
      }

    // Compute the curvatures
    /*
    for(i = 0; i < mlAtom.nVertices; i++) 
      {
      // Only consider atoms affected by the variation
      if(dAtoms[i].user_data)
        {
        // Get the current atom and the derivative (which we are computing)
        MedialAtom &a = xAtoms[i];
        MedialAtom &da = dAtoms[i];

        // Compute the derivative of mean and gauss curvatures
        SMLVec3d Nu = xLoopScheme.Nu(i, xAtoms);
        SMLVec3d Nv = xLoopScheme.Nv(i, xAtoms);
        SMLVec3d dNu = xLoopScheme.Nu(i, dAtoms);
        SMLVec3d dNv = xLoopScheme.Nv(i, dAtoms);

        // TODO: this should be precomputed, figure out a memory-sane solution
        double ee = - dot_product(Nv, a.Xv);
        double ff = -0.5 * (dot_product(Nu, a.Xv) + dot_product(Nv, a.Xu));
        double gg = - dot_product(Nu, a.Xu);

        // Compute the derivatives of 2nd fundamental form
        double dee = - (dot_product(dNv, a.Xv) + dot_product(Nv, da.Xv));
        double dff = -0.5 * 
          ( dot_product(dNv, a.Xu) + dot_product(Nv, da.Xu) + 
            dot_product(dNu, a.Xv) + dot_product(Nu, da.Xv) );
        double dgg = - (dot_product(dNu, a.Xu) + dot_product(Nu, da.Xu));

        // Compute the derivatives of the mean curvature
        da.xMeanCurv = 0.5 * (
          dee * a.G.xContravariantTensor[0][0] + 
          ee * da.G.xContravariantTensor[0][0] +
          2 * (
            dff * a.G.xContravariantTensor[1][0] +
            ff * da.G.xContravariantTensor[1][0] ) +
          dgg * a.G.xContravariantTensor[1][1] + 
          gg * da.G.xContravariantTensor[1][1]);

        // Compute the derivatives of the Gauss curvature
        da.xGaussCurv = 
          ( (dee * gg + ee * dgg - 2 * ff * dff) * a.G.g - 
            (ee * gg - ff * ff) * da.G.g ) / (a.G.g * a.G.g);

        da.Xuu = dNu;
        da.Xvv = dNv;
        }
      }
      */
    }

  delete dt;
}


void
BruteForceSubdivisionMedialModel::
WriteToRegistry(Registry &R)
{
  SubdivisionMedialModel::WriteToRegistry(R);

  // Set the model subtype
  R["Grid.Model.SolverType"] << "BruteForce";
}


void
BruteForceSubdivisionMedialModel::
ReadFromRegistry(Registry &R)
{
}

