:Begin:
:Function:      lpSeparation
:Pattern:       lpSeparation[A_?MatrixQ,B_?MatrixQ]
:Arguments:     {A, B}
:ArgumentTypes: {Manual}
:ReturnType:    Manual
:End:

:Begin:
:Function:      lpNewFeatureSelector
:Pattern:       lpNewFeatureSelector[A_?MatrixQ, B_?MatrixQ, O_?MatrixQ, alpha_Real:5.0]
:Arguments:     {A, B, O, alpha}
:ArgumentTypes: {Manual}
:ReturnType:    Integer
:End:

:Begin:
:Function:      lpDeleteFeatureSelector
:Pattern:       lpDeleteFeatureSelector[handle_Integer]
:Arguments:     {handle}
:ArgumentTypes: {Integer}
:ReturnType:    Integer
:End:


:Begin:
:Function:      lpRunFeatureSelection
:Pattern:       lpRunFeatureSelection[handle_Integer, lambda_Real, eta_Real, nreps_Integer:10]
:Arguments:     {handle, lambda, eta, nreps}
:ArgumentTypes: {Integer, Real, Real, Integer}
:ReturnType:    Manual
:End:
