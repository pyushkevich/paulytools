:Begin:
:Function:      METISPartition
:Pattern:       METISPartition[xAdjIndex_List,xAdj_List,wVtx_List,wAdj_List,wCluster_List]
:Arguments:     {xAdjIndex, xAdj, wVtx, wAdj, wCluster, 0, 0}
:ArgumentTypes: {IntegerList, IntegerList, IntegerList, IntegerList, RealList, Integer, Integer}
:ReturnType:    Manual
:End:

:Begin:
:Function:      BuildDijkstra
:Pattern:       BuildDijkstra[xAdjIndex_List,xAdj_List,wAdj_List]
:Arguments:     {xAdjIndex, xAdj, wAdj}
:ArgumentTypes: {IntegerList, IntegerList, RealList}
:ReturnType:    Integer
:End:

:Begin:
:Function:      ReleaseDijkstra
:Pattern:       ReleaseDijkstra[handle_Integer]
:Arguments:     {handle}
:ArgumentTypes: {Integer}
:ReturnType:    Integer
:End:

:Begin:
:Function:      DistanceDijkstra
:Pattern:       DistanceDijkstra[handle_Integer,iStart_Integer,iEnd_Integer]
:Arguments:     {handle,iStart,iEnd}
:ArgumentTypes: {Integer,Integer,Integer}
:ReturnType:    Real
:End:
