BeginPackage["Pauly`ClutoLink`",
  "DiscreteMath`Combinatorica`"]

Unprotect[ClutoCluster, ClutoGraph, ClutoGraphUndirected]

Options[ClutoCluster] = 
  Sort[{
    Method -> "Direct",
    ClusterCriterion -> Automatic,
    ClusterSelectStrategy -> "SelectLargest",
    VertexPruning -> -1.0,
    EdgePruning -> -1.0,
    MinConnectedComponent -> 0,
    Trials->10,
    Iterations->10,
    RandomSeed->0
    }]

ClutoMethodValues = { "Direct", "RepeatedBisection", "Agglomeration", "GraphCut" }
ClutoClusterCriterionValues = 
  { Automatic, "i1", "i2", "e1", "g1", "g1p", "h1", "h2" }
ClutoSelectStrategyValues = 
  { "SelectBest", "SelectLargest", "SelectLargestSubspace" }

(* ClutoLink = If[$System == "Microsoft Windows",
  Install["Release/clutoML.exe"],
  Install["clutoML"]]; *)
  
    
ClutoCluster::usage = 
"ClutoCluster[graph,nclust,opts] uses CLUTO to partition the graph
into nclust clusters. The graph should be computed using the function
ClutoGraphUndirected"

ClutoGraphUndirected::usage = 
"ClutoGraphUndirected[nverts, edges, weights]
Creates a sparse graph representation of a graph with a given number of 
vertices, and provided edges and edge weights"

Format[ClutoGraph[ai_, a_, w_]] := 
  SequenceForm["\[SkeletonIndicator] ClutoGraph:<", 
    Length[ai]-1, ", ", Length[a], ", ", 
    ">\[SkeletonIndicator]" ]

Begin["`Private`"]
    
ClutoGraphUndirected[nverts_, edges_, weights_] := Module[{y, y1, y2, y3, y4},
  y = Union[Transpose[{edges, weights}], Transpose[{Reverse /@ edges, weights}]];
  y1 = Split[Sort@Join[y[[All, 1, 1]], Range[nverts]]];
  y2 = FoldList[Plus, 1, Length /@ y1 - 1];
  y3 = y[[All, 1, 2]];
  y4 = y[[All, 2]];
  ClutoGraph[y2, y3, y4]];

ClutoCluster[g_ClutoGraph,nclust_Integer,opts___?OptionQ] := 
Module[{meth,crit,strat,vprun,eprun,mincmp,ntrials,niter,seed},

  (* Get the options *)
  {meth,crit,strat,vprun,eprun,mincmp,ntrials,niter,seed} = 
    { Method, ClusterCriterion, ClusterSelectStrategy, VertexPruning, 
      EdgePruning, MinConnectedComponent, 
      Trials, Iterations, RandomSeed } /.  Flatten[{opts,Options[ClutoCluster]}];

      Print[ {meth,crit,strat,vprun,eprun,mincmp,ntrials,niter,seed}//FullForm ]; 

  (* Call the cluto partitioning method *)
  If[meth == "Direct",

    (* Get criterion *)
    crit = If[crit === Automatic, "i1", ToString[crit] ] // Evaluate;
    
    Return @ Global`CLUTOClusterDirect[
      g[[1]], g[[2]], g[[3]], ToString@crit, ntrials, niter, seed, nclust ][[1]]];

  If[meth == "GraphCut",

    (* Call the method *)
    Return @ Global`CLUTOGraphClusterRB[
      g[[1]], g[[2]], g[[3]],40, eprun, vprun, mincmp, ntrials, seed, 
      strat, nclust][[2]]];
]    
  
End[]

Protect[ClutoCluster, ClutoGraph, ClutoGraphUndirected]

EndPackage[]
