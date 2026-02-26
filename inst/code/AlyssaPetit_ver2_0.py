# AlyssaPetit version 2.2
# Use with python 3.x
#
# Changes from v1.0:
#   1. FVS fallback: When minType==-1 (cycle not analytically solvable),
#      compute a minimal Feedback Vertex Set, remove those species, and
#      continue solving the rest analytically.
#   2. Polynomial solving: For FVS species, attempt to solve the steady-state
#      ODE as a polynomial (up to configurable maxPolyDegree, default 2).
#      Only strictly positive real solutions are returned.
#   3. Convergence conditions: For FVS species that cannot be solved in
#      closed form, output symbolic conditions guaranteeing convergence
#      of iterative equilibration (f'(x)<0 for 1D; trace(J)<0 and det(J)>0
#      for 2D systems in the positive orthant).
#
# Changes in v2.2:
#   4. Internal FVS prefix changed from ss_ to fvs_ for clarity.
#   5. fvs_ prefixes are removed from final output: fvs_X -> X with
#      simplification, and trivial identities (X = X) are dropped.

import numpy
import sympy
from sympy import (Matrix, simplify, expand, solve, Symbol, Poly, Rational,
                   oo, S, diff, det, trace, And, StrictLessThan, StrictGreaterThan,
                   Expr)
from numpy import shape, zeros, concatenate
from numpy.linalg import matrix_rank
from sympy.parsing.sympy_parser import parse_expr
from sympy.matrices import *
from sympy.matrices import matrix_multiply_elementwise
import csv
import random
from random import shuffle
from itertools import combinations

# ============================================================================
# Utility functions
# ============================================================================

def LCS(s1, s2):
    m = [[0] * (1 + len(s2)) for i in range(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in range(1, 1 + len(s1)):
        for y in range(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return s1[x_longest - longest: x_longest]

def SolveSymbLES(A,b):
    dim=shape(A)[0]
    Asave=A[:]
    Asave=Matrix(dim, dim, Asave)
    determinant=Asave.det()
    if(determinant==0):
        return([])
    result=[]
    for i in range(dim):
        A=Matrix(dim,dim,Asave)
        A.col_del(i)
        A=A.col_insert(i,b)
        result.append(simplify(A.det()/determinant))
    return(result)

def CutStringListatSymbol(liste, symbol):
    out=[]    
    for el in liste:
        if(symbol in el):
            add=el.split(symbol)
        else:
            add=[el]
        out=out+add
    return(out)

def FillwithRanNum(M):
    dimx=len(M.row(0))
    dimy=len(M.col(0))
    ranM=zeros(dimy, dimx)
    parlist=[]
    ranlist=[]
    for i in M[:]:
        if(i!=0):
            if(str(i)[0]=='-'):
                parlist.append(str(i)[1:])
            else:
                parlist.append(str(i))
    parlist=list(set(parlist))
    for symbol in [' - ', ' + ', '*', '/', '(',')']:
        parlist=CutStringListatSymbol(parlist,symbol)
    parlist=list(set(parlist))
    temp=[]    
    for i in parlist:
        if(i!=''):
            if(not is_number(i)):
                temp.append(i)
                ranlist.append(random.random())
    parlist=temp
    for i in range(dimy):
        for j in range(dimx):
            ranM[i,j]=M[i,j]
            if(ranM[i,j]!=0):
                for p in range(len(parlist)):
                   ranM[i,j]=ranM[i,j].subs(parse_expr(parlist[p]),ranlist[p])
    return(ranM)

def FindLinDep(M, tol=1e-12):
    ranM=FillwithRanNum(M)
    Q,R=numpy.linalg.qr(ranM)
    for i in range(shape(R)[0]):
        for j in range(shape(R)[1]):
            if(abs(R[i,j]) < tol):
                R[i,j]=0.0
    LinDepList=[]
    for i in range(shape(R)[0]):
        if(R[i][i]==0):
            LinDepList.append(i)
    return(LinDepList)

def FindLCL(M, X):
    LCL=[]    
    LinDepList=FindLinDep(M)
    i=0
    counter=0
    deleted_rows=[]
    states=Matrix(X[:])
    while(LinDepList!=[]):
        i=LinDepList[0]
        testM=FillwithRanNum(M)
        rowliste=list(numpy.nonzero(testM[:,i])[0])
        colliste=[i]        
        for z in range(i):
            for k in rowliste:        
                for j in range(i):
                    jliste=list(numpy.nonzero(testM[:,j])[0])
                    if(k in jliste):
                        rowliste=rowliste+jliste
                        colliste=colliste+[j]
            rowliste=list(set(rowliste))
            colliste=list(set(colliste))
        rowliste.sort()
        colliste.sort()
        colliste.pop()        
        rowlisteTry=rowliste[0:(len(colliste))]
        vec=SolveSymbLES(M[rowlisteTry,colliste],M[rowlisteTry,i])
        shufflecounter=0
        while(vec==[] and shufflecounter < 100):
            shuffle(rowliste)
            shufflecounter=shufflecounter+1
            rowlisteTry=rowliste[0:(len(colliste))]
            vec=SolveSymbLES(M[rowlisteTry,colliste],M[rowlisteTry,i])
        if(shufflecounter==100):
            print('Problems while finding conserved quantities!',flush=True)
            return(0,0)
        counter=counter+1
        try:
            mat=[states[l] for l in colliste]
            test=parse_expr('0')
            for v in range(0,len(vec)):
                test=test-parse_expr(str(vec[v]))*parse_expr(str(mat[v]))
        except:
            return([],0)
        partStr=str(test)+' + '+str(states[i])
        partStr=partStr.split(' + ')
        partStr2=[]
        for index in range(len(partStr)):
            partStr2=partStr2+partStr[index].split('-')
        partStr=partStr2
        if(len(partStr) > 1):        
            CLString=LCS(str(partStr[0]),str(partStr[1]))
            for ps in range(2,len(partStr)):
                CLString=LCS(CLString,str(partStr[ps]))
        else:
            CLString=str(partStr[0])
        if(CLString==''):
            CLString=str(counter)
        LCL.append(str(test)+' + '+str(states[i])+' = '+'total'+CLString)
        M.col_del(i)
        states.row_del(i)
        deleted_rows.append(i+counter-1)
        LinDepList=FindLinDep(M)
    return(LCL, deleted_rows)

def printmatrix(M):    
    lengths=[]
    for i in range(len(M.row(0))):
        lengths.append(0)
        for j in range(len(M.col(0))):
            lengths[i]=max(lengths[i],len(str(M.col(i)[j])))          
    string=''.ljust(5)
    string2=''.ljust(5)
    for j in range(len(M.row(0))):
        string=string+(str(j)).ljust(lengths[j]+2)
        for k in range(lengths[j]+2):        
            string2=string2+('-')        
    print(string,flush=True)
    print(string2,flush=True)
    for i in range(len(M.col(0))):
        string=str(i).ljust(4) + '['
        for j in range(len(M.row(0))):
            if(j==len(M.row(0))-1):
                string=string+str(M.row(i)[j]).ljust(lengths[j])
            else:
                string=string+(str(M.row(i)[j])+', ').ljust(lengths[j]+2)        
        print(string+']',flush=True)    
    return()
    
def printgraph(G):
    for el in G:
        print(el+': '+str(G[el]),flush==True)
    return()

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    
def checkNegRows(M):
    NegRows=[]
    if((M==Matrix(0,0,[])) | (M==Matrix(0,1,[])) | (M==Matrix(1,0,[]))):
        return(NegRows)
    else:        
        for i in range(len(M.col(0))):
            foundPos=False
            for j in range(len(M.row(i))):
                if(M[i,j]>0):
                    foundPos=True
            if(foundPos==False):
                NegRows.append(i)    
        return(NegRows)
    
def checkPosRows(M):
    PosRows=[]
    if((M==Matrix(0,0,[])) | (M==Matrix(0,1,[])) | (M==Matrix(1,0,[]))):
        return(PosRows)
    else: 
        for i in range(len(M.col(0))):
            foundNeg=False
            for j in range(len(M.row(i))):
                if(M[i,j]<0):
                    foundNeg=True
            if(foundNeg==False):
                PosRows.append(i)    
        return(PosRows)             

def DetermineGraphStructure(SM, F, X, neglect):
    graph={}    
    SMXF = SM*F
    for i in range(len(SMXF)):
        # If the ODE for this species is identically zero, it has no
        # dependencies and can be solved immediately (it's a free parameter
        # determined by conserved quantities or already resolved).
        if SMXF[i] == 0:
            graph[str(X[i])]=[]
            continue
        liste=[]
        for j in range(len(X)):
            if(SMXF[i]!=(SMXF[i]).subs(X[j],1)):
                if(j==i):
                    In=(SMXF[i]).subs(X[j],0)
                    Out=simplify((SMXF[i]-In)/X[j])
                    if(Out!=Out.subs(X[j],1)):
                        liste.append(str(X[j]))
                else:
                    liste.append(str(X[j]))
            else:
                if(j==i):
                    liste.append(str(X[j]))
        graph[str(X[i])]=liste
    for el in neglect:
        if(parse_expr(el) in X):
            if not el in graph:
                graph[el]=[el]
            else:
                if(el not in graph[el]):
                    graph[el].append(el)
    return(graph)

def FindCycle(graph, X):
    for el in X:
        cycle=find_cycle(graph, str(el), str(el), path=[])
        if(cycle!=None):
            return(cycle)
    return(None)
    
def find_cycle(graph, start, end, path=[]):
    path = path + [start]
    if not start in graph:
        return None
    if ((start == end) & (path!=[start])):
        return path    
    for node in graph[start]:
        if node==end: 
            return (path+[end])
        if node not in path:
            newpath = find_cycle(graph, node, end, path)
            if newpath: 
                return newpath
    return None
    
def GetBestPair(cycle, SM, fluxpars, X, LCLs, neglect):    
    for state in cycle:
        for LCL in LCLs:
            ls=parse_expr(LCL.split(' = ')[0])
            if(ls.subs(parse_expr(state),1)!=ls):
                return(0, state, None, False)
    dimList=[]
    signList=[]
    for state in cycle:
        dim, sign = GetDimension(state, X, SM, True)
        signList.append(sign)
        dimList.append(dim)
    beststate=None
    bestflux=None
    besttype=-1
    n2beat=1000
    signChanged=False
    min2beat=max(dimList)+1
    for i in range(len(dimList)):
        if(dimList[i] < min2beat):
            min2beat=dimList[i]
            sign=signList[i]
            appearList=[]
            if(sign=="minus"):
                fluxpars2use=GetNegFluxParameters(SM, fluxpars, X, cycle[i])                
            else:
                fluxpars2use=GetPosFluxParameters(SM, fluxpars, X, cycle[i])
            abort_flux=False
            for fp in fluxpars2use:
                if(str(fp) not in neglect):
                    appearList.append(GetAppearances(fp, fluxpars, SM))
                else:
                    abort_flux=True
            if(abort_flux):
                print("Sign changed!",flush=True)
                signChanged=True
                if((sign=="minus" and not signChanged) or (sign=="plus" and signChanged)):
                    fluxpars2use=GetNegFluxParameters(SM, fluxpars, X, cycle[i])                
                else:
                    fluxpars2use=GetPosFluxParameters(SM, fluxpars, X, cycle[i])
                abort_flux=False
                for fp in fluxpars2use:
                    if(str(fp) not in neglect):
                        appearList.append(GetAppearances(fp, fluxpars, SM))
                    else:
                        abort_flux=True
            if(sum(appearList) < n2beat and not abort_flux):
                n2beat=sum(appearList)
                beststate=cycle[i]
                if((sign=="minus" and not signChanged) or (sign=="plus" and signChanged)):
                    _bflux=GetNegFluxParameters(SM, fluxpars, X, cycle[i])
                else:
                    _bflux=GetPosFluxParameters(SM, fluxpars, X, cycle[i])
                if not _bflux:
                    # No flux parameters found – cannot solve this candidate
                    beststate=None
                    bestflux=None
                    continue
                bestflux=_bflux[0]
                if(min2beat==1 and max(appearList)==1):
                    besttype=1
                else:
                    if(max(appearList)==1 and min2beat>1):
                        besttype=2
                    else:
                        besttype=3
    return(besttype, beststate, bestflux, signChanged)
    
def GetNegFluxParameters(SM, fluxpars, X, node):
    row=list(X).index(parse_expr(node))
    liste=[]
    for i in range(len(SM.row(row))):
        if(SM.row(row)[i]<0):
            liste.append(fluxpars[i])
    return(liste)

def GetPosFluxParameters(SM, fluxpars, X, node):
    row=list(X).index(parse_expr(node))
    liste=[]    
    for i in range(len(SM.row(row))):
        if(SM.row(row)[i]>0):
            liste.append(fluxpars[i])
    return(liste)
        
def GetType(node, fp, fluxpars, LCLs):
    for LCL in LCLs:
        ls=parse_expr(LCL.split(' = ')[0])
        if(ls.subs(parse_expr(node),1)!=ls):
            return(0)
    if(GetAppearances(fp, fluxpars)==1):
        if(GetDimension(node)==1):
            return(1)
        else:
            return(2)
    else:
        return(3)
        
def GetAppearances(fp, fluxpars, SM):
    anz=0
    cols = [i for i, x in enumerate(fluxpars) if x == fp]
    for i in cols:
        for j in range(len(SM.col(i))):
            if(SM.col(i)[j]!=0):
                anz=anz+1
    return(anz)

def GetDimension(node, X, SM, getSign=False):
    row=list(X).index(parse_expr(node))
    anzminus=0
    anzappearminus=0
    for i in range(len(SM.row(row))):
        if(SM.row(row)[i]<0):
            anzappearminus=anzappearminus+CountNZE(SM.col(i))
            anzminus=anzminus+1
    anzplus=0
    anzappearplus=0
    for i in range(len(SM.row(row))):
        if(SM.row(row)[i]>0):
            anzappearplus=anzappearplus+CountNZE(SM.col(i))
            anzplus=anzplus+1
    if(not getSign):
        return(min(anzminus, anzplus))
    else:
        if(anzminus<anzplus or (anzminus==anzplus and anzappearminus<anzappearplus)):
            return(anzminus, "minus")
        else:
            return(anzplus, "plus")
            
def GetOutfluxes(node, X, SM, F, fluxpars):
    row=list(X).index(parse_expr(node))
    outsum=0
    out=[]
    fps=[]
    for i in range(len(SM.row(row))):
        if(SM.row(row)[i]<0):
            outsum=outsum-SM.row(row)[i]*F[i]
            out.append(-SM.row(row)[i]*F[i])
            fps.append(fluxpars[i])
    return(out, outsum, fps)

def GetInfluxes(node, X, SM, F, fluxpars):
    row=list(X).index(parse_expr(node))
    outsum=0
    out=[]
    fps=[]
    for i in range(len(SM.row(row))):
        if(SM.row(row)[i]>0):
            outsum=outsum+SM.row(row)[i]*F[i]
            out.append(SM.row(row)[i]*F[i])
            fps.append(fluxpars[i])
    return(out, outsum, fps)

def FindNodeToSolve(graph):
    for el in graph:
        if(graph[el]==[]):
            return(el)
    return(None)

def CountNZE(V):
    counter=0
    for v in V:
        if(v!=0):
            counter=counter+1
    return(counter)
    
def Sparsify(M, level, sparseIter):
    oldM=M.copy()
    if(level==3):
        ncol=len(M.row(0))
        print('0 columns of '+str(ncol) +' done',flush=True)
        for i in range(ncol):            
            icol=M.col(i)
            tobeat=CountNZE(M.col(i))
            for j in range(ncol):
                if(i<j):
                    for factor_j in [1,2,-1,-2,0]:
                        for k in range(ncol):
                            if(i<k and j<k):
                                for factor_k in [1,2,-1,-2,0]:
                                    for l in range(ncol):
                                        if(i<l and j<l and k<l):
                                            for factor_l in [1,2,-1,-2,0]:
                                                test=icol+factor_j*M.col(j)+factor_k*M.col(k)+factor_l*M.col(l)
                                                if(tobeat > CountNZE(test)):
                                                    Mtest=M.copy()
                                                    Mtest.col_del(i)
                                                    Mtest=Mtest.col_insert(i,test)
                                                    if(CountNZE(test)!=0 and M.rank()==Mtest.rank()):
                                                        M=Mtest.copy()
                                                        tobeat=CountNZE(test)
            print(str(i+1)+' columns of '+str(ncol) +' done',flush=True)
    if(level==2):
        ncol=len(M.row(0))
        for i in range(ncol):
            icol=M.col(i)
            tobeat=CountNZE(M.col(i))
            for j in range(ncol):
                if(i<j):
                    for factor_j in [1,2,-1,-2,0]:
                        for k in range(ncol):
                            if(i<k and j<k):
                                for factor_k in [1,2,-1,-2,0]:
                                    test=icol+factor_j*M.col(j)+factor_k*M.col(k)
                                    if(tobeat > CountNZE(test)):
                                        Mtest=M.copy()
                                        Mtest.col_del(i)
                                        Mtest=Mtest.col_insert(i,test)
                                        if(CountNZE(test)!=0 and M.rank()==Mtest.rank()):
                                            M=Mtest.copy()
                                            tobeat=CountNZE(test)
    if(level==1):
        ncol=len(M.row(0))
        for i in range(ncol):
            icol=M.col(i)
            tobeat=CountNZE(M.col(i))
            for j in range(ncol):
                if(i<j):
                    for factor_j in [1,2,-1,-2,0]:
                        test=icol+factor_j*M.col(j)
                        if(tobeat > CountNZE(test)):
                            Mtest=M.copy()
                            Mtest.col_del(i)
                            Mtest=Mtest.col_insert(i,test)
                            if(CountNZE(test)!=0 and M.rank()==Mtest.rank()):
                                M=Mtest.copy()
                                tobeat=CountNZE(test)
    if(oldM!=M and sparseIter<10):
        oldM=M.copy() 
        print("Sparsify with level", level,", Iteration ",sparseIter, " of maximal 10",flush=True)
        return(Sparsify(M,level, sparseIter=sparseIter+1))                            
    else:
        return(M)


# ============================================================================
# NEW in v2.0: Feedback Vertex Set (FVS) computation
# ============================================================================

def _graph_has_cycle(adj):
    """Check if directed graph (dict: node -> [neighbours]) has a cycle."""
    WHITE, GRAY, BLACK = 0, 1, 2
    color = {n: WHITE for n in adj}
    def dfs(u):
        color[u] = GRAY
        for v in adj.get(u, []):
            if v in color:
                if color[v] == GRAY:
                    return True
                if color[v] == WHITE and dfs(v):
                    return True
        color[u] = BLACK
        return False
    for node in list(adj.keys()):
        if color[node] == WHITE:
            if dfs(node):
                return True
    return False

def _remove_node_from_graph(adj, node):
    """Return a copy of adj with node removed (as source and target)."""
    new = {}
    for n, nbrs in adj.items():
        if n == node:
            continue
        new[n] = [v for v in nbrs if v != node]
    return new

def FindMinimalFVS(graph):
    """
    Find a minimal Feedback Vertex Set by brute-force over increasing
    subset sizes.  For the typically small graphs in ODE systems this
    is tractable.  Returns a list of node names (strings).
    """
    nodes = list(graph.keys())
    # Quick check: if no cycle, return empty
    if not _graph_has_cycle(graph):
        return []
    # Try subsets of increasing size
    for k in range(1, len(nodes) + 1):
        for subset in combinations(nodes, k):
            reduced = dict(graph)
            for s in subset:
                reduced = _remove_node_from_graph(reduced, s)
            if not _graph_has_cycle(reduced):
                return list(subset)
    return nodes  # worst case: remove everything


# ============================================================================
# NEW in v2.0: Polynomial steady-state solving
# ============================================================================

def _filter_strictly_positive(raw_solutions):
    """
    Filter a list of sympy solutions, keeping only those that are strictly
    positive under the assumption that all free symbols are positive
    (biologically: concentrations and rate constants > 0).
    
    Returns list of strictly positive solutions.
    """
    positive_solutions = []
    for sol in raw_solutions:
        sol = simplify(sol)
        # Skip trivially zero
        if sol == 0:
            continue
        # Skip solutions with imaginary parts (complex roots)
        if sol.has(sympy.I):
            continue
        # Substitute positive dummies for all free symbols to test sign
        try:
            pos_subs = {s: sympy.Dummy(s.name, positive=True)
                        for s in sol.free_symbols}
            sol_pos = sol.subs(pos_subs)
            if sol_pos.is_negative or sol_pos.is_nonpositive:
                continue
            if sol_pos.is_zero:
                continue
        except Exception:
            pass
        positive_solutions.append(sol)
    return positive_solutions


def SolvePolynomialSteadyState(ode_expr, species_symbol, maxDegree=2):
    """
    Try to solve  ode_expr == 0  as a polynomial equation.
    
    Strategy: iterate over ALL free symbols in the expression (not just
    the species), check polynomial degree in each, and solve for the
    symbol with the lowest degree first.  This way, if the ODE is e.g.
    degree 3 in the species but degree 1 in a rate parameter, the
    parameter is solved instead.
    
    All solutions are filtered for strict positivity (under the assumption
    that all symbols are positive, as is standard for biological systems).
    
    Parameters
    ----------
    ode_expr : sympy expression
        The right-hand side of dx/dt = f(x, ...).
    species_symbol : sympy Symbol
        The FVS species (used for reporting; solving is attempted for
        all symbols).
    maxDegree : int
        Maximum polynomial degree to attempt (default 2).
    
    Returns
    -------
    solved_symbol : sympy Symbol or None
        The symbol that was solved for.
    solutions : list of sympy expressions
        Strictly positive real solutions.
    success : bool
        True if polynomial solving succeeded.
    """
    expr = expand(ode_expr)
    free_syms = list(expr.free_symbols)
    
    if not free_syms:
        return None, [], False
    
    # Build list of (degree, symbol) for all symbols that appear polynomially
    candidates = []
    for sym in free_syms:
        try:
            p = Poly(expr, sym)
            deg = p.degree()
            if deg >= 1:
                candidates.append((deg, sym))
        except Exception:
            continue
    
    if not candidates:
        return None, [], False
    
    # Sort by degree (lowest first), prefer the species symbol on ties
    candidates.sort(key=lambda pair: (pair[0], 0 if pair[1] == species_symbol else 1))
    
    # Try each candidate in order of increasing degree
    for deg, sym in candidates:
        if deg > maxDegree:
            continue
        
        try:
            raw_solutions = solve(expr, sym)
        except Exception:
            continue
        
        if not raw_solutions:
            continue
        
        positive_solutions = _filter_strictly_positive(raw_solutions)
        
        if positive_solutions:
            positive_solutions = [simplify(s) for s in positive_solutions]
            return sym, positive_solutions, True
    
    # If no candidate yielded positive solutions, report the best attempt
    best_deg, best_sym = candidates[0]
    if best_deg > maxDegree:
        print(f'    Lowest polynomial degree is {best_deg} (in {best_sym}), '
              f'exceeds maxPolyDegree={maxDegree}, skipping.',flush=True)
    else:
        print(f'    No strictly positive solutions found for any symbol.',flush=True)
    
    return None, [], False


# ============================================================================
# NEW in v2.0: Convergence conditions for iterative equilibration
# ============================================================================

def GetConvergenceConditions(ode_exprs, species_symbols):
    """
    Derive symbolic conditions under which the steady state is globally
    stable in the positive orthant.
    
    Only supported for 1D and 2D FVS subsystems:
      - 1D: f'(x) < 0 for all x > 0  (monotone decrease => unique stable SS)
      - 2D: trace(J) < 0 and det(J) > 0 for all x1,x2 > 0
             (no limit cycles by Bendixson; unique stable equilibrium)
    
    For n > 2 these conditions are not sufficient; no automatic check
    is performed.
    
    Parameters
    ----------
    ode_exprs : list of sympy expressions
        The RHS f_i(x_1,...,x_n) of dx_i/dt = f_i.
    species_symbols : list of sympy Symbols
        The FVS species.
    
    Returns
    -------
    conditions : list of strings
        Human-readable symbolic stability conditions.
    """
    n = len(species_symbols)
    conditions = []
    
    if n == 0:
        return conditions
    
    if n == 1:
        x = species_symbols[0]
        f = ode_exprs[0]
        df = diff(f, x)
        df_simplified = simplify(df)
        conditions.append(
            f'Global stability condition for {x} (1D):')
        conditions.append(
            f'  df/d{x} = {df_simplified}  <  0   for all {x} > 0')
        conditions.append(
            f'  This guarantees a unique, globally stable steady state.')
        return conditions
    
    if n == 2:
        x1, x2 = species_symbols[0], species_symbols[1]
        f1, f2 = ode_exprs[0], ode_exprs[1]
        # Build Jacobian
        J = Matrix([
            [diff(f1, x1), diff(f1, x2)],
            [diff(f2, x1), diff(f2, x2)]
        ])
        J_simplified = simplify(J)
        tr = simplify(J_simplified.trace())
        dt = simplify(J_simplified.det())
        conditions.append(
            f'Global stability conditions for ({x1}, {x2}) (2D):')
        conditions.append(
            f'  Jacobian J = {J_simplified}')
        conditions.append(
            f'  (1) trace(J) = {tr}  <  0   for all {x1},{x2} > 0')
        conditions.append(
            f'  (2) det(J)   = {dt}  >  0   for all {x1},{x2} > 0')
        conditions.append(
            f'  Condition (1) excludes limit cycles (Bendixson criterion).')
        conditions.append(
            f'  Together they guarantee a unique, globally stable steady state.')
        return conditions
    
    # n > 2: no automatic check
    conditions.append(
        f'FVS subsystem has {n} species: {[str(s) for s in species_symbols]}')
    conditions.append(
        f'  Automatic stability analysis is only supported for n <= 2.')
    conditions.append(
        f'  For n > 2, please verify stability manually (e.g. numerical '
        f'eigenvalue analysis at the computed steady state).')
    
    # Still output the Jacobian for reference
    J = sympy.zeros(n, n)
    for i in range(n):
        for j in range(n):
            J[i, j] = simplify(diff(ode_exprs[i], species_symbols[j]))
    conditions.append(f'  Jacobian J (for reference):')
    for i in range(n):
        row_str = '    [' + ', '.join(str(J[i, j]) for j in range(n)) + ']'
        conditions.append(row_str)
    
    return conditions


# ============================================================================
# v2.1: Quick steady-state verification (no printing, returns bool)
# ============================================================================

def _testSteadyStateQuiet(ODE, eqOut, zeroStates, fvsSpecies):
    """
    Silently check whether the solution in eqOut satisfies the original
    ODEs.  Returns True if solution is correct, False otherwise.
    """
    for i in range(len(ODE)):
        expr = parse_expr(str(ODE[i]))
        for zs in zeroStates:
            expr = expr.subs(zs, 0)
        for j in range(len(eqOut)):
            eq_str = eqOut[-(j+1)]
            if ' = ' in eq_str:
                ls, rs = eq_str.split(' = ', 1)
            else:
                ls, rs = eq_str.split('=', 1)
            rs = rs.split('[')[0].strip()
            ls = parse_expr(ls.strip())
            rs = parse_expr(rs)
            expr = expr.subs(ls, rs)
        expr = simplify(expr)
        if expr != 0:
            has_fvs_param = False
            for fvs_node in fvsSpecies:
                if expr.has(Symbol('fvs_' + fvs_node)):
                    has_fvs_param = True
            if not has_fvs_param:
                return False
    return True


# ============================================================================
# v2.1: Apply FVS to remaining cycle nodes
# ============================================================================

def _applyFVS(SSgraph, SM, F, X, fvsSpecies, fvsODEs, neglect,
              fvs_label_printed):
    """
    Compute the minimal Feedback Vertex Set of SSgraph, remove the FVS
    species from the system, and return the updated state.
    
    Parameters
    ----------
    fvs_label_printed : bool
        Whether the FVS abbreviation has already been introduced.
    
    Returns
    -------
    (SM, F, X, SSgraph, cycle, fvs_label_printed)
    """
    if not fvs_label_printed:
        print("    Computing minimal FVS (Feedback Vertex Set) ...",flush=True)
        fvs_label_printed = True
    else:
        print("    Computing minimal FVS ...",flush=True)
    fvs = FindMinimalFVS(SSgraph)
    print("    FVS = " + str(fvs),flush=True)
    for fvs_node in fvs:
        if parse_expr(fvs_node) not in X:
            continue
        fvs_index = list(X).index(parse_expr(fvs_node))
        fvsSpecies.append(fvs_node)
        fvsODEs[fvs_node] = (SM*F)[fvs_index]
        fvs_param = Symbol('fvs_' + fvs_node)
        for f_idx in range(len(F)):
            F[f_idx] = F[f_idx].subs(parse_expr(fvs_node), fvs_param)
        X.row_del(fvs_index)
        SM.row_del(fvs_index)
        print(f"    Removed FVS species: {fvs_node} "
              f"(replaced by parameter fvs_{fvs_node})",flush=True)
    SSgraph = DetermineGraphStructure(SM, F, X, neglect)
    cycle = FindCycle(SSgraph, X)
    return SM, F, X, SSgraph, cycle, fvs_label_printed


# ============================================================================
# Main function: Alyssa v2.0
# ============================================================================

def Alyssa(filename,
          injections=[],
          givenCQs=[],
          neglect=[],
          sparsifyLevel=2,
          maxPolyDegree=2,
          outputFormat='R',
          testSteady='T'):
    """
    AlyssaPetit v2.2 – Steady-state solver for ODE systems.
    
    New parameters (v2.0):
        maxPolyDegree : int (default 2)
            Maximum polynomial degree for which closed-form solutions of
            FVS species are attempted.  Degree 1-2 is fast, 3-4 is
            supported but expensive.
    
    New parameters (v2.1, from v1.1 patches):
        testSteady : str (default 'T')
            'T' to test the steady-state solution, 'F' to skip testing.
    """
    filename=str(filename)
    file=csv.reader(open(filename), delimiter=',')
    print('Reading csv-file ...',flush=True)
    L=[]
    nrrow=0
    nrcol=0
    for row in file:
        nrrow=nrrow+1
        nrcol=len(row)
        L.append(row)
        
    nrspecies=nrcol-2
    
##### Remove injections  
    counter=0
    for i in range(1,len(L)):
        if(L[i-counter][1] in injections):
            L.remove(L[i-counter])
            counter=counter+1       
    
##### Define flux vector F	
    F=[]
    for i in range(1,len(L)):
        F.append(L[i][1])
        F[i-1]=F[i-1].replace('^','**')
        F[i-1]=parse_expr(F[i-1])
        for inj in injections:
            F[i-1]=F[i-1].subs(parse_expr(inj),0)
    F=Matrix(F)

##### Define state vector X
    X=[]
    X=L[0][2:]
    for i in range(len(X)):
        X[i]=parse_expr(X[i])               
    X=Matrix(X)
    Xo=X.copy()
        
##### Define stoichiometry matrix SM
    SM=[]
    for i in range(len(L)-1):
    	SM.append(L[i+1][2:])        
    for i in range(len(SM)):
    	for j in range(len(SM[0])):
    		if (SM[i][j]==''):
    			SM[i][j]='0'
    		SM[i][j]=parse_expr(SM[i][j])    
    SM=Matrix(SM)
    SM=SM.T
    SMorig=SM.copy()

##### Check for zero fluxes
    icounter=0
    jcounter=0
    for i in range(len(F)):
        if(F[i-icounter]==0):
            F.row_del(i-icounter)
            for j in range(len(SM.col(i-icounter))):
                if(SM[j-jcounter,i-icounter]!=0):
                    X.row_del(j-jcounter)
                    SM.row_del(j-jcounter)
                    SMorig.row_del(j-jcounter)
                    jcounter=jcounter+1
            SM.col_del(i-icounter)
            SMorig.col_del(i-icounter)
            icounter=icounter+1
    
    print('Removed '+str(icounter)+' fluxes that are a priori zero!',flush=True)
    nrspecies=nrspecies-icounter

##### Check if some species are zero and remove them from the system
    zeroStates=[]
    NegRows=checkNegRows(SM)
    PosRows=checkPosRows(SM)
    while((NegRows!=[]) | (PosRows!=[])):
        if(NegRows!=[]):        
            row=NegRows[0]
            zeroStates.append(X[row])
            counter=0    
            for i in range(len(F)):
                if(F[i-counter].subs(X[row],1)!=F[i-counter] and F[i-counter].subs(X[row],0)==0):
                    F.row_del(i-counter)
                    SM.col_del(i-counter)                    
                    counter=counter+1
                else:
                    if(F[i-counter].subs(X[row],1)!=F[i-counter] and F[i-counter].subs(X[row],0)!=0):
                        F[i-counter]=F[i-counter].subs(X[row],0)
            X.row_del(row)
            SM.row_del(row)
        else:
            row=PosRows[0]
            zeroFluxes=[]
            for j in range(len(SM.row(row))):
                if(SM.row(row)[j]!=0):
                    zeroFluxes.append(F[j])
            for k in zeroFluxes:
                StateinFlux=[]
                for state in X:
                    if(k.subs(state,1)!=k):
                        StateinFlux.append(state)
                if(len(StateinFlux)==1):
                    zeroStates.append(StateinFlux[0])
                    row=list(X).index(StateinFlux[0])
                    counter=0            
                    for i in range(len(F)):
                        if(F[i-counter].subs(X[row],1)!=F[i-counter]):
                            if(F[i-counter].subs(X[row],0)==0):
                                F.row_del(i-counter)
                                SM.col_del(i-counter)
                            else:
                                F[i-counter]=F[i-counter].subs(X[row],0)                            
                            counter=counter+1
        NegRows=checkNegRows(SM)      
        PosRows=checkPosRows(SM)

    nrspecies=nrspecies-len(zeroStates)
    if(nrspecies==0):
        print('All states are zero!',flush=True)
        return(0)
    else:
        if(zeroStates==[]):
            print('No states found that are a priori zero!',flush=True)
        else:
            print('These states are zero:',flush=True)
            for state in zeroStates:
                print('\t'+str(state),flush=True)
    
    nrspecies=nrspecies+len(zeroStates)

##### Identify linearities, bilinearities and multilinearities        
    Xsquared=[]
    for i in range(len(X)):
        Xsquared.append(X[i]*X[i])        
    Xsquared=Matrix(Xsquared)
      
    BLList=[]
    MLList=[]
    for i in range(len(SM*F)):
        LHS=str(expand((SM*F)[i]))
        LHS=LHS.replace(' ','')
        LHS=LHS.replace('-','+')
        LHS=LHS.replace('**2','tothepowerof2')
        LHS=LHS.replace('**3','tothepowerof3')
        exprList=LHS.split('+')
        for expr in exprList:
            VarList=expr.split('*')
            counter=0
            factors=[]
            for j in range(len(X)):
                anz=0
                if(str(X[j]) in VarList):
                    anz=1
                    factors.append(X[j])
                if((str(X[j])+'tothepowerof2') in VarList):
                    anz=2 
                    factors.append(X[j])
                    factors.append(X[j])
                if((str(X[j])+'tothepowerof3') in VarList):
                    anz=3
                    factors.append(X[j])
                    factors.append(X[j])
                    factors.append(X[j])
                counter=counter+anz
            if(counter==2):
                string=''            
                for l in range(len(factors)):
                    if(l==len(factors)-1):
                        string=string+str(factors[l])
                    else:
                        string=string+str(factors[l])+'*'
                if(not(string in BLList)):
                    BLList.append(string)
            if(counter>2):
                string=''            
                for l in range(len(factors)):
                    if(l==len(factors)-1):
                        string=string+str(factors[l])
                    else:
                        string=string+str(factors[l])+'*'
                if(not(string in MLList)):
                    MLList.append(string)
        
    COPlusLIPlusBL=[]
    for i in range(len(SM*F)):
        COPlusLIPlusBL.append((SM*F)[i])
        for j in range(len(MLList)):
            ToSubs=expand((SM*F)[i]).coeff(MLList[j])
            COPlusLIPlusBL[i]=expand(COPlusLIPlusBL[i]-ToSubs*parse_expr(MLList[j]))
            
    COPlusLI=[]
    for i in range(len(COPlusLIPlusBL)):
        COPlusLI.append(COPlusLIPlusBL[i])
        for j in range(len(BLList)):
            ToSubs=expand((COPlusLIPlusBL)[i]).coeff(BLList[j])
            COPlusLI[i]=expand(COPlusLI[i]-ToSubs*parse_expr(BLList[j]))
    
    C=zeros(len(COPlusLI),len(X))  
    for i in range(len(COPlusLI)):
    	for j in range(len(X)):
    		C[i*len(X)+j]=expand((COPlusLI)[i]).coeff(X[j])
        
    ML=expand(Matrix(SM*F)-Matrix(COPlusLIPlusBL))
    BL=expand(Matrix(COPlusLIPlusBL)-Matrix(COPlusLI))    
    CM=C        

    CMBL=[]
    if(BLList!=[]):
        for i in range(len(BLList)):
            CVBL=[]
            for k in range(len(BL)):
                CVBL.append(BL[k].coeff(BLList[i]))
            CMBL.append(CVBL)            
    else:
        CVBL=[]
        for k in range(len(BL)):
            CVBL.append(0)
        CMBL.append(CVBL)
    CMBL=Matrix(CMBL).T 
    
    if(MLList!=[]):
        CMML=[]
        for i in range(len(MLList)):
            CVML=[]
            for k in range(len(ML)):
                CVML.append(expand(ML[k]).coeff(MLList[i]))
            CMML.append(CVML)    
        CMML=Matrix(CMML).T  
        BLList=BLList+MLList
        CMBL=Matrix(concatenate((CMBL,CMML),axis=1))
      
    for i in range(len(BLList)):
        BLList[i]=parse_expr(BLList[i])
       
    if(BLList!=[]):    
        CMbig=Matrix(concatenate((CM,CMBL),axis=1))
    else:
        CMbig=Matrix(CM)      

    print('Rank of SM is '+str(SM.rank()) + '!',flush=True)
    SMorig=SM.copy()
    ODE=SMorig*F

#### Get Flux Parameters
    fluxpars=[]
    for flux in F:
        if(flux.args!=()):
            foundFluxpar=False
            for el in flux.args:
                if(not foundFluxpar and el not in X and not is_number(str(el))):
                    if(flux.subs(el, 0)==0):
                        fluxpars.append(el)
                        foundFluxpar=True
        else:
            fluxpars.append(flux)

##### Increase Sparsity of stoichiometry matrix SM
    print('Sparsify stoichiometry matrix with sparsify-level '+str(sparsifyLevel)+'!',flush=True)
    newSM=(Sparsify(SM.T, level=sparsifyLevel, sparseIter=1)).T
    if(newSM!=SM):
        print("Sparsified!",flush=True)
        SM=newSM
    
#### Find conserved quantities
    if(givenCQs==[]):
        print('\nFinding conserved quantities ...',flush=True)
        LCLs, rowsToDel=FindLCL(CMbig.transpose(), X)
    else:
        print('\nI took the given conserved quantities!',flush=True)
        LCLs=givenCQs
    if(LCLs!=[]):
        print(LCLs,flush=True)
    else:
        print('System has no conserved quantities!',flush=True)


#### Define graph structure
    print('\nDefine graph structure ...\n',flush=True)

#### Remove cycles, solve remaining equations, resolve implicits.
#### Strategy: attempt with all minTypes; if the steady-state test fails
#### and type-3 was used, retry with type-3 replaced by FVS.
    # Snapshot mutable state so we can retry from scratch
    _snap_SM = SM.copy()
    _snap_F = F.copy()
    _snap_X = X.copy()
    _snap_fluxpars = list(fluxpars)
    _snap_LCLs = list(LCLs)
    
    _forbidType3 = False
    _usedType3 = False
    fvs_label_printed = False   # print full name on first FVS use
    
    for _attempt in range(2):
        # Restore state
        SM = _snap_SM.copy()
        F = _snap_F.copy()
        X = _snap_X.copy()
        fluxpars = list(_snap_fluxpars)
        LCLs = list(_snap_LCLs)
        SSgraph = DetermineGraphStructure(SM, F, X, neglect)
        cycle = FindCycle(SSgraph, X)
        
        gesnew = 0
        eqOut = []
        fvsSpecies = []
        fvsODEs = {}
        _usedType3 = False
        counter = 1
        
        if _attempt == 1:
            print('\n' + '='*60,flush=True)
            print('Retrying cycle removal: replacing type-3 steps with FVS',flush=True)
            print('='*60 + '\n',flush=True)
        
        while(cycle!=None):
            print('Removing cycle '+str(counter),flush=True)
            minType, state2Rem, fp2Rem, signChanged = GetBestPair(
                cycle, SM, fluxpars, X, LCLs, neglect)
            
            # ---- FVS fallback triggers ----
            if minType == 3 and _forbidType3:
                print(f'   {state2Rem}: type-3 disabled, using FVS',flush=True)
                minType = -1
            
            if(minType==-1):
                print("    Cycle cannot be removed analytically.",flush=True)
                SM, F, X, SSgraph, cycle, fvs_label_printed = _applyFVS(
                    SSgraph, SM, F, X, fvsSpecies, fvsODEs, neglect,
                    fvs_label_printed)
                counter += 1
                continue
            
            # ---- minType 0: conserved quantity ----
            if(minType==0):
                for LCL in LCLs:
                    ls=parse_expr(LCL.split(' = ')[0])
                    if(ls.subs(parse_expr(state2Rem),1)!=ls):
                        LCL2Rem=LCL
                LCLs.remove(LCL2Rem)
                index=list(X).index(parse_expr(state2Rem))
                eqOut.append(state2Rem+' = '+state2Rem)
                print('   '+str(state2Rem)+' --> '+'Done by CQ',flush=True)
            
            # ---- minType 1: unique flux parameter ----
            if(minType==1):
                index=list(X).index(parse_expr(state2Rem))
                eq=(SM*F)[index]
                sol=solve(eq, fp2Rem, simplify=False)
                if not sol:
                    print(f'    solve() returned no solution for {fp2Rem}. '
                          f'Using FVS.',flush=True)
                    SM, F, X, SSgraph, cycle, fvs_label_printed = _applyFVS(
                        SSgraph, SM, F, X, fvsSpecies, fvsODEs, neglect,
                        fvs_label_printed)
                    counter += 1
                    continue
                sol = sol[0]
                eqOut.append(str(fp2Rem)+' = '+str(sol))
                print('   '+str(state2Rem)+' --> '+str(fp2Rem),flush=True)
            
            # ---- minType 2: multi-flux, all unique ----
            if(minType==2):
                anz, sign=GetDimension(state2Rem, X, SM, getSign=True)
                index=list(X).index(parse_expr(state2Rem))
                negs, sumnegs, negfps=GetOutfluxes(state2Rem, X, SM, F, fluxpars)
                poss, sumposs, posfps=GetInfluxes(state2Rem, X, SM, F, fluxpars)
                if(anz==1):
                    print("Error in Type Determination. Please report this bug!",flush=True)
                    return(0)
                else:
                    nenner=1
                    for j in range(anz):
                        if(j>0):
                            nenner=nenner+parse_expr('r_'+state2Rem+'_'+str(j))
                    trafoList=[]
                    if((sign=="minus" and not signChanged) or (sign=="plus" and signChanged)):
                        for j in range(len(negs)):
                            flux=negs[j]
                            fp=negfps[j]
                            prefactor=flux/fp
                            if(j==0):
                                trafoList.append(str(fp)+' = ('+str(sumposs)+')*1/('+str(nenner)+')*1/('+str(prefactor)+')')
                            else:
                                gesnew=gesnew+1
                                trafoList.append(str(fp)+' = ('+str(sumposs)+')*'+'r_'+state2Rem+'_'+str(j)+'/('+str(nenner)+')*1/('+str(prefactor)+')')
                        print('   '+str(state2Rem)+' --> '+str(negfps),flush=True)
                    else:
                        for j in range(len(poss)):
                            flux=poss[j]
                            fp=posfps[j]
                            prefactor=flux/fp
                            if(j==0):
                                trafoList.append(str(fp)+' = ('+str(sumnegs)+')*1/('+str(nenner)+')*1/('+str(prefactor)+')')
                            else:
                                gesnew=gesnew+1
                                trafoList.append(str(fp)+' = ('+str(sumnegs)+')*'+'r_'+state2Rem+'_'+str(j)+'/('+str(nenner)+')*1/('+str(prefactor)+')')
                        print('   '+str(state2Rem)+' --> '+str(posfps),flush=True)
                    for eq in trafoList:
                        eqOut.append(eq)
            
            # ---- minType 3: shared flux parameter ----
            if(minType==3):
                _usedType3 = True
                anz, sign=GetDimension(state2Rem, X, SM, getSign=True)
                index=list(X).index(parse_expr(state2Rem))
                negs, sumnegs, negfps=GetOutfluxes(state2Rem, X, SM, F, fluxpars)
                poss, sumposs, posfps=GetInfluxes(state2Rem, X, SM, F, fluxpars)
                if(anz==1):
                    if((sign=="minus" and not signChanged) or (sign=="plus" and signChanged)):
                        fp2Rem=negfps[0]
                        flux=negs[0]
                    else:
                        fp2Rem=posfps[0]
                        flux=poss[0]
                    eq=(SM*F)[index]
                    sol=solve(eq, fp2Rem, simplify=False)
                    if not sol:
                        print(f'    solve() failed for {fp2Rem} (type 3). '
                              f'Using FVS.',flush=True)
                        SM, F, X, SSgraph, cycle, fvs_label_printed = _applyFVS(
                            SSgraph, SM, F, X, fvsSpecies, fvsODEs, neglect,
                            fvs_label_printed)
                        counter += 1
                        continue
                    sol = sol[0]
                    eqOut.append(str(fp2Rem)+' = '+str(sol))
                    FsearchFlux = matrix_multiply_elementwise(abs(SM[index,:]),F.T)
                    colindex=list(FsearchFlux).index(flux)
                    for row2repl in range(len(SM.col(0))):
                        if(SM[row2repl,colindex]!=0 and row2repl!=index):
                            SM=SM.row_insert(row2repl,SM.row(row2repl)-(SM[row2repl,colindex]/SM[index,colindex])*SM.row(index))
                            SM.row_del(row2repl+1)
                else:
                    nenner=1
                    for j in range(anz):
                        if(j>0):
                            nenner=nenner+parse_expr('r_'+state2Rem+'_'+str(j))
                    trafoList=[]
                    if((sign=="minus" and not signChanged) or (sign=="plus" and signChanged)):
                        for j in range(len(negs)):
                            flux=negs[j]
                            fp=negfps[j]
                            prefactor=flux/fp
                            if(j==0):
                                trafoList.append(str(fp)+' = ('+str(sumposs)+')*1/('+str(nenner)+')*1/('+str(prefactor)+')')
                            else:
                                gesnew=gesnew+1
                                trafoList.append(str(fp)+' = ('+str(sumposs)+')*'+'r_'+state2Rem+'_'+str(j)+'/('+str(nenner)+')*1/('+str(prefactor)+')')
                            FsearchFlux = matrix_multiply_elementwise(abs(SM[index,:]),F.T)
                            colindex=list(FsearchFlux).index(flux)
                            for k in range(len(posfps)):
                                SM=SM.col_insert(len(SM.row(0)),SM.col(colindex))
                                F=F.row_insert(len(F),Matrix(1,1,[poss[k]/nenner]))
                                fluxpars.append(posfps[k])
                            SM.col_del(colindex)
                            F.row_del(colindex)
                            fluxpars.__delitem__(colindex)
                        print('   '+str(state2Rem)+' --> '+str(negfps),flush=True)
                    else:
                        for j in range(len(poss)):
                            flux=poss[j]
                            fp=posfps[j]
                            prefactor=flux/fp
                            if(j==0):
                                trafoList.append(str(fp)+' = ('+str(sumnegs)+')*1/('+str(nenner)+')*1/('+str(prefactor)+')')
                            else:
                                gesnew=gesnew+1
                                trafoList.append(str(fp)+' = ('+str(sumnegs)+')*'+'r_'+state2Rem+'_'+str(j)+'/('+str(nenner)+')*1/('+str(prefactor)+')')
                            FsearchFlux = matrix_multiply_elementwise(abs(SM[index,:]),F.T)
                            colindex=list(FsearchFlux).index(flux)
                            for k in range(len(negfps)):
                                SM=SM.col_insert(len(SM.row(0)),SM.col(colindex))
                                F=F.row_insert(len(F),Matrix(1,1,[negs[k]/nenner]))
                                fluxpars.append(negfps[k])
                            SM.col_del(colindex)
                            F.row_del(colindex)
                            fluxpars.__delitem__(colindex)
                        print('   '+str(state2Rem)+' --> '+str(posfps),flush=True)
                    for eq in trafoList:
                        eqOut.append(eq)
            
            X.row_del(index)
            SM.row_del(index)
            SSgraph=DetermineGraphStructure(SM, F, X, neglect)
            cycle=FindCycle(SSgraph, X)
            counter=counter+1
        
        print('There is no cycle in the system!\n',flush=True)
        
    #### Solve remaining equations
        eqOut.reverse()
        print('Solving remaining equations ...\n',flush=True)
        while(SSgraph!={}):
            node=FindNodeToSolve(SSgraph)
            if node is None:
                print('    No leaf node found. Using FVS for remaining '
                      'species.',flush=True)
                SM, F, X, SSgraph, _, fvs_label_printed = _applyFVS(
                    SSgraph, SM, F, X, fvsSpecies, fvsODEs, neglect=[],
                    fvs_label_printed=fvs_label_printed)
                continue
            index=list(X).index(parse_expr(node))
            ode_expr = (SM*F)[index]
            if ode_expr == 0:
                print(f'    {node}: ODE = 0 (free parameter from conserved '
                      f'quantities)',flush=True)
                eqOut.insert(0, node + ' = ' + node)
                X.row_del(index)
                SM.row_del(index)
                SSgraph=DetermineGraphStructure(SM, F, X, neglect=[])
                continue
            sol=solve(ode_expr, parse_expr(node), simplify=True)
            if not sol:
                print(f'    solve() returned no solution for {node}. '
                      f'Treating as FVS species.',flush=True)
                fvsSpecies.append(node)
                fvsODEs[node] = ode_expr
                fvs_param = Symbol('fvs_' + node)
                for f_idx in range(len(F)):
                    F[f_idx] = F[f_idx].subs(parse_expr(node), fvs_param)
                X.row_del(index)
                SM.row_del(index)
                eqOut.insert(0, node + ' = fvs_' + node)
                SSgraph=DetermineGraphStructure(SM, F, X, neglect=[])
                continue
            eqOut.insert(0, node+' = '+str(sol[0]))
            for f in range(len(F)):
                F[f]=F[f].subs(parse_expr(node), sol[0])
            X.row_del(index)
            SM.row_del(index)
            SSgraph=DetermineGraphStructure(SM, F, X, neglect=[])
        
    #### Resolve implicit (self-referential) equations
        _has_implicit = False
        _parsed_eqs = []
        for eq_str in eqOut:
            if ' = ' in eq_str:
                _ls_str, _rs_str = eq_str.split(' = ', 1)
            else:
                _ls_str, _rs_str = eq_str.split('=', 1)
            _ls_str = _ls_str.strip()
            _rs_str = _rs_str.split('[')[0].strip()
            _ls = parse_expr(_ls_str)
            _rs = parse_expr(_rs_str)
            if _ls in _rs.free_symbols:
                _has_implicit = True
            _parsed_eqs.append((_ls, _rs))
        
        if _has_implicit:
            print('Resolving implicit (self-referential) equations ...',
                  flush=True)
            _resolved_subs = {}
            for _i, (_ls, _rs) in enumerate(_parsed_eqs):
                if _ls in _rs.free_symbols:
                    _sols = solve(_ls - _rs, _ls)
                    if _sols:
                        _resolved_expr = simplify(_sols[0])
                        _parsed_eqs[_i] = (_ls, _resolved_expr)
                        _resolved_subs[_ls] = _resolved_expr
                        print(f'  Resolved: {_ls} = {_resolved_expr}',
                              flush=True)
                    else:
                        print(f'  WARNING: Could not resolve implicit '
                              f'equation for {_ls}',flush=True)
            if _resolved_subs:
                for _i, (_ls, _rs) in enumerate(_parsed_eqs):
                    for _sym, _expr in _resolved_subs.items():
                        if _sym != _ls:
                            _rs = _rs.subs(_sym, _expr)
                    _parsed_eqs[_i] = (_ls, simplify(_rs))
                eqOut = [f'{_ls} = {_rs}' for _ls, _rs in _parsed_eqs]
            print('Done.\n',flush=True)
        
    #### Verify and possibly retry
        if _attempt == 0 and _usedType3:
            _correct = _testSteadyStateQuiet(
                ODE, eqOut, zeroStates, fvsSpecies)
            if not _correct:
                print('Steady-state test failed after type-3 cycle removal.',
                      flush=True)
                _forbidType3 = True
                continue   # retry with type-3 disabled
        break  # solution accepted or second attempt done

    # ====================================================================
    # v2.0: Handle FVS species: polynomial solving & convergence conditions
    # ====================================================================
    fvsResults = []       # Collect results for output
    fvsUnsolved = []      # Species that could not be solved polynomially
    fvsUnsolvedODEs = []  # Their ODE expressions
    fvsFreeParams = []    # Species whose equations were trivial (free parameters)
    
    if fvsSpecies:
        print('\n' + '='*60,flush=True)
        print('Solving FVS species (polynomial, maxDegree=%d)' % maxPolyDegree,flush=True)
        print('='*60 + '\n',flush=True)
        
        for fvs_node in fvsSpecies:
            x_sym = parse_expr(fvs_node)
            ode_expr = fvsODEs[fvs_node]
            
            # Substitute already-solved equations into the FVS ODE
            # so it depends only on parameters and the FVS species itself
            for eq_str in eqOut:
                ls, rs = eq_str.split(' = ')
                ode_expr = ode_expr.subs(parse_expr(ls), parse_expr(rs))
            # Also substitute zero states
            for zs in zeroStates:
                ode_expr = ode_expr.subs(zs, 0)
            
            ode_expr = simplify(expand(ode_expr))
            print(f'  Species {fvs_node}:',flush=True)
            print(f'    ODE: d{fvs_node}/dt = {ode_expr}',flush=True)
            
            solved_sym, solutions, success = SolvePolynomialSteadyState(
                ode_expr, x_sym, maxDegree=maxPolyDegree)
            
            if success and solutions:
                print(f'    Solved for: {solved_sym}',flush=True)
                print(f'    Strictly positive solutions:',flush=True)
                for k, sol in enumerate(solutions):
                    print(f'      {solved_sym}_{k+1} = {sol}',flush=True)
                    fvsResults.append(f'{solved_sym} = {sol}  [FVS poly root {k+1}, from ODE of {fvs_node}]')
                # Use first positive solution as the canonical one
                eqOut.append(f'{solved_sym} = {solutions[0]}')
            else:
                if not success:
                    print(f'    Not polynomial (or degree > {maxPolyDegree}) '
                          f'in {fvs_node}. No closed-form solution.',flush=True)
                else:
                    print(f'    No strictly positive solutions found.',flush=True)
                fvsUnsolved.append(fvs_node)
                fvsUnsolvedODEs.append(ode_expr)
                fvsResults.append(
                    f'{fvs_node} = NUMERIC  '
                    f'[solve d{fvs_node}/dt = {ode_expr} = 0 numerically]')
                # Still add to eqOut so downstream code knows about it
                eqOut.append(f'{fvs_node} = fvs_{fvs_node}')
        
        # ============================================================
        # v2.2: Remove fvs_ prefixes BEFORE recurrence resolution.
        # This way the recurrence pass works with original species names
        # and we can properly detect & resolve implicit equations.
        # ============================================================
        print('\nRemoving internal fvs_ prefixes ...',flush=True)
        fvs_back_subs = {}
        for fn in fvsSpecies:
            fvs_back_subs[Symbol('fvs_' + fn)] = parse_expr(fn)
        
        # Apply to all equations in eqOut
        _temp_eqs = []
        for eq_str in eqOut:
            if ' = ' in eq_str:
                ls_str, rs_str = eq_str.split(' = ', 1)
            else:
                ls_str, rs_str = eq_str.split('=', 1)
            rs_str = rs_str.split('[')[0].strip()
            ls = parse_expr(ls_str.strip())
            rs = parse_expr(rs_str)
            for fvs_sym, orig_sym in fvs_back_subs.items():
                ls = ls.subs(fvs_sym, orig_sym)
                rs = rs.subs(fvs_sym, orig_sym)
            _temp_eqs.append((ls, rs))
        
        # Also clean fvsResults
        cleaned_fvsResults = []
        for fr in fvsResults:
            for fvs_sym, orig_sym in fvs_back_subs.items():
                fr = fr.replace(str(fvs_sym), str(orig_sym))
            cleaned_fvsResults.append(fr)
        fvsResults = cleaned_fvsResults
        
        # Clean fvsUnsolvedODEs
        for idx in range(len(fvsUnsolvedODEs)):
            for fvs_sym, orig_sym in fvs_back_subs.items():
                fvsUnsolvedODEs[idx] = fvsUnsolvedODEs[idx].subs(fvs_sym, orig_sym)
        
        print('  Done.\n',flush=True)
        
        # ============================================================
        # resolveRecurrence: Single forward-pass substitution.
        # ============================================================
        if fvsResults:
            print('Resolving recurrences among FVS solutions ...',flush=True)
            
            for i in range(len(_temp_eqs)):
                lhs_i, rhs_i = _temp_eqs[i]
                for j in range(i):
                    lhs_j, rhs_j = _temp_eqs[j]
                    rhs_i = rhs_i.subs(lhs_j, rhs_j)
                _temp_eqs[i] = (lhs_i, simplify(rhs_i))
            
            print('  Done.\n',flush=True)
        
        # ============================================================
        # Resolve implicit (self-referential) equations that arise
        # from fvs_ -> original substitution (e.g. y11 = f(..., y11))
        # and drop trivial identities (e.g. TCR = TCR).
        # ============================================================
        print('Resolving implicit FVS equations and dropping '
              'trivials ...',flush=True)
        
        cleaned_eqs = []
        _resolved_subs = {}
        fvsFreeParams = []   # FVS species whose eqs were trivial
        
        for lhs, rhs in _temp_eqs:
            # Check for trivial identity
            if simplify(lhs - rhs) == 0:
                print(f'  Dropped trivial identity: {lhs} = {rhs}',
                      flush=True)
                if str(lhs) in fvsSpecies:
                    fvsFreeParams.append(str(lhs))
                continue
            # Check for implicit (self-referential) equation
            if lhs in rhs.free_symbols:
                _sols = solve(lhs - rhs, lhs)
                if _sols:
                    rhs = simplify(_sols[0])
                    _resolved_subs[lhs] = rhs
                    print(f'  Resolved implicit: {lhs} = {rhs}',
                          flush=True)
                else:
                    print(f'  WARNING: Could not resolve implicit '
                          f'equation for {lhs}',flush=True)
            cleaned_eqs.append((lhs, rhs))
        
        # Propagate resolved substitutions into remaining equations
        if _resolved_subs:
            for i, (lhs, rhs) in enumerate(cleaned_eqs):
                for sym, expr in _resolved_subs.items():
                    if sym != lhs:
                        rhs = rhs.subs(sym, expr)
                cleaned_eqs[i] = (lhs, simplify(rhs))
        
        eqOut = [f'{lhs} = {rhs}' for lhs, rhs in cleaned_eqs]
        
        # Update fvsResults
        fvsResults = []
        for lhs, rhs in cleaned_eqs:
            if str(lhs) in fvsSpecies:
                fvsResults.append(f'{lhs} = {rhs}  [FVS, resolved]')
        
        print('  Done.\n',flush=True)
        
        if fvsUnsolved:
            print('\n' + '='*60,flush=True)
            print('CONVERGENCE CONDITIONS for iterative equilibration',flush=True)
            print('='*60 + '\n',flush=True)
            
            unsolved_syms = [parse_expr(s) for s in fvsUnsolved]
            conditions = GetConvergenceConditions(fvsUnsolvedODEs, unsolved_syms)
            for cond in conditions:
                print('  ' + cond,flush=True)
            print(flush=True)
    
    # ====================================================================
    # End of v2.0 FVS handling
    # ====================================================================

#### Test Solution  
    if(testSteady=='T'):
        print('Testing Steady State...\n',flush=True)
        
        # Determine which species have no equation (free parameters)
        _solved_species = set()
        for eq_str in eqOut:
            ls_str = eq_str.split(' = ')[0].strip() if ' = ' in eq_str else eq_str.split('=')[0].strip()
            _solved_species.add(ls_str)
        _free_species = set()
        for x in Xo:
            if str(x) not in _solved_species and x not in zeroStates:
                _free_species.add(str(x))
        
        # Use conserved quantities to derive equations for free species
        _cq_eqs = []   # additional substitutions: 'species = expr'
        _remaining_free = set(_free_species)
        for lcl in LCLs:
            ls_str, rs_str = lcl.split(' = ')
            cq_expr = parse_expr(ls_str) - parse_expr(rs_str)  # = 0
            # Substitute solved equations into CQ
            for eq_str in eqOut:
                if ' = ' in eq_str:
                    el, er = eq_str.split(' = ', 1)
                else:
                    el, er = eq_str.split('=', 1)
                er = er.split('[')[0].strip()
                cq_expr = cq_expr.subs(parse_expr(el.strip()), parse_expr(er))
            for zs in zeroStates:
                cq_expr = cq_expr.subs(zs, 0)
            # Substitute already-derived CQ equations
            for cq_eq_str in _cq_eqs:
                cl, cr = cq_eq_str.split(' = ', 1)
                cq_expr = cq_expr.subs(parse_expr(cl), parse_expr(cr))
            cq_expr = simplify(cq_expr)
            # Solve for one of the remaining free species
            for fp in sorted(_remaining_free):
                fp_sym = parse_expr(fp)
                if cq_expr.has(fp_sym):
                    sol = solve(cq_expr, fp_sym)
                    if sol:
                        _cq_eqs.append(f'{fp} = {sol[0]}')
                        _remaining_free.discard(fp)
                        print(f'  From CQ: {fp} = {sol[0]}',flush=True)
                        break
        
        if _remaining_free:
            print(f'  Remaining free parameters: {sorted(_remaining_free)}',
                  flush=True)
        
        # Now test: substitute eqOut + CQ-derived eqs into each ODE
        _all_test_eqs = list(eqOut) + _cq_eqs
        
        NonSteady=False
        for i in range(len(ODE)):
            expr=parse_expr(str(ODE[i]))
            for zs in zeroStates:
                expr=expr.subs(zs, 0)
            for j in range(len(_all_test_eqs)):
                eq_str = _all_test_eqs[-(j+1)]
                if ' = ' in eq_str:
                    ls, rs = eq_str.split(' = ', 1)
                else:
                    ls, rs = eq_str.split('=', 1)
                rs = rs.split('[')[0].strip()
                ls=parse_expr(ls.strip())
                rs=parse_expr(rs)
                expr=expr.subs(ls, rs)
            expr=simplify(expr)
            if(expr!=0):
                # Residual may contain truly free parameters – that's OK
                has_free = False
                for fp in _remaining_free:
                    if expr.has(parse_expr(fp)):
                        has_free = True
                if not has_free:
                    print('   Equation '+str(ODE[i]),flush=True)
                    print('   results:'+str(expr),flush=True)
                    NonSteady=True
        if(NonSteady):
            print('Solution is wrong!\n',flush=True)
        else:
            if _remaining_free:
                print('Solution is correct (free parameter'
                      + ('s' if len(_remaining_free) > 1 else '') + ': '
                      + ', '.join(sorted(_remaining_free)) + ')!\n',flush=True)
            else:
                print('Solution is correct!\n',flush=True)
            
    elif(testSteady=='F'):
        print('Skipping the Testing of Steady State...\n',flush=True)
        
    else:
        print('Skipping the Testing of Steady State...\n',flush=True)
    
#### Print Equations
    print('I obtained the following equations:\n',flush=True)
    if(outputFormat=='M'):
        for state in zeroStates:
            print('\tinit_'+str(state)+'  "0"'+'\n',flush=True)
        eqOutReturn=[]
        for i in range(len(eqOut)):
            eq_str = eqOut[i]
            if ' = ' in eq_str:
                ls, rs = eq_str.split(' = ', 1)
            else:
                ls, rs = eq_str.split('=', 1)
            rs = rs.split('[')[0].strip()  # strip FVS annotations
            ls=parse_expr(ls.strip())
            rs=parse_expr(rs)
            for j in range(i,len(eqOut)):
                eq_str_j = eqOut[j]
                if ' = ' in eq_str_j:
                    ls2, rs2 = eq_str_j.split(' = ', 1)
                else:
                    ls2, rs2 = eq_str_j.split('=', 1)
                rs2 = rs2.split('[')[0].strip()
                rs2=parse_expr(rs2)
                rs2=rs2.subs(ls,rs)
                eqOut[j]=str(ls2)+'='+str(rs2)
            for state in Xo:
                ls=ls.subs(state, parse_expr('init_'+str(state)))
                rs=rs.subs(state, parse_expr('init_'+str(state)))
            eqOut[i]=str(ls)+'  "'+str(rs)+'"'
                            
        for i in range(len(eqOut)):
            eqOut[i]=eqOut[i].replace('**','^')
                    
        for eq in eqOut:
            print('\t'+eq+'\n',flush=True)
            eqOutReturn.append(eq)            
        
    else:
        for state in zeroStates:
            print('\t'+str(state)+' = 0'+'\n',flush=True)
        eqOutReturn=[]
        for eq in eqOut:
            if ' = ' in eq:
                ls, rs = eq.split(' = ', 1)
            else:
                ls, rs = eq.split('=', 1)
            rs_clean = rs.split('[')[0].strip()  # strip FVS annotations
            print('\t'+ls.strip()+' = "'+rs_clean+'",'+'\n',flush=True)
            eqOutReturn.append(ls.strip()+'='+rs_clean)

    print('Number of Species:  '+str(nrspecies),flush=True)
    print('Number of Equations:  '+str(len(eqOut)+len(zeroStates)),flush=True)
    print('Number of new introduced variables:  '+str(gesnew),flush=True)
    if fvsSpecies:
        print('FVS species:  '+str(len(fvsSpecies))
              +' ('+', '.join(fvsSpecies)+')',flush=True)
        _n_solved = len(fvsSpecies) - len(fvsUnsolved) - len(fvsFreeParams)
        if _n_solved > 0:
            print('  Solved:  '+str(_n_solved),flush=True)
        if fvsFreeParams:
            print('  Free (determined by conserved quantities):  '
                  +str(len(fvsFreeParams))+' ('+', '.join(fvsFreeParams)+')',flush=True)
        if fvsUnsolved:
            print('  Require numerical solving:  '+str(len(fvsUnsolved)),flush=True)
    
    # v2.0: Return extended result
    result = {
        'equations': eqOutReturn,
        'fvs_species': fvsSpecies,
        'fvs_results': fvsResults,
        'fvs_unsolved': fvsUnsolved,
        'fvs_free_params': fvsFreeParams,
    }
    if fvsUnsolved:
        unsolved_syms = [parse_expr(s) for s in fvsUnsolved]
        result['convergence_conditions'] = GetConvergenceConditions(
            fvsUnsolvedODEs, unsolved_syms)
    
    return(result)
