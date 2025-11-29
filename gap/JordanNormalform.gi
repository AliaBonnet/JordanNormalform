#Depends on NoFoMa Package by Meinolf Geck 
#https://github.com/geckmf/NoFoMa

#Auxilliary Functions

InstallGlobalFunction(ConvertVecToRowMat, function(vec)
    local mat, n, i;
    n := Length(vec);
    mat := ZeroMatrix(BaseDomain(vec), 1, n);
    for i in [1..n] do 
        mat[1,i] := vec[i];
    od;
    return mat;
end);

#Avoids creating the vector space
InstallGlobalFunction(GenerateRandomVector, function(F, d) #Field d, length d
    local vec, i;
    vec := ZeroVector(F,d);
    for i in [1..d] do
        vec[i] := PseudoRandom(F);
    od;
    return vec;
end);

#Spinning algorithm, checking rank at the end
#Returns basis of <vec>_A
InstallGlobalFunction(SpinningCheck, function(vec, A)
    local F, n, i, res, r;
    F := BaseDomain(A);
    n := NrRows(A);
    res := ZeroMatrix(F,n,n);
    CopySubMatrix(ConvertVecToRowMat(vec), res, [1..1], [1..1], [1..n], [1..n]);
    for i in [2..n] do
        vec := vec * A;
        CopySubMatrix(ConvertVecToRowMat(vec), res, [1..1], [i..i], [1..n], [1..n]);
    od;
    r := RankMat(res);
    return ExtractSubMatrix(res,[1..r],[1..n]);
end);

#Spinning algorithm
#Returns (vec, vec*A, ..., vec*A^(goal-1))
InstallGlobalFunction(SpinUntil, function(vec, A, goal)
    local n, i, res, F;
    F := BaseDomain(A);
    if goal = 0 then 
        return ConvertVecToRowMat(vec); 
    fi;
    n := NrRows(A);
    res := ZeroMatrix(F,goal, n);
    CopySubMatrix(ConvertVecToRowMat(vec), res, [1..1], [1..1], [1..n], [1..n]);
    for i in [2..goal] do
        vec := vec * A;
        CopySubMatrix(ConvertVecToRowMat(vec), res, [1..1], [i..i], [1..n], [1..n]);
    od;
    return res;
end);

#Returns a vector v such that v is not in subspace spanned by gen
#Assumes that gen is already echelonised 
InstallGlobalFunction(FindVectorNotInSubspaceNC, function(gen) #assumes gen is already echelonised
    local w, i, n, F, r, zsf;
    r := NrRows(gen); #dimension of subspace
    F := BaseDomain(gen);
    n := Length(gen[1]);
    if (r = n) then 
        Print("Subspace is already equal to entire space.");
        return Zero(F^n);
    fi;    
    zsf := ZeroMatrix(F,r+1,r+1);
    CopySubMatrix(gen, zsf, [1..r], [1..r], [1..r], [1..r]);
    w := ZeroVector(F,n);
    for i in [1..n] do 
        if zsf[i,i] = Zero(F) then
            w[i] := One(F);
            return w;
        fi;
    od;
    Print("Could not find vector that isn't in subspace!"); 
    return fail;
end);

#TODO: Construct cyclic vector 
#Finds a cyclic vector for A
#Doesn't check if $A$ is cyclic
InstallGlobalFunction(FindCyclicVectorNC, function(A) #Field F, Matrix A, n upper bound of loops;
    local checked, vec, gens, n, i, F;
    n := NrRows(A);
    F := BaseDomain(A);
    checked := []; #Avoid checking double
    vec := GenerateRandomVector(F,n); 
    for i in [1..40] do
        if not (vec in checked) then
            gens := SpinUntil(vec, A, n); #potential basis
            if not (RankMat(gens) < n) then #check if linearly independent, 
                return gens; #returns spinned basis, first element is the cyclic vector
            else
                Add(checked, vec);
            fi;
        fi;
        vec := GenerateRandomVector(F,n);
    od;
    Print("Failed to find cyclic vector!");
    return fail;
end);

InstallGlobalFunction(RemoveZeroRows, function(mat)
    local i,matcopy;
    matcopy := MutableCopyMat(mat);
    for i in Reversed([1..NrRows(matcopy)]) do 
        if IsZero(matcopy[i]) then 
            Remove(matcopy, i);
        fi;
    od;
    return matcopy;
end);

#Returns minimal polynomial of v with regards to A as a 
#coefficient list (ascending powers), along with a basis of <v>_A
InstallGlobalFunction(MinPolVec, function(A,v)
    local i,spun,rel;
    spun := [v];
    for i in [0..NrRows(A)] do
        rel := NullspaceMat(spun);
        if Size(rel) > 0 then 
            return [rel[1],spun];
        fi;
        v := v*A;
        Add(spun,v);
    od;
    Print("Could not compute minimal polynomial of vector!");
    return fail;
end);

#Evaluates v*p(A) VERY VERY efficiently
#Takes (v,vA,...v^n-1A) and polynomial as input
InstallGlobalFunction(PolyEvalFromSpan, function(span,pol)
    local coeffs, i, resu;
    resu := ZeroVector(BaseDomain(span),NrCols(span));
    coeffs := CoefficientsOfUnivariatePolynomial(pol);
    for i in [1..Size(coeffs)] do 
        resu := resu + coeffs[i]*span[i];
    od;
    return resu;
end);

#Primary Decomposition using a modified version of Allan Steel's algorithm 
#Standalone version 
#Returns matrix B such that B*A*B^-1 is in primary decomposition form 
InstallGlobalFunction(PrimaryDecomp, function(A) #returns mat such that mat * A * mat^-1 is primary decomp form 
    local r, rank, F, n, m, f, w, p, j,i, wspan, gens, facs, L_i, qi, k, U_j, v, COB, pot, gs, f2, toAdd;
    rank := 0;
    n := NrRows(A);
    F := DefaultFieldOfMatrix(A);
    v := GenerateRandomVector(F,n);
    gens := []; #Li_s as in Steel paper will go in here (but i will not use separate U)
    gs := []; #distinct factors of minimal polynomial 
    while not rank = n do 
        m := UnivariatePolynomial(F,MinPolVec(A,v)[1]);
        p := One(PolynomialRing(F));
        for i in [1..Size(gens)] do 
            L_i := gens[i];
            if not IsOne(Gcd(m,gs[i])) then
                f := m;
                pot := 0; #this used to be 1... is this ok?
                for j in [1..n] do #n mal ist overkill, ich glaube das ganze geht bestimmt eleganter
                    f2 := Quotient(f,gs[i]);
                    if f2 = fail then 
                        break;
                    fi;
                    pot := pot + 1;
                    f := f2;
                od; 
                w := PolynomialToMatVec(A,CoefficientsOfUnivariatePolynomial(f),v);
                p := p * gs[i]^pot;
                wspan := SpinningCheck(w,A); #TODO: weiß ich noch nicht wie weit ich spinnen werde? (DOCH ICH KENNE GENAU DAS MINIMALPOLYNOM)
                toAdd := EcheloniseMat(Concatenation(wspan,gens[i]));
                if not IsMatrix(toAdd) then 
                    toAdd := [toAdd];  # Convert vector to 1-row matrix
                    toAdd := Matrix(toAdd);
                fi;
                gens[i] := toAdd;
            fi;
        od;
        v :=PolynomialToMatVec(A,CoefficientsOfUnivariatePolynomial(p),v);
        m := Quotient(m,p);
        if not IsOne(m) then 
            facs := Factors(m);
            facs := Collected(facs);
            for i in [1..Size(facs)] do 
                qi := (facs[i][1])^(facs[i][2]);
                w := PolynomialToMatVec(A,CoefficientsOfUnivariatePolynomial(Quotient(m,qi)),v);
                wspan := SpinningCheck(w,A);
                Add(gens,wspan);
                Add(gs, facs[i][1]);
            od;
        fi;
        #alles in eine matrix
        COB := ZeroMatrix(F,n,n);
        k := 0;
        for i in [1..Size(gens)] do 
            L_i := gens[i];
            CopySubMatrix(L_i, COB, [1..NrRows(L_i)], [k+1..k+NrRows(L_i)],[1..n], [1..n]);
            k := k+NrRows(L_i);
        od;
        COB := RemoveZeroRows(COB);
        rank := NrRows(COB);
        if not rank = n then 
            COB := EcheloniseMat(COB);
            v := FindVectorNotInSubspaceNC(COB);
        fi;
    od;
    return COB;
end);

#Primary Decomposition using a modified version of Allan Steel's algorithm for use in the Jordan normal form function
#Takes matrix along with its minimal polynomial as input
#Returns matrix B such that B*A*B^-1 is in primary decomp. form, dimensions of primary subspaces
#and factors of minimal polynomial in correct order
#TODO: STOP COLLECTING MULTIPLICITIES
InstallGlobalFunction(PrimaryDecompositionforJNF, function(A, minpol)
    local r, merge, split1, split2, facOccCorrect, rank, F, n, m, f, w, p, j,i, wspan, gens, facs, L_i, qi, k, U_j, v, COB, pot, gs, f2, toAdd, dims, minpolCollected, pos, minpolFacs, minpolMult;
    rank := 0;
    n := NrRows(A);
    F := DefaultFieldOfMatrix(A);
    v := GenerateRandomVector(F,n);
    gens := []; #Li_s as in Steel paper will go in here (but i will not use separate U)
    gs := []; #distinct factors of minimal polynomial 
    dims := [];
    minpolFacs := [];
    minpolMult := [];
    while not rank = n do 
        m := UnivariatePolynomial(F,MinPolVec(A,v)[1]);
        p := One(PolynomialRing(F));
        for i in [1..Size(gens)] do 
            L_i := gens[i];
            if not IsOne(Gcd(m,gs[i])) then
                f := m;
                pot := 0; #this used to be 1... is this ok?
                for j in [1..n] do #n mal ist overkill
                    f2 := Quotient(f,gs[i]);
                    if f2 = fail then 
                        break;
                    fi;
                    pot := pot + 1;
                    f := f2;
                od; 
                w := PolynomialToMatVec(A,CoefficientsOfUnivariatePolynomial(f),v);
                p := p * gs[i]^pot;
                wspan := SpinningCheck(w,A); #TODO: SPINUNTIL, ICH KENNE DAS MINPOL
                toAdd := EcheloniseMat(Concatenation(wspan,gens[i]));
                if not IsMatrix(toAdd) then 
                    toAdd := [toAdd];  # Convert vector to 1-row matrix
                    toAdd := Matrix(toAdd);
                fi;
                gens[i] := toAdd;
            fi;
        od;
        v :=PolynomialToMatVec(A,CoefficientsOfUnivariatePolynomial(p),v);
        m := Quotient(m,p);
        if not IsOne(m) then 
            facs := Factors(m);
            facs := Collected(facs);
            for i in [1..Size(facs)] do 
                pos := Position(minpolFacs, facs[i][1]);
                if pos = fail then #TODO: einmal mit index abspeichern
                    Add(minpolFacs, facs[i][1]);
                    Add(minpolMult, facs[i][2]);
                else
                    minpolMult[pos] := minpolMult[pos] + facs[i][2];
                fi;
                qi := (facs[i][1])^(facs[i][2]);
                w := PolynomialToMatVec(A,CoefficientsOfUnivariatePolynomial(Quotient(m,qi)),v);
                wspan := SpinningCheck(w,A);
                Add(gens,wspan);
                Add(gs, facs[i][1]);
            od;
        fi;
        #alles in eine matrix
        COB := ZeroMatrix(F,n,n);
        k := 0;
        for i in [1..Size(gens)] do 
            L_i := gens[i];
            CopySubMatrix(L_i, COB, [1..NrRows(L_i)], [k+1..k+NrRows(L_i)],[1..n], [1..n]);
            k := k+NrRows(L_i);
        od;
        COB := RemoveZeroRows(COB);
        rank := NrRows(COB);
        if not rank = n then 
            COB := EcheloniseMat(COB);
            v := FindVectorNotInSubspaceNC(COB);
        fi;
    od;
    #collect dimensions
    for i in [1..Size(gens)] do
        Add(dims, Size(gens[i]));
    od;
    #nicht gut
    facOccCorrect := Collected(Factors(minpol));
    split1 := [];
    split2 := [];
    minpolMult := [];
    for j in [1..Size(facOccCorrect)] do 
        Add(split1, facOccCorrect[j][1]);
        Add(split2, facOccCorrect[j][2]);
    od;
    for i in [1..Size(minpolFacs)] do 
        pos := Position(split1, minpolFacs[i]);
        Add(minpolMult, split2[pos]);
    od;
    merge := [];
    for i in [1..Size(minpolFacs)] do 
        Add(merge, [minpolFacs[i],minpolMult[i]]);
    od;
    return [COB, dims, merge];
end);

#Primary Decomposition for cyclic matrices 
#Jordan normal form will call this function if a cyclic matrix is detected
InstallGlobalFunction(PrimaryDecompositionforJNFCyclic, function(A, minpol)
    local r, vspan,rank, F, n, m, f, w, j,i, wspan, gens, facs, L_i, qi, k, v, COB, dims;
    rank := 0;
    n := NrRows(A);
    F := DefaultFieldOfMatrix(A);
    A := Matrix(F,A);
    vspan := FindCyclicVectorNC(A);
    v := vspan[1];
    dims := [];
    m := minpol;
    facs := Factors(m);
    facs := Collected(facs);
    COB := ZeroMatrix(F,n,n);
    k := 0;
    for i in [1..Size(facs)] do 
        qi := (facs[i][1])^(facs[i][2]);
        w := PolyEvalFromSpan(vspan,Quotient(m,qi));
        wspan := SpinUntil(w,A,Degree(facs[i][1])*facs[i][2]);
        CopySubMatrix(wspan, COB, [1..NrRows(wspan)], [k+1..k+NrRows(wspan)], [1..n], [1..n]);
        Add(dims, NrRows(wspan));
        k := k + NrRows(wspan);
    od;
    return [COB, dims, facs];
end);

#Input: matrix A with minimal polynomial p^m, vector v and p(A)
#Returns: v, length r of v, vp^(r-1)(A)
InstallGlobalFunction(GetMinPolPowerWithVec, function(A,p,m,v,Ainp) 
    #this doesnt even need A and p
    local j, veccopy, lastcopy;
    if IsZero(v) then
        return [v,0,v];
    fi;
    lastcopy := ShallowCopy(v);
    veccopy := ShallowCopy(v);
    for j in [1..m] do 
        veccopy := veccopy * Ainp; 
        if IsZero(veccopy) then #check if vector got nulled
            return [v, j, lastcopy]; #why not just return v? GLEICH MAL GUCKEN 
        fi;
        lastcopy := veccopy;
    od;
    Print("Failed to find minimal polynomial of ", v);
    return fail; #shouldn't happen
end);

#Returns linear dependence q_1,...,q_k as described in paper for cyclic decomposition
InstallGlobalFunction(FindLinearDependenceNC, function(vecs, A, d) #vecs, A, degree of p, returns coeffs of qis (ascending degree) #THIS ONLY WORKS FOR SETTING IN THEOREM 
    local n, F, i, rel, tosolve, currdim, qis;
    F := DefaultFieldOfMatrix(A);
    n := NrRows(A);
    tosolve := MutableCopyMatrix(ZeroMatrix(F,Length(vecs)*d,n));
    currdim := 1; 
    for i in [1..Length(vecs)] do
        CopySubMatrix(SpinUntil(vecs[i], A, d), tosolve, [1..d], [currdim..currdim +  d-1], [1..n], [1..n]);
        currdim := currdim +  d;
    od;
    rel := NullspaceMat(tosolve);
    if rel = fail or rel = [] then #shouldn't happen
        Print("Couldn't find F[A]-linear dependence"); 
    fi;
    rel := rel[1];
    #turn into usable polynomial coeffs
    qis := [];
    currdim := 1; #ja das muss so #aber kann man auch einfach arithmetisch machen
    for i in [1..Length(vecs)] do 
        Add(qis, ExtractSubVector(rel, [currdim..currdim+d-1]));
        currdim := currdim + d;
    od;
    return qis;
end);

#Input: matrix A with minimalpolynomial p^m
#Returns matrix B such that A^Inverse(B) is in cyclic decomposition form, dimensions of cyclic subspaces
InstallGlobalFunction(CyclicDecompositionOfPrimarySubspace, function (A, p, m) 
    local F, n, d, Ainp, ws, allspun, wspun, k, tomult, minpolpowers, dims, wtrip, i, conj, sumdim, wdim, qis, j, r, wstrich, currdim, w, vecs;
    F := DefaultFieldOfMatrix(A);
    n := NrRows(A);
    d := Degree(p);
    if m * d = n then #return if it's already cyclic
        return [One(GL(n,F)), [n]]; 
    fi;
    Ainp := p(A);
    ws := [];
    w := ZeroVector(F,n);
    while IsZero(w) do #make sure we aren't spinning zero vector
        w := GenerateRandomVector(F,n);
    od;
    Add(ws,w);
    allspun := EcheloniseMat(SpinUntil(w,A,n));
    while NrRows(allspun) < n do 
        w := FindVectorNotInSubspaceNC(allspun);
        wspun := SpinUntil(w,A,n);
        Append(allspun, wspun);
        Add(ws, w);
        allspun := EcheloniseMat(allspun); #DOING THIS SHOULD
    od;
    minpolpowers := ZeroVector(F,NrRows(ws));
    for i in [1..NrRows(ws)] do 
        minpolpowers[i] := GetMinPolPowerWithVec(A,p,m,ws[i],Ainp);#for all ws: v, r such that p^r(A)(v)=0 and p^(r-1)(A)(v)
    od;
    dims := []; #dimensions of my cyclic subspaces
    sumdim := 0; #dimension of my direct sum
    conj := ZeroMatrix(F,n,n);
    while not Sum(minpolpowers, function(v) return v[2]*d; end) = n do #generally bigger than n, working our way down  
        SortBy(minpolpowers, function(v) return v[2]; end); #Sort my ws by power of their minimal polynomial (ascending)
        vecs := ZeroVector(F,Length(minpolpowers));
        for i in [1..Length(minpolpowers)] do
            vecs[i] := minpolpowers[i][3]; #p(A)^(r-1)(w)
        od;
        qis := FindLinearDependenceNC(vecs,A,d); #coefficients of qis as described in theorem
        j := 0;
        for i in [1..Length(qis)] do #find j as described in theorem (works because we sorted this list beforehand)
            if not IsZero(qis[i]) then
                j := i;
                break;
            fi;
        od;
        r := minpolpowers[j][2]; #r as described in theorem
        wstrich := ZeroVector(F,n);
        for i in [1..Length(qis)] do #calculate new reduced w_j
            if not IsZero(qis[i]) then #avoid unnecessary computations
                tomult := minpolpowers[i][1] * qis[i](A); #TODO: GANZ DRINGEND DURCH NOFOMA MATPOLVEC ERSETZEN
                for k in [1..minpolpowers[i][2] - r] do
                    tomult := tomult * Ainp;
                od;
                wstrich := wstrich + tomult;
                #wstrich := wstrich + minpolpowers[i][1] * PolToMat(A,qis[i]) * Ainp^(minpolpowers[i][2] - r); #Alles an den Vektor dranmultiplizieren statt Matrixmultiplikation
            fi;
        od;
        if IsZero(wstrich) then 
            Remove(minpolpowers,j); #if its zero vec we get no more information from it
        else 
            minpolpowers[j] := GetMinPolPowerWithVec(A,p,m,wstrich,Ainp);
        fi;
    od;
    currdim := 1;
    SortBy(minpolpowers, function(v) return v[2]; end); #sort by length
    for i in [1..Length(minpolpowers)] do
        wtrip := minpolpowers[i];
        conj{[currdim..currdim+(wtrip[2]*d)-1]}{[1..n]} := SpinUntil(wtrip[1],A,wtrip[2]*d);
        Add(dims, wtrip[2]*d);
        currdim := currdim + wtrip[2]*d;
    od;
    return [conj,dims];
end);

#Input: Cyclic matrix A with minimal polynomial p^m
#Returns matrix B such that A^Inverse(B) is in Jordan block form
InstallGlobalFunction(JordanBlock, function(A, p, m) #For JordanNormalform
    local F, i, spun, n, b, basis,d,r,v, Ainp, newA, newAinp, Ainps;
    n := NrRows(A);
    d := Degree(p);
    F := DefaultFieldOfMatrix(A);
    spun := FindCyclicVectorNC(A);
    basis := ZeroMatrix(F,n,n);
    CopySubMatrix(spun,basis,[1..d],[1..d],[1..n],[1..n]);
    for r in [1..d] do
        #b := spun[r];
        for i in [1..m-1] do 
            b := PolyEvalFromSpan(ExtractSubMatrix(spun,[r..n],[1..n]),p^i);
            basis[d*i+r]:=b; #Das steht falschrum im skript aber kann ich sicherlich eleganter lösen
        od;
    od;
    return basis;
end);

#Input: Matrix A with irreducible Minimal polynomial
#Returns matrix B such that A^Inverse(B) is in Jordan normal form
InstallGlobalFunction(JordanNormalformIrred, function(A)
    local F,n, cobrank, COB, blockdim, spun, v, w;
    n := Size(A);
    F := DefaultFieldOfMatrix(A); # get underlying field
    v := GenerateRandomVector(F,n);
    spun := SpinUntil(v, A, n);
    blockdim := Rank(spun);  #works bc minpol of vector divides minpol of matrix
    COB := ZeroMatrix(F,n,n);
    CopySubMatrix(spun, COB, [1..blockdim],[1..blockdim],[1..n],[1..n]);
    cobrank := blockdim;
    while not cobrank = n do 
        w := FindVectorNotInSubspaceNC(COB{[1..cobrank]}{[1..n]});
        spun := SpinUntil(w,A,blockdim);
        CopySubMatrix(spun, COB, [1..blockdim],[cobrank+1..cobrank+blockdim],[1..n],[1..n]);
        cobrank := cobrank + blockdim;
    od;
    return COB;
end);

#TODO: CHANGE VARIABLE NAMES
#Input: Matrix A
#Returns matrix B such that A^Inverse(B) is in Jordan normal form 
InstallGlobalFunction(JordanNormalform, function(A) # JordanNormalform with Jordanblocks
    local n, F, d, pol, minpol, collected, sortedfacs, split1, pos, split2, facOccCorrect, facs, hauptraum, hauptraumdims, crhr, cyclicdims, subsubCOB, COB, subA, cy, i, j, facOcc, subCOB, crcy, subsubA, prepreCOB, preCOB;
    F := DefaultFieldOfMatrix(A);
    A := Matrix(F,A);
    n := NrRows(A);
    minpol := MinimalPolynomial(F,A);
    if IsIrreducible(minpol) then 
        return JordanNormalformIrred(A);
    fi;
    if Degree(minpol) = n then 
        hauptraum := PrimaryDecompositionforJNFCyclic(A, minpol); 
    else
        hauptraum := PrimaryDecompositionforJNF(A, minpol);
    fi;
    hauptraumdims := hauptraum[2]; #dimensions of generalized eigenspaces
    facOcc := hauptraum[3]; #factors of minimalpolynomial and their multiplicity
    COB := hauptraum[1]; #Change of basis matrix, this will be the final COB from A to JNF
    A := COB*A*Inverse(COB); #A in hauptraumform
    crhr := 1; #current row (hauptraum)
    preCOB := ZeroMatrix(F,n,n); #COB matrix from A in hauptraumform to haupträume in cyclic form
    for i in [1..Size(hauptraumdims)] do 
        if hauptraumdims[i] = 1 then 
            preCOB[crhr,crhr] := One(F); #could be avoided by directly intializing identity matrix, but i guess it doesnt matter
            crhr := crhr + 1;
            continue;
        fi;
        pol := facOcc[i][1]; 
        d := Degree(pol); #d as in...
        subA := ExtractSubMatrix(A,[crhr..crhr+hauptraumdims[i]-1],[crhr..crhr+hauptraumdims[i]-1]);
        cy := CyclicDecompositionOfPrimarySubspace(subA, facOcc[i][1], facOcc[i][2]); #zerlege haupträume in zyklische unterräume
        cyclicdims := cy[2]; #dimensions of cyclic subspaces
        subCOB := cy[1]; #subCOB to be assembled 
        subA := subCOB*subA*Inverse(subCOB); #subA in cyclic decomposition form
        crcy := 1; #current row (cyclic)
        prepreCOB := ZeroMatrix(F,hauptraumdims[i], hauptraumdims[i]);
        for j in [1..Size(cyclicdims)] do 
            if cyclicdims[j] = 1 then
                prepreCOB[crcy,crcy] := One(F);
                crcy := crcy + 1;
                continue;
            fi;
            subsubA := ExtractSubMatrix(subA, [crcy..crcy+cyclicdims[j]-1], [crcy..crcy+cyclicdims[j]-1]);
            subsubCOB := JordanBlock(subsubA,pol,cyclicdims[j]/d);
            CopySubMatrix(subsubCOB, prepreCOB, [1..cyclicdims[j]], [crcy..crcy+cyclicdims[j]-1], [1..cyclicdims[j]], [crcy..crcy+cyclicdims[j]-1]);
            crcy := crcy + cyclicdims[j];
        od;
        subCOB := prepreCOB * subCOB;
        CopySubMatrix(subCOB, preCOB, [1..hauptraumdims[i]], [crhr..crhr+hauptraumdims[i]-1], [1..hauptraumdims[i]], [crhr..crhr+hauptraumdims[i]-1]);
        crhr := crhr + hauptraumdims[i];
    od;
    return preCOB*COB;
end);


