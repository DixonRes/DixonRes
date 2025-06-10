// Dixon Resultant Package for Magma
// Author: Haohai Suo
// E-mail: 1362172421@qq.com
// Released under GPL Version 3  
// Copyright (C) Haohai Suo, 2025";

load "huang.txt";
load "Interpolation.txt";

function GetLinearIndex(exponent, dims)
    index := 0;
    stride := 1;
    for k in [#dims..1 by -1] do
        index +:= exponent[k] * stride;
        stride *:= dims[k];
    end for;
    return index + 1; 
end function;

function BuildCoefficientMatrix(DM, d0, d1, nv, np, gens)
    S := Parent(DM);
    K := BaseRing(BaseRing(S));
    R := CoefficientRing(K);
    nrows := &*d0;
    ncols := &*d1;
    dix := ZeroMatrix(K, nrows, ncols);
    dixR := ZeroMatrix(R, nrows, ncols);
    terms := Monomials(DM);
    coeffs := Coefficients(DM);

    p := 3;
    primes := [];
    for i in [1..np] do
        Append(~primes, R!p);
        p := NextPrime(p);
    end for;
    psi := hom< K -> R | primes >;
    rows_set := {};
    cols_set := {};
    for i in [1..#terms] do
        term := terms[i];
        coeff := coeffs[i];
        ncoeff := K!Numerator(coeff);
        exponents := Exponents(term);
        e0 := [exponents[k] : k in [1..nv]];
        e1 := [exponents[k] : k in [nv + 1..nv + nv]];
        
        valid := true;
        for k in [1..nv] do
            if e0[k] ge d0[k] or e1[k] ge d1[k] then
                valid := false;
                break;
            end if;
        end for;

        if valid then
            row := GetLinearIndex(e0, d0);
            col := GetLinearIndex(e1, d1);
            Include(~rows_set, row);
            Include(~cols_set, col);
            dix[row, col] := ncoeff;
            dixR[row, col] := psi(ncoeff);
        end if;
    end for;
    //print rows,cols;
    rows := Sort(SetToSequence(rows_set));
    cols := Sort(SetToSequence(cols_set));
    dix2 := Submatrix(dix, rows, cols);
    dixR2 := Submatrix(dixR, rows, cols);
//print dixR,dixR2;
    return dix2, dixR2;
end function;

function SGCD(M)
    R := BaseRing(M);
    nrows := NumberOfRows(M);
    ncols := NumberOfColumns(M);

    gcds := [];

    for i in [1..nrows] do
        row := [M[i][j] : j in [1..ncols]];
        g := row[1];
        for j in [2..#row] do
            g := GCD(g, row[j]);
            if IsOne(g) then
                break;
    end if;
        end for;

        if not IsOne(g) then
            Append(~gcds, g);
            for j in [1..ncols] do
                ;//M[i][j] := M[i][j] div g;
            end for;
        end if;
    end for;

    for j in [1..ncols] do
        col := [M[i][j] : i in [1..nrows]];
        g := col[1];
        for i in [2..#col] do
            g := GCD(g, col[i]);
            if IsOne(g) then
                break;
    end if;
        end for;
        
        if not IsOne(g) then
            Append(~gcds, g);
            for i in [1..nrows] do
                ;//M[i][j] := M[i][j] div g;
            end for;
        end if;
    end for;

    return M, gcds;
end function;

function ExtractParameters(ps, vars)
    R := Parent(ps[1]);
    n := Rank(R);
    basis := [R.i : i in [1..n]];
    used := [false : i in [1..n]];

    for f in ps do
        mons := Monomials(f); 
        for m in mons do
            exps := Exponents(m); 
            for i in [1..n] do
                if exps[i] ne 0 then
                    used[i] := true;
                end if;
            end for;
        end for;
    end for;

    vars_set := Seqset(vars);
    pars := [basis[i] : i in [1..n] | used[i] and not basis[i] in vars_set];
    return pars;
end function;

function ExportMatrix(M, filename)
    fp := Open(filename, "w");
    p := Characteristic(BaseRing(BaseRing(M)));
    m := NumberOfRows(M);
    n := NumberOfColumns(M);
    degs := [Degree(M[i][j]) : i in [1..m], j in [1..n]];
    maxd := n * Max(degs);
    fprintf fp, "PMATRIX\n%o %o %o\n", p, m, maxd;
    
    for i := 1 to m do
        for j := 1 to n do
            poly := M[i,j];
            terms := Terms(poly);
            fprintf fp, "P %o %o %o\n", i-1, j-1, #terms;
            
            for term in terms do
                fprintf fp, "T %o %o\n", Degree(term), Integers()!LeadingCoefficient(term);
            end for;
        end for;
    end for;
    delete fp;
    return Sprintf("PML: %o", filename);
end function;

function ExportPivots(M)
    E := EchelonForm(M);
    rowPivots := [];
    for i in [1..Nrows(E)] do
        row := Eltseq(E[i]);
        for j in [1..#row] do
            if row[j] ne 0 then
                Append(~rowPivots, j);
                break;
            end if;
        end for;
    end for;
    rowPivots := [pos : pos in rowPivots | pos ne 0];
    
    ET := EchelonForm(Transpose(M));
    colPivots := [];
    for i in [1..Nrows(ET)] do
        row := Eltseq(ET[i]);
        for j in [1..#row] do
            if row[j] ne 0 then
                Append(~colPivots, j);
                break;
            end if;
        end for;
    end for;
    colPivots := [pos : pos in colPivots | pos ne 0];
    return rowPivots, colPivots;
end function;

function dixon(ps, vars, is_interpolation)

print("\nStep0: Preprocessing\n");

    if #vars ne #ps - 1 then
        error "Number of variables must be equal to the number of equations minus one.";
    end if;

    pars := ExtractParameters(ps, vars);
//print pars;
    print "Moving original variables and polynomials into a new polynomial ring...";
    R0 := Universe(ps);
    K := CoefficientRing(R0);
    nd := &+[Degree(p) : p in ps];
    ds := &*[Degree(p) : p in ps];
    //print ds, #K;
    if ds gt #K-2 then
        ds := #K-2;
    end if;

    allGens := [ R0.i : i in [1..Ngens(R0)] ];

    varIdx := [ Type(v) eq RngMPolElt select Index(allGens, v) else v : v in vars ];
    parIdx := [ Type(p) eq RngMPolElt select Index(allGens, p) else p : p in pars ];

    nv := #varIdx;
    np := #parIdx;

    R0names := Names(R0);
    origNames := [ R0names[i] : i in varIdx ];
    dualNames := [ "~" cat s : s in origNames ];
    paraNames := [ R0names[i] : i in parIdx ];

if false then //true false
    R1  := PolynomialRing(K, 2*nv + np);
    gensR1 := [ R1.i : i in [1..2*nv+np] ];
    x := gensR1[1..nv];
    p := gensR1[2*nv+1..2*nv+np];

    newNames := origNames cat dualNames cat paraNames;
    AssignNames(~R1, newNames);
    images := [];
    for i in [1..Ngens(R0)] do
        if i in varIdx then
            k := Index(varIdx, i);
            Append(~images, x[k]);
        elif i in parIdx then
            k := Index(parIdx, i);
            Append(~images, p[k]);
        else
            Append(~images, R1!0);
        end if;
    end for;
    phi := hom< R0 -> R1 | images >;
    psd := [ phi(f) : f in ps ];

print("\nStep1: Construct the Cancellation Matrix and Dixon Polynomial\n");

    M := ZeroMatrix(R1, nv+1, nv+1);
    for i in [1..nv+1] do
        subst := [ k le i-1 select gensR1[nv+k] else gensR1[k] : k in [1..2*nv+np] ];
        for j in [1..nv+1] do
            M[i,j] := Evaluate(psd[j], subst);
        end for;
    end for;
print "Size of Cancellation Matrix:", nv+1, nv+1;
//print("Cancellation Matrix:");
//print M;
print("Computing Dixon Polynomial...");
    deg1 := Degree(ps[1]);
    if false then //false true is_interpolation eq 1
        time DM := Determinant(M);
        factors := [ gensR1[i] - gensR1[nv+i] : i in [1..nv] ];
        time for f in factors do
            Q := DM div f;
            DM := Q;
        end for;

    else
        rows := [ M[i] : i in [1..nv+1] ];  // nv+1 rows
        newrows := [rows[1]];  // First row unchanged
        time for i in [1..nv] do
            delta := rows[i+1] - rows[i];
            denom := gensR1[i] - gensR1[nv+i];
            newrow := delta div denom;  // Symbolic division (exact)
            Append(~newrows, newrow);
        end for;
        newM := Matrix(newrows);
        time DM := Determinant(newM);
        print("Degree of each variable: ");
        print("Number of terms: ");
        print #Terms(DM);
            
        degs := [(deg1-1)*i-1 : i in [1..nv]] cat [(deg1-1)*(nv+1-i)-1 : i in [1..nv]] cat [(nv+np)*(Degree(ps[1])-1)];
        print  (nv+np)*(Degree(ps[1])-1)+1;
        print Parent(newM[1][1]);
        print Degree(DM);
        //time DM1 := MultiDeterminant(newM,degs);
        //time DM2 := MultiDeterminant_TotalDegreeBound(newM,(nv+np)*(Degree(ps[1])-1)+1);
        //time DM3 := ComputePolyMatrixDet(newM, 2*nv+np, (nv+np)*(Degree(ps[1])-1)+1);
        //print "\n",DM,"\n",DM2,"\n",DM3,"\n";
        //print [Degree(DM, gensR1[i]): i in [1..2*nv+np]];
    end if;

    if np eq 1 then
        P<x> := PolynomialRing(K);
    else
        P := PolynomialRing(K, np);
        AssignNames(~P, paraNames);
    end if;
    FR := FunctionField(P);
    R := PolynomialRing(FR, 2*nv);

    newNames2 := origNames cat dualNames;
    AssignNames(~R, newNames2);

    gensR := [ R.i : i in [1..2*nv] ];
    gensP := [ P.i : i in [1..np] ];
    x := gensR[1..2*nv];
    p := gensP[1..np];

    images := [];
    for i in [1..2*nv] do
        Append(~images, x[i]);
    end for;
    for i in [1..np] do
        Append(~images, R!(p[k]));
    end for;
    phi2 := hom< R1 -> R | images >;

    DM := phi2(DM);

else
    if np eq 1 then
        P<x> := PolynomialRing(K);
    else
        P := PolynomialRing(K, np);
        paraNames := [ R0names[i] : i in parIdx ];
        AssignNames(~P, paraNames);
    end if;
    FR := FunctionField(P);
    R := PolynomialRing(FR, 2*nv);

    origNames := [ R0names[i] : i in varIdx ];
    dualNames := [ "~" cat s : s in origNames ];
    newNames := origNames cat dualNames;
    AssignNames(~R, newNames);

    gensR := [ R.i : i in [1..2*nv] ];
    gensP := [ P.i : i in [1..np] ];
    x := gensR[1..2*nv];
    p := gensP[1..np];

    images := [];
    for i in [1..Ngens(R0)] do
        if i in varIdx then
            k := Index(varIdx, i);
            Append(~images, x[k]);
        elif i in parIdx then
            k := Index(parIdx, i);
            Append(~images, R!(p[k]));
        else
            Append(~images, R!0);
        end if;
    end for;
    phi := hom< R0 -> R | images >;

    psd := [ phi(f) : f in ps ];

print("\nStep1: Construct the Cancellation Matrix and Dixon Polynomial\n");

    M := ZeroMatrix(R, nv+1, nv+1);
    for i in [1..nv+1] do
        subst := [ k le i-1 select gensR[nv+k] else gensR[k] : k in [1..2*nv] ];
        for j in [1..nv+1] do
            M[i,j] := Evaluate(psd[j], subst);
        end for;
    end for;
print "Size of Cancellation Matrix:", nv+1, nv+1;
//print("Cancellation Matrix:");
//print M;
deg1 := Degree(ps[1]);
print("Computing Dixon Polynomial...");
    if false then //false true
        time DM := Determinant(M);
        factors := [ gensR[i] - gensR[nv+i] : i in [1..nv] ];
        time for f in factors do
            Q := DM div f;
            DM := Q;
        end for;
    else
        rows := [ M[i] : i in [1..nv+1] ];  // nv+1 rows
        newrows := [rows[1]];  // First row unchanged
        time for i in [1..nv] do
            delta := rows[i+1] - rows[i];
            denom := gensR[i] - gensR[nv+i];
            newrow := delta div denom;  // Symbolic division (exact)
            Append(~newrows, newrow);
        end for;
        newM := Matrix(newrows);
        time DM := Determinant(newM);
        //print DM;
        print("Degree of each variable: ");
        print  [Degree(DM, gensR[i]): i in [1..2*nv]];
        print("Number of terms: ");
        print #Terms(DM);
            
            degs := [(deg1-1)*i-1 : i in [1..nv]] cat [(deg1-1)*(nv+1-i)-1 : i in [1..nv]];
            //time DM1 := MultiDeterminant(newM,degs);
            //time DM2 := MultiDeterminant_TotalDegreeBound(newM,(nv+np)*(Degree(ps[1])-1)+1);
        //print "\n",DM,"\n",DM1,"\n",DM2,"\n";
    end if;
end if;
    
//print("Dixon Polynomial:");
//print DM;


print("\nStep2: Build the Dixon Matrix\n");

    d0 := [Degree(DM, gensR[i]) + 1 : i in [1..nv]];
    d1 := [Degree(DM, gensR[nv + i]) + 1 : i in [1..nv]];
    //d0 := [Degree(DM, gensR[i]) + 1 : i in [1..nv-1]] cat [Degree(DM, gensR[nv-1]) + 1];
    //d1 := [Degree(DM, gensR[nv + 2]) + 1] cat [Degree(DM, gensR[nv + i]) + 1 : i in [2..nv]];

//print d0,d1; 
    time dix, dixM := BuildCoefficientMatrix(DM, d0, d1, nv, np, gensR);

print "Size of Dixon Matrix: ", NumberOfRows(dixM), NumberOfColumns(dixM);
//print("Dixon Matrix: ");
//print dix;


print("\nStep3: Extract Maximal Rank Submatrix\n");

    time rowPivots, colPivots := ExportPivots(dixM);
    //print rowPivots, colPivots;
    dix2 := Submatrix(dix, colPivots, rowPivots);

    dix2, gcds := SGCD(dix2);

print "Size of Maximal Rank Submatrix: ", NumberOfRows(dix2), NumberOfColumns(dix2);
//print("Maximal Rank Submatrix: ");
//print(dix2);

print("\nStep4: Compute the Dixon Resultant\n");
if is_interpolation eq 1 then
    print("Using Interpolation");
elif is_interpolation eq 0 then
    print("Using Gauss elimination");
else
    print("Using PML");
end if;
extf := [];
time ExportMatrix(dix2, "dixon_matrix");
    if np eq 1 then
        if is_interpolation eq 0 then
            time d1 := Determinant(dix2);
        elif is_interpolation eq 1 then
            time d1 := UniDeterminant(dix2, ds);
        else
            time ExportMatrix(dix2, "dixon_matrix_pml");
            return 1;
        end if;
        d2 := d1/LeadingCoefficient(d1);
        // inv_images := R0.parIdx[1];
        // d2 := Evaluate(d2, inv_images);
//print "Dixon Resultant:", d2;
    else
        if is_interpolation eq 0 then
            time d1 := Determinant(dix2);
        else
            dss := [ds : _ in [1..np]];
            time d1 := MultiDeterminant(dix2, dss);
            //time d1 := ComputePolyMatrixDet(dix2, np, ds+1);
//print d1;
print [Degree(d1, gensP[i]): i in [1..np]];
        end if;
        d2 := d1/LeadingCoefficient(d1);
//print("Dixon Resultant:");
//print d2;
print("\nStep5: Remove Extraneous Factors\n");
        inv_images := [R0.parIdx[k]: k in [1..np]];
        d2 := Evaluate(d2, inv_images);
        time fac := Factorization(d2);
        ps_eval := [];
        evals := [];
        if IsFinite(K) then
            pe := PrimitiveElement(K);
        else
            pe := 2;
        end if;
        for i in [1..Ngens(R0)] do
            if i in varIdx then
                Append(~evals, allGens[i]);
            elif i in parIdx then
                if i eq parIdx[1] then
                    Append(~evals, allGens[i]);
                else
                    Append(~evals, R0!pe^i);
                end if;
            else
                Append(~evals, R0!0);
            end if;
        end for;
        for i in [1..#ps] do
            Append(~ps_eval, Evaluate(ps[i], evals));
        end for;
        I := ideal<R0|ps_eval>;
        G := GroebnerBasis(I);
        g := G[#G];
        d2f := [];

        for i in [1..#fac] do
            fac_eval := Evaluate(fac[i][1], evals);
            if not (Resultant(fac_eval, g, pars[1]) eq 0) then
                Append(~extf, fac[i][1]);
            else
                Append(~d2f, fac[i][1]);
            end if;
        end for;
        if #d2f eq 0 then
            d2 := &*extf;
        else
            d2 := &*d2f;
        end if;
print("Extraneous Factors:");
print #extf;
    end if;
    return d2, gcds, extf;
end function;
