// Dixon Resultant Package for Magma
// Author: Haohai Suo
// E-mail: 1362172421@qq.com
// Released under GPL Version 3  
// Copyright (C) Haohai Suo, 2025";

function Lagrange(nodes, values, X)
    k := #nodes;
    L := Parent(values[1])!0;
    for j in [1..k] do
        num := Parent(values[1])!1;
        den := Parent(nodes[1])!1;
        for m in [1..k] do
            if m ne j then
                num *:= (X - nodes[m]);
                den *:= (nodes[j] - nodes[m]);
            end if;
        end for;
        L +:= values[j] * num * den^-1;
    end for;
    return L;
end function;

function TensorInterpolation(Vars, Grids, Vals, P)
    d := #Vars;
    if d eq 1 then
        return Lagrange(Grids[1], Vals, Vars[1]); 
    end if;
    blockSize := &* [ #Grids[i] : i in [1..d-1] ];
    Hd := [];
    
    s := [];
    for j in [1..#Grids[d]] do
        for _ in [1..#s] do
            printf "\b";
        end for;
        s := Sprint(j) cat "/" cat Sprint(#Grids[d]) cat " ";
        printf s;
        subVals := Vals[(j-1)*blockSize + 1 .. j*blockSize];
        H_j := TensorInterpolation(
                   Vars[1..d-1],
                   Grids[1..d-1],
                   subVals,
                   P
               );
        Append(~Hd, H_j);
    end for;
    Xd := Vars[d];
    deg_d := Max([ Degree(h, Xd) : h in Hd ]);
    C := [];
    for k in [0..deg_d] do
        seq := [ Coefficient(Hd[j], Xd, k) : j in [1..#Grids[d]] ];
        Ck := Lagrange(Grids[d], seq, Xd);
        Append(~C, Ck);
    end for;

    R := P!0;
    for k in [0..deg_d] do
        R +:= C[k+1] * Xd^k;
    end for;
    return R;
end function;

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
            dix[row, col] := ncoeff;
            dixR[row, col] := psi(ncoeff);
        end if;
    end for;
    
    return dix, dixR;
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
print("Computing Dixon Polynomial...");
    if false then
        time DM := Determinant(M);
        factors := [ gensR[i] - gensR[nv+i] : i in [1..nv] ];
        for f in factors do
            Q := DM div f;
            DM := Q;
        end for;
    else
        rows := [ M[i] : i in [1..nv+1] ];  // nv+1 rows
        newrows := [rows[1]];  // First row unchanged
        for i in [1..nv] do
            delta := rows[i+1] - rows[i];
            denom := gensR[i] - gensR[nv+i];
            newrow := delta div denom;  // Symbolic division (exact)
            Append(~newrows, newrow);
        end for;
        newM := Matrix(newrows);
        time DM := Determinant(newM);
    end if;
//print M,newM;

//print("Dixon Polynomial:");
//print DM;

print("\nStep2: Build the Dixon Matrix\n");

    d0 := [Degree(DM, gensR[i]) + 1 : i in [1..nv]];
    d1 := [Degree(DM, gensR[nv + i]) + 1 : i in [1..nv]];
print("Degree of each variable: ");
print d0,d1; 
    time dix, dixM := BuildCoefficientMatrix(DM, d0, d1, nv, np, gensR);
print "Size of Dixon Matrix: ", NumberOfRows(dixM), NumberOfColumns(dixM);
//print("Dixon Matrix: ");
//print dix;


print("\nStep3: Extract Maximal Rank Submatrix\n");

    E := EchelonForm(dixM);
    rowPivots := [];
    time for i in [1..Nrows(E)] do
        row := Eltseq(E[i]);
        for j in [1..#row] do
            if row[j] ne 0 then
                Append(~rowPivots, j);
                break;
            end if;
        end for;
    end for;
    rowPivots := [pos : pos in rowPivots | pos ne 0];
    
    ET := EchelonForm(Transpose(dixM));
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

    dix2 := Submatrix(dix, colPivots, rowPivots);

    dix2, gcds := SGCD(dix2);

print "Size of Maximal Rank Submatrix: ", NumberOfRows(dix2), NumberOfColumns(dix2);
//print("Maximal Rank Submatrix: ");
//print(dix2);

print("\nStep4: Compute the Dixon Resultant\n");
if is_interpolation eq 1 then
    print("Using Interpolation");
else
    print("Using Gauss elimination");
end if;
    if IsFinite(K) then
        pe := PrimitiveElement(K);
    else
        pe := 2;
    end if;
    extf := [];
    if np eq 1 then
        if is_interpolation eq 0 then
            time d1 := Determinant(dix2);
        else
            points := [ pe^i : i in [0..ds] ] cat [K!0];
            values := [ ];
            npoint := ds+2;
            ipoint := 0;
            s := [];
print "Computing Interpolation Points...";
            time for a in points do
                ipoint := ipoint + 1;
                for j in [1..#s] do
                    printf "\b";
                end for;
                s := Sprint(ipoint) cat "/" cat Sprint(npoint) cat " ";
                printf s;
                N := Evaluate(dix2, a); 
                Append(~values, Determinant(N)); 
            end for;
print "Interpolating...";
            d1 := Interpolation(points, values);
        end if;
        d2 := d1/LeadingCoefficient(d1);
        // inv_images := R0.parIdx[1];
        // d2 := Evaluate(d2, inv_images);
//print "Dixon Resultant:", d2;
    else
        if is_interpolation eq 0 then
            time d1 := Determinant(dix2);
        else
            uppers := [ds : _ in [1..np]];
            grids := [];

            for u in uppers do
                Append(~grids, [ pe^i : i in [0..u]] cat [K!0] );
            end for;
//print grids;
            vals := [];
            npoint := (ds+2)^np;
            ipoint := 0;
            s := [];
print "Computing Interpolation Network...";
            time for point in CartesianProduct([ grids[k] : k in [1..np] ]) do
                reversed_point := [ point[np - i + 1] : i in [1..np] ];
                 M_eval := Matrix(K, Nrows(dix2), Ncols(dix2), 
                          [Evaluate(elt, reversed_point) : elt in Eltseq(dix2)]);
                Append(~vals, Determinant(M_eval)); 
                ipoint := ipoint + 1;
                for j in [1..#s] do
                    printf "\b";
                end for;
                s := Sprint(ipoint) cat "/" cat Sprint(npoint) cat " ";
                printf s;
            end for;
            vars :=  [ P.i : i in [1..Ngens(P)] ];
print "Interpolating...";
            time d1 := TensorInterpolation(vars, grids, vals, P);
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
print extf;
    end if;
    return d2, gcds, extf;
end function;

//SetColumns(170);

K := GF(65537);
R< x0, x1, x2, x3>  := PolynomialRing(K, 4);
ps := [3*x0^2*x2^2*x3 + 5*x0*x1*x2*x3^2 - 2*x2^4 - 5*x0^2*x2 + 2*x1^2*x2 + 4*x0*x1 - x0 - 5,
 5*x0^3*x1*x2 - 4*x1^4*x2 + 3*x0^2*x1*x2^2 + 5*x1*x2^2*x3^2 - x2^3*x3^2 - 4*x0*x1*x2^2 + 2*x0*x1*x2 - 3*x0*x3 + 5*x2*x3,
 -x0^5 + 4*x0^2*x1^3 + 5*x0^3*x1*x3 - 5*x0^3*x2*x3 + x0*x2^3*x3 + 3*x0*x1*x3^3 + 3*x1*x3^2 + 4*x0^2 + 4*x2];
vars := [x0,x1];
is_interpolation := 1;
time d := dixon(ps, vars, is_interpolation);
