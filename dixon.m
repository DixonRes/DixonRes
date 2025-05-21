function Lagrange1D(nodes, values, X)
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
        return Lagrange1D(Grids[1], Vals, Vars[1]); 
    end if;
    blockSize := &* [ #Grids[i] : i in [1..d-1] ];
    Hd := [];
    for j in [1..#Grids[d]] do
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
print Grids;
        Ck := Lagrange1D(Grids[d], seq, Xd);
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

function dixon(ps, vars, pars, is_interpolation)

print(1);

    R0 := Universe(ps);
    K := CoefficientRing(R0);
    ds := &*[Degree(p) : p in ps];
    print ds, #K;
    if ds gt #K-1 then
        ds := #K-1;
    end if;

    allGens := [ R0.i : i in [1..Ngens(R0)] ];

    varIdx := [ Type(v) eq RngMPolElt select Index(allGens, v) else v : v in vars ];
    parIdx := [ Type(p) eq RngMPolElt select Index(allGens, p) else p : p in pars ];

    nv := #varIdx;
    np := #parIdx;
    if np eq 1 then
        P<x> := PolynomialRing(K);
    else
        P := PolynomialRing(K, np);
    end if;
    FR := FunctionField(P);
    R := PolynomialRing(FR, 2*nv);

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

    M := ZeroMatrix(R, nv+1, nv+1);
    for i in [1..nv+1] do
        subst := [ k le i-1 select gensR[nv+k] else gensR[k] : k in [1..2*nv] ];
        for j in [1..nv+1] do
            M[i,j] := Evaluate(psd[j], subst);
        end for;
    end for;

print(2);

    time DM := Determinant(M);
    factors := [ gensR[i] - gensR[nv+i] : i in [1..nv] ];
    
    time for f in factors do
        Q := DM div f;
        DM := Q;
    end for;

print(3);

    d0 := [Degree(DM, gensR[i]) + 1 : i in [1..nv]];
    d1 := [Degree(DM, gensR[nv + i]) + 1 : i in [1..nv]];
    
    time dix, dixM := BuildCoefficientMatrix(DM, d0, d1, nv, np, gensR);
    print Parent(dixM[1,1]);
    print NumberOfRows(dixM), NumberOfColumns(dixM);
print(4);

    time E := EchelonForm(dixM);
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
    
    time ET := EchelonForm(Transpose(dixM));
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

    print NumberOfRows(dix2), NumberOfColumns(dix2);
print(6);

    if is_interpolation eq 0 then
        time d1 := Determinant(dix2);
        if np eq 1 then
            return d1/LeadingCoefficient(d1);
            inv_images := R0.parIdx[1];
        else
            inv_images := [R0.parIdx[k]: k in [1..np]];
        end if;
        d2 := Evaluate(d1, inv_images);
        d2 := d2/LeadingCoefficient(d2);

        return d2;
    else
        pe := PrimitiveElement(K);
        if np eq 1 then 
            points := [ pe^i : i in [1..ds] ] cat [K!0];
            values := [ ];
            ii := 0;
            time for a in points do
                ii := ii + 1;
                print ii, a;
                N := Evaluate(dix2, a); 
                Append(~values, Determinant(N)); 
            end for;
            time d2 := Interpolation(points, values);
            d2 := d2/LeadingCoefficient(d2);
            return d2;
        else
            uppers := [ds : _ in [1..np]];
            grids := [];
            for u in uppers do
                Append(~grids, [ pe^i : i in [1..u]] cat [K!0] );
            end for;
            vals := [];
            npoint := (ds+1)^np;
            ipoint := 0;
            for point in CartesianProduct([ grids[k] : k in [1..np] ]) do
                reversed_point := [ point[np - i + 1] : i in [1..np] ];
                 M_eval := Matrix(K, Nrows(dix2), Ncols(dix2), 
                          [Evaluate(elt, reversed_point) : elt in Eltseq(dix2)]);
                Append(~vals, Determinant(M_eval)); 
                ipoint := ipoint + 1;
                print ipoint,npoint;
            end for;

            vars :=  [ P.i : i in [1..Ngens(P)] ];
            time d1 := TensorInterpolation(vars, grids, vals, P);
            // time d1 := Determinant(dix2);
            if np eq 1 then
                return d1;
                inv_images := R0.parIdx[1];
            else
                inv_images := [R0.parIdx[k]: k in [1..np]];
            end if;
            d2 := Evaluate(d1, inv_images);
            d2 := d2/LeadingCoefficient(d2);
            return d2;
        end if;
    end if;
print(7);
end function;

K := GF(65537);
R< x0, x1, x2, x3 >  := PolynomialRing(K, 4);
ps := [4748*x0^2 + 11101*x0*x1 - 6109*x1^2 - 31021*x0*x2 + 9318*x1*x2 + 29439*x2^2 - 21354*x0*x3 + 24482*x1*x3 + 21080*x2*x3 - 16286*x3^2 - 28333*x0 + 18483*x1 - 31088*x2 - 12606*x3 - 14615,
 1124*x0^2 + 24743*x0*x1 + 4039*x1^2 - 7804*x0*x2 + 26702*x1*x2 + 15844*x2^2 + 3073*x0*x3 + 275*x1*x3 - 15498*x2*x3 - 4836*x3^2 + 22767*x0 + 4408*x1 - 26504*x2 + 32131*x3 + 16732,
 6706*x0^2 + 1956*x0*x1 - 13851*x1^2 - 7861*x0*x2 - 8872*x1*x2 + 27248*x2^2 + 28478*x0*x3 - 13449*x1*x3 - 12087*x2*x3 + 31514*x3^2 + 29552*x0 + 25623*x1 - 12628*x2 - 14744*x3 - 30708];
vars := [x0,x1];
pars := [x3,x2];
is_interpolation := 1;
time d := dixon(ps, vars, pars, is_interpolation);

//rs := Roots(d);

I := ideal<R|ps>;
//SetVerbose("Groebner", 1);
IsZeroDimensional(I);
Dimension(I);
time G := GroebnerBasis(I);
g := G[#G];
// v := Variety(I);
time r := Resultant(ps[2],ps[1],x0);
r := r/LeadingCoefficient(r);
