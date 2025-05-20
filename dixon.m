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
        Ck := Lagrange1D(Grids[d], seq, Xd);
        Append(~C, Ck);
    end for;

    R := P!0;
    for k in [0..deg_d] do
        R +:= C[k+1] * Xd^k;
    end for;
    return R;
end function;


function dixon(ps, vars, pars, is_interpolation)

print(1);

    R0 := Universe(ps);
    K := CoefficientRing(R0);
    ds := &*[Degree(p) : p in ps];
    print ds, #K;
    if ds gt #K then
        ds := #K-2;
    end if;


    allGens := [ R0.i : i in [1..Ngens(R0)] ];

    varIdx := [ Type(v) eq RngMPolElt select Index(allGens, v) else v : v in vars ];
    parIdx := [ Type(p) eq RngMPolElt select Index(allGens, p) else p : p in pars ];

    nv := #varIdx;
    np := #parIdx;

    R  := PolynomialRing(K, 2*nv + np);
    gens := [ R.i : i in [1..2*nv+np] ];
    x := gens[1..nv];
    p := gens[nv+1..nv+np];

    images := [];
    for i in [1..Ngens(R0)] do
        if i in varIdx then
            k := Index(varIdx, i);
            Append(~images, x[k]);
        elif i in parIdx then
            k := Index(parIdx, i);
            Append(~images, p[k]);
        else
            Append(~images, R!0);
        end if;
    end for;
    phi := hom< R0 -> R | images >;

    psd := [ phi(f) : f in ps ];

    M := ZeroMatrix(R, nv+1, nv+1);
    for i in [1..nv+1] do
        subst := [ k le i-1 select gens[nv+np+k] else gens[k] : k in [1..2*nv+np] ];
        for j in [1..nv+1] do
            M[i,j] := Evaluate(psd[j], subst);
        end for;
    end for;

print(2);

    DM := Determinant(M);
    factors := [ gens[i] - gens[nv+np+i] : i in [1..nv] ];
    
    for f in factors do
        Q := DM div f;
        DM := Q;
    end for;

print(3);

    d0 := [ Degree(DM, gens[i])    + 1 : i in [1..nv] ];
    d1 := [ Degree(DM, gens[nv+np+i]) + 1 : i in [1..nv] ];

    cs := [];
    for ex1 in CartesianProduct([ [0..d1[k]-1] : k in [1..nv] ]) do
        poly := DM;
        for k in [1..nv] do
            poly := Coefficient(poly, gens[nv+np+k], ex1[k]);
        end for;
        Append(~cs, poly);
    end for;

    Drows := [];
    for ex0 in CartesianProduct([ [0..d0[k]-1] : k in [1..nv] ]) do
        row := [];
        for poly in cs do
            tmp := poly;
            for k in [1..nv] do
                tmp := Coefficient(tmp, gens[k], ex0[k]);
            end for;
            Append(~row, tmp);
        end for;
        Append(~Drows, row);
    end for;
    dix := Matrix(Drows);

print(4);

    if np eq 1 then
        S<x> := PolynomialRing(K);
        images := [ i eq nv+1 select S.1 else K!0 : i in [1..2*nv+np] ];
    else
        S := PolynomialRing(K, np);
        images := [
        i ge nv+1 and i le nv+np
          select S.(i - nv)
          else K!0
          : i in [1..2*nv+np]
        ];
    end if;

    phi := hom< R -> S | images >;
    dix1 := Matrix(S,
        [ [ phi(dix[i,j]) : j in [1..Ncols(dix)] ]
          : i in [1..Nrows(dix)] ]);

    p := 47;
    primes := [];
    for i in [1..np] do
        Append(~primes, K!p);
        p := NextPrime(p);
    end for;
    psi := hom< S -> K | primes >;
    dixM := Matrix(K,
        [ [ psi(dix1[i,j]) : j in [1..Ncols(dix1)] ]
          : i in [1..Nrows(dix1)] ]);

print(5);

    E := EchelonForm(dixM);
    rowPivots := [];
    for i in [1..Nrows(E)] do
        for j in [1..Ncols(E)] do
            if E[i,j] ne 0 then
                Append(~rowPivots, j);
                break;
            end if;
        end for;
    end for;
    
    ET := EchelonForm(Transpose(dixM));
    colPivots := [];
    for i in [1..Nrows(ET)] do
        for j in [1..Ncols(ET)] do
            if ET[i,j] ne 0 then
                Append(~colPivots, j);
                break;
            end if;
        end for;
    end for;
    
    dix2 := Submatrix(dix1, colPivots, rowPivots);

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
            points := [ pe^i : i in [0..ds] ];
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
                Append(~grids, [ pe^i : i in [0..u]]);
            end for;
            vals := [];
            for point in CartesianProduct([ grids[k] : k in [1..np] ]) do
                reversed_point := [ point[np - i + 1] : i in [1..np] ];
                 M_eval := Matrix(K, Nrows(dix2), Ncols(dix2), 
                          [Evaluate(elt, reversed_point) : elt in Eltseq(dix2)]);
                Append(~vals, Determinant(M_eval)); 
            end for;

            vars :=  [ S.i : i in [1..Ngens(S)] ];
            time d1 := TensorInterpolation(vars, grids, vals, S);
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
R< x0, x1, x2 >  := PolynomialRing(K, 3);
ps := [5772*x0^3 - 7175*x0^2*x1 - 10624*x0*x1^2 + 11290*x1^3 - 19348*x0^2*x2 - 17179*x0*x1*x2 + 26711*x1^2*x2 - 29117*x0*x2^2 + 2909*x1*x2^2 - 31941*x2^3 + 6176*x0^2 - 25584*x0*x1 + 445*x1^2 + 19873*x0*x2 - 1617*x1*x2 + 29544*x2^2 + 11167*x0 - 25730*x1 + 19977*x2 + 22688,
 -9270*x0^3 - 10033*x0^2*x1 + 20435*x0*x1^2 - 4556*x1^3 - 23678*x0^2*x2 - 20624*x0*x1*x2 + 7510*x1^2*x2 + 12862*x0*x2^2 - 18101*x1*x2^2 - 24472*x2^3 + 15656*x0^2 + 7872*x0*x1 + 20425*x1^2 - 19358*x0*x2 - 30017*x1*x2 + 18744*x2^2 - 19*x0 - 13006*x1 + 14675*x2 - 23126];
vars := [x0];
pars := [x1,x2];

time d := dixon(ps,vars,pars,1);

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
