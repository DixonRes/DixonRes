function dixon(ps, vars, pars)

print(1);

    R0 := Universe(ps);
    K := CoefficientRing(R0);
    
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
        S := PolynomialRing(K);
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

    time d1 := Determinant(dix2);

print(7);

    if np eq 1 then
        //return d1;
        inv_images := R0.parIdx[1];
    else
        inv_images := [R0.parIdx[k]: k in [1..np]];
    end if;
    d2 := Evaluate(d1, inv_images);
    d2 := d2/LeadingCoefficient(d2);
    return d2;
end function;


GF2:=FiniteField(2);
GF2_polyring<z8>:=PolynomialRing(GF2);
irr_poly:=z8^8 + z8^4 + z8^3 + z8 + 1;
GF2_n<z8>:=ext<GF2|irr_poly>;
GF2_n_poly_ring<x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10>:=PolynomialRing(GF2_n,11);

p0 := x1^4*x3^4*z8^7 + x2^4*x3^4*z8^7 + x2^4*x3^4*z8^6 + x2^4*x3^2*z8^7 + x2^2*x3^4*z8^7 + x1^4*x3^4*z8^4 + x1^4*x3^2*z8^6 + x1^4*x3*z8^7 + x1^2*x3^4*z8^6 + x1*x3^4*z8^7 + x2^4*x3^4*z8^4 + x2^4*x3^2*z8^6 + x2^4*x3*z8^7 + x2^2*x3^4*z8^6 + x2*x3^4*z8^7 + x1^4*x3^4*z8^3 + x1^4*x3*z8^6 + x1^2*x3^2*z8^7 + x1*x3^4*z8^6 + x2^4*x3^2*z8^5 + x2^2*x3^4*z8^5 + x3^4*z8^7 + x1^2*x3*z8^7 + x1*x3^2*z8^7 + x2^4*x3^4*z8^2 + x2^4*x3*z8^5 + x2^2*x3^2*z8^6 + x2^2*x3*z8^7 + x2*x3^4*z8^5 + x2*x3^2*z8^7 + x1^4*x3^2*z8^3 + x1^4*x3*z8^4 + x1^2*x3^4*z8^3 + x1^2*x3^2*z8^5 + x1^2*x3*z8^6 + x1*x3^4*z8^4 + x1*x3^2*z8^6 + x1*x3*z8^7 + x2^4*x3*z8^4 + x2^2*x3^2*z8^5 + x2*x3^4*z8^4 + x2*x3*z8^7 + x1^2*x3^2*z8^4 + x1^2*x3*z8^5 + x1*x3^2*z8^5 + x2^4*x3*z8^3 + x2*x3^4*z8^3 + x2*x3*z8^6 + x3^4*z8^4 + x3*z8^7 + x1^4*x3^2*z8 + x1^4*x3*z8^2 + x1^2*x3^4*z8 + x1*x3^4*z8^2 + x1*x3*z8^5 + x2^4*x3*z8^2 + x2^2*x3*z8^4 + x2*x3^4*z8^2 + x2*x3^2*z8^4 + x2*x3*z8^5 + x3*z8^6 + x1^4*x3^2 + x1^2*x3^4 + x1*x3*z8^4 + x2^4*x3*z8 + x2^2*x3^2*z8^2 + x2*x3^4*z8 + x3^4*z8^2 + x3^2*z8^4 + x1*x3*z8^3 + x2^2*x3^2*z8 + x3^2*z8^3 + x3*z8^4 + x1^2*x3^2 + x1*x3*z8^2 + x2^2*x3^2 + x3^2*z8^2 + x3*z8^3 + x1*x3*z8 + x3*z8^2 + x2*x3 + x3*z8 + 1;
p1 :=  x1^4*x4^4*z8^7 + x2^4*x4^4*z8^7 + x2^4*x4^4*z8^6 + x1^4*x4^4*z8^5 + x1^4*x4^2*z8^7 + x1^2*x4^4*z8^7 + x2^4*x4^4*z8^5 + x1^4*x4^4*z8^4 + x1^4*x4^2*z8^6 + x1^2*x4^4*z8^6 + x2^4*x4^4*z8^4 + x2^4*x4^2*z8^6 + x2^2*x4^4*z8^6 + x1^4*x4*z8^6 + x1^2*x4^2*z8^7 + x1*x4^4*z8^6 + x2^4*x4^4*z8^3 + x2^4*x4^2*z8^5 + x2^2*x4^4*z8^5 + x1^4*x4^4*z8^2 + x1^4*x4^2*z8^4 + x1^4*x4*z8^5 + x1^2*x4^4*z8^4 + x1^2*x4^2*z8^6 + x1*x4^4*z8^5 + x2^4*x4^2*z8^4 + x2^2*x4^4*z8^4 + x1^4*x4^2*z8^3 + x1^2*x4^4*z8^3 + x1*x4*z8^7 + x2^2*x4*z8^6 + x2*x4^2*z8^6 + x2*x4*z8^7 + x4^2*z8^7 + x1^4*x4^4 + x1^4*x4^2*z8^2 + x1^2*x4^4*z8^2 + x1^2*x4*z8^5 + x1*x4^2*z8^5 + x1*x4*z8^6 + x2^4*x4^4 + x2^4*x4^2*z8^2 + x2^4*x4*z8^3 + x2^2*x4^4*z8^2 + x2^2*x4^2*z8^4 + x2*x4^4*z8^3 + x1^2*x4^2*z8^3 + x1^2*x4*z8^4 + x1*x4^2*z8^4 + x2^4*x4^2*z8 + x2^2*x4^4*z8 + x2^2*x4^2*z8^3 + x4^2*z8^5 + x1^4*x4^2 + x1^2*x4^4 + x1^2*x4^2*z8^2 + x1^2*x4*z8^3 + x1*x4^2*z8^3 + x1*x4*z8^4 + x2^4*x4*z8 + x2^2*x4*z8^3 + x2*x4^4*z8 + x2*x4^2*z8^3 + x4^4*z8^2 + x1^4*x4 + x1^2*x4^2*z8 + x1^2*x4*z8^2 + x1*x4^4 + x1*x4^2*z8^2 + x1*x4*z8^3 + x2^4*x4 + x2^2*x4*z8^2 + x2*x4^4 + x2*x4^2*z8^2 + x4^4*z8 + x1*x4*z8^2 + x4^2*z8^2 + x1^2*x4 + x1*x4^2 + x1*x4*z8 + x2^2*x4 + x2*x4^2 + x4*z8^2 + x1*x4 + x4^2 + x4 + 1;
p2 :=  x5^4*x7^4*z8^7 + x6^4*x7^4*z8^7 + x6^4*x7^4*z8^6 + x6^4*x7^2*z8^7 + x6^2*x7^4*z8^7 + x5^4*x7^4*z8^4 + x5^4*x7^2*z8^6 + x5^4*x7*z8^7 + x5^2*x7^4*z8^6 + x5*x7^4*z8^7 + x6^4*x7^4*z8^4 + x6^4*x7^2*z8^6 + x6^4*x7*z8^7 + x6^2*x7^4*z8^6 + x6*x7^4*z8^7 + x5^4*x7^4*z8^3 + x5^4*x7*z8^6 + x5^2*x7^2*z8^7 + x5*x7^4*z8^6 + x6^4*x7^2*z8^5 + x6^2*x7^4*z8^5 + x5^2*x7*z8^7 + x5*x7^2*z8^7 + x6^4*x7^4*z8^2 + x6^4*x7*z8^5 + x6^2*x7^2*z8^6 + x6^2*x7*z8^7 + x6*x7^4*z8^5 + x6*x7^2*z8^7 + x5^4*x7^2*z8^3 + x5^4*x7*z8^4 + x5^2*x7^4*z8^3 + x5^2*x7^2*z8^5 + x5^2*x7*z8^6 + x5*x7^4*z8^4 + x5*x7^2*z8^6 + x5*x7*z8^7 + x6^4*x7*z8^4 + x6^2*x7^2*z8^5 + x6*x7^4*z8^4 + x6*x7*z8^7 + x7^4*z8^5 + x7^2*z8^7 + x5^2*x7^2*z8^4 + x5^2*x7*z8^5 + x5*x7^2*z8^5 + x6^4*x7*z8^3 + x6*x7^4*z8^3 + x6*x7*z8^6 + x7^4*z8^4 + x7^2*z8^6 + x5^4*x7^2*z8 + x5^4*x7*z8^2 + x5^2*x7^4*z8 + x5*x7^4*z8^2 + x5*x7*z8^5 + x6^4*x7*z8^2 + x6^2*x7*z8^4 + x6*x7^4*z8^2 + x6*x7^2*z8^4 + x6*x7*z8^5 + x7^4*z8^3 + x5^4*x7^2 + x5^2*x7^4 + x5*x7*z8^4 + x6^4*x7*z8 + x6^2*x7^2*z8^2 + x6*x7^4*z8 + x7^2*z8^4 + x7*z8^5 + x5*x7*z8^3 + x6^2*x7^2*z8 + x5^2*x7^2 + x5*x7*z8^2 + x6^2*x7^2 + x5*x7*z8 + x7*z8^2 + x6*x7 + 1;
p3 :=  x5^4*x8^4*z8^7 + x6^4*x8^4*z8^7 + x6^4*x8^4*z8^6 + x5^4*x8^4*z8^5 + x5^4*x8^2*z8^7 + x5^2*x8^4*z8^7 + x6^4*x8^4*z8^5 + x5^4*x8^4*z8^4 + x5^4*x8^2*z8^6 + x5^2*x8^4*z8^6 + x6^4*x8^4*z8^4 + x6^4*x8^2*z8^6 + x6^2*x8^4*z8^6 + x5^4*x8*z8^6 + x5^2*x8^2*z8^7 + x5*x8^4*z8^6 + x6^4*x8^4*z8^3 + x6^4*x8^2*z8^5 + x6^2*x8^4*z8^5 + x5^4*x8^4*z8^2 + x5^4*x8^2*z8^4 + x5^4*x8*z8^5 + x5^2*x8^4*z8^4 + x5^2*x8^2*z8^6 + x5*x8^4*z8^5 + x6^4*x8^2*z8^4 + x6^2*x8^4*z8^4 + x5^4*x8^2*z8^3 + x5^2*x8^4*z8^3 + x5*x8*z8^7 + x6^2*x8*z8^6 + x6*x8^2*z8^6 + x6*x8*z8^7 + x8^2*z8^7 + x5^4*x8^4 + x5^4*x8^2*z8^2 + x5^2*x8^4*z8^2 + x5^2*x8*z8^5 + x5*x8^2*z8^5 + x5*x8*z8^6 + x6^4*x8^4 + x6^4*x8^2*z8^2 + x6^4*x8*z8^3 + x6^2*x8^4*z8^2 + x6^2*x8^2*z8^4 + x6*x8^4*z8^3 + x5^2*x8^2*z8^3 + x5^2*x8*z8^4 + x5*x8^2*z8^4 + x6^4*x8^2*z8 + x6^2*x8^4*z8 + x6^2*x8^2*z8^3 + x8^2*z8^5 + x5^4*x8^2 + x5^2*x8^4 + x5^2*x8^2*z8^2 + x5^2*x8*z8^3 + x5*x8^2*z8^3 + x5*x8*z8^4 + x6^4*x8*z8 + x6^2*x8*z8^3 + x6*x8^4*z8 + x6*x8^2*z8^3 + x8^4*z8^2 + x5^4*x8 + x5^2*x8^2*z8 + x5^2*x8*z8^2 + x5*x8^4 + x5*x8^2*z8^2 + x5*x8*z8^3 + x6^4*x8 + x6^2*x8*z8^2 + x6*x8^4 + x6*x8^2*z8^2 + x8^4*z8 + x5*x8*z8^2 + x8^2*z8^2 + x5^2*x8 + x5*x8^2 + x5*x8*z8 + x6^2*x8 + x6*x8^2 + x8*z8^2 + x5*x8 + x8^2 + x8 + 1;
p4 :=  x1*z8^7 + x1*z8^4 + x1*z8^3 + x0*x1*z8 + x1*z8^2 + 1;
p5 :=  x2*z8^6 + x0*x2*z8^2 + x0*x2*z8 + x2*z8 + 1;
p6 :=  x5*z8^7 + x5*z8^4 + x3*x5*z8 + x4*x5*z8 + x4*x5 + x5*z8 + x5 + 1;
p7 :=  x6*z8^5 + x3*x6*z8^2 + x4*x6*z8^2 + x6*z8^3 + x3*x6*z8 + x4*x6*z8 + x4*x6 + x6 + 1;
p8 :=  x9*z8^6 + x9*z8^3 + x7*x9*z8 + x8*x9*z8 + x9*z8^2 + x8*x9 + x9 + 1;
p9 :=  x10*z8^7 + x10*z8^5 + x7*x10*z8^2 + x8*x10*z8^2 + x10*z8^3 + x7*x10*z8 + x8*x10*z8 + x8*x10 + x10 + 1;
p10 :=  x9^4*z8^5 + x9^2*z8^7 + x10^4*z8^5 + x9^4*z8^4 + x10^2*z8^6 + x10*z8^7 + x9^4*z8^3 + x9^2*z8^5 + x10^2*z8^5 + z8^7 + x9^2*z8^4 + x9*z8^5 + x10^4*z8^2 + x10*z8^5 + z8^6 + x9^4*z8 + x10^4*z8 + x10*z8^4 + x10^4 + x10^2*z8^2 + x10*z8^3 + z8^4 + x9^2*z8 + x9*z8^2 + x9^2 + x9*z8 + x10*z8 + x9;

time d1 := dixon([p4,p5],[x0],[x1,x2]);
time r1 := Resultant(d1,p0,x1);
time r2 := Resultant(d1,p1,x1);
time r3 := Resultant(r1,r2,x2);
time r4 := Resultant(r3,p6,x3);
time r5 := Resultant(r3,p7,x3);
//time r6 := Resultant(r4,r5,x4);
time d2 := dixon([r3,p6,p7],[x3,x4],[x5,x6]);
time r7 := Resultant(p10,p9,x10);
time r8 := Resultant(r7,p8,x9);
time r9 := Resultant(r8,p3,x8);
time r10 := Resultant(r9,p2,x7);
//time r11 := Resultant(d2,r10,x5);
time d3 := dixon([d2,r10],[x5],[x6]);