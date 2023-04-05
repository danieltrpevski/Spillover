TITLE Glutamate-induced glutamate release (GIGR)

UNITS {
    (molar) = (1/liter)
    (mM) = (millimolar)
    (uM) = (micromolar)
}

NEURON {
    POINT_PROCESS gigr
    RANGE K2, K5, KX, KY, KZ, kc, kf, eps, V0, V1, VM2, VM3, V4, VM5
    RANGE n, r1, r2, r3, ru, rb, V, B, Rmax, VmaxA, VmaxN, KN, KA, abg
    RANGE x0, y0, z0, R0, Rin
}

PARAMETER {
    K2 = 0.1 (uM)
    K5 = 1.0 (uM)
    Kd = 0.6 (uM)
    KX = 0.3 (uM)
    KY = 0.2 (uM)
    KZ = 0.1 (uM)
    kc = 0.166e-3 (/ms)
    kf = 0.0166e-3 (/ms)
    eps = 0.0166e-3 (/ms)
    V0 = 0.033e-3 (uM/ms)
    V1 = 0.033e-3 (uM/ms)
    VM2 = 0.1e-3 (uM/ms)
    VM3 = 0.333e-3 (uM/ms)
    V4 = 0.0416e-3 (uM/ms)
    VM5 = 0.5e-3 (uM/ms)

    n = 10000
    r1 = 1000e-3 (/uM /ms)
    r2 = 100e-3 (/ms)
    r3 = 4000e-3 (/ms)
    ru = 100e-3 (/ms)
    rb = 1e-10 (/uM^4 /ms)
    V = 100000 (uM)
    B = 1 (uM)
    Rmax = 13287 (uM)
    VmaxA = 510e-3 (uM/ms)
    VmaxN = 1530e-3 (uM/ms)
    KN = 2 (uM)
    KA = 200 (uM)

    x0 = 0
    y0 = 0
    z0 = 0
    R0 = 0
    Rin = 0
} 

ASSIGNED { 
    v2
    v3
    v5
    alpha
    beta 
    gamma
    abg
}

STATE { x y z R }

BREAKPOINT {
    SOLVE states METHOD derivimplicit
}

INITIAL {  
    x = x0
    y = y0
    z = z0
    R = R0
}

DERIVATIVE states { 
    LOCAL Rv
    fluxes()
    x' = V0 + V1*(R/Rmax) - v2 + v3 + kf*y -kc*x
    y' = v2 - v3 - kf*y
    z' = V4*(R/Rmax) - v5 - eps*z
    abg = alpha*x^4/(beta + gamma*x^4) 
    Rv =  abg - VmaxN*R/(R + KN) - VmaxA*R/(R + KA)
    R' = Rv + Rin
}

PROCEDURE fluxes() {
    UNITSOFF
    v2 = VM2 * x^2/(K2^2 + x^2)
    v3 = VM3 * x^4/(KX^4 + x^4) * y^2/(KY^2 + y^2) * z^4/(KZ^4 + z^4)
    v5 = VM5 * z/(K5 + z) * x^2/(Kd^2 + x^2)

    alpha = n*r1*r3*B*V
    beta = (r2 + r3)*ru/rb
    gamma = r2 + r3 + r1*B 
    UNITSON
}

NET_RECEIVE( rin ) {
    R = R + rin
}


COMMENT
Glutamate-induced dlutamate release mechanism taken from [1].

[1] Larter R. and Glendening Craig M., Glutamate-induced glutamate release: A proposed mechanismfor calcium bursting in astrocytes, : Chaos: An Interdisciplinary Journal of Nonlinear Science 15, 047511 (2005); doi: 10.1063/1.2102467

ENDCOMMENT
