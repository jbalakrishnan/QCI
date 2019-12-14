################################################################################
# Local heights code and helper functions
################################################################################


def int_omega(C,P):
    x,y = C.local_coord(P)
    b = P[1]
    a = P[0]
    omega = (x.derivative()*b/(y*(x-a)))
    return (omega- omega.parent().gens()[0]^-1).integral() - log(2*b)

def prune(D):
    newlist=[]
    p = D[-1][0].parent().order()
    GFp = GF(p)
    for x in D:
        if (x[0],-x[1]) not in newlist:
            newlist.append(x)
    return newlist

def find_Q(C,P,i=0):
    try:
        p = P[1].parent().prime()
    except AttributeError:
        p = P[1].parent().order()
    while True:
        try:
            Q = C.lift_x(ZZ(P[0]) + i*p)
            break
        except ValueError:
            i += 1
    try:
        if ZZ(Q[1][0]) != ZZ(P[1][0]):
            Q = C(Q[0],-Q[1])
    except TypeError:
        if ZZ(Q[1][0]) != ZZ(P[1]):
            Q = C(Q[0],-Q[1])
    return Q


def height_minus(C, P, prec):
    """
    this gives the normalized local height between -P and P
    """
    negP = C(P[0],-P[1])
    div1 = [(1,P),(-1,negP)]
    int_eta = C.init_height(div1,div1,prec)
    int_eta = C.eta_integral(div1,div1)
    Q = find_Q(C,P,1)
    negQ = C(Q[0],-Q[1])
    int_omega_at_P = int_omega(C,P)
    int_omega_P_to_Q = int_omega_at_P(Q[0]-P[0]) + log(Q[0]-P[0],0)
    C.init_height([(1,P),(-1,negP)], [(1,Q),(-1,negQ)])
    int_omega_wQ_to_Q = C.omega_integral([(1,P),(-1,negP)], [(1,Q),(-1,negQ)])
    return -2*int_omega_P_to_Q + int_omega_wQ_to_Q - int_eta

def X_to_E1map(Xpt,E1):
    if Xpt[2] == 0:
        pass
    else:
        x,y,_ = Xpt
        return E1(x^2,y)

def X_to_E2map(Xpt,E2):
    if Xpt[2] == 0:
        return E2(0,1)
    else:
        x,y,_ = Xpt
        if x !=0:
            return E2(x^-2,y*x^-3)
        else:
            return E2(0,1,0)

def X_to_E2map_plus(Xpt,E2):
    if Xpt[2] == 0:
        pass
    else:
        x,y,_ = Xpt
        if x !=0:
            return E2(x^-2,y*x^-3) + E2(0,1)
        else:
            return E2(0,1,0) +E2(0,1)

def X_to_E2map_minus(Xpt,E2):
    if Xpt[2] == 0:
        pass
    else:
        x,y,_ = Xpt
        if x !=0:
            return E2(x^-2,y*x^-3) - E2(0,1)
        else:
            return E2(0,1,0) - E2(0,1)

def X_to_E1map_plus(Xpt,E1):
    if Xpt[2] == 0:
        return E1((0,1))
    else:
        x,y,_ = Xpt
        return E1(x^2,y) + E1.lift_x((0,1))

def X_to_E1map_minus(Xpt,E1):
    if Xpt[2] == 0:
        return  E1((0,-1))
    else:
        x,y,_ = Xpt
        return  E1(x^2,y) + E1.lift_x((0,-1))

def local_height_at_p(EK,P,prec):
    #this assumes EK is over the padics
    hmin = height_minus(EK,P,prec)
    b = P[1]
    tau_pt = 1/4*log(4*b^2) + 1/4*hmin
    return tau_pt

def param_f2Xplus(x,y,E2K):
    """
    gives the x-coordinate of  local parameter on E2K at the point f2(z) + (0,1)
    uses the local coordinate (x,y) at z
    """
    x1 = x^-2
    x2 = 0
    y1 = x^-3*y
    y2 = 1
    a = E2K.a2()
    lam = (y1-y2)/(x1-x2)
    return lam^2 - a - x1 - x2

def param_f2Xmin(x,y,E2K):
    """
    gives the x-coordinate of the local parameter on E2K at hte point f2(z) - (0,1) = f2(z) + (0,-1)
    uses the local coordinate (x,y) at z
    """
    x1 = x^-2
    x2 = 0
    y1 = x^-3*y
    y2 = -1
    a = E2K.a2()
    lam = (y1-y2)/(x1-x2)
    return lam^2 - a - x1 - x2

def param_f1(x,y):
    """
    gives the x-coordinate of the local parameter on E1K at the point f1(z)
    uses the local coordinate (x,y) at z
    """
    return x^2

def param_f2(x,y):
    """
    gives the x-coordinate of the local parameter on E2K at the point f2(z)
    uses the local coordinate (x,y) at z
    """
    return x^-2

def param_f1Xplus(x,y,E1K):
    """
    gives the x-coordinate of  local parameter on E1K at the point f1(z) + (0,1)
    uses the local coordinate (x,y) at z
    """
    x1 = x^2
    y1 = y
    x2 = 0
    y2 = 1
    if y1.valuation() < 0:
        Kt.<t> = LaurentSeriesRing(Qp(p,prec))
    else:
        pass
    a = E1K.a2()
    try:
        lam = (y1-y2)/(x1-x2)
    except TypeError:
        y1 = Kt(y1)
        x1 = Kt(x1)
        lam = Kt((y1-y2)/(x1-x2))
        a = Kt(a)
    return lam^2 - a - x1 - x2

def param_f1Xmin(x,y,E1K):
    """
    gives the x-coordinate of the local parameter on E1K at hte point f1(z) + (0,-1)
    uses the local coordinate (x,y) at z
    """
    x1 = x^2
    y1 = y
    x2 = 0
    y2 = -1
    if y1.valuation() < 0:
        Kt.<t> = LaurentSeriesRing(Qp(p,prec))
    else:
        pass
    a = E1K.a2()
    try:
        lam = (y1-y2)/(x1-x2)
    except TypeError:
        y1 = Kt(y1)
        x1 = Kt(x1)
        a = Kt(a)
        lam = Kt((y1-y2)/(x1-x2))
    return lam^2 - a - x1 - x2

def solver(X, XK, pt, xx, yy, rho_z, T):
    points=[]
    D2zunshift = rho_z
    for c in T:
        print('beta = ', c)
        D2z = D2zunshift - c
        val = min([x.valuation() for x in D2z.list()])
        D2z = D2z*p**(-val)
        if D2z.valuation() > 0:
            t = D2z.parent().gen()
            try:
                D2z = (D2z/t**D2z.valuation()).power_series().polynomial()
            except AttributeError:
                D2z = (D2z/t**D2z.valuation()).polynomial()
            roots = gp.polrootspadic(D2z,p,prec-2)
            roots_new = [(sage_eval('%s'%roots[i+1])).add_bigoh(prec-5) for i in range(len(roots))]
            roots_new = [0] + [x for x in roots_new if x.valuation() > 0]
        else:
            try:
                roots = gp.polrootspadic(D2z,p,prec-2)
            except TypeError:
                try:
                    D2z = D2z.power_series().polynomial()
                except AttributeError:
                    D2z = D2z.polynomial()
                roots = gp.polrootspadic(D2z,p,prec-2)
            roots_new = [(sage_eval('%s'%roots[i+1])).add_bigoh(prec-5) for i in range(len(roots))]
            roots_new = [x for x in roots_new if x.valuation() > 0]
        for r in roots_new:
            #first handle the non-Weierstrass case
            if pt[1] != 0 and pt[2] == 1:
                r = XK(xx(r),yy(r))
                try:
                    Xpt = X(r[0],r[1])
                    ('Found point on X: ', Xpt)
                    points = points + [Xpt]
                except TypeError:
                    print('Found point on X: ', r)
                    points = points + [r]
            elif pt[1] == 0:
                print('Finite Weierstrass case! Note that the example curve and prime did not run into this case.')
            else:
                #this is the disk at infinity
                if r != 0:
                    try:
                        Xpt = X(xx(r),yy(r))
                        print('Found point on X: ', Xpt)
                        points = points + [Xpt]
                    except (TypeError, ValueError):
                        Xpt = XK(xx(r),yy(r))
                        print('Found point on X: ', Xpt)
                        points = points + [Xpt]
                else:
                    print('Found point on X: (1: 1: 0)')
                    points = points + ['infty']
    return points

################################################################################
# Quadratic Chabauty for bielliptic genus 2 curves over QQ
################################################################################

def qcb(X, E1, E2, H1K, H2K, T0list, Tnot0list, p, prec):
    """
    The main routine, which first handles the residue disks where x = 0 and then the disks where x != 0.
    """
    assert E1.rank() == 1 and E2.rank() == 1
    prechere = prec
    K = Qp(p,prec)
    XK = X.change_ring(K)
    E1K = E1.change_ring(K)
    E2K = E2.change_ring(K)

    F = GF(p)
    E1p = E1.change_ring(F)
    E2p = E2.change_ring(F)
    E1m = E1.minimal_model()
    E2m = E2.minimal_model()
    E1mp = E1m.change_ring(F)
    E2mp = E2m.change_ring(F)

    h1 = E1m.padic_height(p)
    h2 = E2m.padic_height(p)
    z1 = E1m.gens()[0]
    z2 = E2m.gens()[0]
    gE1mtoE1 = E1m.isomorphism_to(E1)
    gE2mtoE2 = E2m.isomorphism_to(E2)
    g2 = gE2mtoE2(z2)
    g2 = H2K(g2[0],g2[1])
    g1 = gE1mtoE1(z1)
    g1 = H1K(g1[0],g1[1])
    inf1 = H1K(0,1,0)
    inf2 = H2K(0,1,0)

    alpha2 =  h2(z2)/(H2K.coleman_integrals_on_basis(inf2,g2)[0])^2
    alpha1 =  h1(z1)/(H1K.coleman_integrals_on_basis(inf1, g1)[0])^2

    Xp = X.change_ring(F)

    from sage.schemes.hyperelliptic_curves.monsky_washnitzer import matrix_of_frobenius_hyperelliptic
    N1 = H1K.cup_product_matrix()
    M1,_ = matrix_of_frobenius_hyperelliptic(H1K)
    M1prec = M1^prec
    M1prec1 = M1prec.columns()[1]
    A1 = matrix(2,2,[1,M1prec1[0],0,M1prec1[1]])
    A1inv = A1^-1
    a = A1inv.list()[1]
    w0dualH1 = vector([a,-1])

    N2 = H2K.cup_product_matrix()
    M2,_ = matrix_of_frobenius_hyperelliptic(H2K)
    M2prec = M2^prec
    M2prec1 = M2prec.columns()[1]
    A2 = matrix(2,2,[1,M2prec1[0],0,M2prec1[1]])
    A2inv = A2^-1
    a = A2inv.list()[1]
    w0dualH2 = vector([a,-1])
    points=[]

    #This is just one residue disk with x = 0.
    D0 = [Xp.lift_x(0)]
    for P in D0:
        Ito01 = H2K.coleman_integrals_on_basis(H2K(0,1,0),H2K(0,1))[0]
        if P[1] == 0:
            print('Note: the example curve did not have any Weierstrass disks! If you end up here, add Weierstrass local coordinates!')
        Q = find_Q(XK,P)
        xx,yy = XK.local_coord(Q)
        f2plus = X_to_E2map_plus(Q,E2K)
        x_in_disk_of_f2plus = param_f2Xplus(xx,yy,E2K)
        f2min = X_to_E2map_minus(Q,E2K)
        x_in_disk_of_f2min = param_f2Xmin(xx,yy,E2K)
        f2plus = H2K(f2plus[0],f2plus[1])
        f2min = H2K(f2min[0],f2min[1])
        f2 = X_to_E2map(Q,E2K)
        x_in_disk_of_f2 = param_f2(xx,yy)
        f2 = H2K(f2[0],f2[1],f2[2])
        f1 = X_to_E1map(Q,E1K)
        x_in_disk_of_f1 = param_f1(xx,yy)
        f1 = H1K(f1[0],f1[1])
        tauf2plus = local_height_at_p(H2K, f2plus, prec)
        tauf2min = local_height_at_p(H2K, f2min, prec)
        tauf1 = local_height_at_p(H1K, f1,prec)
        tauzlist = []
        for pt in [f2plus, f2min, f1]:
            if pt == f2plus or pt == f2min:
                Y = H2K
            elif pt == f1:
                Y = H1K
            PQ_sing = vector(Y.tiny_integrals_on_basis_to_z(pt))
            PQ_doub = Y.tiny_double_integrals_on_basis_to_z(pt)
            if pt[2] != 0:
                bP_sing = Y.coleman_integrals_on_basis(Y(0,1,0),pt)
            else:
                bP_sing = PQ_sing
            if pt[1] != 0 and pt[2] != 0:
                bz_sing = Y.coleman_integrals_on_basis_to_z(Y(0,1,0),pt)
            elif pt[1] == 0 or pt[2] == 0:
                bz_sing = PQ_sing
            if pt == f2plus or pt == f2min:
                w0bar = w0dualH2
            elif pt == f1:
                w0bar = w0dualH1
            diff = -2*(vector(PQ_doub[0:2])*w0bar + PQ_sing[0]*(bP_sing*w0bar))
            if pt == f2plus:
                ##assumes f2plus is non-weier
                tnewfplus = x_in_disk_of_f2plus - f2plus[0]
                tau_z = tauf2plus + diff(tnewfplus)
            elif pt == f2min:
                ##assumes f2min is non-weier
                tnewfmin = x_in_disk_of_f2min - f2min[0]
                tau_z = tauf2min + diff(tnewfmin)
            elif pt == f1:
                ##assumes f1 is non-weier
                tnewf1 = x_in_disk_of_f1 - f1[0]
                tau_z = tauf1 + diff(tnewf1)
            tauzlist.append(tau_z)
        tauzf2plus = tauzlist[0]
        tauzf2min = tauzlist[1]
        tauzf1 = tauzlist[2]
        if P[0] == 0:  #to handle disk of x = 0
            xhere = xx^-2
            yhere = xx^-3*yy
            n31 = (xhere.derivative()/(2*yhere)).integral().power_series()
        else:
            x_in_disk_of_f2 = param_f2(xx,yy)
            tnewf2 = x_in_disk_of_f2 - f2[0]
            n31 = H2K.coleman_integrals_on_basis_to_z(H2K(0,1,0),f2)[0](tnewf2)
        V = n31.parent()
        alpha2term = 2*(n31^2 +Ito01^2)*V(alpha2)
        n51 = H1K.coleman_integrals_on_basis_to_z(H1K(0,1,0),f1)[0](tnewf1)
        V = n51.parent()
        alpha1term = 2*n51^2*V(alpha1)
        rho_z = 2*tauzf1 - tauzf2plus - tauzf2min + alpha2term  - alpha1term
        newpts = solver(X, XK, Q, xx, yy, rho_z, T0list)
        points = points + newpts

    #Now the other residue disks
    D = Xp.rational_points()
    D.remove(Xp(0,1))
    D.remove(Xp(0,-1))
    D = prune(D)
    for P in D:
        Ito01 = H1K.coleman_integrals_on_basis(H1K(0,1,0),H1K(0,1))[0]
        if P[1] == 0:
            print('Note: the example curve did not have any finite Weierstrass disks! If you end up here, add Weierstrass local coordinates!')
        if P[2] != 0:
            Q = find_Q(XK,P)
            xx,yy = XK.local_coord(Q)
        else:
            Q = XK(0,1,0)
            xx,yy = XK.local_coordinates_at_infinity()
            t = xx.parent().gens()[0]
        f1plus = X_to_E1map_plus(Q,E1K)
        f1plus = H1K(f1plus[0],f1plus[1])
        x_in_disk_of_f1_plus = param_f1Xplus(xx,yy, E1K)
        tauf1plus = local_height_at_p(H1K,f1plus,prec)

        f1min = X_to_E1map_minus(Q,E1K)
        f1min = H1K(f1min[0],f1min[1])
        x_in_disk_of_f1_minus = param_f1Xmin(xx,yy, E1K)
        tauf1min = local_height_at_p(H1K,f1min,prec)

        f2 = X_to_E2map(Q,E2K)
        f2 = H2K(f2[0],f2[1])
        x_in_disk_of_f2 = param_f2(xx,yy)
        tauf2 = t4 = local_height_at_p(H2K, f2,prec)
        if Q!= XK(0,1,0):
            f1 = X_to_E1map(Q,E1K)
            x_in_disk_of_f1 = param_f1(xx,yy)
            f1 = H1K(f1[0],f1[1])
        try:
            x_in_disk_of_f2 = (x_in_disk_of_f2).power_series()
        except AttributeError:
            pass
        tauzlist = []
        for pt in [f1plus, f1min, f2]:
            if pt == f1plus or pt == f1min:
                Y = H1K
            elif pt == f2:
                Y = H2K
            PQ_sing = vector(Y.tiny_integrals_on_basis_to_z(pt))
            PQ_doub = Y.tiny_double_integrals_on_basis_to_z(pt)
            if pt[2] != 0:
                bP_sing = Y.coleman_integrals_on_basis(Y(0,1,0),pt)
            else:
                bP_sing = PQ_sing
            if pt[1] != 0 and pt[2] != 0:
                bz_sing = Y.coleman_integrals_on_basis_to_z(Y(0,1,0),pt)
            elif pt[1] == 0 or pt[2] == 0:
                bz_sing = PQ_sing
            if pt == f1plus or pt == f1min:
                w0bar = w0dualH1
            elif pt == f2:
                w0bar = w0dualH2
            diff = -2*(vector(PQ_doub[0:2])*w0bar + PQ_sing[0]*(bP_sing*w0bar))
            if pt == f1plus:
                try:
                    tnewf1plus = (x_in_disk_of_f1_plus - f1plus[0]).power_series()
                except (TypeError, AttributeError):
                    tnewf1plus = x_in_disk_of_f1_plus - f1plus[0]
                tau_z = tauf1plus + diff(tnewf1plus)
            elif pt == f1min:
                try:
                    tnewf1minus = (x_in_disk_of_f1_minus - f1min[0]).power_series()
                except (TypeError, AttributeError):
                    tnewf1minus = x_in_disk_of_f1_minus - f1min[0]
                tau_z = tauf1min + diff(tnewf1minus)
            elif pt == f2:
                tnewf2 = x_in_disk_of_f2 - f2[0]
                tau_z = t4 + diff(tnewf2)
            tauzlist.append(tau_z)
        n51 = H2K.coleman_integrals_on_basis_to_z(H2K(0,1,0),f2)[0](tnewf2)
        V = n51.parent()
        alpha2term = 2*n51^2*V(alpha2)
        if Q[2] != 0:
            tnewf1 = x_in_disk_of_f1 - f1[0]
            n31 = H1K.coleman_integrals_on_basis_to_z(H1K(0,1,0),f1)[0](tnewf1)
        else:
            xhere = xx^2
            yhere = yy
            n31 = (xhere.derivative()/(2*yhere)).integral().power_series()
        V = n31.parent()
        alpha1term = 2*(n31^2 +Ito01^2)*V(alpha1)
        tauzf1plus = tauzlist[0]
        tauzf1min = tauzlist[1]
        tauzf2 = tauzlist[2]
        rho_z = 2*tauzf2 -tauzf1plus - tauzf1min + alpha1term - alpha2term
        newpts = solver(X, XK, Q, xx, yy, rho_z, Tnot0list)
        points = points + newpts
    return points

################################################################################
# Example 1 (Section 8.3)
################################################################################

R.<x> = QQ['x']
a = -2
b = -1
X = HyperellipticCurve(x^6 + a*x^4 + b*x^2 + 1)
Xpts = [ X(0,-1),X(0,1),X(-3/2,-1/8),X(-3/2,1/8),X(3/2,-1/8),X(3/2,1/8)]
p = 3
prec = 15
K = Qp(p,prec)
E1 = EllipticCurve([0,a,0,b,1])
E2 = EllipticCurve([0,b,0,a,1])
E1K = HyperellipticCurve(x^3 + a*x^2 + b*x + 1).change_ring(K)
E2K = HyperellipticCurve(x^3 + b*x^2 + a*x + 1).change_ring(K)
points = qcb(X,E1,E2,E1K,E2K,[-8/3*log(K(2)), -4/3*log(K(2))],[0, 4/3*log(K(2))], p, prec)
print("These are all of the found points (up to sign): ", points)