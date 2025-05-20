import sys
from sage.modular.dims import dimension_cusp_forms
from itertools import combinations


def get_sigma_k_1_at_0(k, chi):
    return (-1/(2*k)) * chi.bernoulli(k)


def get_sigma_k_1(k, chi1, chi2, n):
    D1 = chi1.conductor()
    D2 = chi2.conductor()
    if n == 0 and D2 != 1:
        return 0
    if n == 0 and D2 == 1:
        return get_sigma_k_1_at_0(k, chi1)
    ret = 0
    for d1 in divisors(n):
        d2 = n // d1
        ret += chi1(d1) * chi2(d2) * d1^(k-1)
    return ret





def get_a(K, l, chi, n):
    D = chi.level()
    k = K - l
    assert K % 2 == 0
    assert chi.is_primitive()
    assert chi(-1) == (-1)^l
    assert D % 2 == 1
    assert is_squarefree(D)
    chi_decomp = {}
    for chi_p in chi.decomposition():
        chi_decomp[chi_p.conductor()] = chi_p.extend(D)
    def get_chi_decomp(D_pm):
        ret_ = DirichletGroup(1)[0].extend(D)
        for p,_ in factor(D_pm):
            ret_ *= chi_decomp[p]
        return ret_.primitive_character()

    ret = 0
    for D1 in divisors(D):
        D2 = D // D1
        chi1 = get_chi_decomp(D1)
        chi2 = get_chi_decomp(D2)
        chi1_bar = chi1.bar()
        chi2_bar = chi2.bar()
        in_sum = 0
        for a1 in range(n*D2+1):
            a2 = n*D2 - a1
            in_sum += get_sigma_k_1(l,chi1,chi2_bar,a1) * get_sigma_k_1(k,chi1_bar,chi2,a2)
        ret += chi2_bar(-1) * in_sum
    sig0sig0 = get_sigma_k_1_at_0(l, chi) * get_sigma_k_1_at_0(k, chi.bar())
    ret -= 2 * sig0sig0 * sigma(n, K-1) / zeta(1 - K)
    return ret



def get_b(K, l, chi, n):
    D = chi.level()
    k = (K-2) - l
    assert K % 2 == 0
    assert chi.is_primitive()
    assert chi(-1) == (-1)^l
    assert D % 2 == 1
    assert is_squarefree(D)
    chi_decomp = {}
    for chi_p in chi.decomposition():
        chi_decomp[chi_p.conductor()] = chi_p.extend(D)
    def get_chi_decomp(D_pm):
        ret_ = DirichletGroup(1)[0].extend(D)
        for p,_ in factor(D_pm):
            ret_ *= chi_decomp[p]
        return ret_.primitive_character()
    # Now compute the summation
    ret = 0
    for D1 in divisors(D):
        D2 = D // D1
        chi1 = get_chi_decomp(D1)
        chi2 = get_chi_decomp(D2)
        chi1_bar = chi1.bar()
        chi2_bar = chi2.bar()
        in_sum = 0
        for a1 in range(n*D2+1):
            a2 = n*D2 - a1
            in_sum += get_sigma_k_1(l,chi1,chi2_bar,a1) * get_sigma_k_1(k,chi1_bar,chi2,a2) * (l*a2-k*a1)
        ret += (chi2_bar(-1) / D2) * in_sum
    return ret




def verify_form(coeffs, K):
    col = vector(coeffs)
    if K <= 10 or K == 14:
        return all(aa == 0 for aa in coeffs)
    prec = len(coeffs)
    SK_basis = CuspForms(1, K).basis()
    mat_tp = []
    for frm in SK_basis:
        frm_qexp = frm.qexp(prec+1)
        mat_tp.append( [QQ(frm_qexp[i]) for i in range(1,prec+1)] )
    mat = matrix(mat_tp).transpose() 
    try:
        mat.solve_right(col)
        return True
    except:
        print('ERROR:')
        print(col, 'not in span of')
        print(mat)
        return False



def compute_coeff_formulas(D, K, typ):
    if   typ == 'a':
        _0_4 = 0
        get_ab = get_a
    elif typ == 'b':
        _0_4 = 4
        get_ab = get_b

    prec = K//12 + 10
    prim_chars = [chi for chi in DirichletGroup(D) if chi.is_primitive()]

    for l in range(3, (K-_0_4)//2+1):
        for chi in prim_chars:
            if chi(-1) != (-1)^l: 
                continue
            coeffs = [get_ab(K,l,chi,n) for n in range(1,prec+1)]
            print(f'[{typ}] K={K}, l={l}, chi={chi}')
            print('coeffs =', coeffs)
            assert verify_form(coeffs, K)
    


def compute_coeffs(K_UB, D_UB, typ):
    for K in range(8, K_UB+1, 2):
        for D in range(1, D_UB+1, 2):
            if not is_squarefree(D):
                continue
            print(f'COEFFICIENTS FOR: K={K}, D={D} ##########################')
            compute_coeff_formulas(D, K, typ)





def compute_dets_l(K_UB, D_UB, typ):
    if   typ == 'a': 
        get_ab = get_a
        _0_4 = 0
    elif typ == 'b': 
        _0_4 = 4
        get_ab = get_b

    for D in range(1, D_UB+1):
        if D % 2 == 0 or not is_squarefree(D):
            continue
        for K in range(8, K_UB+1, 2):
            print(f'[{typ}] CHECKING: D={D}, K={K}')
            for chi in DirichletGroup(D):
                if not chi.is_primitive():
                    continue
                if chi(-1) == -1:
                    l_vals = list(range(3, (K-_0_4)//2+1, 2)) 
                else:
                    l_vals = list(range(4, (K-_0_4)//2+1, 2))
                dim = dimension_cusp_forms(Gamma0(1),K)
                n = min(len(l_vals), dim)
                if n < dim: 
                    print(f'WARNING: n={n} < dim={dim}; K={K}, l_vals={l_vals} ')
                for li in combinations(l_vals, n):
                    mat = []
                    for i in range(n):
                        l = li[i]
                        rw = [get_ab(K,l,chi,t) for t in range(1,n+1)]
                        mat.append(rw)
                    if matrix(mat).determinant() == 0:
                        print(f'[{typ}] FOUND!!!!: K={K}, li={li}')
                        print(chi)
                        print(mat.matrix())
                        exit()
                        



def compute_dets_chi(K_UB, D_UB, typ):
    if   typ == 'a': 
        get_ab = get_a
        _0_4 = 2
    elif typ == 'b': 
        _0_4 = 4
        get_ab = get_b

    for K in range(8, K_UB+1, 2):
        for D in range(1, D_UB+1):
            if D % 2 == 0 or not is_squarefree(D):
                continue
            print(f'[{typ}] CHECKING: K={K}, D={D}')
            for l in range(3, (K-_0_4)//2+1):
                chi_vals = []
                for chi in DirichletGroup(D):
                    if chi.is_primitive() and chi(-1) == (-1)^l:    
                        chi_vals.append(chi)
                dim = dimension_cusp_forms(Gamma0(1),K)
                n = min(len(chi_vals), dim)
                if n < dim: 
                    print(f'WARNING: n={n} < dim={dim}; K={K}, l={l}')
                for chii in combinations(chi_vals, n):
                    mat = []
                    for i in range(n):
                        chi = chii[i]
                        rw = [get_ab(K,l,chi,t) for t in range(1,n+1)]
                        mat.append(rw)
                    if matrix(mat).determinant() == 0:
                        print(f'[{typ}] FOUND!!!!: K={K}, l={l}')
                        [print(chi) for chi in chii]
                        print(mat.matrix())
                        exit()
                        








PARAM = sys.argv[1]
K_UB, D_UB = int(sys.argv[2]), int(sys.argv[3])
if PARAM in ['coeffs_a', 'coeffs_b']:
    compute_coeffs(K_UB, D_UB, typ=PARAM[-1])
    # RUN WITH: 'coeffs_a 16 15'
elif PARAM in ['dets_l_a', 'dets_l_b']: 
    compute_dets_l(K_UB, D_UB, typ=PARAM[-1])
    # RUN_WITH: 'dets_l_a 40 40'
elif PARAM in ['dets_chi_a', 'dets_chi_b']: 
    compute_dets_chi(K_UB, D_UB, typ=PARAM[-1])
    # RUN_WITH: 'dets_chi_a 40 40'
else: 
    assert False, 'invalid PARAM'






