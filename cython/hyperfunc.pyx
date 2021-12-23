# distutils: language = c++

cdef extern from "stdlib.h":
    ctypedef void const_void "const void"
    void qsort(void *base, int nmemb, int size,
            int(*compar)(const void *, const void *)) nogil

from libcpp.vector cimport vector
from libc.stdlib cimport malloc, free, rand, RAND_MAX
from libcpp cimport bool

import cython

import numpy as np
cimport numpy as np

np.import_array()

from cpython cimport array
import array

DTYPE = np.intc   

cpdef bool clipping_condition_A_numpy(np.ndarray[double, ndim=2] A, 
                                      np.ndarray[double, ndim=1] R):
    cdef int n = A.shape[0]
    cdef int m = A.shape[1]
    
    cdef double **A_c, *B_c
    
    A_c = convert_double_matrix_numpy_c(A)
    R_c = convert_double_vector_numpy_c(R)
    
    return(clipping_condition_A(n, m, A_c, R_c))

cpdef bool clipping_condition_B_numpy(np.ndarray[double, ndim=2] A, 
                                      np.ndarray[double, ndim=1] R):
    cdef int n = A.shape[0]
    cdef int m = A.shape[1]
    
    cdef double **A_c, *B_c
    
    A_c = convert_double_matrix_numpy_c(A)
    R_c = convert_double_vector_numpy_c(R)
    
    return(clipping_condition_B(n, m, A_c, R_c))

cdef bool clipping_condition_A(int n, 
                               int m, 
                               double **A, 
                               double  *R):
    cdef Py_ssize_t card_I, q_I, q_K, j_K, j_K_bar, r_J
    
    cdef int bc
    
    cdef int ***list_I, *I, **list_K, ***list_J_indices, *K, *K_bar, k, card_J, q_J, bc_J, *J, card_K
        
    cdef double * v 
    
    list_I = get_list_I(m - 1)
    
    for card_I in range( (m-1) + 1):
        bc = binomial_coefficient(m-1, card_I)
        for q_I in range(bc):
            I = list_I[card_I][q_I]
            # print('---s', 'card_I', card_I, 'q_I', q_I, 'bc=',bc)
            # print('I', [I[p] for p in range(card_I)])
            
            card_K = card_I - 1
            list_K = get_list_K(n, card_K)
            
            list_J_indices = combination_indices_full(n - card_K, n - card_K)
                
            for q_K in range(binomial_coefficient(n, card_K)):
                K = list_K[q_K]
                K_bar = bar(n, K, card_I-1)
                
                # print('K', [K[r_K] for r_K in range(card_K)])
                
                # print('K_bar', [K_bar[r_K_bar] for r_K_bar in range(n - card_K)])
                
                for card_J in range( n - card_K + 1):
                    for q_J in range(binomial_coefficient(n - card_K, card_J)):
                        J = <int *> malloc(card_J * sizeof(int *))
                        for r_J in range(card_J):
                            J[r_J] = K_bar[list_J_indices[card_J][q_J][r_J]]
                        # print('J :', [J[r_J] for r_J in range(card_J)])
                        
                        v = compute_vertex(n, A, R, I, card_I, J, card_J, K, card_K)
                        # print('v :', [v[i_v] for i_v in range(n)])w
                        
                        if test_vertex(n, m, A, R, I, card_I, v):
                            # one vertex is solution !
                            # according to the clipping condition (A) definition
                            # the condition is then unsatisfied
                            # False is therefore returned
                            return(False)
            # print('---e')
    return(True)

cdef bool clipping_condition_B(int n, 
                                int m, 
                                double **A, 
                                double *R):
    cdef Py_ssize_t card_I, q_I, q_K, j_K, j_K_bar, r_J
    
    cdef int bc
    
    cdef int ***list_I, *I, **list_K, ***list_J_indices, *K, *K_bar, k, card_J, q_J, bc_J, *J, card_K
        
    cdef double * v 
    cdef double random
    
    list_I = get_list_I(m - 1)
    
    for card_I in range( (m-1) + 1):
        bc = binomial_coefficient(m-1, card_I)
        for q_I in range(bc):
            I = list_I[card_I][q_I]
            # print('---s', 'card_I', card_I, 'q_I', q_I, 'bc=',bc)
            # print('I', [I[p] for p in range(card_I)])
            
            card_K = card_I
            list_K = get_list_K(n, card_K)
            
            list_J_indices = combination_indices_full(n - card_K, n - card_K)
                
            for q_K in range(binomial_coefficient(n, card_K)):
                K = list_K[q_K]
                K_bar = bar(n, K, card_I-1)
                
                # print('K', [K[r_K] for r_K in range(card_K)])
                
                # print('K_bar', [K_bar[r_K_bar] for r_K_bar in range(n - card_K)])
                
                for card_J in range( n - card_K + 1):
                    for q_J in range(binomial_coefficient(n - card_K, card_J)):
                        J = <int *> malloc(card_J * sizeof(int *))
                        for r_J in range(card_J):
                            J[r_J] = K_bar[list_J_indices[card_J][q_J][r_J]]
                        # print('J :', [J[r_J] for r_J in range(card_J)])
                        # print('===')
                        v = compute_vertex(n, A, R, I, card_I, J, card_J, K, card_K)
                        
                        
                        # print(random)
                        # print('I', [I[id_I] for id_I in range(card_I)])
                        # print('J', [J[r_J] for r_J in range(card_J)])
                        # print('K', [K[r_K] for r_K in range(card_K)])
                        # print('v', [v[i] for i in range(n)])
                        
                        if test_vertex(n, m, A, R, I, card_I, v):
                            print('accepted')
                            v_star = get_vertex_star(n, v)
                            print('v :', [v[i_v] for i_v in range(n)])
                            print('v_star :', [v_star[i_v] for i_v in range(1, v_star[0]+1)])
                            # v_01 = np.delete(np.arange(self.n), v_star-1) + 1
                            
                            # prod = 1
                            # for t in I:
                            #     I_union_m_remove_t = list(I.copy())
                            #     I_union_m_remove_t.append(self.m)
                            #     I_union_m_remove_t.remove(t)
                            #     I_union_m_remove_t = np.sort(I_union_m_remove_t)
                                
                            #     prod *= self._get_det_sub_A(v_star, I_union_m_remove_t)
                            
                            # for t in v_01:
                            #     I_union_m = list(I.copy())
                            #     I_union_m.append(self.m)
                            #     I_union_m = np.sort(I_union_m)
                                
                            #     v_star_union_t = list(v_star.copy())
                            #     v_star_union_t.append(t)
                            #     v_star_union_t = np.sort(v_star_union_t)
                                
                            #     prod *= self._get_det_sub_A(v_star_union_t, I_union_m)
            # print('---e')
    return(True)

cdef int * get_vertex_star(int n,
                           double *v):
    cdef Py_ssize_t i
    cdef int n_v_star = 0
    cdef int *v_star_too_long, *v_star
    
    v_star_too_long = <int *> malloc(n * sizeof(int *))
    
    for i in range(n):
        if v[i] != 0 and v[i] != 1:
            v_star_too_long[n_v_star] = i
            n_v_star = n_v_star + 1
    
    v_star = <int *> malloc((n_v_star + 1) * sizeof(int *))
    
    v_star[0] = n_v_star
    for i in range(1, n_v_star + 1):
        v_star[i] = v_star_too_long[i-1]
        
    return(v_star)

cdef double ** convert_double_matrix_numpy_c(np.ndarray[double, ndim=2] A):
    cdef Py_ssize_t i, j
    cdef int n = A.shape[0]
    cdef int m = A.shape[1]
    cdef double **A_c
    
    A_c = <double **> malloc(n * sizeof(double **))
        
    for i in range(n):
        A_c[i] = <double *> malloc(n * sizeof(double *))
        for j in range(m):
            A_c[i][j] = A[i][j]
    
    return(A_c)

cdef double * convert_double_vector_numpy_c(np.ndarray[double, ndim=1] V):
    cdef Py_ssize_t i
    cdef int n = V.shape[0]
    cdef double *V_c
    
    V_c = <double *> malloc(n * sizeof(double *))
        
    for i in range(n):
        V_c[i] = V[i]
    
    return(V_c)

cdef double * compute_vertex(int n, 
                             double **A, 
                             double *R, 
                             int *I,
                             int card_I,
                             int *J,
                             int card_J,
                             int *K,
                             int card_K):
    
    cdef Py_ssize_t i, id_I, id_J, id_K
    cdef double *v
    v = <double *> malloc(n * sizeof(double *))
    # generate system
    # A_sys is directly transposed
    # in order to correspond to the classical gaussian solving A_sys X = B_sys.
    # we should have card_I >= card_K
    # to prevent card_I > card_K which lead to unsolvable system
    # we limit to card_I equations to have (card_K, card_K) system.
    
    cdef double **A_sys, *B_sys
    
    A_sys = <double **> malloc(card_I * sizeof(double **))
    B_sys = <double *> malloc(card_I * sizeof(double *))
    
    for id_I in range(card_I):
        A_sys[id_I] = <double *> malloc(card_K * sizeof(double *))
        for id_K in range(card_K):
            A_sys[id_I][id_K] = A[K[id_K]][I[id_I]]
        
        B_sys[id_I] = - R[I[id_I]]
        # print('!', B_sys[id_I])
        for id_J in range(card_J):
            B_sys[id_I] = B_sys[id_I] - A[J[id_J]][I[id_I]]
        # print('!!', B_sys[id_I])
    # we have now A_sys a (card_K, card_K) matrix and B_sys a (card_K,) vector
    
    # solve
    v_K = gauss(card_K, A_sys, B_sys)
    
    for i in range(n):
        v[i] = 0
    
    for id_K in range(card_K):
        v[K[id_K]] = v_K[id_K]
        
    for id_J in range(card_J):
        v[J[id_J]] = 1    
    
    return(v)

cdef bool test_vertex(int n,
                      int m,
                      double **A,
                      double *R,
                      int *I,
                      int card_I,
                      double *v):
    cdef Py_ssize_t id_I, i
    cdef int *I_bar
    cdef double g
    
    
    I_bar = bar(m, I, card_I)
        
    # H = 
    for id_I in range(card_I):
        g = R[I[id_I]]
        for i in range(n):
            g = g + A[i][I[id_I]] * v[i]
        # print('g=', g)
        if abs(g) > 0.0000001: # 10**-7
            return(False)
        
    # H +
    for id_I in range(n - card_I - 1):
        # print('I_bar[id_I]', I_bar[id_I])
        g = R[I_bar[id_I]]
        for i in range(n):
            g = g + A[i][I_bar[id_I]] * v[i]
        
        if g < 0:
            return(False)
    
    # # # Hypercube
    for i in range(n):
        if v[i] < 0 or v[i] > 1:
            return(False)
    
    return(True)
    

@cython.boundscheck(False)  # Deactivate bounds checking.
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True) # Deactivate zero division checking.
cdef double * gauss(int n, 
                    double **A, 
                    double *B):
    # inspired by https://martin-thoma.com/solving-linear-equations-with-gaussian-elimination/
    cdef Py_ssize_t i, j, k
    cdef double **AB, *x, c
    cdef double maxEl, tmp
    cdef int maxRow
    
    AB = generate_AB(n, A, B)
    
    for i in range(n):
        # Search for maximum in this column
        maxEl = abs(AB[i][i])
        maxRow = i
        for k in range(i+1,n):
            if abs(AB[k][i]) > maxEl:
                maxEl = abs(AB[k][i])
                maxRow = k

        # Swap maximum row with current row (column by column)
        for k in range(i, n+1):
            tmp = AB[maxRow][k]
            AB[maxRow][k] = AB[i][k]
            AB[i][k] = tmp

        # Make all rows below this one 0 in current column
        for k in range(i+1, n):
            if AB[i][i] == 0.0:
                AB[i][i] = 0.00000001
                # print('1', i, AB[k][i], AB[i][i])
                
            c = -AB[k][i]/AB[i][i]
            for j in range(i, n+1):
                if i==j:
                    AB[k][j] = 0;
                else:
                    AB[k][j] += c * AB[i][j]

    # Solve equation Ax=b for an upper triangular matrix A
    x = <double *> malloc(n * sizeof(double *))
    for i in range(n):
        x[i] = 0.0
        
    for i in range(n - 1, -1, -1):
        if AB[i][i] == 0.0:
            # print('2', i, AB[i][n], AB[i][i])
            AB[i][i] = 0.00000001
            
            
        x[i] = AB[i][n]/AB[i][i]
        for k in range(i - 1, -1, -1):
            AB[k][n] -= AB[k][i] * x[i]
    return x

cdef double ** random_perturbation(int n, 
                                   double **A):
    cdef Py_ssize_t i, j
    cdef double **A_pert
    A_pert = <double **> malloc(n * sizeof(double **))
    
    for i in range(n):
        A_pert[i] = <double *> malloc(n * sizeof(double *))
        for j in range(n):
            A_pert[i][j] = A[i][j] + (rand() / (RAND_MAX + 1.0) * 2 - 1) * 0.000001
    
    return(A_pert)

cdef double ** generate_AB(int n,
                           double **A,
                           double *B):
    cdef Py_ssize_t i, j
    cdef double **AB
    AB = <double **> malloc(n * sizeof(double **))
    
    for i in range(n):
        AB[i] = <double *> malloc((n + 1) * sizeof(double **))
        for j in range(n):
            AB[i][j] = A[i][j]
    
    for i in range(n):
        AB[i][n] = B[i]
    
    return(AB)

cdef int * bar(int n, int *K, int card_K):
    cdef Py_ssize_t i
    cdef int * K_bar
    cdef int j_K, j_K_bar
    cdef bool trigger
    K_bar = <int *> malloc((n - card_K) * sizeof(int *))
    
    j_K = 0
    j_K_bar = 0
    for i in range(n):
        trigger = True
        for j_K in range(card_K):
            if K[j_K] == i:
                trigger = False
                break
        if trigger :
            K_bar[j_K_bar] = i
            j_K_bar = j_K_bar + 1
    
    return(K_bar)

cdef int *** get_list_I(int m):
    cdef int ***list_I
    
    list_I = combination_indices_full(m, m)
    
    return(list_I)

cdef int ** get_list_K(int n, int card):
    
    cdef int **list_K
    list_K = combinations_indices(n, card)
    
    return(list_K)

# cdef int ** get_list_J(n, list_K, card_K):
#     cdef int *** list_J
    
#     # cdef int *
    
#     for K in list_K:
#         pool = np.arange(self.n) + 1
        
#         if len(K) > 0:
#             pool = np.delete(pool, np.array(K) - 1)
        
#         list_J__K = []
        
#         for n_ones in range(self.n - len(K) + 1):
#             list_J__K += list(itertools.combinations(pool, n_ones))
        
#         list_J__K = [np.array(list(J__K)).astype(int) for J__K in list_J__K]
#         list_J.append(list_J__K)
            
#     return(list_J)

@cython.boundscheck(False)  # Deactivate bounds checking.
@cython.wraparound(False)   # Deactivate negative indexing.
cdef int *** combination_indices_full(int n, int k):
    cdef Py_ssize_t i, j, p
    cdef int ***ci_full, **ci, bc
    
    ci_full = <int ***> malloc((k + 1) * sizeof(int ***))
    
    for i in range(k+1):
        ci = combinations_indices(n, i)
        
        bc = binomial_coefficient(n, i)
        
        ci_full[i] = <int **> malloc(bc * sizeof(int **))
        
        for j in range(bc):
            ci_full[i][j] = <int *> malloc(i * sizeof(int *))
            
            for p in range(i):
                ci_full[i][j][p] = ci[j][p]
            
    
    return(ci_full)
            
@cython.boundscheck(False)  # Deactivate bounds checking.
@cython.wraparound(False)   # Deactivate negative indexing.
cdef int ** combinations_indices(int n, int k):
    cdef Py_ssize_t i, j
    cdef int **ci, *pool, *indices, id_ci, bc
    
    # compute the binomial coefficient
    # i.e. the number of combinations
    bc = binomial_coefficient(n, k)
    
    ci = <int **> malloc(bc * sizeof(int*))
    
    pool = list_range(n)
    
    indices = list_range(k)
    
    id_ci = 0
    ci[id_ci] = <int *> malloc(k * sizeof(int))
    for i in range(k):
        ci[id_ci][i] = indices[i]
    id_ci = id_ci + 1

    while True:
        for i in reversed(range(k)):
            if indices[i] != i + n - k:
                break
        else:
            return ci
        indices[i] += 1
        for j in range(i+1, k):
            indices[j] = indices[j-1] + 1
            
        ci[id_ci] = <int *> malloc(k * sizeof(int))
        for i in range(k):
            ci[id_ci][i] = indices[i]
        id_ci = id_ci + 1

cdef int factorial(int n):
    cdef int f = 1
    cdef Py_ssize_t i
    
    for i in range(1,n+1):
        f = f * i
    
    return(f)

@cython.cdivision(True)
cdef int binomial_coefficient(int n, int k):
    return(factorial(n) / (factorial(k) * factorial(n-k)))

cdef int * list_range(int n):
    cdef Py_ssize_t i
    cdef int *lst
    
    lst = <int *> malloc(n * sizeof(int))
    
    for i in range(n):
        lst[i] = i
    
    return(lst)