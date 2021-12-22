# distutils: language = c++

cdef extern from "stdlib.h":
    ctypedef void const_void "const void"
    void qsort(void *base, int nmemb, int size,
            int(*compar)(const void *, const void *)) nogil

from libcpp.vector cimport vector
from libc.stdlib cimport malloc, free
from libcpp cimport bool

import cython

import numpy as np
cimport numpy as np

np.import_array()

from cpython cimport array
import array

DTYPE = np.intc

cpdef int clipping_condition_A(int n, 
                               int m, 
                               np.ndarray[double, ndim=2] A, 
                               np.ndarray[double, ndim=1] R):
    cdef Py_ssize_t p_I, q_I, q_K, j_K, j_K_bar, r_J
    
    cdef int bc
    
    cdef int ***list_I, *I, **list_K, ***list_J_indices, *K, *K_bar, k, p_J, q_J, bc_J, *J, card_K
    
    cdef np.ndarray[double, ndim=1] v 
    
    list_I = get_list_I(m - 1)
    
    for p_I in range( (m-1) + 1):
        bc = binomial_coefficient(m-1, p_I)
        for q_I in range(bc):
            I = list_I[p_I][q_I]
            print('---s', 'p_I', p_I, 'q_I', q_I, 'bc=',bc)
            print('I', [I[p] for p in range(p_I)])
            
            card_K = p_I - 1
            list_K = get_list_K(n, card_K)
            
            # if n - (p_I-1) <= n:
            # print(n, n - (p_I - 1))
            list_J_indices = combination_indices_full(n - card_K, n - card_K)
            
            
                
            for q_K in range(binomial_coefficient(n, card_K)):
                K = list_K[q_K]
                K_bar = bar(n, K, p_I-1)
                
                print('K', [K[r_K] for r_K in range(card_K)])
                
                print('K_bar', [K_bar[r_K_bar] for r_K_bar in range(n - card_K)])
                
                for p_J in range( n - card_K + 1):
                    for q_J in range(binomial_coefficient(n - card_K, p_J)):
                        J = <int *> malloc(p_J * sizeof(int *))
                        for r_J in range(p_J):
                            J[r_J] = K_bar[list_J_indices[p_J][q_J][r_J]]
                        print('J :', [J[r_J] for r_J in range(p_J)])
                        
                        v = compute_vertex(n, m, A, R, I, p_I, J, p_J, K, card_K)
                        print('v :', v)
            print('---e')
    return(n)

cdef np.ndarray[double, ndim=1] compute_vertex(int n, 
                             int m, 
                             np.ndarray[double, ndim=2] A, 
                             np.ndarray[double, ndim=1] R, 
                             int *I,
                             int card_I,
                             int *J,
                             int card_J,
                             int *K,
                             int card_K):
    
    cdef Py_ssize_t i, j
    
    # generate system
    A_sys = np.zeros((card_K, card_I), dtype=np.double)
    cdef np.ndarray[double, ndim=2] A_sys_view = A_sys
    
    R_sys = np.zeros(card_I, dtype=np.double)
    cdef np.ndarray[double, ndim=1] R_sys_view = R_sys
    
    for id_I in range(card_I):
        
        for id_K in range(card_K):
            A_sys_view[id_K, id_I] = A[K[id_K], I[id_I]]
        
        R_sys_view[id_I] = R[I[id_I]]
        for id_J in range(card_J):
            R_sys_view[id_I] = R_sys_view[id_I] + A[J[id_J], I[id_I]]
    
    # solve
    # cdef np.ndarray[double, ndim=1] v = np.linalg.solve(A_sys, -R_sys)
    cdef np.ndarray[double, ndim=1] v = np.zeros(n, dtype=np.double)
    
    return(v)
    

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