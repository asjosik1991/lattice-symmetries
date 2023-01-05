#pragma once

#if defined(__cplusplus)
#include <complex>
#include <cstddef>
#include <cstdint>
#else
#include <stdatomic.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#endif

#if defined __has_include
#if __has_include("chpl-external-array.h")
#include "chpl-external-array.h"
#else
typedef struct {
  void *elts;
  uint64_t num_elts;

  void *freer;
} chpl_external_array;
#endif
#endif

#if defined(__cplusplus)
extern "C" {
#endif

#if defined(__cplusplus)
using ls_hs_scalar = std::complex<double>;
#else
#if defined(LS_NO_STD_COMPLEX)
typedef struct ls_hs_scalar {
  double real;
  double imag;
} ls_hs_scalar;
#else
typedef _Complex double ls_hs_scalar;
#endif
#endif

// {{{ Misc
void ls_hs_init(void);
void ls_hs_exit(void);
__attribute__((noreturn)) void ls_hs_fatal_error(char const *func, int line,
                                                 char const *message);

void ls_hs_set_exception_handler(void (*handler)(char const *message));
void ls_hs_error(char const *message);

#define LS_FATAL_ERROR(msg) ls_hs_fatal_error(__func__, __LINE__, msg)

#define LS_CHECK(cond, msg)                                                    \
  ((cond) ? ((void)0) : ls_hs_fatal_error(__func__, __LINE__, msg))

// void ls_hs_internal_set_free_stable_ptr(void (*f)(void *));

// void ls_hs_free_stable_ptr(void *);
// }}}

void ls_hs_destroy_external_array(chpl_external_array *arr);

// {{{ Symmetry
typedef struct ls_hs_symmetry ls_hs_symmetry;
ls_hs_symmetry *ls_hs_symmetry_from_json(const char *json_string);
void ls_hs_destroy_symmetry(ls_hs_symmetry *);

int ls_hs_symmetry_sector(ls_hs_symmetry const *);
int ls_hs_symmetry_length(ls_hs_symmetry const *);
int *ls_hs_symmetry_permutation(ls_hs_symmetry const *);
void ls_hs_destroy_permutation(int *);
// }}}

// {{{ Symmetries
typedef struct ls_hs_symmetries ls_hs_symmetries;
ls_hs_symmetries *ls_hs_symmetries_from_json(const char *json_string);
void ls_hs_destroy_symmetries(ls_hs_symmetries *);
// }}}

// {{{ Basis
typedef enum ls_hs_particle_type {
  LS_HS_SPIN,
  LS_HS_SPINFUL_FERMION,
  LS_HS_SPINLESS_FERMION
} ls_hs_particle_type;

typedef void (*ls_hs_internal_state_index_kernel_type)(
    ptrdiff_t batch_size, uint64_t const *alphas, ptrdiff_t alphas_stride,
    ptrdiff_t *indices, ptrdiff_t indices_stride, void const *private_data);

typedef void (*ls_hs_internal_is_representative_kernel_type)(
    ptrdiff_t batch_size, uint64_t const *alphas, ptrdiff_t alphas_stride,
    uint8_t *are_representatives, double *norms, void const *private_data);

typedef void (*ls_hs_internal_state_info_kernel_type)(
    ptrdiff_t batch_size, uint64_t const *alphas, ptrdiff_t alphas_stride,
    uint64_t *betas, ptrdiff_t betas_stride, ls_hs_scalar *characters,
    double *norms, void const *private_data);

typedef struct ls_hs_basis_kernels {
  ls_hs_internal_state_info_kernel_type state_info_kernel;
  void *state_info_data;
  ls_hs_internal_is_representative_kernel_type is_representative_kernel;
  void *is_representative_data;
  ls_hs_internal_state_index_kernel_type state_index_kernel;
  void *state_index_data;
} ls_hs_basis_kernels;

typedef struct ls_hs_permutation_group {
  _Atomic int refcount;
  int number_bits;
  int number_shifts;
  int number_masks;
  uint64_t *masks;
  uint64_t *shifts;
  double *eigvals_re;
  double *eigvals_im;
  void *haskell_payload;
} ls_hs_permutation_group;

typedef struct ls_hs_basis {
  _Atomic int refcount;
  int number_sites;
  int number_particles;
  int number_up;
  ls_hs_particle_type particle_type;
  bool state_index_is_identity;
  bool requires_projection;
  ls_hs_basis_kernels *kernels;
  chpl_external_array representatives;
  void *haskell_payload;
} ls_hs_basis;

ls_hs_basis *ls_hs_basis_from_json(char const *json_string);
ls_hs_basis *ls_hs_clone_basis(ls_hs_basis const *);
void ls_hs_destroy_basis(ls_hs_basis *);
uint64_t ls_hs_max_state_estimate(ls_hs_basis const *);
uint64_t ls_hs_min_state_estimate(ls_hs_basis const *);
int ls_hs_basis_number_bits(ls_hs_basis const *basis);
int ls_hs_basis_number_words(ls_hs_basis const *basis);
bool ls_hs_basis_has_fixed_hamming_weight(ls_hs_basis const *);

ptrdiff_t ls_hs_fixed_hamming_state_to_index(uint64_t basis_state);
uint64_t ls_hs_fixed_hamming_index_to_state(ptrdiff_t index,
                                            int hamming_weight);

char const *ls_hs_basis_to_json(ls_hs_basis const *);
void ls_hs_destroy_string(char const *);

void ls_hs_basis_build(ls_hs_basis *basis);
bool ls_hs_basis_is_built(ls_hs_basis const *basis);

char const *ls_hs_basis_state_to_string(ls_hs_basis const *,
                                        uint64_t const *state);

void ls_hs_state_index(ls_hs_basis const *basis, ptrdiff_t batch_size,
                       uint64_t const *spins, ptrdiff_t spins_stride,
                       ptrdiff_t *indices, ptrdiff_t indices_stride);

void ls_hs_is_representative(ls_hs_basis const *basis, ptrdiff_t batch_size,
                             uint64_t const *alphas, ptrdiff_t alphas_stride,
                             uint8_t *are_representatives, double *norms);

void ls_hs_state_info(ls_hs_basis const *basis, ptrdiff_t batch_size,
                      uint64_t const *alphas, ptrdiff_t alphas_stride,
                      uint64_t *betas, ptrdiff_t betas_stride,
                      ls_hs_scalar *characters, double *norms);

void ls_hs_build_representatives(ls_hs_basis *basis, uint64_t lower,
                                 uint64_t upper);

void ls_hs_unchecked_set_representatives(ls_hs_basis *basis,
                                         chpl_external_array const *states);
// }}}

// {{{ Expr

typedef struct ls_hs_expr ls_hs_expr;
// typedef struct ls_hs_expr {
//   _Atomic int refcount;
//   void *haskell_payload;
// } ls_hs_expr;

// ls_hs_expr *ls_hs_create_expr(char const *expression);
ls_hs_expr *ls_hs_expr_from_json(char const *json_string);
void ls_hs_destroy_expr(ls_hs_expr *expr);

char const *ls_hs_expr_to_json(ls_hs_expr const *expr);

typedef void (*ls_hs_index_replacement_type)(int spin, int site, int *new_spin,
                                             int *new_site);
ls_hs_expr *ls_hs_replace_indices(ls_hs_expr const *expr,
                                  ls_hs_index_replacement_type callback);

char const *ls_hs_expr_to_string(ls_hs_expr const *expr);

ls_hs_expr *ls_hs_expr_plus(ls_hs_expr const *a, ls_hs_expr const *b);
ls_hs_expr *ls_hs_expr_minus(ls_hs_expr const *a, ls_hs_expr const *b);
ls_hs_expr *ls_hs_expr_times(ls_hs_expr const *a, ls_hs_expr const *b);
ls_hs_expr *ls_hs_expr_scale(ls_hs_scalar const *z, ls_hs_expr const *a);

// }}}

typedef struct ls_internal_operator_kernel_data
    ls_internal_operator_kernel_data;

typedef struct ls_hs_nonbranching_terms {
  int number_terms;
  int number_bits;
  // number_words = ceil(number_bits / 64)
  ls_hs_scalar const *v; // array of shape [number_terms]
  uint64_t const *m;     // array of shape [number_terms, number_words]
  uint64_t const *l;     // array of shape [number_terms, number_words]
  uint64_t const *r;     // array of shape [number_terms, number_words]
  uint64_t const *x;     // array of shape [number_terms, number_words]
  uint64_t const *s;     // array of shape [number_terms, number_words]
  // all arrays are contiguous in row-major order
} ls_hs_nonbranching_terms;

typedef struct ls_hs_operator {
  _Atomic int refcount;
  ls_hs_basis const *basis;
  ls_hs_nonbranching_terms const *off_diag_terms;
  ls_hs_nonbranching_terms const *diag_terms;
  // ls_internal_operator_kernel_data const *apply_off_diag_cxt;
  // ls_internal_operator_kernel_data const *apply_diag_cxt;
  void *haskell_payload;
} ls_hs_operator;

ls_hs_operator *ls_hs_create_operator(ls_hs_basis const *basis,
                                      ls_hs_expr const *expr);
ls_hs_operator *ls_hs_clone_operator(ls_hs_operator const *);
void ls_hs_destroy_operator(ls_hs_operator *);
int ls_hs_operator_max_number_off_diag(ls_hs_operator const *);
ls_hs_expr const *ls_hs_operator_get_expr(ls_hs_operator const *);
ls_hs_basis const *ls_hs_operator_get_basis(ls_hs_operator const *);

// ls_hs_operator *ls_hs_create_operator(ls_hs_basis const *basis,
//                                       char const *expression, int
//                                       number_tuples, int tuple_size, int
//                                       const *indices);
// char const *ls_hs_operator_to_json(ls_hs_operator const *, bool
// include_basis); ls_hs_operator *ls_hs_operator_from_json(char const
// *json_string,
//                                          ls_hs_basis const *basis);

// ls_hs_operator *ls_hs_operator_plus(ls_hs_operator const *,
//                                     ls_hs_operator const *);
// ls_hs_operator *ls_hs_operator_minus(ls_hs_operator const *,
//                                      ls_hs_operator const *);
// ls_hs_operator *ls_hs_operator_times(ls_hs_operator const *,
//                                      ls_hs_operator const *);
// ls_hs_operator *ls_hs_operator_scale(ls_hs_scalar const *,
//                                      ls_hs_operator const *);
// ls_hs_operator *ls_hs_operator_hermitian_conjugate(ls_hs_operator const *);

// bool ls_hs_operator_is_hermitian(ls_hs_operator const *);
// bool ls_hs_operator_is_identity(ls_hs_operator const *);
// bool ls_hs_operator_is_real(ls_hs_operator const *);

// char const *ls_hs_operator_pretty_terms(ls_hs_operator const *);

// ls_hs_operator *ls_hs_load_hamiltonian_from_yaml(char const *);

typedef struct ls_hs_yaml_config {
  ls_hs_basis const *basis;
  ls_hs_operator const *hamiltonian;
  int number_observables;
  ls_hs_operator const *const *observables;
} ls_hs_yaml_config;

ls_hs_yaml_config *ls_hs_load_yaml_config(char const *);
void ls_hs_destroy_yaml_config(ls_hs_yaml_config *);

// ls_internal_operator_kernel_data *ls_internal_create_apply_diag_kernel_data(
//     ls_hs_nonbranching_terms const *diag_terms);
// ls_internal_operator_kernel_data
// *ls_internal_create_apply_off_diag_kernel_data(
//     ls_hs_nonbranching_terms const *off_diag_terms);
// void ls_internal_destroy_operator_kernel_data(
//     ls_internal_operator_kernel_data *p);

// void ls_internal_operator_apply_off_diag(
//     ptrdiff_t batch_size, uint64_t const *alphas, ptrdiff_t alphas_stride,
//     uint64_t *betas, ptrdiff_t betas_stride, ls_hs_scalar *coeffs,
//     ls_internal_operator_kernel_data const *cxt);
// void ls_internal_operator_apply_diag(
//     ptrdiff_t batch_size, uint64_t const *alphas, ptrdiff_t alphas_stride,
//     ls_hs_scalar *coeffs, ls_internal_operator_kernel_data const *cxt);

// void ls_hs_operator_apply_diag_kernel(ls_hs_operator const *op,
//                                       ptrdiff_t batch_size,
//                                       uint64_t const *alphas,
//                                       ptrdiff_t alphas_stride,
//                                       ls_hs_scalar *coeffs);

// void ls_hs_operator_apply_off_diag_kernel(
//     ls_hs_operator const *op, ptrdiff_t batch_size, uint64_t const *alphas,
//     ptrdiff_t alphas_stride, uint64_t *betas, ptrdiff_t betas_stride,
//     ls_hs_scalar *coeffs);

void ls_internal_operator_apply_off_diag_x1(
    ls_hs_operator const *op, ptrdiff_t batch_size, uint64_t const *alphas,
    uint64_t *betas, ls_hs_scalar *coeffs, ptrdiff_t *offsets,
    double const *xs);
void ls_internal_operator_apply_diag_x1(ls_hs_operator const *op,
                                        ptrdiff_t batch_size,
                                        uint64_t const *alphas, double *ys,
                                        double const *xs);

// {{{ Binomials
// typedef struct {
//   int dimension;
//   bool is_per_sector;
//   uint64_t *coefficients; // array of shape [dimension, dimension]
//                           // stored in row-major order
//   // coefficients[n, k] corresponds to binomial(n, k)
// } ls_hs_combinadics_kernel_data;
//
// ls_hs_combinadics_kernel_data *
// ls_hs_internal_create_combinadics_kernel_data(int number_bits,
//                                               bool is_per_sector);
// void ls_hs_internal_destroy_combinadics_kernel_data(
//     ls_hs_combinadics_kernel_data *p);

// void ls_hs_state_index_combinadics_kernel(ptrdiff_t batch_size,
//                                           uint64_t const *spins,
//                                           ptrdiff_t spins_stride,
//                                           ptrdiff_t *indices,
//                                           ptrdiff_t indices_stride,
//                                           void const *private_kernel_data);
// }}}

// void ls_hs_state_index_identity_kernel(ptrdiff_t batch_size,
//                                        uint64_t const *spins,
//                                        ptrdiff_t spins_stride,
//                                        ptrdiff_t *indices,
//                                        ptrdiff_t indices_stride,
//                                        void const *private_kernel_data);
//
// {{{ Indexing

// typedef struct ls_hs_internal_lookup_cache {
//   int number_prefix_bits;
//   int bit_shift;
//   ptrdiff_t number_representatives;
//   ptrdiff_t number_ranges;
//   uint64_t *representatives;
//   ptrdiff_t *ranges;
// } ls_hs_internal_lookup_cache;

void ls_hs_state_index_binary_search_kernel(ptrdiff_t batch_size,
                                            uint64_t const *spins,
                                            ptrdiff_t spins_stride,
                                            ptrdiff_t *indices,
                                            ptrdiff_t indices_stride,
                                            void const *private_kernel_data);

// }}}

typedef struct ls_internal_halide_kernel_data ls_internal_halide_kernel_data;

ls_internal_halide_kernel_data *
ls_internal_create_halide_kernel_data(ls_hs_permutation_group const *g,
                                      int spin_inversion);

void ls_internal_destroy_halide_kernel_data(ls_internal_halide_kernel_data *p);

void ls_hs_is_representative_halide_kernel(
    ptrdiff_t batch_size, uint64_t const *alphas, ptrdiff_t alphas_stride,
    uint8_t *are_representatives, double *norms, void const *private_data);

void ls_hs_state_info_halide_kernel(ptrdiff_t batch_size,
                                    uint64_t const *alphas,
                                    ptrdiff_t alphas_stride, uint64_t *betas,
                                    ptrdiff_t betas_stride,
                                    ls_hs_scalar *characters, double *norms,
                                    void const *private_data);

// void ls_hs_evaluate_wavefunction_via_statevector(
//     ls_hs_basis const *basis, ptrdiff_t batch_size, uint64_t const *alphas,
//     ptrdiff_t alphas_stride, void const *state_vector, size_t element_size,
//     void *coeffs, );

typedef struct ls_hs_state_index_binary_search_data
    ls_hs_state_index_binary_search_data;

ls_hs_state_index_binary_search_data *
ls_hs_create_state_index_binary_search_kernel_data(
    chpl_external_array const *representatives, int number_bits);

void ls_hs_destroy_state_index_binary_search_kernel_data(
    ls_hs_state_index_binary_search_data *cache);

void ls_hs_state_index_binary_search_kernel(ptrdiff_t batch_size,
                                            uint64_t const *spins,
                                            ptrdiff_t spins_stride,
                                            ptrdiff_t *indices,
                                            ptrdiff_t indices_stride,
                                            void const *private_kernel_data);

// void ls_hs_internal_block_binary_search(ptrdiff_t const block_size,
//                                         uint64_t const haystack[],
//                                         ptrdiff_t const haystack_size,
//                                         uint64_t const needles[],
//                                         ptrdiff_t indices[]);

// void ls_hs_internal_binary_search_x1(uint64_t const haystack[],
//                                      ptrdiff_t haystack_size,
//                                      uint64_t needle,
//                                      ptrdiff_t *index);

typedef struct ls_chpl_kernels {
  void (*enumerate_states)(ls_hs_basis const *, uint64_t, uint64_t,
                           chpl_external_array *);
  void (*operator_apply_off_diag)(ls_hs_operator *, int64_t, uint64_t *,
                                  chpl_external_array *, chpl_external_array *,
                                  chpl_external_array *, int64_t);
  void (*operator_apply_diag)(ls_hs_operator *, int64_t, uint64_t *,
                              chpl_external_array *, int64_t);
  void (*matrix_vector_product)(ls_hs_operator *, int, double const *,
                                double *);
} ls_chpl_kernels;

ls_chpl_kernels const *ls_hs_internal_get_chpl_kernels(void);

void ls_chpl_init(void);
void ls_chpl_finalize(void);

void ls_hs_internal_set_chpl_kernels(ls_chpl_kernels const *kernels);

#ifdef __cplusplus
} // extern "C"
#endif
