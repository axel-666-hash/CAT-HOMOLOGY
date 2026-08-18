"""
Computational Algebraic Topology
Blocks 1-6: simplicial complexes, boundary matrices over F_2,
Smith reduction, Betti numbers, elementary collapses, comparison pipeline.

Implements the reduction algorithm of Zomorodian & Carlsson (SoCG 2004, §2.5)
and elementary collapses, following Nanda's lecture notes.
"""

import itertools
from collections import defaultdict


# =============================================================
# BLOCK 1: simplicial complexes
# =============================================================

def get_all_faces(simplex):
    """
    All nonempty faces of a simplex, i.e. all nonempty subsets.
    A k-simplex has 2^(k+1) - 1 faces (Nanda, closure of a simplex).
    """
    faces = set()
    s_tuple = tuple(sorted(simplex))
    for r in range(1, len(s_tuple) + 1):
        for combo in itertools.combinations(s_tuple, r):
            faces.add(combo)
    return faces


class SimplicialComplex:
    """
    A simplicial complex stored as a set of sorted tuples.
    Downward closure (Definition 1.1, rule 2) is enforced at insertion,
    so an invalid complex cannot be constructed.
    """

    def __init__(self, simplices):
        self.simplices = set()
        for s in simplices:
            s_tuple = tuple(sorted(s))
            for face in get_all_faces(s_tuple):
                self.simplices.add(face)

    def k_simplices(self, k):
        """The set K_k of k-dimensional simplices, sorted for determinism."""
        return sorted(s for s in self.simplices if len(s) == k + 1)

    def dimension(self):
        """dim K = max over simplices of (#vertices - 1)."""
        if not self.simplices:
            return -1
        return max(len(s) - 1 for s in self.simplices)

    def num_simplices(self):
        return len(self.simplices)

    def f_vector(self):
        """(#K_0, #K_1, ..., #K_dim)."""
        d = self.dimension()
        return tuple(len(self.k_simplices(k)) for k in range(d + 1))

    def euler_characteristic(self):
        f = self.f_vector()
        return sum((-1) ** k * f[k] for k in range(len(f)))

    def copy(self):
        K = SimplicialComplex([])
        K.simplices = set(self.simplices)
        return K

    def remove_pair(self, sigma, tau):
        """Delete a free pair during an elementary collapse."""
        self.simplices.discard(tuple(sorted(sigma)))
        self.simplices.discard(tuple(sorted(tau)))

    def __repr__(self):
        return f"SimplicialComplex(dim={self.dimension()}, f={self.f_vector()})"


# =============================================================
# BLOCK 2: boundary matrices over F_2
# =============================================================

def boundary_matrix_f2(K, k):
    """
    The boundary matrix d_k : C_k -> C_{k-1} over F_2.

    Rows are (k-1)-simplices, columns are k-simplices. Over F_2 the signs
    (-1)^i vanish, so each column is stored as the SET of row indices where
    the entry is 1. Every column has exactly k+1 entries.

    Returns (matrix, rows, cols).
    """
    cols = K.k_simplices(k)

    if k == 0:
        return [], [], cols          # d_0 = 0 by convention

    rows = K.k_simplices(k - 1)
    rows_dict = {s: i for i, s in enumerate(rows)}

    matrix_cols = []
    for k_simplex in cols:
        col = set()
        for i in range(len(k_simplex)):
            face = k_simplex[:i] + k_simplex[i+1:]     # sigma_{-i}
            col.add(rows_dict[face])
        matrix_cols.append(col)

    return matrix_cols, rows, cols


# =============================================================
# BLOCK 3: Smith reduction over F_2, with basis tracking
# =============================================================

def smith_f2(matrix, n_rows, n_cols):
    """
    Reduction algorithm (ZC §2.5) over F_2.

    Columns are processed left to right. A surviving column claims a pivot
    row, taken to be the least index in its support; a later column containing
    a claimed row is reduced against the claiming column. Over F_2 the only
    elementary operation needed is C_i <- C_i + C_j, realised as symmetric
    difference of the index sets.

    Q is initialised to the identity and updated alongside every column
    operation, so Q[j] records which combination of ORIGINAL simplices the
    current column j represents. A column reduced to the empty set is a
    kernel vector, and Q[j] is then the corresponding cycle.

    Returns (rank, Q, non_pivot_columns).
    """
    A = [set(col) for col in matrix]
    Q = [{j} for j in range(n_cols)]
    pivot_col = {}                      # row -> column that claimed it
    rank = 0

    for j in range(n_cols):
        col = A[j]

        changed = True
        while changed:
            changed = False
            for r in list(col):
                if r in pivot_col:
                    pc = pivot_col[r]
                    A[j] = A[j].symmetric_difference(A[pc])
                    Q[j] = Q[j].symmetric_difference(Q[pc])
                    col = A[j]
                    changed = True
                    break

        if not col:
            continue                    # dead column -> kernel vector

        pivot_col[min(col)] = j
        rank += 1

    pivots = set(pivot_col.values())
    non_pivot = [j for j in range(n_cols) if j not in pivots]
    return rank, Q, non_pivot


# =============================================================
# BLOCK 4: Betti numbers
# =============================================================

def compute_homology(K, max_dim=None):
    """
    Betti numbers over F_2 via beta_k = #K_k - (r_k + r_{k+1})
    (Nanda, Proposition 3.13).

    Returns {k: beta_k}.
    """
    if max_dim is None:
        max_dim = K.dimension()

    ranks = {}
    n_simplices = {}

    for k in range(0, max_dim + 2):     # +2 so that r_{k+1} exists at the top
        mat, rows, cols = boundary_matrix_f2(K, k)
        n_simplices[k] = len(cols)

        if k == 0 or len(cols) == 0:
            ranks[k] = 0                # d_0 = 0, or no simplices in this dim
            continue

        ranks[k] = smith_f2(mat, len(rows), len(cols))[0]

    return {k: n_simplices[k] - ranks[k] - ranks[k + 1]
            for k in range(0, max_dim + 1)}


# =============================================================
# BLOCK 5: elementary collapses
# =============================================================

def find_free_face(K):
    """
    A free face is a codimension-one face contained in EXACTLY one simplex.
    Searches from the top dimension downward and returns (sigma, tau),
    or None if the complex is collapsed.
    """
    for d in range(K.dimension(), 0, -1):
        face_count = defaultdict(list)

        for tau in K.k_simplices(d):
            for i in range(len(tau)):
                face = tau[:i] + tau[i+1:]
                face_count[face].append(tau)

        for sigma, cofaces in face_count.items():
            if len(cofaces) == 1:
                return sigma, cofaces[0]

    return None


def collapse(K):
    """
    Apply elementary collapses until no free face remains.
    Each collapse is a deformation retraction, so homology is preserved.
    Returns (collapsed complex, number of collapses).
    """
    K_collapsed = K.copy()
    count = 0
    while True:
        pair = find_free_face(K_collapsed)
        if pair is None:
            break
        K_collapsed.remove_pair(*pair)
        count += 1
    return K_collapsed, count


def min_coface_count(K):
    """
    Smallest number of cofaces of any codimension-one face.
    Equals 1 exactly when a free face exists; used to diagnose why a
    complex admits no collapse.
    """
    best = None
    for d in range(K.dimension(), 0, -1):
        counts = defaultdict(int)
        for tau in K.k_simplices(d):
            for i in range(len(tau)):
                counts[tau[:i] + tau[i+1:]] += 1
        if counts:
            m = min(counts.values())
            best = m if best is None else min(best, m)
    return best


# =============================================================
# BLOCK 6: the comparison pipeline
# =============================================================

def pipeline(K, verbose=True):
    """
    Compute homology, collapse to exhaustion, recompute, compare.
    Betti numbers must agree, since collapses preserve homotopy type.
    """
    h_original = compute_homology(K)
    n_before = K.num_simplices()

    K_collapsed, n_collapses = collapse(K)
    n_after = K_collapsed.num_simplices()

    h_collapsed = compute_homology(K_collapsed)

    all_k = set(h_original) | set(h_collapsed)
    match = all(h_original.get(k, 0) == h_collapsed.get(k, 0) for k in all_k)

    if verbose:
        print(f"Original:  {n_before:6} simplices, beta = {h_original}")
        print(f"Collapsed: {n_after:6} simplices, beta = {h_collapsed}")
        print(f"Collapses: {n_collapses}")
        print(f"Homology preserved: {match}")

    return h_original, h_collapsed, match
