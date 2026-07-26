# ============================================================
#  SO(5) reduced Clebsch--Gordon coefficients
#  C2 + C4 composite diagonalization version, v5 high-performance version
#
#  Main improvement:
#
#      C2_SO5 alone cannot distinguish some SO(5) irreps.
#      In particular, in the code convention,
#
#          C2(6,0) = C2(5,3) = 30.
#
#      Therefore this version constructs the independent quartic
#      Casimir C4_SO5 from
#
#          W_a = (1/16) epsilon_{abcde} {L_bc, L_de},
#          C4  = sum_a W_a W_a,
#
#      and selects target states by residuals of C2, C4, J^2, K^2.
#
#  Compatibility:
#
#      The main function is still
#
#          compute_reduced_so5_cg_composite(...)
#
#      and the returned NamedTuple still contains:
#
#          reduced_keyed, reduced_pretty, selected_vectors,
#          diagnostics, H, HL, HR, O, basis, product_states.
#
#      Extra fields are added:
#
#          H4, epsC4, target_C4, c4_info.
#
#      Your previous Colab interface can call this code unchanged.
#
# ============================================================
#
#  v5 高性能版说明（数学定义与主接口保持不变）
#
#  关键优化：
#
#  1. 不再构造完整两粒子空间上的 C4_full。
#     对固定磁扇区 S，直接计算
#
#         P_S C4 P_S = sum_a (W_a P_S)' * (W_a P_S),
#
#     只构造 W_a 在目标扇区对应的列。该等式来自 W_a 的 Hermiticity，
#     与原始 C4 = sum_a W_a W_a 完全相同，但避免最昂贵的完整 W_a^2。
#
#  2. 两粒子总生成元只保存 a<b 的 10 个独立 L_ab，避免重复生成下三角。
#
#  3. 缓存组合算符 O 的完整本征分解。同一 (p,M1,M2,eps...) 扫描多个
#     (P,Q,J,K) 时只对角化一次。
#
#  4. 交换子诊断可缓存，且主函数默认关闭；需要验证时显式打开。
#
#  5. 目标选择先用 O 的目标本征值预筛选，再计算四种严格残差；
#     每个矩阵只做一次 matrix-vector product，同时得到期望值和残差。
#
#  6. 共享 SU(2) CG cache；约化系数提取时先预计算所有 SO(4) CG。
#
#  注意：首次 include 会有 Julia JIT 编译开销。性能比较应先预热一次。
#
# ============================================================

using LinearAlgebra
using SparseArrays
using Statistics

# ------------------------------------------------------------
# 0. Basic type aliases and utilities
# ------------------------------------------------------------

const SO4State = NTuple{4,Int}       # (tj,tm,tk,tn) = 2*(j,m,k,n)
const ProductState = Tuple{Int,Int}  # pair of single-particle basis indices

# ------------------------------------------------------------
# Global caches
# ------------------------------------------------------------
# 这些 cache 的内容都只依赖明确写入 key 的参数。
# 修改生成元 convention 后必须执行 clear_SO5_caches!()。

const __SO5_BASIS_CACHE__         = Dict{Int,Any}()
const __SO5_SINGLE_L_CACHE__      = Dict{Int,Any}()
const __SO5_SECTOR_MATRIX_CACHE__ = Dict{Any,Any}()

# 兼容旧接口的完整 5x5 total-L cache（主流程不再使用）。
const __SO5_TOTAL_L_CACHE__       = Dict{Int,Any}()

# 主流程使用：只保存 10 个 a<b 的 total L_ab。
const __SO5_TOTAL_L_UPPER_CACHE__ = Dict{Int,Any}()
const __SO5_C4_CALIB_CACHE__      = Dict{Any,Any}()

# 同一扇区扫描多个目标时，避免重复全谱对角化及交换子计算。
const __SO5_COMPOSITE_EIG_CACHE__ = Dict{Any,Any}()
const __SO5_COMMUTATOR_CACHE__    = Dict{Any,Any}()

# CGCache 类型在后文定义，因此这里用 Ref{Any} 保存共享对象。
const __SO5_SHARED_CG_CACHE__ = Ref{Any}(nothing)

function clear_SO5_caches!()
    empty!(__SO5_BASIS_CACHE__)
    empty!(__SO5_SINGLE_L_CACHE__)
    empty!(__SO5_SECTOR_MATRIX_CACHE__)
    empty!(__SO5_TOTAL_L_CACHE__)
    empty!(__SO5_TOTAL_L_UPPER_CACHE__)
    empty!(__SO5_C4_CALIB_CACHE__)
    empty!(__SO5_COMPOSITE_EIG_CACHE__)
    empty!(__SO5_COMMUTATOR_CACHE__)
    __SO5_SHARED_CG_CACHE__[] = nothing
    println("SO(5) caches cleared.")
    return nothing
end

function SO5_cache_status()
    return (
        basis_cache_entries          = length(__SO5_BASIS_CACHE__),
        single_L_cache_entries       = length(__SO5_SINGLE_L_CACHE__),
        sector_matrix_cache_entries  = length(__SO5_SECTOR_MATRIX_CACHE__),
        legacy_total_L_entries       = length(__SO5_TOTAL_L_CACHE__),
        upper_total_L_entries        = length(__SO5_TOTAL_L_UPPER_CACHE__),
        c4_calibration_cache_entries = length(__SO5_C4_CALIB_CACHE__),
        composite_eigen_cache_entries = length(__SO5_COMPOSITE_EIG_CACHE__),
        commutator_cache_entries      = length(__SO5_COMMUTATOR_CACHE__),
        shared_cg_entries = (__SO5_SHARED_CG_CACHE__[] === nothing ? 0 :
                             length(__SO5_SHARED_CG_CACHE__[].table))
    )
end

qval(t::Int) = 0.5 * t

function to2(x; tol::Float64=1e-10)
    y = 2.0 * Float64(x)
    r = round(Int, y)

    if abs(y - r) > tol
        error("Quantum number $x is not an integer or half-integer.")
    end

    return r
end

function half_string(t::Int)
    if iseven(t)
        return string(div(t,2))
    else
        return string(t, "/2")
    end
end

function state_string(s::SO4State)
    tj,tm,tk,tn = s

    return "(j=$(half_string(tj)), m=$(half_string(tm)); " *
           "k=$(half_string(tk)), n=$(half_string(tn)))"
end

# ------------------------------------------------------------
# 1. SO(5) representation data
# ------------------------------------------------------------

"""
    so5_casimir(P,Q)

Quadratic Casimir of SO(5) irrep (P,Q), in the convention of this code:

    C2(P,Q) = (P^2 + Q^2)/2 + 2P + Q.
"""
so5_casimir(P::Integer, Q::Integer) =
    (Float64(P)^2 + Float64(Q)^2)/2 + 2Float64(P) + Float64(Q)

"""
    so5_casimir4(P,Q)

Quartic Casimir eigenvalue in the same (P,Q) convention.

This corresponds to

    C4 = sum_a W_a W_a,
    W_a = (1/16) epsilon_{abcde} {L_bc,L_de},

up to the normalization fixed by single-particle calibration.

Formula:

    C4(P,Q) =
        (P-Q)(P-Q+2)(P+Q+2)(P+Q+4) / 16.

It distinguishes, for example,

    C2(6,0) = C2(5,3) = 30,
    C4(6,0) = 240,
    C4(5,3) = 60.
"""
so5_casimir4(P::Integer, Q::Integer) =
    Float64((P-Q)*(P-Q+2)*(P+Q+2)*(P+Q+4)) / 16.0

"""
    so5_dim(P,Q)

Dimension of SO(5) irrep (P,Q):

    dim(P,Q) = (1+Q)(1+P-Q)(P+2)(P+Q+3)/6.
"""
function so5_dim(P::Integer, Q::Integer)
    return Int((1+Q)*(1+P-Q)*(P+2)*(P+Q+3) ÷ 6)
end

"""
    so5_p0_dimension(p)

Dimension of the SO(5) irrep (p,0):

    dim(p,0) = (p+1)(p+2)(p+3)/6.
"""
function so5_p0_dimension(p::Integer)
    return Int((p+1)*(p+2)*(p+3) ÷ 6)
end

# ------------------------------------------------------------
# 2. Log-factorial cache and SU(2) CG coefficients
# ------------------------------------------------------------

mutable struct CGCache
    table::Dict{NTuple{6,Int},Float64}
    logfact::Vector{Float64}
end

function logfactorials(nmax::Int)
    lf = zeros(Float64, nmax + 1)

    for n in 1:nmax
        lf[n+1] = lf[n] + log(Float64(n))
    end

    return lf
end

function CGCache(nmax::Int=128)
    return CGCache(Dict{NTuple{6,Int},Float64}(), logfactorials(nmax))
end


"""
    get_shared_CGCache(nmax)

返回全局共享的 SU(2) CG cache。若需要更大的 log-factorial 表，只原地扩展，
从而避免扫描多个目标 channel 时重复计算相同的 SU(2) CG 系数。
"""
function get_shared_CGCache(nmax::Int=128)
    obj = __SO5_SHARED_CG_CACHE__[]

    if obj === nothing
        cache = CGCache(nmax)
        __SO5_SHARED_CG_CACHE__[] = cache
        return cache
    end

    cache = obj::CGCache
    ensure_logfact!(cache, nmax)
    return cache
end

function ensure_logfact!(cache::CGCache, nmax::Int)
    oldmax = length(cache.logfact) - 1

    if nmax <= oldmax
        return nothing
    end

    resize!(cache.logfact, nmax + 1)

    for n in (oldmax+1):nmax
        cache.logfact[n+1] = cache.logfact[n] + log(Float64(n))
    end

    return nothing
end

@inline function lfact(cache::CGCache, n::Int)
    if n < 0
        return Inf
    end

    ensure_logfact!(cache, n)
    return cache.logfact[n+1]
end

@inline function can_couple2(tj1::Int, tj2::Int, tJ::Int)
    return abs(tj1 - tj2) <= tJ <= tj1 + tj2 &&
           iseven(tj1 + tj2 - tJ)
end

"""
    su2_cg2(cache,tj1,tm1,tj2,tm2,tJ,tM)

Compute ordinary SU(2) Clebsch--Gordon coefficient

    <j1 m1, j2 m2 | J M>.

All quantum numbers are doubled integers.
"""
function su2_cg2(cache::CGCache,
                 tj1::Int, tm1::Int,
                 tj2::Int, tm2::Int,
                 tJ::Int,  tM::Int)

    key = (tj1,tm1,tj2,tm2,tJ,tM)

    if haskey(cache.table, key)
        return cache.table[key]
    end

    if tm1 + tm2 != tM
        cache.table[key] = 0.0
        return 0.0
    end

    if abs(tm1) > tj1 || abs(tm2) > tj2 || abs(tM) > tJ
        cache.table[key] = 0.0
        return 0.0
    end

    if !can_couple2(tj1,tj2,tJ)
        cache.table[key] = 0.0
        return 0.0
    end

    a1   = div(tj1 + tj2 - tJ, 2)
    a2   = div(tj1 - tj2 + tJ, 2)
    a3   = div(-tj1 + tj2 + tJ, 2)
    aden = div(tj1 + tj2 + tJ, 2) + 1

    c1 = div(tj1 + tm1, 2)
    c2 = div(tj1 - tm1, 2)
    c3 = div(tj2 + tm2, 2)
    c4 = div(tj2 - tm2, 2)
    c5 = div(tJ + tM, 2)
    c6 = div(tJ - tM, 2)

    if minimum((a1,a2,a3,aden,c1,c2,c3,c4,c5,c6)) < 0
        cache.table[key] = 0.0
        return 0.0
    end

    maxarg = maximum((a1,a2,a3,aden,c1,c2,c3,c4,c5,c6))
    ensure_logfact!(cache, maxarg)

    log_pref =
        0.5 * (
            log(Float64(tJ + 1)) +
            lfact(cache,a1) + lfact(cache,a2) + lfact(cache,a3) -
            lfact(cache,aden) +
            lfact(cache,c1) + lfact(cache,c2) +
            lfact(cache,c3) + lfact(cache,c4) +
            lfact(cache,c5) + lfact(cache,c6)
        )

    zmin = max(
        0,
        cld(tj2 - tJ - tm1, 2),
        cld(tj1 - tJ + tm2, 2)
    )

    zmax = min(
        div(tj1 + tj2 - tJ, 2),
        div(tj1 - tm1, 2),
        div(tj2 + tm2, 2)
    )

    if zmin > zmax
        cache.table[key] = 0.0
        return 0.0
    end

    b2 = div(tj1 + tj2 - tJ, 2)
    b3 = div(tj1 - tm1, 2)
    b4 = div(tj2 + tm2, 2)
    b5 = div(tJ - tj2 + tm1, 2)
    b6 = div(tJ - tj1 - tm2, 2)

    s = 0.0

    for z in zmin:zmax
        d1 = z
        d2 = b2 - z
        d3 = b3 - z
        d4 = b4 - z
        d5 = b5 + z
        d6 = b6 + z

        if minimum((d1,d2,d3,d4,d5,d6)) < 0
            continue
        end

        log_den =
            lfact(cache,d1) + lfact(cache,d2) + lfact(cache,d3) +
            lfact(cache,d4) + lfact(cache,d5) + lfact(cache,d6)

        s += (isodd(z) ? -1.0 : 1.0) * exp(-log_den)
    end

    val = exp(log_pref) * s
    cache.table[key] = val

    return val
end

# ------------------------------------------------------------
# 3. SO(4)-adapted basis for (p,0)
# ------------------------------------------------------------

"""
    generate_so4_basis_p0(p)

Generate SO(4)-adapted basis of the SO(5) irrep (p,0).

State convention:

    state = (tj,tm,tk,tn) = 2*(j,m,k,n),

where

    tk = p - tj.
"""
function generate_so4_basis_p0(p::Int)
    p < 0 && error("p must be nonnegative.")

    if haskey(__SO5_BASIS_CACHE__, p)
        cached = __SO5_BASIS_CACHE__[p]
        return cached.basis, cached.label
    end

    basis = SO4State[]
    label = Dict{SO4State,Int}()

    for tj in 0:p
        tk = p - tj

        for tm in -tj:2:tj
            for tn in -tk:2:tk
                s = (tj,tm,tk,tn)
                push!(basis, s)
                label[s] = length(basis)
            end
        end
    end

    expected = so5_p0_dimension(p)

    if length(basis) != expected
        error("Generated basis dimension $(length(basis)) != expected $expected.")
    end

    __SO5_BASIS_CACHE__[p] = (basis = basis, label = label)

    return basis, label
end

# ------------------------------------------------------------
# 4. Product sector at fixed M1,M2
# ------------------------------------------------------------

"""
    generate_product_sector_fixed_M(p,tM1,tM2)

Generate all ordered two-particle product states in

    (p,0) ⊗ (p,0)

with fixed magnetic quantum numbers

    m1 + m2 = M1,
    n1 + n2 = M2.

This function does not restrict by J,K.
"""
function generate_product_sector_fixed_M(p::Int, tM1::Int, tM2::Int)
    basis, label = generate_so4_basis_p0(p)

    by_mag = Dict{Tuple{Int,Int},Vector{Int}}()

    for i in eachindex(basis)
        tj,tm,tk,tn = basis[i]
        key = (tm,tn)

        if !haskey(by_mag,key)
            by_mag[key] = Int[]
        end

        push!(by_mag[key], i)
    end

    product_states = ProductState[]

    for i in eachindex(basis)
        tj1,tm1,tk1,tn1 = basis[i]
        needed = (tM1 - tm1, tM2 - tn1)

        for j in get(by_mag, needed, Int[])
            push!(product_states, (i,j))
        end
    end

    return product_states, basis, label
end

function build_product_lookup(product_states::Vector{ProductState})
    lookup = Dict{ProductState,Int}()

    for i in eachindex(product_states)
        lookup[product_states[i]] = i
    end

    return lookup
end

# ------------------------------------------------------------
# 5. SU(2) ladder coefficients
# ------------------------------------------------------------

@inline function su2_casimir_from_tj(tj::Int)
    j = qval(tj)
    return j * (j + 1.0)
end

@inline function jplus_coeff2(tj::Int, tm::Int)
    if tm >= tj
        return 0.0
    end

    j = qval(tj)
    m = qval(tm)

    return sqrt(max(0.0, j*(j+1.0) - m*(m+1.0)))
end

@inline function jminus_coeff2(tj::Int, tm::Int)
    if tm <= -tj
        return 0.0
    end

    j = qval(tj)
    m = qval(tm)

    return sqrt(max(0.0, j*(j+1.0) - m*(m-1.0)))
end

@inline function sqrt_nonneg(x::Real)
    return sqrt(max(0.0, Float64(x)))
end

@inline function sqrtprod(a::Real, b::Real)
    return sqrt_nonneg(Float64(a) * Float64(b))
end

# ------------------------------------------------------------
# 6. Matrix insertion helpers
# ------------------------------------------------------------

@inline function add_pair!(A::Matrix{Float64},
                           row::Int,
                           label::Dict{SO4State,Int},
                           lookup::Dict{ProductState,Int},
                           snew1::SO4State,
                           snew2::SO4State,
                           coeff::Float64)

    coeff == 0.0 && return nothing

    i1 = get(label, snew1, 0)
    i1 == 0 && return nothing

    i2 = get(label, snew2, 0)
    i2 == 0 && return nothing

    col = get(lookup, (i1,i2), 0)
    col == 0 && return nothing

    A[row,col] += coeff

    return nothing
end

@inline function add_single!(A::SparseMatrixCSC{ComplexF64,Int},
                             col::Int,
                             label::Dict{SO4State,Int},
                             snew::SO4State,
                             coeff)

    coeff == 0 && return nothing

    row = get(label, snew, 0)
    row == 0 && return nothing

    A[row,col] += ComplexF64(coeff)

    return nothing
end

# ------------------------------------------------------------
# 7. SO(5)/SO(4) T-generator cross terms for C2
# ------------------------------------------------------------

function add_T_terms!(H::Matrix{Float64},
                      row::Int,
                      s1::SO4State,
                      s2::SO4State,
                      label::Dict{SO4State,Int},
                      lookup::Dict{ProductState,Int})

    tj1,tm1,tk1,tn1 = s1
    tj2,tm2,tk2,tn2 = s2

    j1 = qval(tj1); m1 = qval(tm1); k1 = qval(tk1); n1 = qval(tn1)
    j2 = qval(tj2); m2 = qval(tm2); k2 = qval(tk2); n2 = qval(tn2)

    # T_{++}(1) T_{--}(2)
    add_pair!(H,row,label,lookup,
        (tj1+1,tm1+1,tk1-1,tn1+1),
        (tj2-1,tm2-1,tk2+1,tn2-1),
        sqrtprod(j1+m1+1, k1-n1) *
        sqrtprod(j2+m2, k2-n2+1))

    add_pair!(H,row,label,lookup,
        (tj1+1,tm1+1,tk1-1,tn1+1),
        (tj2+1,tm2-1,tk2-1,tn2-1),
        sqrtprod(j1+m1+1, k1-n1) *
        sqrtprod(k2+n2, j2-m2+1))

    add_pair!(H,row,label,lookup,
        (tj1-1,tm1+1,tk1+1,tn1+1),
        (tj2-1,tm2-1,tk2+1,tn2-1),
        sqrtprod(j1-m1, k1+n1+1) *
        sqrtprod(k2-n2+1, j2+m2))

    add_pair!(H,row,label,lookup,
        (tj1-1,tm1+1,tk1+1,tn1+1),
        (tj2+1,tm2-1,tk2-1,tn2-1),
        sqrtprod(j1-m1, k1+n1+1) *
        sqrtprod(k2+n2, j2-m2+1))

    # T_{--}(1) T_{++}(2)
    add_pair!(H,row,label,lookup,
        (tj1+1,tm1-1,tk1-1,tn1-1),
        (tj2-1,tm2+1,tk2+1,tn2+1),
        sqrtprod(j1-m1+1, k1+n1) *
        sqrtprod(j2-m2, k2+n2+1))

    add_pair!(H,row,label,lookup,
        (tj1+1,tm1-1,tk1-1,tn1-1),
        (tj2+1,tm2+1,tk2-1,tn2+1),
        sqrtprod(j1-m1+1, k1+n1) *
        sqrtprod(k2-n2, j2+m2+1))

    add_pair!(H,row,label,lookup,
        (tj1-1,tm1-1,tk1+1,tn1-1),
        (tj2-1,tm2+1,tk2+1,tn2+1),
        sqrtprod(j1+m1, k1-n1+1) *
        sqrtprod(k2+n2+1, j2-m2))

    add_pair!(H,row,label,lookup,
        (tj1-1,tm1-1,tk1+1,tn1-1),
        (tj2+1,tm2+1,tk2-1,tn2+1),
        sqrtprod(j1+m1, k1-n1+1) *
        sqrtprod(k2-n2, j2+m2+1))

    # T_{-+}(1) T_{+-}(2)
    add_pair!(H,row,label,lookup,
        (tj1-1,tm1-1,tk1+1,tn1+1),
        (tj2+1,tm2+1,tk2-1,tn2-1),
        sqrtprod(j1+m1, k1+n1+1) *
        sqrtprod(k2+n2, j2+m2+1))

    add_pair!(H,row,label,lookup,
        (tj1-1,tm1-1,tk1+1,tn1+1),
        (tj2-1,tm2+1,tk2+1,tn2-1),
        -sqrtprod(j1+m1, k1+n1+1) *
        sqrtprod(k2-n2+1, j2-m2))

    add_pair!(H,row,label,lookup,
        (tj1+1,tm1-1,tk1-1,tn1+1),
        (tj2+1,tm2+1,tk2-1,tn2-1),
        -sqrtprod(j1-m1+1, k1-n1) *
        sqrtprod(k2+n2, j2+m2+1))

    add_pair!(H,row,label,lookup,
        (tj1+1,tm1-1,tk1-1,tn1+1),
        (tj2-1,tm2+1,tk2+1,tn2-1),
        sqrtprod(j1-m1+1, k1-n1) *
        sqrtprod(k2-n2+1, j2-m2))

    # T_{+-}(1) T_{-+}(2)
    add_pair!(H,row,label,lookup,
        (tj1-1,tm1+1,tk1+1,tn1-1),
        (tj2+1,tm2-1,tk2-1,tn2+1),
        sqrtprod(j1-m1, k1-n1+1) *
        sqrtprod(k2-n2, j2-m2+1))

    add_pair!(H,row,label,lookup,
        (tj1-1,tm1+1,tk1+1,tn1-1),
        (tj2-1,tm2-1,tk2+1,tn2+1),
        -sqrtprod(j1-m1, k1-n1+1) *
        sqrtprod(k2+n2+1, j2+m2))

    add_pair!(H,row,label,lookup,
        (tj1+1,tm1+1,tk1-1,tn1-1),
        (tj2+1,tm2-1,tk2-1,tn2+1),
        -sqrtprod(j1+m1+1, k1+n1) *
        sqrtprod(k2-n2, j2-m2+1))

    add_pair!(H,row,label,lookup,
        (tj1+1,tm1+1,tk1-1,tn1-1),
        (tj2-1,tm2-1,tk2+1,tn2+1),
        sqrtprod(j1+m1+1, k1+n1) *
        sqrtprod(k2+n2+1, j2+m2))

    return nothing
end

# ------------------------------------------------------------
# 8. Build C2_SO5, J^2, K^2 in fixed sector
# ------------------------------------------------------------

function build_casimir_matrices(p::Int,
                                product_states::Vector{ProductState},
                                basis::Vector{SO4State},
                                label::Dict{SO4State,Int},
                                lookup::Dict{ProductState,Int})

    d = length(product_states)

    H  = zeros(Float64, d, d)
    HL = zeros(Float64, d, d)
    HR = zeros(Float64, d, d)

    for row in eachindex(product_states)
        idx1, idx2 = product_states[row]

        s1 = basis[idx1]
        s2 = basis[idx2]

        tj1,tm1,tk1,tn1 = s1
        tj2,tm2,tk2,tn2 = s2

        j1 = qval(tj1); m1 = qval(tm1); k1 = qval(tk1); n1 = qval(tn1)
        j2 = qval(tj2); m2 = qval(tm2); k2 = qval(tk2); n2 = qval(tn2)

        H[row,row]  += p^2 + 4p + 4m1*m2 + 4n1*n2
        HL[row,row] += j1*(j1+1.0) + j2*(j2+1.0) + 2m1*m2
        HR[row,row] += k1*(k1+1.0) + k2*(k2+1.0) + 2n1*n2

        # Left SU(2): J_+(1)J_-(2)
        if tm1 < tj1 && tm2 > -tj2
            c = jplus_coeff2(tj1,tm1) * jminus_coeff2(tj2,tm2)

            snew1 = (tj1,tm1+2,tk1,tn1)
            snew2 = (tj2,tm2-2,tk2,tn2)

            add_pair!(H,  row, label, lookup, snew1, snew2, 2c)
            add_pair!(HL, row, label, lookup, snew1, snew2, c)
        end

        # Left SU(2): J_-(1)J_+(2)
        if tm1 > -tj1 && tm2 < tj2
            c = jminus_coeff2(tj1,tm1) * jplus_coeff2(tj2,tm2)

            snew1 = (tj1,tm1-2,tk1,tn1)
            snew2 = (tj2,tm2+2,tk2,tn2)

            add_pair!(H,  row, label, lookup, snew1, snew2, 2c)
            add_pair!(HL, row, label, lookup, snew1, snew2, c)
        end

        # Right SU(2): K_+(1)K_-(2)
        if tn1 < tk1 && tn2 > -tk2
            c = jplus_coeff2(tk1,tn1) * jminus_coeff2(tk2,tn2)

            snew1 = (tj1,tm1,tk1,tn1+2)
            snew2 = (tj2,tm2,tk2,tn2-2)

            add_pair!(H,  row, label, lookup, snew1, snew2, 2c)
            add_pair!(HR, row, label, lookup, snew1, snew2, c)
        end

        # Right SU(2): K_-(1)K_+(2)
        if tn1 > -tk1 && tn2 < tk2
            c = jminus_coeff2(tk1,tn1) * jplus_coeff2(tk2,tn2)

            snew1 = (tj1,tm1,tk1,tn1-2)
            snew2 = (tj2,tm2,tk2,tn2+2)

            add_pair!(H,  row, label, lookup, snew1, snew2, 2c)
            add_pair!(HR, row, label, lookup, snew1, snew2, c)
        end

        add_T_terms!(H, row, s1, s2, label, lookup)
    end

    H  = 0.5 .* (H  .+ transpose(H))
    HL = 0.5 .* (HL .+ transpose(HL))
    HR = 0.5 .* (HR .+ transpose(HR))

    return H, HL, HR
end

# ------------------------------------------------------------
# 9. Single-particle J,K,T generators and L_ab construction
# ------------------------------------------------------------

"""
    build_single_particle_JKT(p,basis,label)

Construct single-particle generators in the SO(4)-adapted basis.

J_i and K_i are Hermitian SU(2) generators.
Tpp,Tpm,Tmp,Tmm are non-Hermitian spherical components of the SO(5)/SO(4)
vector, chosen to reproduce the T cross terms in the C2 matrix.
"""
function build_single_particle_JKT(p::Int,
                                   basis::Vector{SO4State},
                                   label::Dict{SO4State,Int})

    d = length(basis)

    Jp = spzeros(ComplexF64,d,d)
    Jm = spzeros(ComplexF64,d,d)
    Jz = spzeros(ComplexF64,d,d)

    Kp = spzeros(ComplexF64,d,d)
    Km = spzeros(ComplexF64,d,d)
    Kz = spzeros(ComplexF64,d,d)

    Tpp = spzeros(ComplexF64,d,d)
    Tpm = spzeros(ComplexF64,d,d)  # T_{+-}
    Tmp = spzeros(ComplexF64,d,d)  # T_{-+}
    Tmm = spzeros(ComplexF64,d,d)

    for col in eachindex(basis)
        s = basis[col]
        tj,tm,tk,tn = s

        j = qval(tj); m = qval(tm)
        k = qval(tk); n = qval(tn)

        # J generators.
        Jz[col,col] += ComplexF64(m)
        if tm < tj
            add_single!(Jp,col,label,(tj,tm+2,tk,tn),jplus_coeff2(tj,tm))
        end
        if tm > -tj
            add_single!(Jm,col,label,(tj,tm-2,tk,tn),jminus_coeff2(tj,tm))
        end

        # K generators.
        Kz[col,col] += ComplexF64(n)
        if tn < tk
            add_single!(Kp,col,label,(tj,tm,tk,tn+2),jplus_coeff2(tk,tn))
        end
        if tn > -tk
            add_single!(Km,col,label,(tj,tm,tk,tn-2),jminus_coeff2(tk,tn))
        end

        # T_{++}
        add_single!(Tpp,col,label,
            (tj+1,tm+1,tk-1,tn+1),
            sqrtprod(j+m+1, k-n))
        add_single!(Tpp,col,label,
            (tj-1,tm+1,tk+1,tn+1),
            sqrtprod(j-m, k+n+1))

        # T_{--}
        add_single!(Tmm,col,label,
            (tj+1,tm-1,tk-1,tn-1),
            sqrtprod(j-m+1, k+n))
        add_single!(Tmm,col,label,
            (tj-1,tm-1,tk+1,tn-1),
            sqrtprod(j+m, k-n+1))

        # T_{+-}
        # IMPORTANT SIGN CONVENTION:
        # These two signs are fixed by the original add_T_terms! bilinear
        # C2 cross terms.  With this convention adjoint(T_{+-}) = T_{-+}.
        add_single!(Tpm,col,label,
            (tj+1,tm+1,tk-1,tn-1),
            sqrtprod(j+m+1, k+n))
        add_single!(Tpm,col,label,
            (tj-1,tm+1,tk+1,tn-1),
            -sqrtprod(j-m, k-n+1))

        # T_{-+}
        add_single!(Tmp,col,label,
            (tj-1,tm-1,tk+1,tn+1),
            sqrtprod(j+m, k+n+1))
        add_single!(Tmp,col,label,
            (tj+1,tm-1,tk-1,tn+1),
            -sqrtprod(j-m+1, k-n))
    end

    J1 = 0.5 .* (Jp .+ Jm)
    J2 = (Jp .- Jm) ./ (2im)
    J3 = Jz

    K1 = 0.5 .* (Kp .+ Km)
    K2 = (Kp .- Km) ./ (2im)
    K3 = Kz

    return (
        J = [J1,J2,J3],
        K = [K1,K2,K3],
        Tpp = Tpp,
        Tpm = Tpm,
        Tmp = Tmp,
        Tmm = Tmm
    )
end

"""
    build_so5_L_single(p,basis,label)

Construct the ten single-particle SO(5) generators L_ab, a,b=1,...,5.

Convention:

    A_i = J_i + K_i,
    B_i = J_i - K_i,

    L_ij = epsilon_ijk A_k,
    L_i4 = B_i,

and the coset vector V_mu = L_{mu,5} is built from T components by the
Hermitian convention compatible with the original C2 cross terms

    V1 = (T_{+-} + T_{-+})/(2i),
    V2 = (T_{+-} - T_{-+})/2,
    V3 = (T_{++} - T_{--})/(2i),
    V4 = (T_{++} + T_{--})/2.

With the sign convention used above, Tpm_adj = -Tmp and Tpp_adj = Tmm,
so V_mu are Hermitian.
"""
function build_so5_L_single(p::Int,
                            basis::Vector{SO4State},
                            label::Dict{SO4State,Int})

    if haskey(__SO5_SINGLE_L_CACHE__, p)
        cached = __SO5_SINGLE_L_CACHE__[p]
        return cached.L, cached.gens
    end

    gens = build_single_particle_JKT(p,basis,label)

    J = gens.J
    K = gens.K

    d = length(basis)
    zero = spzeros(ComplexF64,d,d)

    L = Array{SparseMatrixCSC{ComplexF64,Int},2}(undef,5,5)

    for a in 1:5, b in 1:5
        L[a,b] = copy(zero)
    end

    A = [J[i] .+ K[i] for i in 1:3]
    B = [J[i] .- K[i] for i in 1:3]

    # SO(4) subalgebra.
    L[1,2] = A[3]
    L[2,1] = -A[3]

    L[2,3] = A[1]
    L[3,2] = -A[1]

    L[3,1] = A[2]
    L[1,3] = -A[2]

    for i in 1:3
        L[i,4] = B[i]
        L[4,i] = -B[i]
    end

    # Coset SO(5)/SO(4).
    Tpp = gens.Tpp
    Tpm = gens.Tpm
    Tmp = gens.Tmp
    Tmm = gens.Tmm

    # Correct Cartesian components of the SO(5)/SO(4) vector.
    #
    # This is the crucial convention.  It is fixed by requiring that the
    # ten matrices L_ab satisfy the SO(5) commutation relations and that
    # the single-particle C2 and C4 are scalar on the irrep (p,0).
    #
    # With the T signs above:
    #
    #     adjoint(Tpp) = Tmm,
    #     adjoint(Tpm) = Tmp.
    #
    # The Cartesian vector components are
    #
    #     V1 = (Tpp + Tmm)/2,
    #     V2 = (Tpp - Tmm)/(2i),
    #     V3 = (Tpm + Tmp)/2,
    #     V4 = (Tpm - Tmp)/(2i).
    #
    # This ordering is different from the naive guess
    # (V1,V2,V3,V4) = functions of (Tpm,Tmp,Tpp,Tmm).  The diagnostics
    # single_particle_scalar_residual and the commutators with J^2,K^2
    # are very sensitive to this convention.
    V1 = 0.5 .* (Tpp .+ Tmm)
    V2 = (Tpp .- Tmm) ./ (2im)
    V3 = 0.5 .* (Tpm .+ Tmp)
    V4 = (Tpm .- Tmp) ./ (2im)

    V = [V1,V2,V3,V4]

    for μ in 1:4
        L[μ,5] = V[μ]
        L[5,μ] = -V[μ]
    end

    __SO5_SINGLE_L_CACHE__[p] = (L = L, gens = gens)

    return L, gens
end

# ------------------------------------------------------------
# 10. Quartic Casimir C4 from L_ab -- fast fixed-sector version
# ------------------------------------------------------------

function levi_civita5(a::Int,b::Int,c::Int,d::Int,e::Int)
    x = (a,b,c,d,e)

    if length(unique(x)) < 5
        return 0
    end

    inv = 0
    y = collect(x)

    for i in 1:4
        for j in (i+1):5
            if y[i] > y[j]
                inv += 1
            end
        end
    end

    return iseven(inv) ? 1 : -1
end

function sparse_identity_complex(d::Int)
    return spdiagm(0 => ones(ComplexF64,d))
end

# SO(5) 反对称生成元只有 10 个独立的上三角分量。
const __SO5_UPPER_PAIRS__ = (
    (1,2), (1,3), (1,4), (1,5),
    (2,3), (2,4), (2,5),
    (3,4), (3,5),
    (4,5)
)

# 每个 tuple 是 (effective_epsilon,b,c,d,e)，并且始终满足 b<c、d<e。
# effective_epsilon 已经吸收了将下三角 L_de=-L_ed 改写到上三角时的符号。
const __SO5_C4_UPPER_PAIRINGS__ = (
    (( 1,2,3,4,5), (-1,2,4,3,5), ( 1,2,5,3,4)),
    ((-1,1,3,4,5), ( 1,1,4,3,5), (-1,1,5,3,4)),
    (( 1,1,2,4,5), (-1,1,4,2,5), ( 1,1,5,2,4)),
    ((-1,1,2,3,5), ( 1,1,3,2,5), (-1,1,5,2,3)),
    (( 1,1,2,3,4), (-1,1,3,2,4), ( 1,1,4,2,3))
)

@inline function so5_upper_pair_index(a::Int, b::Int)
    @boundscheck 1 <= a < b <= 5 || throw(ArgumentError("require 1 <= a < b <= 5"))

    if a == 1
        return b - 1               # (1,2)..(1,5) -> 1..4
    elseif a == 2
        return b + 2               # (2,3)..(2,5) -> 5..7
    elseif a == 3
        return b + 4               # (3,4)..(3,5) -> 8..9
    else
        return 10                  # (4,5)
    end
end

"""
只保存 a<b 的十个独立 SO(5) 生成元。

主流程中所有 C4 pairing 已经统一改写为上三角分量，因此不需要保存
L_ba=-L_ab 的第二份稀疏矩阵。
"""
struct SO5UpperGenerators
    mats::NTuple{10,SparseMatrixCSC{ComplexF64,Int}}
    dim::Int
end

@inline function upper_L(G::SO5UpperGenerators, a::Int, b::Int)
    @boundscheck a < b || throw(ArgumentError("upper_L requires a < b"))
    return G.mats[so5_upper_pair_index(a,b)]
end

function upper_generators_from_full(
    L::Array{SparseMatrixCSC{ComplexF64,Int},2}
)
    mats = ntuple(i -> begin
        a,b = __SO5_UPPER_PAIRS__[i]
        L[a,b]
    end, 10)

    return SO5UpperGenerators(mats, size(L[1,1],1))
end

"""
    build_total_L_two_particle(Lsp)

兼容旧代码的完整 5x5 total-L 构造。高性能主流程不会调用这个函数。
"""
function build_total_L_two_particle(
    Lsp::Array{SparseMatrixCSC{ComplexF64,Int},2}
)
    d1 = size(Lsp[1,1],1)
    I1 = sparse_identity_complex(d1)
    Z  = spzeros(ComplexF64,d1*d1,d1*d1)

    Ltot = Array{SparseMatrixCSC{ComplexF64,Int},2}(undef,5,5)

    for a in 1:5, b in 1:5
        if nnz(Lsp[a,b]) == 0
            Ltot[a,b] = copy(Z)
        else
            Ltot[a,b] = kron(Lsp[a,b],I1) .+ kron(I1,Lsp[a,b])
            dropzeros!(Ltot[a,b])
        end
    end

    return Ltot
end

function build_total_L_two_particle_cached(
    p::Int,
    Lsp::Array{SparseMatrixCSC{ComplexF64,Int},2};
    use_cache::Bool=true
)
    if use_cache && haskey(__SO5_TOTAL_L_CACHE__, p)
        return __SO5_TOTAL_L_CACHE__[p].Ltot
    end

    Ltot = build_total_L_two_particle(Lsp)

    if use_cache
        __SO5_TOTAL_L_CACHE__[p] = (Ltot = Ltot,)
    end

    return Ltot
end

"""
    build_total_L_two_particle_upper(Lsp)

高性能版本：只构造十个独立的

    L_ab^tot = L_ab kron I + I kron L_ab,    a<b.

相比完整 5x5 数组，避免了下三角生成元的重复 Kronecker 构造和缓存。
"""
function build_total_L_two_particle_upper(
    Lsp::Array{SparseMatrixCSC{ComplexF64,Int},2}
)
    d1 = size(Lsp[1,1],1)
    I1 = sparse_identity_complex(d1)

    mats = ntuple(idx -> begin
        a,b = __SO5_UPPER_PAIRS__[idx]
        Lab = Lsp[a,b]
        T = kron(Lab,I1) .+ kron(I1,Lab)
        dropzeros!(T)
        T
    end, 10)

    return SO5UpperGenerators(mats, d1*d1)
end

function build_total_L_two_particle_upper_cached(
    p::Int,
    Lsp::Array{SparseMatrixCSC{ComplexF64,Int},2};
    use_cache::Bool=true
)
    if use_cache && haskey(__SO5_TOTAL_L_UPPER_CACHE__, p)
        return __SO5_TOTAL_L_UPPER_CACHE__[p].Ltot_upper
    end

    Ltot_upper = build_total_L_two_particle_upper(Lsp)

    if use_cache
        __SO5_TOTAL_L_UPPER_CACHE__[p] = (Ltot_upper = Ltot_upper,)
    end

    return Ltot_upper
end

"""
    build_C4_from_upper(L; return_W=false)

在一个完整空间上构造 C4。这个函数主要用于尺寸很小的单粒子 C4 calibration。
内部逐个构造 W_a，并立即累加 W_a^2；默认不同时保存五个 W_a。
"""
function build_C4_from_upper(
    L::SO5UpperGenerators;
    verbose::Bool=false,
    return_W::Bool=false
)
    dmat = L.dim
    H4 = spzeros(ComplexF64,dmat,dmat)

    Wsave = return_W ?
        Vector{SparseMatrixCSC{ComplexF64,Int}}(undef,5) : nothing

    for a in 1:5
        pairings = __SO5_C4_UPPER_PAIRINGS__[a]

        eps,b,c,d,e = pairings[1]
        Lbc = upper_L(L,b,c)
        Lde = upper_L(L,d,e)
        Wa = (eps/2.0) .* (Lbc*Lde .+ Lde*Lbc)

        for q in 2:3
            eps,b,c,d,e = pairings[q]
            Lbc = upper_L(L,b,c)
            Lde = upper_L(L,d,e)
            Wa = Wa .+ (eps/2.0) .* (Lbc*Lde .+ Lde*Lbc)
        end

        # 数学上 Wa 已经 Hermitian；保留原实现的显式数值对称化。
        Wa = 0.5 .* (Wa .+ Wa')
        dropzeros!(Wa)

        H4 = H4 .+ Wa*Wa

        if return_W
            Wsave[a] = Wa
        end
    end

    H4 = 0.5 .* (H4 .+ H4')
    dropzeros!(H4)

    if verbose
        println("Built full C4 on dimension $dmat with streaming W_a accumulation.")
    end

    return H4, Wsave
end

"""
兼容原接口。主流程只在单粒子 calibration 中调用。
"""
function build_C4_from_L(
    L::Array{SparseMatrixCSC{ComplexF64,Int},2};
    verbose::Bool=false,
    return_W::Bool=true
)
    return build_C4_from_upper(
        upper_generators_from_full(L);
        verbose=verbose,
        return_W=return_W
    )
end

"""
    _build_Wa_sector_gram(a,Ltot,Lcols)

直接构造固定输入扇区 S 上的 W_a P_S，并返回 Gram matrix

    (W_a P_S)' * (W_a P_S) = P_S W_a^2 P_S.

因为每个上三角 L_ab 都是 Hermitian，反对易组合 W_a 也是 Hermitian，
所以不需要形成完整 W_a，更不需要形成完整 W_a^2。
"""
function _build_Wa_sector_gram(
    a::Int,
    Ltot::SO5UpperGenerators,
    Lcols::NTuple{10,SparseMatrixCSC{ComplexF64,Int}}
)
    pairings = __SO5_C4_UPPER_PAIRINGS__[a]

    eps,b,c,d,e = pairings[1]
    ibc = so5_upper_pair_index(b,c)
    ide = so5_upper_pair_index(d,e)
    Lbc = Ltot.mats[ibc]
    Lde = Ltot.mats[ide]

    Wa_cols = (eps/2.0) .* (
        Lbc * Lcols[ide] .+ Lde * Lcols[ibc]
    )

    for q in 2:3
        eps,b,c,d,e = pairings[q]
        ibc = so5_upper_pair_index(b,c)
        ide = so5_upper_pair_index(d,e)
        Lbc = Ltot.mats[ibc]
        Lde = Ltot.mats[ide]

        Wa_cols = Wa_cols .+ (eps/2.0) .* (
            Lbc * Lcols[ide] .+ Lde * Lcols[ibc]
        )
    end

    dropzeros!(Wa_cols)
    gram = Matrix(adjoint(Wa_cols) * Wa_cols)
    return gram, nnz(Wa_cols)
end

"""
    build_C4_sector_from_upper(Ltot,sector_indices; parallel=false)

这是 v5 的核心提速函数。原算法为：

    构造完整 W_a (N x N) -> 构造完整 C4 (N x N) -> 截取 d x d 扇区。

这里改为：

    只取 L_ab 的扇区列 -> 构造 W_a P_S (N x d)
    -> 计算 (W_a P_S)'(W_a P_S) (d x d)。

当 d 远小于完整乘积维数 N 时，时间和内存都会大幅下降。
"""
function build_C4_sector_from_upper(
    Ltot::SO5UpperGenerators,
    sector_indices::Vector{Int};
    parallel::Bool=false,
    verbose::Bool=false
)
    dsec = length(sector_indices)

    if dsec == 0
        return zeros(ComplexF64,0,0), zeros(Int,5)
    end

    # 每个独立生成元只截取一次目标扇区列，后续 15 个 pairing 重复使用。
    Lcols = ntuple(k -> Ltot.mats[k][:,sector_indices], 10)

    wa_nnz = zeros(Int,5)
    H4 = zeros(ComplexF64,dsec,dsec)

    if parallel && Threads.nthreads() > 1
        # 并行路径：五个 a 独立计算，最后按固定顺序求和。
        gram_parts = Vector{Matrix{ComplexF64}}(undef,5)
        Threads.@threads for a in 1:5
            gram_parts[a], wa_nnz[a] = _build_Wa_sector_gram(a,Ltot,Lcols)
        end
        for a in 1:5
            H4 .+= gram_parts[a]
        end
    else
        # 串行路径：每个 Gram matrix 立即累加，降低峰值内存。
        for a in 1:5
            gram, wa_nnz[a] = _build_Wa_sector_gram(a,Ltot,Lcols)
            H4 .+= gram
        end
    end

    # 固定顺序累加后再显式 Hermitize，避免并行路径中的微小舍入不对称。
    @inbounds for j in 1:dsec
        H4[j,j] = ComplexF64(real(H4[j,j]),0.0)
        for i in 1:(j-1)
            z = 0.5 * (H4[i,j] + conj(H4[j,i]))
            H4[i,j] = z
            H4[j,i] = conj(z)
        end
    end

    if verbose
        println("Built fixed-sector C4 directly: full dim=$(Ltot.dim), sector dim=$dsec")
        println("nnz(W_a P_sector) = ", wa_nnz)
        println("parallel C4 construction = ", parallel && Threads.nthreads() > 1)
    end

    return H4, wa_nnz
end

function dense_real_checked(A;
                            imag_tol::Float64=1e-8,
                            name::String="matrix",
                            allow_warning::Bool=true)

    # 只分配一个复数 dense 输入和一个最终实矩阵，避免 imag/abs/real 广播临时数组。
    Ad = A isa Matrix{ComplexF64} ? A : Matrix(A)
    R = Matrix{Float64}(undef,size(Ad))
    imax = 0.0

    @inbounds for idx in eachindex(Ad)
        z = Ad[idx]
        R[idx] = real(z)
        imax = max(imax,abs(imag(z)))
    end

    if imax > imag_tol
        msg = "Imaginary part of $name is not negligible: max |Im| = $imax."
        if allow_warning
            @warn msg
        else
            error(msg)
        end
    end

    return R, imax
end

"""
    build_C4_sector_matrix(...; c4_method=:sector_columns)

默认使用 :sector_columns 直接构造固定扇区 C4。
保留 :full_legacy 仅用于小 p 下与原算法交叉验证。
"""
function build_C4_sector_matrix(p::Int,
                                product_states::Vector{ProductState},
                                basis::Vector{SO4State},
                                label::Dict{SO4State,Int};
                                calibrate::Bool=true,
                                imag_tol::Float64=1e-8,
                                use_cache::Bool=true,
                                c4_method::Symbol=:sector_columns,
                                parallel_C4::Bool=false,
                                verbose::Bool=false)

    d1 = length(basis)

    # Single-particle L_ab.
    Lsp, gens = build_so5_L_single(p,basis,label)

    # Single-particle C4 calibration. This depends only on p and calibration mode.
    calib_key = (p,calibrate)

    if use_cache && haskey(__SO5_C4_CALIB_CACHE__,calib_key)
        calib = __SO5_C4_CALIB_CACHE__[calib_key]
        target_sp = calib.target_sp
        raw_mean = calib.raw_mean
        sp_scalar_residual = calib.sp_scalar_residual
        imax_sp = calib.imax_sp
        scale = calib.scale
    else
        # 单粒子空间很小；这里使用 streaming full-C4，且不保存五个 W_a。
        H4_sp_raw, _ = build_C4_from_L(
            Lsp;
            verbose=false,
            return_W=false
        )
        H4_sp_dense, imax_sp = dense_real_checked(
            H4_sp_raw;
            imag_tol=imag_tol,
            name="single-particle C4",
            allow_warning=true
        )

        target_sp = so5_casimir4(p,0)
        raw_mean = tr(H4_sp_dense) / d1

        # 不构造 dense identity；仅修改对角元计算标量残差。
        residual_matrix = copy(H4_sp_dense)
        @inbounds for i in 1:d1
            residual_matrix[i,i] -= raw_mean
        end
        sp_scalar_residual = norm(residual_matrix) / max(norm(H4_sp_dense),1e-14)

        scale = 1.0
        if calibrate
            if abs(raw_mean) < 1e-12
                @warn "Single-particle raw C4 mean is too small; no C4 normalization calibration applied."
            else
                scale = target_sp / raw_mean
            end
        end

        if use_cache
            __SO5_C4_CALIB_CACHE__[calib_key] = (
                target_sp=target_sp,
                raw_mean=raw_mean,
                sp_scalar_residual=sp_scalar_residual,
                imax_sp=imax_sp,
                scale=scale
            )
        end
    end

    if verbose
        println("--------- C4 single-particle calibration ---------")
        println("target C4(p,0)         = ",target_sp)
        println("raw single-particle C4 = ",raw_mean)
        println("C4 scale               = ",scale)
        println("single-particle C4 scalar residual = ",sp_scalar_residual)
        println("max imaginary part in C4_sp = ",imax_sp)
        println("--------------------------------------------------")
    end

    sector_indices = Vector{Int}(undef,length(product_states))
    @inbounds for a in eachindex(product_states)
        i,j = product_states[a]
        sector_indices[a] = (i-1)*d1 + j
    end

    wa_nnz = nothing

    if c4_method === :sector_columns
        Ltot_upper = build_total_L_two_particle_upper_cached(
            p,Lsp; use_cache=use_cache
        )
        H4_sector_complex, wa_nnz = build_C4_sector_from_upper(
            Ltot_upper,
            sector_indices;
            parallel=parallel_C4,
            verbose=verbose
        )
        H4_sector_dense, imax_sec = dense_real_checked(
            H4_sector_complex;
            imag_tol=imag_tol,
            name="sector C4",
            allow_warning=true
        )

    elseif c4_method === :full_legacy
        # 仅用于小 p 验证。这个分支保留原来的 full-space 算法。
        Ltot = build_total_L_two_particle_cached(p,Lsp;use_cache=use_cache)
        H4_full_raw, _ = build_C4_from_L(
            Ltot;
            verbose=verbose,
            return_W=false
        )
        H4_sector_raw = H4_full_raw[sector_indices,sector_indices]
        H4_sector_dense, imax_sec = dense_real_checked(
            H4_sector_raw;
            imag_tol=imag_tol,
            name="sector C4 legacy",
            allow_warning=true
        )
    else
        error("Unknown c4_method=$c4_method. Use :sector_columns or :full_legacy.")
    end

    H4_sector = scale .* H4_sector_dense

    # 原地对称化，避免 H4 + transpose(H4) 的额外大矩阵分配。
    dsec = size(H4_sector,1)
    @inbounds for j in 1:dsec
        for i in 1:(j-1)
            x = 0.5 * (H4_sector[i,j] + H4_sector[j,i])
            H4_sector[i,j] = x
            H4_sector[j,i] = x
        end
    end

    info = (
        scale=scale,
        target_single_particle_C4=target_sp,
        raw_single_particle_C4_mean=raw_mean,
        single_particle_scalar_residual=sp_scalar_residual,
        single_particle_imag_max=imax_sp,
        sector_imag_max=imax_sec,
        calibrated=calibrate,
        method=c4_method,
        wa_sector_nnz=wa_nnz,
        parallel_C4=(parallel_C4 && Threads.nthreads() > 1)
    )

    return H4_sector,info
end

# ------------------------------------------------------------
# 11. Expectation, residual, and commutator diagnostics
# ------------------------------------------------------------

expectation(v::AbstractVector, A::AbstractMatrix) =
    real(dot(v, A*v))

angular_momentum_from_casimir(x::Real) =
    sqrt(max(0.0, Float64(x) + 0.25)) - 0.5

function fix_phase(v::Vector{Float64})
    isempty(v) && return v

    i = argmax(abs.(v))

    if v[i] < 0
        return -v
    else
        return v
    end
end

function residual_norm(A::AbstractMatrix, target::Real, v::AbstractVector)
    nv = norm(v)

    if nv < 1e-14
        return Inf
    end

    return norm(A*v .- Float64(target).*v) / nv
end

function relative_commutator_norm(A::AbstractMatrix, B::AbstractMatrix)
    nA = max(norm(A), 1e-14)
    nB = max(norm(B), 1e-14)
    return norm(A*B - B*A) / (nA*nB)
end

function commutator_diagnostics(H::AbstractMatrix,
                                HL::AbstractMatrix,
                                HR::AbstractMatrix;
                                H4=nothing,
                                verbose::Bool=true)

    c_H_HL  = relative_commutator_norm(H,HL)
    c_H_HR  = relative_commutator_norm(H,HR)
    c_HL_HR = relative_commutator_norm(HL,HR)

    if H4 === nothing
        if verbose
            println("--------- commutator diagnostics ---------")
            println("relative ||[H,HL]||  = ", c_H_HL)
            println("relative ||[H,HR]||  = ", c_H_HR)
            println("relative ||[HL,HR]|| = ", c_HL_HR)
            println("------------------------------------------")
        end

        return (
            H_HL  = c_H_HL,
            H_HR  = c_H_HR,
            HL_HR = c_HL_HR
        )
    else
        c_H_H4   = relative_commutator_norm(H,H4)
        c_H4_HL  = relative_commutator_norm(H4,HL)
        c_H4_HR  = relative_commutator_norm(H4,HR)

        if verbose
            println("--------- commutator diagnostics ---------")
            println("relative ||[H,HL]||   = ", c_H_HL)
            println("relative ||[H,HR]||   = ", c_H_HR)
            println("relative ||[HL,HR]||  = ", c_HL_HR)
            println("relative ||[H,H4]||   = ", c_H_H4)
            println("relative ||[H4,HL]||  = ", c_H4_HL)
            println("relative ||[H4,HR]||  = ", c_H4_HR)
            println("------------------------------------------")
        end

        return (
            H_HL  = c_H_HL,
            H_HR  = c_H_HR,
            HL_HR = c_HL_HR,
            H_H4  = c_H_H4,
            H4_HL = c_H4_HL,
            H4_HR = c_H4_HR
        )
    end
end


"""
    commutator_diagnostics_cached(cache_key,...)

交换子只依赖扇区矩阵，不依赖目标 (P,Q,J,K)。同一扇区最多计算一次。
"""
function commutator_diagnostics_cached(
    cache_key,
    H::AbstractMatrix,
    HL::AbstractMatrix,
    HR::AbstractMatrix;
    H4=nothing,
    use_cache::Bool=true,
    verbose::Bool=true
)
    key = (cache_key,H4 === nothing ? :without_H4 : :with_H4)

    if use_cache && haskey(__SO5_COMMUTATOR_CACHE__,key)
        verbose && println("Using cached commutator diagnostics for key = ",key)
        return __SO5_COMMUTATOR_CACHE__[key]
    end

    result = commutator_diagnostics(H,HL,HR;H4=H4,verbose=verbose)
    if use_cache
        __SO5_COMMUTATOR_CACHE__[key] = result
    end
    return result
end

# ------------------------------------------------------------
# 12. Composite diagonalization and residual selection
# ------------------------------------------------------------

@inline function symmetrize_real_inplace!(A::Matrix{Float64})
    d = size(A,1)
    @boundscheck size(A,2) == d || throw(DimensionMismatch("matrix must be square"))

    @inbounds for j in 1:d
        for i in 1:(j-1)
            x = 0.5 * (A[i,j] + A[j,i])
            A[i,j] = x
            A[j,i] = x
        end
    end
    return A
end

function diagonalize_composite(H::Matrix{Float64},
                               H4::Matrix{Float64},
                               HL::Matrix{Float64},
                               HR::Matrix{Float64};
                               epsC4::Float64=1e-5,
                               epsJ::Float64=1e-3,
                               epsK::Float64=1e-6)

    # 单次 fused broadcast 构造 O，避免四个大矩阵表达式产生临时数组。
    O = similar(H)
    @. O = H + epsC4*H4 + epsJ*HL + epsK*HR
    symmetrize_real_inplace!(O)

    eigO = eigen(Symmetric(O))
    return eigO,O
end

"""
    diagonalize_composite_cached(...)

同一扇区及同一组 eps 只进行一次完整对角化。
当扫描多个 (P,Q,J,K) 时，这是热缓存阶段最重要的加速。
"""
function diagonalize_composite_cached(
    sector_cache_key,
    H::Matrix{Float64},
    H4::Matrix{Float64},
    HL::Matrix{Float64},
    HR::Matrix{Float64};
    epsC4::Float64=1e-5,
    epsJ::Float64=1e-3,
    epsK::Float64=1e-6,
    use_cache::Bool=true,
    verbose::Bool=false
)
    key = (sector_cache_key,epsC4,epsJ,epsK)

    if use_cache && haskey(__SO5_COMPOSITE_EIG_CACHE__,key)
        verbose && println("Using cached composite eigensystem for key = ",key)
        return __SO5_COMPOSITE_EIG_CACHE__[key]
    end

    eigO,O = diagonalize_composite(
        H,H4,HL,HR;
        epsC4=epsC4,
        epsJ=epsJ,
        epsK=epsK
    )
    result = (eigO=eigO,O=O)

    if use_cache
        __SO5_COMPOSITE_EIG_CACHE__[key] = result
    end
    return result
end

"""
一次 matrix-vector product 同时计算期望值与目标本征方程残差。
原代码分别调用 expectation 和 residual_norm，会对同一矩阵重复乘 v。
"""
function expectation_and_residual(
    A::AbstractMatrix,
    target::Real,
    v::AbstractVector
)
    nv = norm(v)
    nv < 1e-14 && return (NaN,Inf)

    Av = A*v
    expval = real(dot(v,Av)) / (nv*nv)

    @inbounds for i in eachindex(Av,v)
        Av[i] -= Float64(target)*v[i]
    end
    residual = norm(Av)/nv
    return expval,residual
end

"""
    select_target_vectors_from_composite(...)

先根据 O 的目标本征值筛掉绝大多数不可能的 eigenvectors，再用
C2、C4、J2、K2 四种残差做最终严格判定。
"""
function select_target_vectors_from_composite(eigO,
                                              H::Matrix{Float64},
                                              H4::Matrix{Float64},
                                              HL::Matrix{Float64},
                                              HR::Matrix{Float64},
                                              P::Int,
                                              Q::Int,
                                              tJ::Int,
                                              tK::Int;
                                              epsC4::Float64=1e-5,
                                              epsJ::Float64=1e-3,
                                              epsK::Float64=1e-6,
                                              tolC::Float64=1e-8,
                                              tolC4::Float64=1e-6,
                                              tolJ::Float64=1e-8,
                                              tolK::Float64=1e-8,
                                              prefilter::Bool=true,
                                              prefilter_safety::Float64=10.0)

    targetC  = so5_casimir(P,Q)
    targetC4 = so5_casimir4(P,Q)

    Jval = qval(tJ)
    Kval = qval(tK)
    targetJ2 = Jval*(Jval+1.0)
    targetK2 = Kval*(Kval+1.0)

    targetO = targetC + epsC4*targetC4 + epsJ*targetJ2 + epsK*targetK2

    # 若一个向量满足四个残差容差，则它对 O 的残差至多为该加权和。
    base_window = tolC + abs(epsC4)*tolC4 +
                  abs(epsJ)*tolJ + abs(epsK)*tolK
    candidate_window = max(
        prefilter_safety*base_window,
        100.0*eps(Float64)*(1.0+abs(targetO))
    )

    candidate_indices = prefilter ?
        findall(x -> abs(x-targetO) <= candidate_window,eigO.values) :
        collect(eachindex(eigO.values))

    d = size(H,1)
    selected = Vector{Vector{Float64}}()
    diagnostics = Vector{Dict{String,Float64}}()

    for a in candidate_indices
        v = fix_phase(Vector(@view eigO.vectors[:,a]))

        c_exp,resC   = expectation_and_residual(H, targetC, v)
        c4_exp,resC4 = expectation_and_residual(H4,targetC4,v)
        j_exp,resJ   = expectation_and_residual(HL,targetJ2,v)
        k_exp,resK   = expectation_and_residual(HR,targetK2,v)

        if resC < tolC && resC4 < tolC4 && resJ < tolJ && resK < tolK
            push!(selected,v)
            push!(diagnostics,Dict(
                "O_eigenvalue" => Float64(eigO.values[a]),
                "C2_SO5" => Float64(c_exp),
                "C4_SO5" => Float64(c4_exp),
                "J2" => Float64(j_exp),
                "K2" => Float64(k_exp),
                "J_measured" => Float64(angular_momentum_from_casimir(j_exp)),
                "K_measured" => Float64(angular_momentum_from_casimir(k_exp)),
                "resC" => Float64(resC),
                "resC4" => Float64(resC4),
                "resJ" => Float64(resJ),
                "resK" => Float64(resK)
            ))
        end
    end

    if isempty(selected)
        return zeros(Float64,d,0),diagnostics
    end
    return hcat(selected...),diagnostics
end

function remove_linearly_dependent_vectors_with_diagnostics(
    V::Matrix{Float64},
    diagnostics;
    tol::Float64=1e-10
)
    d,n = size(V)
    n == 0 && return V,diagnostics

    # 预分配 Q，避免原实现每加入一个向量都重新 hcat 全部已保留向量。
    Q = Matrix{Float64}(undef,d,n)
    kept_diag = typeof(diagnostics)()
    nkeep = 0

    for a in 1:n
        r = copy(@view V[:,a])

        # 两遍 modified Gram-Schmidt，提高接近简并情形的数值稳定性。
        for _ in 1:2
            for k in 1:nkeep
                q = @view Q[:,k]
                alpha = dot(q,r)
                @inbounds @simd for i in 1:d
                    r[i] -= alpha*q[i]
                end
            end
        end

        nr = norm(r)
        if nr > tol
            nkeep += 1
            @inbounds @simd for i in 1:d
                Q[i,nkeep] = r[i]/nr
            end
            push!(kept_diag,diagnostics[a])
        end
    end

    if nkeep == 0
        return zeros(Float64,d,0),kept_diag
    end
    return Q[:,1:nkeep],kept_diag
end

# ------------------------------------------------------------
# 13. SO(4) CG and reduced coefficient extraction
# ------------------------------------------------------------

function so4_cg2(cache::CGCache,
                 s1::SO4State,
                 s2::SO4State,
                 tJ::Int,
                 tM1::Int,
                 tK::Int,
                 tM2::Int)

    tj1,tm1,tk1,tn1 = s1
    tj2,tm2,tk2,tn2 = s2

    cgL = su2_cg2(cache, tj1,tm1,tj2,tm2,tJ,tM1)
    cgR = su2_cg2(cache, tk1,tn1,tk2,tn2,tK,tM2)

    return cgL * cgR
end

function reduced_key_string(channel::Int,
                            P::Int,
                            Q::Int,
                            tJ::Int,
                            tK::Int,
                            tj1::Int,
                            tk1::Int,
                            tj2::Int,
                            tk2::Int)

    return "channel=$channel : " *
           "SO5=($P,$Q), " *
           "SO4=(J=$(half_string(tJ)),K=$(half_string(tK))), " *
           "(j1=$(half_string(tj1)),k1=$(half_string(tk1))) ⊗ " *
           "(j2=$(half_string(tj2)),k2=$(half_string(tk2)))"
end

function extract_reduced_coefficients(vectors::Matrix{Float64},
                                      product_states::Vector{ProductState},
                                      basis::Vector{SO4State},
                                      cache::CGCache,
                                      P::Int,
                                      Q::Int,
                                      tJ::Int,
                                      tK::Int,
                                      tM1::Int,
                                      tM2::Int;
                                      cg_tol::Float64=1e-12)

    sums = Dict{NTuple{9,Int},Float64}()
    counts = Dict{NTuple{9,Int},Int}()
    nch = size(vectors,2)

    # SO(4) CG 与 outer channel 无关，先一次性预计算。
    valid_positions = Int[]
    cg_values = Float64[]
    channel_independent_labels = Vector{NTuple{4,Int}}()

    for a in eachindex(product_states)
        idx1,idx2 = product_states[a]
        s1 = basis[idx1]
        s2 = basis[idx2]
        c_so4 = so4_cg2(cache,s1,s2,tJ,tM1,tK,tM2)

        if abs(c_so4) >= cg_tol
            tj1,tm1,tk1,tn1 = s1
            tj2,tm2,tk2,tn2 = s2
            push!(valid_positions,a)
            push!(cg_values,c_so4)
            push!(channel_independent_labels,(tj1,tk1,tj2,tk2))
        end
    end

    for ch in 1:nch
        v = fix_phase(Vector(@view vectors[:,ch]))

        @inbounds for q in eachindex(valid_positions)
            a = valid_positions[q]
            reduced = v[a]/cg_values[q]
            tj1,tk1,tj2,tk2 = channel_independent_labels[q]
            key = (ch,P,Q,tJ,tK,tj1,tk1,tj2,tk2)
            sums[key] = get(sums,key,0.0) + reduced
            counts[key] = get(counts,key,0) + 1
        end
    end

    reduced_keyed = Dict{NTuple{9,Int},Float64}()
    reduced_pretty = Dict{String,Float64}()

    for (key,total) in sums
        val = total/counts[key]
        reduced_keyed[key] = val

        ch,P0,Q0,tJ0,tK0,tj1,tk1,tj2,tk2 = key
        pretty_key = reduced_key_string(ch,P0,Q0,tJ0,tK0,tj1,tk1,tj2,tk2)
        reduced_pretty[pretty_key] = val
    end

    return reduced_keyed,reduced_pretty
end

# ------------------------------------------------------------
# 14. Cached sector matrix construction
# ------------------------------------------------------------

function build_sector_matrices_cached(p::Int,
                                      tM1::Int,
                                      tM2::Int;
                                      calibrate_C4::Bool=true,
                                      c4_imag_tol::Float64=1e-8,
                                      c4_method::Symbol=:sector_columns,
                                      parallel_C4::Bool=false,
                                      use_cache::Bool=true,
                                      verbose::Bool=false)

    # H, H4, HL, HR depend only on the fixed magnetic sector and on the
    # C4 normalization convention.  They do NOT depend on the target
    # (P,Q,J,K).  This is the main cache for scanning many target channels.
    key = (p,tM1,tM2,calibrate_C4,c4_imag_tol,c4_method)

    if use_cache && haskey(__SO5_SECTOR_MATRIX_CACHE__, key)
        data = __SO5_SECTOR_MATRIX_CACHE__[key]
        if verbose
            println("Using cached sector matrices for key = ", key)
        end
        return data
    end

    product_states, basis, label = generate_product_sector_fixed_M(p,tM1,tM2)
    lookup = build_product_lookup(product_states)

    if isempty(product_states)
        data = (
            H = zeros(Float64,0,0),
            H4 = zeros(Float64,0,0),
            HL = zeros(Float64,0,0),
            HR = zeros(Float64,0,0),
            basis = basis,
            product_states = product_states,
            c4_info = nothing,
            cache_key = key,
            from_cache = false
        )
        if use_cache
            __SO5_SECTOR_MATRIX_CACHE__[key] = data
        end
        return data
    end

    H, HL, HR = build_casimir_matrices(p, product_states, basis, label, lookup)

    H4, c4_info = build_C4_sector_matrix(
        p,
        product_states,
        basis,
        label;
        calibrate = calibrate_C4,
        imag_tol = c4_imag_tol,
        use_cache = use_cache,
        c4_method = c4_method,
        parallel_C4 = parallel_C4,
        verbose = verbose
    )

    data = (
        H = H,
        H4 = H4,
        HL = HL,
        HR = HR,
        basis = basis,
        product_states = product_states,
        c4_info = c4_info,
        cache_key = key,
        from_cache = false
    )

    if use_cache
        __SO5_SECTOR_MATRIX_CACHE__[key] = data
    end

    return data
end

# ------------------------------------------------------------
# 14. Main calling function
# ------------------------------------------------------------

function compute_reduced_so5_cg_composite(p::Int,
                                           P::Int,
                                           Q::Int,
                                           J,
                                           M1,
                                           K,
                                           M2;
                                           epsC4::Float64=1e-5,
                                           epsJ::Float64=1e-3,
                                           epsK::Float64=1e-6,
                                           tolC::Float64=1e-8,
                                           tolC4::Float64=1e-6,
                                           tolJ::Float64=1e-8,
                                           tolK::Float64=1e-8,
                                           cg_tol::Float64=1e-12,
                                           check_commutators::Bool=false,
                                           cache_commutators::Bool=true,
                                           cache_eigensystem::Bool=true,
                                           prefilter_candidates::Bool=true,
                                           independence_tol::Float64=1e-10,
                                           calibrate_C4::Bool=true,
                                           c4_imag_tol::Float64=1e-8,
                                           c4_method::Symbol=:sector_columns,
                                           parallel_C4::Bool=false,
                                           use_cache::Bool=true,
                                           verbose::Bool=true)

    tJ  = to2(J)
    tM1 = to2(M1)
    tK  = to2(K)
    tM2 = to2(M2)

    if abs(tM1) > tJ
        error("Invalid quantum numbers: |M1| must be <= J.")
    end

    if abs(tM2) > tK
        error("Invalid quantum numbers: |M2| must be <= K.")
    end

    cache = get_shared_CGCache(max(128,4p + 2P + 2Q + 32))

    sector_data = build_sector_matrices_cached(
        p,
        tM1,
        tM2;
        calibrate_C4 = calibrate_C4,
        c4_imag_tol = c4_imag_tol,
        c4_method = c4_method,
        parallel_C4 = parallel_C4,
        use_cache = use_cache,
        verbose = verbose
    )

    H  = sector_data.H
    H4 = sector_data.H4
    HL = sector_data.HL
    HR = sector_data.HR
    basis = sector_data.basis
    product_states = sector_data.product_states
    c4_info = sector_data.c4_info

    if verbose
        println("==============================================")
        println("SO(5) reduced CG by C2+C4 composite diagonalization")
        println("Optimizations active: direct-sector C4 columns + upper-L cache + eigensystem cache + CG cache")
        println("C4 method: $c4_method, parallel_C4=$(parallel_C4 && Threads.nthreads()>1)")
        println("Single-particle irrep: (p,0) = ($p,0)")
        println("Target SO(5):          (P,Q) = ($P,$Q)")
        println("Target SO(4):          J=$(qval(tJ)), M1=$(qval(tM1)); K=$(qval(tK)), M2=$(qval(tM2))")
        println("Single-particle dim:   ", length(basis))
        println("Product sector dim:    ", length(product_states))
        println("epsC4 = $epsC4, epsJ = $epsJ, epsK = $epsK")
        println("Residual tolerances: tolC=$tolC, tolC4=$tolC4, tolJ=$tolJ, tolK=$tolK")
        println("Target C2(P,Q) = ", so5_casimir(P,Q))
        println("Target C4(P,Q) = ", so5_casimir4(P,Q))
        println("Sector cache key = ", sector_data.cache_key)
        println("==============================================")
    end

    if isempty(product_states)
        return (
            reduced_keyed = Dict{NTuple{9,Int},Float64}(),
            reduced_pretty = Dict{String,Float64}(),
            selected_vectors = zeros(Float64,0,0),
            diagnostics = Vector{Dict{String,Float64}}(),
            commutators = nothing,
            c4_info = c4_info,
            H = H,
            H4 = H4,
            HL = HL,
            HR = HR,
            O = zeros(Float64,0,0),
            basis = basis,
            product_states = product_states,
            target_so5 = (P,Q),
            target_so4_doubled = (tJ,tM1,tK,tM2),
            target_C4 = so5_casimir4(P,Q),
            epsC4 = epsC4,
            epsJ = epsJ,
            epsK = epsK,
            cache_key = sector_data.cache_key,
            c4_method = c4_method,
            parallel_C4 = parallel_C4,
            cache_eigensystem = cache_eigensystem,
            use_cache = use_cache,
            selection_rule = "residual_C2_C4_J_K"
        )
    end

    comms = nothing

    if check_commutators
        comms = commutator_diagnostics_cached(
            sector_data.cache_key,
            H,HL,HR;
            H4=H4,
            use_cache=(use_cache && cache_commutators),
            verbose=verbose
        )
    end

    eig_data = diagonalize_composite_cached(
        sector_data.cache_key,
        H,H4,HL,HR;
        epsC4=epsC4,
        epsJ=epsJ,
        epsK=epsK,
        use_cache=(use_cache && cache_eigensystem),
        verbose=verbose
    )
    eigO = eig_data.eigO
    O = eig_data.O

    selected_vectors, diagnostics = select_target_vectors_from_composite(
        eigO,H,H4,HL,HR,P,Q,tJ,tK;
        epsC4=epsC4,
        epsJ=epsJ,
        epsK=epsK,
        tolC=tolC,
        tolC4=tolC4,
        tolJ=tolJ,
        tolK=tolK,
        prefilter=prefilter_candidates
    )

    selected_vectors, diagnostics = remove_linearly_dependent_vectors_with_diagnostics(
        selected_vectors,
        diagnostics;
        tol = independence_tol
    )

    if verbose
        println("Selected target-vector number: ", size(selected_vectors,2))
        println("Target C2(P,Q) = ", so5_casimir(P,Q))
        println("Target C4(P,Q) = ", so5_casimir4(P,Q))
        println("Target J(J+1)  = ", qval(tJ)*(qval(tJ)+1.0))
        println("Target K(K+1)  = ", qval(tK)*(qval(tK)+1.0))

        for (i,diag) in enumerate(diagnostics)
            println("---- selected state $i ----")
            println("O eigenvalue = ", diag["O_eigenvalue"])
            println("<C2_SO5>     = ", diag["C2_SO5"])
            println("<C4_SO5>     = ", diag["C4_SO5"])
            println("<J^2>        = ", diag["J2"], "   J = ", diag["J_measured"])
            println("<K^2>        = ", diag["K2"], "   K = ", diag["K_measured"])
            println("resC         = ", diag["resC"])
            println("resC4        = ", diag["resC4"])
            println("resJ         = ", diag["resJ"])
            println("resK         = ", diag["resK"])
        end
    end

    reduced_keyed, reduced_pretty = extract_reduced_coefficients(
        selected_vectors,
        product_states,
        basis,
        cache,
        P,Q,tJ,tK,tM1,tM2;
        cg_tol=cg_tol
    )

    return (
        reduced_keyed = reduced_keyed,
        reduced_pretty = reduced_pretty,
        selected_vectors = selected_vectors,
        diagnostics = diagnostics,
        commutators = comms,
        c4_info = c4_info,
        H = H,
        H4 = H4,
        HL = HL,
        HR = HR,
        O = O,
        basis = basis,
        product_states = product_states,
        target_so5 = (P,Q),
        target_so4_doubled = (tJ,tM1,tK,tM2),
        target_C4 = so5_casimir4(P,Q),
        epsC4 = epsC4,
        epsJ = epsJ,
        epsK = epsK,
        cache_key = sector_data.cache_key,
        c4_method = c4_method,
        parallel_C4 = parallel_C4,
        cache_eigensystem = cache_eigensystem,
        use_cache = use_cache,
        selection_rule = "residual_C2_C4_J_K"
    )
end

# ------------------------------------------------------------
# 15. Convenience wrapper and diagnostics
# ------------------------------------------------------------

function SO5CG_composite(p::Int,
                         P::Int,
                         Q::Int,
                         J,
                         M1,
                         K,
                         M2;
                         kwargs...)

    result = compute_reduced_so5_cg_composite(
        p,P,Q,J,M1,K,M2;
        kwargs...
    )

    return result.reduced_pretty, result
end

function print_reduced(result; digits::Int=12)
    keys_sorted = sort(collect(keys(result.reduced_pretty)))

    for k in keys_sorted
        val = result.reduced_pretty[k]
        println(rpad(k, 95), " => ", round(val; digits=digits))
    end
end

function summarize_spectrum(A::Matrix{Float64}; tol::Float64=1e-6)
    vals = eigen(Symmetric(A)).values

    used = falses(length(vals))
    groups = Vector{Tuple{Float64,Int}}()

    for i in eachindex(vals)
        used[i] && continue

        close = findall(j -> !used[j] && abs(vals[j]-vals[i]) < tol,
                        eachindex(vals))

        for j in close
            used[j] = true
        end

        push!(groups, (mean(vals[close]), length(close)))
    end

    sort!(groups, by=x -> x[1])

    return groups
end

function check_basis_dimension(p::Int)
    basis, label = generate_so4_basis_p0(p)
    return length(basis) == so5_p0_dimension(p)
end

function check_C2_C4_collision_example()
    println("C2(6,0) = ", so5_casimir(6,0))
    println("C2(5,3) = ", so5_casimir(5,3))
    println("C4(6,0) = ", so5_casimir4(6,0))
    println("C4(5,3) = ", so5_casimir4(5,3))
end


"""
    verify_fast_C4_against_legacy(p,tM1,tM2; ...)

建议只对小 p 使用。比较新 fixed-sector columns 算法与原 full-space 算法，
返回相对 Frobenius 误差和最大元素误差。
"""
function verify_fast_C4_against_legacy(
    p::Int,
    tM1::Int,
    tM2::Int;
    calibrate_C4::Bool=true,
    c4_imag_tol::Float64=1e-8,
    verbose::Bool=true
)
    product_states,basis,label = generate_product_sector_fixed_M(p,tM1,tM2)

    H4_fast,_ = build_C4_sector_matrix(
        p,product_states,basis,label;
        calibrate=calibrate_C4,
        imag_tol=c4_imag_tol,
        use_cache=false,
        c4_method=:sector_columns,
        parallel_C4=false,
        verbose=verbose
    )

    H4_legacy,_ = build_C4_sector_matrix(
        p,product_states,basis,label;
        calibrate=calibrate_C4,
        imag_tol=c4_imag_tol,
        use_cache=false,
        c4_method=:full_legacy,
        parallel_C4=false,
        verbose=verbose
    )

    Delta = H4_fast-H4_legacy
    return (
        relative_error=norm(Delta)/max(norm(H4_legacy),1e-14),
        max_abs_error=maximum(abs,Delta),
        sector_dimension=length(product_states)
    )
end

"""
    recommended_SO5_runtime_setup!(; blas_threads=...)

设置 dense eigen 使用的 BLAS 线程数。请用 JULIA_NUM_THREADS 启动 Julia，
才能让 parallel_C4=true 的五个 W_a 任务并行。
"""
function recommended_SO5_runtime_setup!(;
    blas_threads::Int=max(1,min(Sys.CPU_THREADS,8))
)
    BLAS.set_num_threads(blas_threads)
    println("Julia threads = ",Threads.nthreads())
    println("BLAS threads  = ",BLAS.get_num_threads())
    return nothing
end

println("SO(5) core loaded: v5 fast fixed-sector C4 + cached composite eigensystem.")


