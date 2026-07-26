using HTTP
using JSON3
using Dates
using LinearAlgebra

include(joinpath(@__DIR__, "SO5_CG_v5_fast.jl"))

# ============================================================
# SO(5) CG Web API
#
# Endpoints:
#   GET  /health   健康检查
#   POST /api/cg   计算约化 SO(5) Clebsch--Gordon 系数
#
# 说明：
# - 核心计算缓存会常驻进程，因此同一 (p,M1,M2) 扇区的后续请求会更快。
# - 核心缓存中的 Dict 不是并发写安全的，所以使用全局锁串行化重计算请求。
# - 对公开服务限制最大 p，避免单个请求耗尽服务器内存。
# ============================================================

const STARTED_AT = now(UTC)
const COMPUTE_LOCK = ReentrantLock()
const MAX_P = parse(Int, get(ENV, "SO5_MAX_P", "10"))
const MAX_ROWS = parse(Int, get(ENV, "SO5_MAX_ROWS", "2000"))
const MAX_BODY_BYTES = parse(Int, get(ENV, "SO5_MAX_BODY_BYTES", "16384"))
const API_VERSION = "1.0.0"

# 使用逗号分隔，例如：
#   SO5_ALLOWED_ORIGINS=https://maggiexheuw.github.io,http://localhost:4000
# 默认 * 便于初次部署；正式上线建议设置为你的博客域名。
const ALLOWED_ORIGINS_RAW = strip(get(ENV, "SO5_ALLOWED_ORIGINS", "*"))
const ALLOWED_ORIGINS = ALLOWED_ORIGINS_RAW == "*" ? String["*"] :
    filter(!isempty, strip.(split(ALLOWED_ORIGINS_RAW, ',')))

function request_origin(req::HTTP.Request)
    return get(req.headers, "Origin", "")
end

function cors_origin(req::HTTP.Request)
    origin = request_origin(req)
    if "*" in ALLOWED_ORIGINS
        return "*"
    end
    return origin in ALLOWED_ORIGINS ? origin : ""
end

function response_headers(req::HTTP.Request; json::Bool=true)
    headers = Pair{String,String}[
        "Cache-Control" => "no-store",
        "X-Content-Type-Options" => "nosniff",
        "Vary" => "Origin",
    ]

    if json
        push!(headers, "Content-Type" => "application/json; charset=utf-8")
    end

    allowed = cors_origin(req)
    if !isempty(allowed)
        push!(headers, "Access-Control-Allow-Origin" => allowed)
        push!(headers, "Access-Control-Allow-Methods" => "GET, POST, OPTIONS")
        push!(headers, "Access-Control-Allow-Headers" => "Content-Type")
        push!(headers, "Access-Control-Max-Age" => "86400")
    end

    return headers
end

function json_response(req::HTTP.Request, status::Integer, payload)
    return HTTP.Response(
        Int(status);
        headers = response_headers(req),
        body = JSON3.write(payload),
    )
end

function empty_response(req::HTTP.Request, status::Integer=204)
    return HTTP.Response(
        Int(status);
        headers = response_headers(req; json=false),
        body = "",
    )
end

# JSON3.Object 在不同版本中允许 Symbol/String 键访问。这个辅助函数兼容两者。
function json_get(obj, name::AbstractString, default)
    skey = Symbol(name)
    if haskey(obj, skey)
        return obj[skey]
    elseif haskey(obj, name)
        return obj[name]
    end
    return default
end

function as_int(x, name::AbstractString)
    if x isa Integer
        return Int(x)
    elseif x isa Real && isfinite(Float64(x)) && isinteger(Float64(x))
        return Int(round(Float64(x)))
    end
    throw(ArgumentError("$name 必须是整数。"))
end

function as_float(x, name::AbstractString)
    x isa Real || throw(ArgumentError("$name 必须是数值。"))
    y = Float64(x)
    isfinite(y) || throw(ArgumentError("$name 必须是有限数值。"))
    return y
end

function as_bool(x, name::AbstractString)
    x isa Bool || throw(ArgumentError("$name 必须是布尔值。"))
    return x
end

is_half_integer(x::Float64; tol::Float64=1e-10) = abs(2x - round(2x)) <= tol

function bounded_positive(x::Float64, name::AbstractString; max_value::Float64=1.0)
    0.0 < x <= max_value || throw(ArgumentError("$name 必须满足 0 < $name <= $max_value。"))
    return x
end

function parse_request(obj)
    p = as_int(json_get(obj, "p", 10), "p")
    P = as_int(json_get(obj, "P", p), "P")
    Q = as_int(json_get(obj, "Q", 0), "Q")

    J  = as_float(json_get(obj, "J", 0.0), "J")
    M1 = as_float(json_get(obj, "M1", 0.0), "M1")
    K  = as_float(json_get(obj, "K", 0.0), "K")
    M2 = as_float(json_get(obj, "M2", 0.0), "M2")

    epsC4 = bounded_positive(as_float(json_get(obj, "epsC4", 1e-5), "epsC4"), "epsC4")
    epsJ  = bounded_positive(as_float(json_get(obj, "epsJ", 1e-4), "epsJ"), "epsJ")
    epsK  = bounded_positive(as_float(json_get(obj, "epsK", 1e-6), "epsK"), "epsK")

    tolC  = bounded_positive(as_float(json_get(obj, "tolC", 1e-6), "tolC"), "tolC")
    tolC4 = bounded_positive(as_float(json_get(obj, "tolC4", 1e-6), "tolC4"), "tolC4")
    tolJ  = bounded_positive(as_float(json_get(obj, "tolJ", 1e-6), "tolJ"), "tolJ")
    tolK  = bounded_positive(as_float(json_get(obj, "tolK", 1e-6), "tolK"), "tolK")
    cg_tol = bounded_positive(as_float(json_get(obj, "cg_tol", 1e-12), "cg_tol"), "cg_tol")

    calibrate_C4 = as_bool(json_get(obj, "calibrate_C4", true), "calibrate_C4")
    show_diagnostics = as_bool(json_get(obj, "show_diagnostics", true), "show_diagnostics")
    recognize_radicals = as_bool(json_get(obj, "recognize_radicals", true), "recognize_radicals")

    max_rows = as_int(json_get(obj, "max_rows", 200), "max_rows")
    digits = as_int(json_get(obj, "digits", 10), "digits")

    0 <= p <= MAX_P || throw(ArgumentError("公开服务要求 0 <= p <= $MAX_P。"))
    0 <= Q <= P <= 2p || throw(ArgumentError("需要满足 0 <= Q <= P <= 2p。"))

    for (name, value) in (("J", J), ("M1", M1), ("K", K), ("M2", M2))
        is_half_integer(value) || throw(ArgumentError("$name 必须是整数或半整数。"))
    end

    0 <= J <= p || throw(ArgumentError("需要满足 0 <= J <= p。"))
    0 <= K <= p || throw(ArgumentError("需要满足 0 <= K <= p。"))
    abs(M1) <= J || throw(ArgumentError("需要满足 |M1| <= J。"))
    abs(M2) <= K || throw(ArgumentError("需要满足 |M2| <= K。"))

    1 <= max_rows <= MAX_ROWS || throw(ArgumentError("max_rows 必须在 1 到 $MAX_ROWS 之间。"))
    4 <= digits <= 16 || throw(ArgumentError("digits 必须在 4 到 16 之间。"))

    return (
        p=p, P=P, Q=Q, J=J, M1=M1, K=K, M2=M2,
        epsC4=epsC4, epsJ=epsJ, epsK=epsK,
        tolC=tolC, tolC4=tolC4, tolJ=tolJ, tolK=tolK,
        cg_tol=cg_tol, calibrate_C4=calibrate_C4,
        show_diagnostics=show_diagnostics,
        recognize_radicals=recognize_radicals,
        max_rows=max_rows, digits=digits,
    )
end

function half_label(t::Int)
    return iseven(t) ? string(div(t, 2)) : string(t, "/2")
end

# 将 sqrt(n/d) 写成 a*sqrt(r)/b。CG 系数常见为有理数平方根，
# 这个识别算法只对 x^2 做一次 rationalize，比旧界面逐个扫描根式快得多。
function square_split(n::Int)
    n >= 0 || throw(ArgumentError("square_split requires n >= 0"))
    outside = 1
    inside = n
    d = 2
    while d*d <= inside
        dd = d*d
        while inside % dd == 0
            outside *= d
            inside = div(inside, dd)
        end
        d += 1
    end
    return outside, inside
end

function radical_label(x::Real; tol::Float64=1e-10, max_den::Int=100_000)
    y = Float64(x)
    abs(y) <= tol && return "0"

    r = rationalize(Int, y*y; tol=tol)
    n = numerator(r)
    d = denominator(r)

    d <= max_den || return nothing
    abs(Float64(r) - y*y) <= tol * max(1.0, y*y) || return nothing
    n >= 0 || return nothing

    out_n, in_n = square_split(n)
    out_d, in_d = square_split(d)

    # sqrt(n/d) = out_n * sqrt(in_n*in_d) / (out_d*in_d)
    a = out_n
    b = out_d * in_d
    rad = in_n * in_d

    g = gcd(a, b)
    a = div(a, g)
    b = div(b, g)

    sign = y < 0 ? "-" : ""

    if rad == 1
        return b == 1 ? string(sign, a) : string(sign, a, "/", b)
    end

    root = a == 1 ? "√$(rad)" : "$(a)√$(rad)"
    return b == 1 ? string(sign, root) : string(sign, root, "/", b)
end

function coefficient_payload(result, cfg)
    items = collect(result.reduced_keyed)
    sort!(items, by = first)

    total = length(items)
    nshow = min(total, cfg.max_rows)
    rows = Vector{NamedTuple}(undef, nshow)

    for idx in 1:nshow
        key, value = items[idx]
        channel, P0, Q0, tJ, tK, tj1, tk1, tj2, tk2 = key
        exact = cfg.recognize_radicals ? radical_label(value; tol=1e-10) : nothing

        rows[idx] = (
            channel=channel,
            P=P0,
            Q=Q0,
            J=0.5tJ,
            K=0.5tK,
            J_label=half_label(tJ),
            K_label=half_label(tK),
            j1=0.5tj1,
            k1=0.5tk1,
            j2=0.5tj2,
            k2=0.5tk2,
            j1_label=half_label(tj1),
            k1_label=half_label(tk1),
            j2_label=half_label(tj2),
            k2_label=half_label(tk2),
            coefficient=Float64(value),
            exact=exact,
        )
    end

    return rows, total
end

function calculate(cfg)
    t0 = time_ns()

    result = compute_reduced_so5_cg_composite(
        cfg.p, cfg.P, cfg.Q,
        cfg.J, cfg.M1, cfg.K, cfg.M2;
        epsC4=cfg.epsC4,
        epsJ=cfg.epsJ,
        epsK=cfg.epsK,
        tolC=cfg.tolC,
        tolC4=cfg.tolC4,
        tolJ=cfg.tolJ,
        tolK=cfg.tolK,
        cg_tol=cfg.cg_tol,
        check_commutators=false,
        cache_commutators=true,
        cache_eigensystem=true,
        prefilter_candidates=true,
        calibrate_C4=cfg.calibrate_C4,
        c4_method=:sector_columns,
        parallel_C4=Threads.nthreads() > 1,
        use_cache=true,
        verbose=false,
    )

    elapsed = (time_ns() - t0) / 1.0e9
    coefficients, coefficient_count = coefficient_payload(result, cfg)

    diagnostics = cfg.show_diagnostics ? result.diagnostics : Any[]

    return (
        ok=true,
        api_version=API_VERSION,
        elapsed_seconds=elapsed,
        request=(
            p=cfg.p, P=cfg.P, Q=cfg.Q,
            J=cfg.J, M1=cfg.M1, K=cfg.K, M2=cfg.M2,
            epsC4=cfg.epsC4, epsJ=cfg.epsJ, epsK=cfg.epsK,
        ),
        targets=(
            C2=so5_casimir(cfg.P, cfg.Q),
            C4=so5_casimir4(cfg.P, cfg.Q),
            J2=cfg.J*(cfg.J+1.0),
            K2=cfg.K*(cfg.K+1.0),
        ),
        dimensions=(
            single_particle=length(result.basis),
            product_sector=length(result.product_states),
            selected_vectors=size(result.selected_vectors, 2),
            coefficient_count=coefficient_count,
            returned_coefficients=length(coefficients),
        ),
        coefficients=coefficients,
        diagnostics=diagnostics,
        c4_info=result.c4_info,
        cache=SO5_cache_status(),
    )
end

function handle_health(req::HTTP.Request)
    return json_response(req, 200, (
        ok=true,
        service="SO(5) reduced CG calculator",
        api_version=API_VERSION,
        started_at=string(STARTED_AT),
        julia_version=string(VERSION),
        julia_threads=Threads.nthreads(),
        blas_threads=BLAS.get_num_threads(),
        max_p=MAX_P,
        cache=SO5_cache_status(),
    ))
end

function handle_calculation(req::HTTP.Request)
    body_bytes = req.body
    length(body_bytes) <= MAX_BODY_BYTES ||
        return json_response(req, 413, (ok=false, error="请求体过大。"))

    obj = try
        JSON3.read(String(copy(body_bytes)))
    catch err
        return json_response(req, 400, (ok=false, error="JSON 格式错误。", detail=sprint(showerror, err)))
    end

    cfg = try
        parse_request(obj)
    catch err
        return json_response(req, 422, (ok=false, error=sprint(showerror, err)))
    end

    # 核心的全局缓存会在构造阶段写 Dict。为保证公开 API 稳定，计算请求串行执行。
    lock(COMPUTE_LOCK)
    try
        payload = calculate(cfg)
        return json_response(req, 200, payload)
    catch err
        bt = catch_backtrace()
        @error "SO(5) calculation failed" exception=(err, bt)
        return json_response(req, 500, (
            ok=false,
            error="计算失败。请检查量子数是否对应实际存在的分解通道。",
            detail=sprint(showerror, err),
        ))
    finally
        unlock(COMPUTE_LOCK)
    end
end

function app(req::HTTP.Request)
    method = uppercase(String(req.method))
    path = HTTP.URI(req.target).path

    if method == "OPTIONS"
        return empty_response(req)
    elseif method == "GET" && path == "/health"
        return handle_health(req)
    elseif method == "POST" && path == "/api/cg"
        return handle_calculation(req)
    elseif method == "GET" && path == "/"
        return json_response(req, 200, (
            ok=true,
            message="SO(5) CG API is running.",
            endpoints=["GET /health", "POST /api/cg"],
        ))
    else
        return json_response(req, 404, (ok=false, error="Not found"))
    end
end

function warmup!()
    get(ENV, "SO5_WARMUP", "1") == "1" || return
    @info "Starting a small Julia/JIT warm-up calculation"
    try
        compute_reduced_so5_cg_composite(
            1, 2, 0, 0, 0, 0, 0;
            epsC4=1e-5,
            epsJ=1e-4,
            epsK=1e-6,
            tolC=1e-6,
            tolC4=1e-6,
            tolJ=1e-6,
            tolK=1e-6,
            check_commutators=false,
            c4_method=:sector_columns,
            parallel_C4=false,
            use_cache=true,
            verbose=false,
        )
        @info "Warm-up finished"
    catch err
        @warn "Warm-up failed; the server will still start" exception=(err, catch_backtrace())
    end
end

function main()
    host = get(ENV, "HOST", "0.0.0.0")
    port = parse(Int, get(ENV, "PORT", "8080"))
    blas_threads = parse(Int, get(ENV, "SO5_BLAS_THREADS", string(max(1, min(Sys.CPU_THREADS, 8)))))
    BLAS.set_num_threads(blas_threads)

    @info "SO(5) API configuration" host port max_p=MAX_P julia_threads=Threads.nthreads() blas_threads=BLAS.get_num_threads()
    warmup!()

    # HTTP.serve 是阻塞式服务入口，适合容器进程。
    HTTP.serve(app, host, port; verbose=true, max_body_bytes=MAX_BODY_BYTES)
end

abspath(PROGRAM_FILE) == @__FILE__ && main()
