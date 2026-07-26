FROM julia:1.12.6-bookworm

WORKDIR /app

# 先复制依赖清单，让 Docker 可以复用包安装层。
COPY Project.toml /app/Project.toml
RUN julia --project=/app -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

COPY SO5_CG_v5_fast.jl /app/SO5_CG_v5_fast.jl
COPY server.jl /app/server.jl

ENV HOST=0.0.0.0 \
    PORT=8080 \
    JULIA_NUM_THREADS=auto \
    SO5_MAX_P=10 \
    SO5_MAX_ROWS=2000 \
    SO5_WARMUP=1

EXPOSE 8080

CMD ["julia", "--project=/app", "/app/server.jl"]
