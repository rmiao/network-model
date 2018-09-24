# @time using DelimitedFiles
# include("ecn_mark.jl")
#
#
#
# N = 10000
# M = 1000
# # nomark = Vector{Float64}(undef, N)
# nomark = ones(Float64, N)
# # queue = rand(M)
# const out = rand(3*N+M)
# queue = view(out, 3*N+1:3*N+M)
#
# const matrix = rand(N, M)
#
# @time ecn_mark!(nomark, queue)

#
# using DelayDiffEq
# const p0 = 0.2; const q0 = 0.3; const v0 = 1; const d0 = 5
# const p1 = 0.2; const q1 = 0.3; const v1 = 1; const d1 = 1
# const d2 = 1; const beta0 = 1; const beta1 = 1; const tau = 1
# const out = zeros(3)
# function bc_model(du,u,h,p,t)
#   h(out, p, t)
#   # newout = ones(3)
#   du[1] = (v0/(1+beta0*(out[3]^2))) * (p0 - q0)*u[1] - d0*u[1]
#   du[2] = (v0/(1+beta0*(out[3]^2))) * (1 - p0 + q0)*u[1] +
#           (v1/(1+beta1*(out[3]^2))) * (p1 - q1)*u[2] - d1*u[2]
#   du[3] = (v1/(1+beta1*(out[3]^2))) * (1 - p1 + q1)*u[2] - d2*u[3]
# end
# lags = [tau]
# h(out, p, t) = (out.=1.0)
#
# tspan = (0.0,10.0)
# u0 = [1.0,1.0,1.0]
# prob = DDEProblem(bc_model,u0,h,tspan,constant_lags = lags)
# alg = MethodOfSteps(Tsit5())
#
# sol = solve(prob,alg)
# println(sol)
# using Plots; plot(sol)

@with_kw struct DCQCN_Param
    g::Float64
    t1::UInt
    alpha_init::Float64
    rate_first_cnp::Float64
    Rmin::Float64
    t0::UInt
    t2::UInt
    timer::UInt
    F::UInt
    Rai::Float64
    Rhai::Float64
    nic_rate_cap::Float64
    qp_per_flow::Float64
    N::UInt
    M::UInt
end
