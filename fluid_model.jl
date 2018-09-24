
# [slack] <chrisrackauckas> p is for parameters. You can pass anything through there.
# [slack] <chrisrackauckas> or you can use whatever Julia method you like: closures, call-overloaded types, etc. to encapsulate data locally.
# [slack] <chrisrackauckas> And yes globals can be used if you wish, though they should be const if you want full type-stability and performance.



# const Kmin = 100
# const Kmax = 200
# const pmax = 0.2


function ecn_mark( q::SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{UInt64}},true},
                    N::UInt)
    # nomark = prod(ones(N,M) - matrix .* mark, dims=(2))
    # flow_queue_mark = similar(matrix, Float64)

    nomark = ones(Float64, N)
    mark = 0.0
    for i = 1:length(q)
        if (q[i] <= Kmin[i])
            mark = 0.0
        elseif (q[i] > Kmax[i])
            mark = 1.0
        else
            mark = pmax[i] * (q[i] - Kmin[i]) / (Kmax[i] - Kmin[i]);
        end
        nomark .*= (ones(N) - MATRIX[:,i] .* mark)
    end

    return nomark
end

const EGRESS_RATE = Vector{Float64}(undef, NUMNIC)
function extrapolate!(lag::SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{UInt64}},true},
                        num_qp::Float64,
                        nic_rate_cap::Float64)

    # println(lag)
    @. lag *= num_qp
    # println(lag)
    EGRESS_RATE[:] = lag' * HOST
    # egress_rate[egress_rate.>nic_rate_cap] = nic_rate_cap/
    for i = 1:NUMNIC
        if EGRESS_RATE[i] > nic_rate_cap
            EGRESS_RATE[i] = (nic_rate_cap / EGRESS_RATE[i])
        else
            EGRESS_RATE[i] = 1.0
        end
    end
    # println(EGRESS_RATE)
    # testb = HOST * EGRESS_RATE
    # println(testb)
    lag .= lag .* (HOST * EGRESS_RATE)
    # println(lag)

    # overflow .= nic_rate_cap ./ egress_rate[egress_rate .> nic_rate_cap]
    # reduction = host[:, egress_rate > nic_rate_cap] * transpose(overflow[overflow < 1.0])
    # reduction[reduction.==0] .= 1
    # @. lag *= reduction
end



const out = zeros(3 * NUMFLOW + NUMLINK)
const global_t = zeros(1)

function fluid_model(du, u, h, p, t)
    # @unpack (timer, Rai, Rhai, g, t2, t1, t0, F, Rmin, qp_per_flow, nic_rate_cap, N, M, C, Kmin, Kmax, pmax) = p
    @unpack (timer, Rai, Rhai, g, t2, t1, t0, F, Rmin, qp_per_flow, nic_rate_cap, N, M) = p

    h(out, p, t-t2)

    # if (t / 1000 > global_t[1])
    #     global_t[1] += 1
    #     println(t)
    # end
    println(t)

    # println("newrun $qp_per_flow")
    # println(out)
    @views extrapolate!(out[1:N], qp_per_flow, nic_rate_cap)
    # println(out)
    # nomark = Vector{Float64}(undef, N)

    # queue = view(out, 3*N+1:3*N+M)
    # @views nomark = ecn_mark(out[3*N+1:3*N+M], Kmin, Kmax, pmax, N)
    @views nomark = ecn_mark(out[3*N+1:3*N+M], N)


    for i = 1:N
        if MASK[i] == 0
            continue
        end

        if nomark[i]==1
            du[i] = (u[N+i]-u[i])/2/timer
            du[N+i] = Rai/timer
            du[2*N+i] = -g/t1 * u[2*N+i]

        elseif nomark[i]==0
            # min_dec
            # rgp_gd
            # TODO: how to incoporate min_dec and rpg_gd

            # du[i] = max(-u[i] * u[2*N+i]/2/t0, -0.5 * u[i]/t0)
            du[i] = -u[i] * u[2*N+i]/2/t0
            du[N+i] = (u[i]-u[N+i]) / t0
            du[2*N+i] = g/t1 * (1 - u[2*N+i])
        else
            a = 1-nomark[i]^(t0*out[i])
            c = 1-nomark[i]^(t1*out[i])
            d = (1-nomark[i]) / (nomark[i]^(-timer*out[i])-1)
            e = d * (nomark[i]^((F * timer * out[i])))
            f = d * (nomark[i]^((2 * F * timer * out[i])))

            # du[i] = max(-u[i] * u[2*N+i] / 2 / t0 * a, -0.5 * u[i] / t0) +(u[N+i]-u[i]) * out[i] * d / 2;
            du[i] = -u[i] * u[2*N+i] / 2 / t0 * a +(u[N+i]-u[i]) * out[i] * d / 2;
            du[N+i] = (u[i]-u[N+i]) * a / t0 + Rai * out[i] * e + Rhai * out[i] * f;
            du[2*N+i] = g/t1 * (c - u[2*N+i]);
        end

        if u[i] <= Rmin
            du[i] = max(du[i], 0)
        end
    end

    # println(out)
    traf_sum = out[1:N]' * MATRIX;

    for i=1:M
        if u[3*N+i] > 0
            du[3*N+i] = (traf_sum[i] - C[i])
        else
            du[3*N+i] = max(traf_sum[i] - C[i], 0)
        end
    end
end
