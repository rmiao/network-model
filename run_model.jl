function get_real_rate!(y_dcqcn, y_real, y, dcqcn_param::DCQCN_Param)
    @unpack (nic_rate_cap) = dcqcn_param

    egress_t = y' * HOST

    for i = 1:NUMNIC
        if egress_t[i] > nic_rate_cap
            egress_t[i] = (nic_rate_cap / egress_t[i])
        else
            egress_t[i] = 1.0
        end
    end

    y_real .= (HOST * egress_t') .* y
    y_dcqcn .= y

    # println(y_real)
    # println(y_dcqcn)

end

function get_cnp_rate!(cnp, rate, delta, queue, dcqcn_param::DCQCN_Param)
    @unpack (N) = dcqcn_param

    nomark = ecn_mark(queue, N)
    # println(nomark)
    # println(pkts)
    cnp .+= (ones(N) - nomark) .* rate .* delta
end

function get_queue_size!(queue_avg, queue)
    queue_avg .= queue
end

function get_pause!()

end


function averaging(chunk, t, u, dcqcn_param::DCQCN_Param)
    @unpack (nic_rate_cap, qp_per_flow, M, N) = dcqcn_param
    # extrapolate

    # println(size(u,1))
    # println(u[1])
    # println(u[1])
    for i=1:size(u,1)
        u[i][1:N] = u[i][1:N] * qp_per_flow
    end
    # println(qp_per_flow)
    # println(u[1])

    # println(u[1])

    # u[:][1:N] = u[:][1:N] .* qp_per_flow

    # println("unreached here")

    t_prev = t[1]

    time = zeros(Float64, 0)
    y_real = zeros(Float64, N, 0)
    y = zeros(Float64, N, 0)
    cnp = zeros(Float64, N, 0)
    queue = zeros(Float64, M, 0)
    # pause_count = zeros(Float64, N, 0)
    # pause_duration = zeros(Float64, N, 0)

    inds = 1.0
    delta = 0.0


    temp_y_real = zeros(Float64, N)
    temp_y = zeros(Float64, N)
    temp_cnp = zeros(Float64, N)
    temp_queue = zeros(Float64, M)
    # temp_pause_count = zeros(Float64, N)
    # temp_pause_duration = zeros(Float64, N)

    index_count = 0
    accu_count = 0.0
    for i = 2:length(t)
        delta = t[i] - t_prev
        t_prev = t[i]
        accu_count += delta

        # 1. Calcuate the real flow rate

        @views get_real_rate!(temp_y, temp_y_real, u[i-1][1:N], dcqcn_param)

        # 2. Calcuate the queue size
        @views get_queue_size!(temp_queue, u[i-1][3*N+1:3*N+M])

        # 3. Calcuate the CNP rate

        @views get_cnp_rate!(temp_cnp, temp_y_real, delta, u[i-1][3*N+1:3*N+M], dcqcn_param)


        index_count += 1

        if (floor(t[i] / chunk) >= inds)
            # println("new line: $(inds*chunk) $accu_count")
            # println(temp_cnp)
            append!(time, inds * chunk)
            y_real = hcat(y_real, temp_y_real)

            # println(temp_y_real / chunk)

            # println((temp_y_real / accu_count)' * HOST)
            # println((temp_y_real / accu_count)' * EGRESS)

            y = hcat(y, temp_y)
            cnp = hcat(cnp, temp_cnp)
            queue = hcat(queue, temp_queue)

            fill!(temp_y_real, 0.0)
            fill!(temp_y, 0.0)
            fill!(temp_cnp, 0.0)
            fill!(temp_queue, 0.0)

            inds += 1.0
            index_count = 0
            accu_count = 0.0
        end
    end


    return time, y, y_real, queue, cnp

end

function run_model(env_param::Env_Param, dcqcn_param::DCQCN_Param, sim_len)
    #
    # incidents = importdata(PROJ_PATH + "incidents.config");
    # incident_index = 1;
    # matrix = importdata(PROJ_PATH + sprintf("switch_%d.config", incidents(incident_index)));
    # host = importdata(PROJ_PATH + sprintf("host_%d.config", incidents(incident_index)));
    # mask = importdata(PROJ_PATH + sprintf("mask_%d.config", incidents(incident_index)));

    # @unpack (nic_rate) = env_param
    # @unpack (alpha_init, t2, nic_rate_cap, N, M, Kmin, Kmax, pmax) = dcqcn_param
    @unpack (alpha_init, t2, nic_rate_cap, N, M) = dcqcn_param

    # terminal = [ones(M,1); zeros(M,1)];
    # shared_buffer=importdata(proj_folder + "device_resource.config")./1024;
    # device_port=importdata(proj_folder + "device_to_port.config");
    # threshold_offset = zeros(M,1);
    # under_pfc = false(N,1);
    # ingress_offset = 12 * 256 / 1024; % 12 cells threshold_offset, 3KB

    # rc0 =  ones(N) .* LINERATE
    rc0 = zeros(Float64, N)
    rt0 = copy(rc0)
    q0 =  zeros(Float64, M)
    alpha0 =  ones(Float64, N) .* alpha_init


    rc1=rc0
    rt1=rt0
    alpha1=alpha0
    q1=q0


    start1 = vcat(rc1, rt1, alpha1, q1)
    # start1 = rand(N)
    fill!(MASK, false)

    # d = (EGRESS[:,126].==1)
    # @views start1[1:NUMFLOW][d] .= LINERATE
    # @views start1[NUMFLOW+1:2NUMFLOW][d] .= LINERATE
    # @views MASK[d] .= 1


    lags = [t2]

    h(out, p, t) = (out.= start1)
    # h(p, t) = ones(3*N+M)

    tspan = (0.0, sim_len)
    prob = DDEProblem(fluid_model, start1, h, tspan, dcqcn_param; constant_lags=lags)

    # alg = MethodOfSteps(Tsit5())
    alg = MethodOfSteps(BS3())
    # alg = MethodOfSteps(Rosenbrock23(autodiff=false))
    # alg = MethodOfSteps(Rodas4(autodiff=false))

    # @time sol = solve(prob, alg, reltol=1e-6, )
    cb1 = ContinuousCallback(pfc, nothing, pause!)
    cb2 = ContinuousCallback(unpfc, nothing, unpause!)
    cb3 = ContinuousCallback(arrival, nothing, inject!)
    cb4 = ContinuousCallback(departure, nothing, remove!)
    # cb5 = ContinuousCallback(firstecn, nothing, rate_at_first_ecn!)

    cbset = CallbackSet(cb1, cb2, cb3, cb4)


    @time sol = solve(prob, alg, callback=cbset)

    # @time sol = solve(prob, alg)

    println(sol.retcode)

    (t, y, y_real, queue, cnp) = averaging(CHUNK, sol.t, sol.u, dcqcn_param)



    pause_duration = PAUSE_DURATION'
    pause_num = PAUSE_COUNT'

    return (t, y, y_real, queue, cnp, pause_duration, pause_num)
end
