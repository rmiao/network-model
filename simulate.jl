# @time using DifferentialEquations
@time using DelayDiffEq
@time using DelimitedFiles
@time using Parameters
@time import JSON
# @time using Plots

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
    # C::Vector{Float64}
    # Kmin::Vector{Float64}
    # Kmax::Vector{Float64}
    # pmax::Vector{Float64}
end



@with_kw struct Env_Param
    min_dec::Float64
    rgp_gd::Float64
    clamp_tgt_rate::UInt
    clamp_tgt_rate_after_time_inc::UInt
    byte::UInt
    # nic_rate::Float64
    num_of_qp::UInt
    pktsize::UInt
end

# const PROJ_PATH = "/tmp/dump/"
# const DATA_PATH = "/tmp/dump/data/"
const PROJ_PATH = "./data/"
const DATA_PATH = "./data/result/"
const SIM_LEN = 0.01
const CHUNK = 1000


# const INCIDENT = readdlm(PROJ_PATH * "incidents.config", ',', Int, use_mmap=true)
incident_index = 0
const MATRIX = readdlm(PROJ_PATH * "switch_$incident_index.config", ',', Int, use_mmap=true)
const HOST = readdlm(PROJ_PATH * "host_$incident_index.config", ',', Int, use_mmap=true)
const MASK = readdlm(PROJ_PATH * "mask_$incident_index.config", ',', Int, use_mmap=true)
const EGRESS = readdlm(PROJ_PATH * "egress.config", ',', Int, use_mmap=true)
const INGRESS = readdlm(PROJ_PATH * "ingress.config", ',', Int, use_mmap=true)
const LINK_STATUS = readdlm(PROJ_PATH * "link_status.config", ',', Int, use_mmap=true)
const DEVICE_PORT = readdlm(PROJ_PATH * "device_to_port.config", ',', Int, use_mmap=true)
const SHARED_BUFFER = readdlm(PROJ_PATH * "device_resource.config", ',', Int, use_mmap=true) ./ 1024

const (NUMFLOW, NUMLINK) = size(MATRIX)
const NUMNIC = size(HOST, 2)
const NUMHOST = size(EGRESS, 2)


const NOT_ARRIVED = trues(NUMHOST)
const ARRIVAL_TIME = readdlm(PROJ_PATH * "arrival.config", ',', Int, use_mmap=true)
const DEPARTURE_TIME = readdlm(PROJ_PATH * "departure.config", ',', Int, use_mmap=true)


const DCQCN_SETTINGS = JSON.parse(read(PROJ_PATH * "dcqcn_settings.json", String))
const LINERATE = DCQCN_SETTINGS["nic_rate"] * 1024 * 1024 / 8 / 1000000
const C = LINK_STATUS[:,1] * 1024 * 1024 / 8 / 1000000 # KB/us
const Kmin = LINK_STATUS[:,2] .* LINK_STATUS[:,3] / 1024  # KB
const Kmax = LINK_STATUS[:,2] .* LINK_STATUS[:,4] / 1024  # KB
const pmax = LINK_STATUS[:,5] ./100  # probability


# include("ecn_mark.jl")
include("fluid_model.jl")
include("run_model.jl")
include("pfc.jl")


function write_to_file(t, egress, dcqcn_egress, queue, cnp, pause_duration, pause_num)
    writedlm(DATA_PATH * "t.dat", t)
    writedlm(DATA_PATH * "egress.dat", egress)
    writedlm(DATA_PATH * "dcqcn_egress.dat", dcqcn_egress)
    writedlm(DATA_PATH * "queue.dat", queue)
    writedlm(DATA_PATH * "cnp.dat", cnp)
    writedlm(DATA_PATH * "pause_count.dat", pause_num)
    writedlm(DATA_PATH * "pause_duration.dat", pause_duration)
end


function aggregate(t, y, y_real, q, cnp, _pause_duration, _pause_num)
    time = t./1000000
    (numflow, numhost) = size(EGRESS)
    numtime = size(y, 2)
    egress_traffic = Array{Float64}(undef, numhost, numtime)
    dcqcn_egress_traffic = Array{Float64}(undef, numhost, numtime)

    for i = 1:numtime
        egress_traffic[:,i] = EGRESS' * y_real[:,i] .* 8
        dcqcn_egress_traffic[:,i] = EGRESS' * y[:,i] .* 8
    end
    # egress_traffic = [x'*EGRESS for x in y_real] .* 8
    # dcqcn_egress_traffic = [x'*EGRESS for x in y] .* 8

    cnp_sent = cnp' * INGRESS
    cnp_handled = cnp' * EGRESS

    queue = q

    pause_duration = _pause_duration * INGRESS
    pause_num = _pause_num * INGRESS

    return time, egress_traffic', dcqcn_egress_traffic', queue', cnp_sent, cnp_handled, pause_duration, pause_num
end




function simulate()
    dcqcn_settings = DCQCN_SETTINGS

    env_param = Env_Param(
        min_dec = dcqcn_settings["rpg_min_dec_fac"] / 100,
        rgp_gd = 1023 / 2^dcqcn_settings["rpg_gd"],
        clamp_tgt_rate = dcqcn_settings["clamp_tgt_rate"],  # default
        clamp_tgt_rate_after_time_inc = dcqcn_settings["clamp_tgt_rate_after_time_inc"], # default
        byte = dcqcn_settings["rpg_byte_reset"],  # not used
        # nic_rate = dcqcn_settings["nic_rate"] * 1024 * 1024 / 8 / 1000000,
        pktsize = dcqcn_settings["pkt_size"], # packet size (KB)
        num_of_qp = dcqcn_settings["num_of_qp"],  # num of qp per host
        # N = NUMFLOW,
        # M = NUMLINK,
        # nH = size(HOST, 2)
    )

    dcqcn_param = DCQCN_Param(
        g = (1024.0 - dcqcn_settings["dce_tcp_g"])/1024.0,
        t1 = dcqcn_settings["dce_tcp_rtt"],
        t2 = dcqcn_settings["rtt_estimate"], # control loop, us
        alpha_init = dcqcn_settings["initial_alpha_value"] / 1024.0,
        rate_first_cnp = dcqcn_settings["rate_to_set_on_first_cnp"] * 1024 / 8 / 1000000,  # Mbit/s
        Rmin = dcqcn_settings["rpg_min_rate"] * 1024 / 8 / 1000000,
        t0 = dcqcn_settings["rate_reduce_monitor_period"], # Rate update period upon CNP
        timer = dcqcn_settings["rpg_time_reset"],
        F = dcqcn_settings["rpg_threshold"],
        Rai = dcqcn_settings["rpg_ai_rate"] * 1024 / 8 / 1000000,
        Rhai = dcqcn_settings["rpg_hai_rate"] * 1024 / 8 / 1000000,  # KB/us = pkt/us
        nic_rate_cap = dcqcn_settings["nic_rate_cap"] * 1024 * 1024 / 8 / 1000000,
        # num_of_qp = dcqcn_settings["num_of_qp"],  # num of qp per host
        qp_per_flow = dcqcn_settings["num_of_qp"] / maximum(sum(EGRESS, dims=(1))),  # TODO: works for equal number of QP for each host only
        N = NUMFLOW,
        M = NUMLINK
        # C = LINK_STATUS[:,1] * 1024 * 1024 / 8 / 1000000, # KB/us
        # Kmin = LINK_STATUS[:,2] .* LINK_STATUS[:,3] / 1024,  # KB
        # Kmax = LINK_STATUS[:,2] .* LINK_STATUS[:,4] / 1024,  # KB
        # pmax = LINK_STATUS[:,5] ./100  # probability
    )


    # testfunc(dcqcn_param)
    (t, y, y_real, q, cnp, _pause_duration, _pause_num) = run_model(env_param, dcqcn_param, SIM_LEN * 1000000)
    (time, egress_traffic, dcqcn_egress_traffic, queue, cnp_sent, cnp_handled, pause_duration, pause_num) = aggregate(t, y, y_real, q, cnp, _pause_duration, _pause_num)
    write_to_file(time, egress_traffic, dcqcn_egress_traffic, queue, cnp_sent, pause_duration, pause_num)

end
