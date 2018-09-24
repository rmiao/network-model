const BUFFER_DYNAMIC = 1.0/16
const STATE = zeros(Float64, 2 * NUMFLOW)
const EVENTMASK = trues(NUMLINK)
const PFCRMIN = 1*1024/8.0/1000000
const OFFSET = Array{Float64}(undef, NUMLINK)
const INGRESS_OFFSET = 12.0 * 256 / 1024
const DUMMY_EVENT = 10.0
const PAUSE_COUNT = zeros(UInt,  NUMFLOW, floor(Int, SIM_LEN * 1000000/CHUNK))
const PAUSE_DURATION = zeros(Float64, NUMFLOW, floor(Int, SIM_LEN * 1000000/CHUNK))
const PAUSE_TIME = Array{Float64}(undef, NUMLINK)
# const FIRST_ECN_ENABLE = trues(NUMLINK)


function arrival(u, t, integrator)
    minimum(vcat((ARRIVAL_TIME .- t)[NOT_ARRIVED], DUMMY_EVENT))
end

function inject!(integrator)
    u = integrator.u

    time_to_arrival = ARRIVAL_TIME .- integrator.t
    minvalue = minimum(time_to_arrival[NOT_ARRIVED])

    host = 0
    for i = 1:NUMHOST
        if(time_to_arrival[i] == minvalue && NOT_ARRIVED[i])
            host = i
            break
        end
    end

    d = (EGRESS[:,host].==1)

    @views u[1:NUMFLOW][d] .= LINERATE
    @views u[NUMFLOW+1:2NUMFLOW][d] .= LINERATE
    @views MASK[d] .= 1

    # println("start host: $host")
    # println(u)
    NOT_ARRIVED[host] = false

end

function departure(u, t, integrator)
    minimum(vcat(((SIM_LEN * 1000000 - t) .- DEPARTURE_TIME)[.~NOT_ARRIVED], DUMMY_EVENT))
end

function remove!(integrator)
    u = integrator.u

    time_to_departure = (SIM_LEN * 1000000 - integrator.t) .- DEPARTURE_TIME
    minvalue = minimum(time_to_departure[.~NOT_ARRIVED])

    host = 0
    for i = 1:NUMHOST
        if(time_to_departure[i] == minvalue && ~NOT_ARRIVED[i])
            host = i
            break
        end
    end

    d = (EGRESS[:,host].==1)

    @views u[1:NUMFLOW][d] .= 0.0
    @views u[NUMFLOW+1:2NUMFLOW][d] .= 0.0
    @views MASK[d] .= 0

    NOT_ARRIVED[host] = true

end

function pfc(u, t, integrator)
    # println(u[end-NUMLINK+1:end])
    # println((DEVICE_PORT' * (SHARED_BUFFER - DEVICE_PORT * u[end-NUMLINK+1:end]) *BUFFER_DYNAMIC - u[end-NUMLINK+1:end])[EVENTMASK])
    minimum(vcat((DEVICE_PORT' * (SHARED_BUFFER - DEVICE_PORT * u[end-NUMLINK+1:end]) *BUFFER_DYNAMIC - u[end-NUMLINK+1:end])[EVENTMASK], DUMMY_EVENT))
    # println(minimum(vcat((DEVICE_PORT' * (SHARED_BUFFER - DEVICE_PORT * u[end-NUMLINK+1:end]) *BUFFER_DYNAMIC - u[end-NUMLINK+1:end])[EVENTMASK], DUMMY_EVENT)))
end


function pause!(integrator)
    # println("Pause")

    u = integrator.u

    buffer_leftover = DEVICE_PORT' * (SHARED_BUFFER - DEVICE_PORT * u[end-NUMLINK+1:end]) * BUFFER_DYNAMIC - u[end-NUMLINK+1:end]
    minvalue = minimum(buffer_leftover[EVENTMASK])

    pfcqueue = 0
    for i = 1:NUMLINK
        if(buffer_leftover[i] == minvalue && EVENTMASK[i])
            pfcqueue = i
            break
        end
    end

    d = (MATRIX[:,pfcqueue].==1)

    @views STATE[1:NUMFLOW][d] = u[1:NUMFLOW][d]
    @views STATE[NUMFLOW+1:2NUMFLOW][d] = u[NUMFLOW+1:2NUMFLOW][d]

    @views u[1:NUMFLOW][d] .= 0.0
    @views u[NUMFLOW+1:2NUMFLOW][d] .= 0.0


    # we count the pause separately per each queue
    # since we assume the pauses happen sparsely in the network

    # count pfc pause
    @. @views PAUSE_COUNT[:,floor(Int, integrator.t/CHUNK)][d] += 1
    PAUSE_TIME[pfcqueue] = integrator.t

    OFFSET[pfcqueue] = u[3*NUMFLOW+pfcqueue] - INGRESS_OFFSET
    EVENTMASK[pfcqueue] = false
end


function unpfc(u, t, integrator)
    minimum(vcat((u[end-NUMLINK+1:end] - OFFSET)[.~EVENTMASK], DUMMY_EVENT))
end

function unpause!(integrator)
    # println("Unpause")

    u = integrator.u
    buffer_offset = u[end-NUMLINK+1:end] - OFFSET
    minvalue = minimum(buffer_offset[.~EVENTMASK])

    pfcqueue = 0
    for i = 1:NUMLINK
        if(buffer_offset[i] == minvalue && ~EVENTMASK[i])
            pfcqueue = i
            break
        end
    end

    d = (MATRIX[:, pfcqueue].==1)

    @. @views PAUSE_DURATION[:, floor(Int, integrator.t/CHUNK)][d] += integrator.t - PAUSE_TIME[pfcqueue]

    @views u[1:NUMFLOW][d] = STATE[1:NUMFLOW][d]
    @views u[NUMFLOW+1:2NUMFLOW][d] = STATE[NUMFLOW+1:2NUMFLOW][d]

    EVENTMASK[pfcqueue] = true
end




# depreciated
function firstecn(u, t, integerator)
    minimum(vcat((Kmin - u[end-NUMLINK+1:end])[FIRST_ECN_ENABLE], DUMMY_EVENT))
end

function rate_at_first_ecn!(integrator)
    u = integrator.u
    # we have to loop the links twice, any better idea?
    buffer_leftover = Kmin - u[end-NUMLINK+1:end]
    minvalue = minimum(buffer_leftover[FIRST_ECN_ENABLE])

    firstecnqueue = 0
    for i = 1:NUMLINK
        if(buffer_leftover[i] == minvalue && FIRST_ECN_ENABLE[i])
            firstecnqueue = i
            break
        end
    end

    # println("firstecn on link $firstecnqueue")

    d = (MATRIX[:,firstecnqueue].==1)

    @views u[1:NUMFLOW][d] .= PFCRMIN
    @views u[NUMFLOW+1:2NUMFLOW][d] .= PFCRMIN
    FIRST_ECN_ENABLE[firstecnqueue] = false
end
