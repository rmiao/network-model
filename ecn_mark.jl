

    # drop_rate = zeros(1,M);
    #
    # drop_rate(lag(3*N+1: 3*N+M) > Kmax) = 1;
    #
    # ecnmark = (lag(3*N+1: 3*N+M) - Kmin)./(Kmax-Kmin).*pmax;
    # watermark = lag(3*N+1: 3*N+M)<=Kmax & lag(3*N+1: 3*N+M)>Kmin;
    # drop_rate(watermark) = ecnmark(watermark);
    #
    # nodrop = prod(ones(N,M) - matrix .* drop_rate, 2);
