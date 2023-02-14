module AiryJU
# Module for calculation of current-voltage curve using 
# Airy functions 
using SpecialFunctions


function An(x,n)
    p=1.0/(n+2)
    EPS = 1/400
    w = (2.0/(n+2))*(x.^((n+2)/2))
    return p*sqrt.(x) *( besseli.(-p,w) - besseli.(p,w) )
end

function Bn(x,n)
    p=1.0/(n+2)
    w = (2.0/(n+2))* x.^( (n+2)/2 )
    return sqrt.(p*x)*( besseli.(-p,w) + besseli.(p,w) )
end


function JU(n = 1, eps = 0.08,  Umax=0.5, Npts=300)
    U = LinRange(1/Npts,Umax,Npts)
    p = 1/(n+2) 
    beta =(2*U/(eps*eps)).^(1.0/(2+n))
    ai = An.(beta,n)
    bi = Bn.(beta,n)
    ai0 = p^(1-p) / gamma(1-p)
    bi0 = p^(0.5-p) / gamma(1-p)
    sinu = sin.(U/eps)
    V = 2*sin(pi/(n+2)).*beta*eps.*sinu/(pi*sqrt(n+2))
    J = V ./ (ai0*bi-bi0*ai)
    return ([0,U],[0,J])
end

# JU-curve for n=1. Implementation is more direct than using 
# current(1,...)
function JU_n0(eps = 0.08, Umax = 0.5, M = 200)
    U = LinRange(0.0, Umax, M+1)
    beta = (2*U/(eps*eps)).^(1.0/3.0)
    ai = airyai.(beta)
    dai = airyaiprime.(beta)
    bi = airybi.(beta)
    dbi = airybiprime.(beta)

    sinu = sin.(U/eps)
    V = beta .*sinu*eps/pi
    divisor = ai[1]*bi-bi[1]*ai
    divisor[1] = 1 # J[1] mora biti nula
    J = V ./ divisor
    return (U,J)
end

end # module