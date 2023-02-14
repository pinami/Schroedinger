import Statistics
using Plots
using NLsolve

include("airyJU.jl")
import .AiryJU

N = 300    # Globalna varijabla - broj elemenata mreže 
M = 2*N+2  # Globalna varijabla - dimenzija sustava
h = 1.0/N  # Globalna varijabla - finoća mreže
ε = 0.08   # Globalna varijabla - skalirani Planck
τ = 2      # Globalna varijabla - relaksacijski faktor
γ = 1.5    # Globalna varijabla - isentropni koeficijent
Umax = 0.5 # Globalna varijabla - maksimalni napon
U = 0.0    # Globalna varijabla - trenutni napon (za računanje J-U krivulje)
relax = 0  # Globalna varijabla - uključi isključi relaksaciju
press = 2  # Globalna varijabla - uključi isključi tlak
T = 0.01    # Globalna varijabla - skalirana temperatura


function isothermal(n)
    @assert n > 0
    T*log(n)
end

function d_isothermal(n)
    @assert n > 0
    T*1.0/n
end

function isentropic(n)
    @assert n > 0
    @assert γ > 1
    T*γ*n^(γ-1)/(γ-1)
end

function d_isentropic(n)
    @assert n > 0
    @assert γ > 1
    T*γ*n^(γ-2)
end

# Ovdje se bira tip tlaka. 
f = isentropic
df = d_isentropic
if press == 2 
    f = isothermal
    df = d_isothermal
end
# Funkcija Z bez relaksacijskog člana 
function Z(x,u,v)
    (h/ε)^2 * (U*x + press *f(u^2+v^2) 
                   + relax * (ε/τ) * atan(v,u))
end

function dZ_du(u,v)
    (h/ε)^2 * ( press * 2*u*df(u^2+v^2) 
              - relax * (ε/τ) * v /(u^2+v^2) )
end

function dZ_dv(u,v)
    (h/ε)^2 * (press * 2*v*df(u^2+v^2) 
             + relax * (ε/τ) * u /(u^2+v^2) )
end

# Y = je trenutno rješenje = [u,v], U je potencijal
function residual!(Res, Y)
    M = length(Y)
    N::Int64 = (M-2)/2
    @assert 2*N+2 == M
    x = LinRange(0.0,1.0, N+1) # broj čvorova je N+1
    u = Y[1:N+1]
    v = Y[N+2:M]
    @assert size(Res) == size(Y)
    Res[1] = u[1] - 1
    for i=2:N
        Res[i] = -u[i-1] + 2*(1 + Z(x[i],u[i],v[i]))*u[i] - u[i+1]
    end
    Res[N+1] = u[N+1] - cos(U/ε) 
    K = N+1
    Res[K+1] = v[1]
    for i=2:N
        Res[K+i] = -v[i-1] + 2*(1 + Z(x[i],u[i],v[i]))*v[i] - v[i+1]
    end
    Res[K+N+1] = v[N+1] - sin(U/ε)
    return Res
end

# Ima grešku. Za sada ga ne koristim
function jacobian!(A, Y)
    M = length(Y)
    N::Int64 = (M-2)/2
    @assert 2*N+2 == M
    x = LinRange(0.0,1.0,N+1) # broj čvorova je N+1
    u = Y[1:N+1]
    v = Y[N+2:M]
    K=N+1  # indeks baze za drugi dio matrice
    
    A[1,1] = 1.0
    for i=2:N
        A[i,i] = 2*(1 + Z(x[i],u[i],v[i]) ) + 2*dZ_du(u[i],v[i])*u[i]
        A[i,i+1] = -1 
        A[i,i-1] = -1
        A[i,K+i] = 2*dZ_dv(u[i],v[i])*u[i]
    end
    A[N+1,N+1] = 1.0
    # Donji dio matrice
    A[K+1,K+1] = 1.0
    for i=2:N
        A[K+i,K+i] = 2*(1 + Z(x[i],u[i],v[i]) ) + 2*dZ_dv(u[i],v[i])*v[i]
        A[K+i,K+i+1] = -1 
        A[K+i,K+i-1] = -1
        A[K+i,i] = 2*dZ_du(u[i],v[i])*v[i]
    end
    A[K+N+1,K+N+1] = 1.0
end

# Current J calculated form the solution Y
function J_current(Y)
    M = length(Y)
    N::Int64 = (M-2)/2
    @assert 2*N+2 == M
    x = LinRange(0.0,1.0,N+1) # broj čvorova je N+1
    u_re = Y[1:N+1]
    u_im = Y[N+2:M] 
    J = Array{Float64}(undef,N+1)
    J .= zero(Float64)
    for i=2:N
        J[i] = ε*(u_re[i]*(u_im[i+1]-u_im[i-1]) 
                - u_im[i]*(u_re[i+1]-u_re[i-1]) )/(2*h)
    end
    J[1] = J[2] 
    J[N+1] = J[N]
    return J
end

function plotdensity()                          
    u_re = Sol[1:N+1]
    u_im = Sol[N+2:M]                          
    x = LinRange(0.0,1.0,N+1) 
    plot(x, u_re.*u_re.+u_im.*u_im, label="n") 
end



println("Početak petlje:")
UPtsMax = 100
Sol  = Array{Float64}(undef,M)
JUcurve = Array{Float64}(undef,UPtsMax+1)
JUcurve[1]=0.0

function solve()
    println("solving with ε = $ε, press=$press, relax=$relax")
for i=1:UPtsMax
    global U = i*Umax/UPtsMax
    print("step $i, U = $U\n")
    Sol .= one(Float64)
#   local res =nlsolve(residual!, jacobian!, Sol, factor=0.1)
    local res =nlsolve(residual!, Sol, autodiff = :forward)
#    local res =nlsolve(residual!, Sol) #finite diff
    if !converged(res)
        println("Newton did not converged!")
    end
    global Sol = copy(res.zero)
    J = J_current(Sol)  
    #local x = LinRange(0.0,1.0, N+1) 
    #plot(x, J, label="J")
    global JUcurve[i+1] = Statistics.mean(J)
    if JUcurve[i+1] > 1E-8
        if Statistics.std(J) > JUcurve[i+1]*1E-3
            println("Warning: J may not be constant")
        end
    end
#    println("srednja vrijednost = ", JUcurve[i])
#    println("standardna devijacija = ", Statistics.std(J))
end
end

solve()
if press == 0
    titleString = "ε = $ε"
    if relax == 1
        titleString = "ε = $ε, τ = $τ"
    end
elseif press == 1
    titleString = "ε = $ε, pressure isentropic"
    if relax == 1
        titleString = "ε = $ε, τ = $τ, pressure isentropic"
    end
elseif press == 2
    titleString = "ε = $ε, pressure isothermal"
    if relax == 1
        titleString = "ε = $ε, τ = $τ, pressure isothermal"
    end
end
plot(LinRange(0.0,U,UPtsMax+1), JUcurve,  
     xlabel="Potential U", ylabel="Current J",
     title = titleString,  label="JU curve", linewidth=2)

 # Usporedba s rješenjem bez tlaka i relaksacije
 U,J1 = AiryJU.JU(1,eps, Umax, UPtsMax+1)
 plot!(U, J1,   label="simple", linewidth=2)    
 savefig("isothermal_0.1.png") 
#plot!(x, J, label="J")
