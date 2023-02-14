using Plots
include("airyJU.jl")
import .AiryJU

eps = 0.08

U,J1 = AiryJU.JU(1,eps)
_,J2 = AiryJU.JU(2,eps)
_,J3 = AiryJU.JU(3,eps)

plot(U,J1, label="n=1", xlabel="voltage U", ylabel="current J", 
           title="JU curve, eps = $eps", linewidth=2)
plot!(U,J2, label="n=2",linewidth=2)
plot!(U,J3, label="n=3", linewidth=2)
#savefig("JU-2.png")


U,J1 = AiryJU.JU(1,0.08) 
_,J2 = AiryJU.JU(1,0.06)
_,J3 = AiryJU.JU(1,0.04)
plot(U,J1,  label="eps=0.08", xlabel="voltage U", 
            ylabel="current J", title="JU curve, n=1", linewidth=2)
plot!(U,J2, label="eps=0.06", linewidth=2) 
plot!(U,J3, label="eps=0.04", linewidth=2)
savefig("JU-3.png")

U,J1 = AiryJU.JU(2,0.08) 
_,J2 = AiryJU.JU(2,0.06)
_,J3 = AiryJU.JU(2,0.04)
plot(U,J1,  label="eps=0.08", xlabel="voltage U", 
            ylabel="current J", title="JU curve, n=2", linewidth=2)
plot!(U,J2, label="eps=0.06", linewidth=2) 
plot!(U,J3, label="eps=0.04", linewidth=2)
savefig("JU-4.png")

