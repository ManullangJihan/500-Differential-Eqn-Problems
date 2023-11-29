using DifferentialEquations
using LinearAlgebra
using Plots
using Plots: plot, plot!
using ModelingToolkit


# Example 1.1
f(m, k, t) = -k * m
m0 = 10
k = -log(0.5) / 6.0
t = (0.0, 10.0)
problem = ODEProblem(f, m0, t, k)
sol = solve(problem)
plot(
    sol,
    lw=2.0, 
    label="Numeric Solution", 
    xaxis="Time", 
    yaxis="Harm Pesticide in Pounds",
    title="Pesticide Decomposition"
)

plot!(
    sol.t, t -> 10exp(-k*t),
    lw=1.5,
    label="Analytic solution"
)

savefig("Pesticide_decomposition.png")


# Example 1.2 
f(p, k, t) = k * p
p0 = 120
k = (1/3) * log(5/3)
t = (0.0, 10.0)
problem = ODEProblem(f, p0, t, k)
sol = solve(problem)
plot(sol,
     label="Numeric Solution",
     xlabel="Time",
     ylabel="Population",
     title="Bacteria growth")
plot!(sol.t,
      t -> 120exp(log(5/3)/3 * t),
      label="Analytic Solution")

savefig("ex_12_Bacteria_growth.png")

# Example 1.3
f(I, a, x) = -a * I
I₀ = 100.0 
a = 0.5
x = (0.0, 20.0)
problem = ODEProblem(f, I₀, x, a)
sol = solve(problem)
plot(sol,
     label="Numeric Solution",
     xlabel="Distance",
     ylabel="Intensity",
     title="Sound's waves Intensity travel through the air")
plot!(sol.t,
      t -> I₀*exp(-a*t),
      label="Analytic Solution")

savefig("ex_13_Sound_intensity.png")

# Example 1.4
f(I, D, r) = -D * I 
I₀ = 100.0
D = 0.7
r = (0.0, 30.0)
problem = ODEProblem(f, I₀, r, D)
sol = solve(problem)
plot(sol,
     label="Numeric Solution",
     xlabel="Thickness",
     ylabel="Intensity",
     title="X-ray absorption")
plot!(sol.t,
      t -> I₀*exp(-D*t),
      label="Analytic Solution")

savefig("ex_14_Xray_absorption.png")

# Example 1.5
f(P, r, t) = r * P
r = 0.4
P₀ = 100.0
t = (0.0, 3.0)
problem = ODEProblem(f, P₀, t, r)
sol = solve(problem)
plot(sol,
     label="Numeric Solution",
     xlabel="Time",
     ylabel="Balance",
     title="Bank Deposits")
plot!(sol.t,
      t -> P₀*exp(r*t),
      label="Analytic Solution")

savefig("ex_15_bank_deposit.png")

# Example 1.6 
# Radioactive decay

f(N, k, t) = k * N

half_life_dict = Dict(
    "Alumunium" => 7.4e5,
    "Beryllium" => 1.51e6,
    "Carbon" => 5730,
    "Chlorine" => 3.01e5,
    "Iodine" => 8.05,
    "Potassium" => 1.2e9,
    "Polonium209" => 100,
    "Polonium210" => 0.3780821917808219,
    "Radon" => 0.010465753424657534,
    "Radium226" => 1700,
    "Thorium" => 75000,
    "Uranium" => 4.51e9
)

N₀ = 1.0 
k = -(1/half_life_dict["Alumunium"]) * log(2)
t = (0.0, 10000.0)
problem = ODEProblem(f, N₀, t, k)
sol = solve(problem)
plot(sol, label="Numeric Solution", lw=1)

N(t, element) = N₀ * exp(-log(2)*t/half_life_dict[element])
plot!(sol.t, t->N(t, "Alumunium"), lw=4, label="Analytic Solution", opacity=0.5)
savefig("ex_16_radioactive_decay.png")

# Example 1.7
function number_scientist(ds, p, t)
    b, d = p
    ds = (b - d) * ds
end

p = [4, 3]
s₀ = 10
t = (0.0, 5.0)
problem = ODEProblem(number_scientist, s₀, t, p)
sol = solve(problem)
plot(sol, xlabel="Time in years", ylabel="Number of scientist", 
     title="Scientist population", label="Numeric solution")

S(t) = s₀ * exp((p[1] - p[2]) * t)
plot!(sol.t, t -> S(t), label="Analytic Solution")
savefig("ex_17_number_of_scientist.png")

# Example 1.8
function thomas_robert_population(dp, p, t)
    k, m = p
    dp = k * dp + m
end

p₀ = 1000
p = [0.016, 5000]
t = (0.0, 5.0)
problem = ODEProblem(thomas_robert_population, p₀, t, p)
sol = solve(problem)

P(t) = exp(p[1]*t) * (p₀ + p[2]/p[1]) - p[2]/p[1]

plot(sol,
     label="Numeric Solution",
     xlabel="Time",
     ylabel="Population",
     title="Thomas Robert Population law")
plot!(sol.t,
      t -> P(t),
      label="Analyic Solution")
savefig("ex_18_thomas_population.png")

rate(t) = exp(p[1]*t) * (p[1] * p₀ + p[2])
PP(t) = (1/p[1]) * (rate(t) - p[2])
plot!(sol.t, t->PP(t), label="x")

# Example 1.9
function law_of_cooling(dtheta, p, t)
    k, T = p
    dtheta = - k * (dtheta - T)
end

θ₀ = 70
p = [1.2, 50]
t = (0.0, 10.0)
problem = ODEProblem(law_of_cooling, θ₀, t, p)
sol = solve(problem)
plot(
    sol, 
    xlabel="Time in second", 
    ylabel="Temperature", 
    title="Law of cooling", 
    label="Numeric Solution"
)

T(t) = p[2] + exp(-p[1] * t) * (θ₀ - p[2])
plot!(
    sol.t,
    t -> T(t),
    label="Analytic Solution",
)
savefig("ex_19_law_of_cooling.png")

# Example 1.10
function commodity_price(dp, p, t)
    k, α, β = p
    dp = k * (α - β*dp)
end

p = [0.1, 2.0, 2.0]
p₀ = 100.0
t = (0.0, 5.0)
problem = ODEProblem(commodity_price, p₀, t, p)
sol = solve(problem, Tsit5())
plot(
    sol,
    xlabel="Time",
    ylabel="Commodity Price",
    title="Price of Commodity over time",
    label="Numeric Solution"
)

k, β, α = p
C = -β*p₀ + α
P(t) = -1/β * (C*exp(-β*k*t) - α)
plot!(
    sol.t,
    t -> P(t),
    label="Analytic Solution"
)
savefig("ex_20_commodity_price")


# Example 1.11
# Equation of continuity

function conservation_equation(dy, p, t)
    dy = 10 - (4 / (500 + t)) * dy
end

y₀ = 50.0
t = (0.0, 10.0)
problem = ODEProblem(conservation_equation, y₀, t)
sol = solve(problem)
plot(
    sol, xlabel="Time", ylabel="Salt in grams",  title="Equation of continuity",
    label="Numeric Solution", color=:blue, linewidth=2.0)

y(t) = 2*(500 + t) - 950*(1 + t / 500)^-4
plot!(
    sol.t, t -> y(t), label="Analytic Solution", color=:red, linewidth=1.5)
savefig("ex_111_equation_of_continuity.png")

# Example 1.12
function glucose_volume(dg, p, t)
    I, k, V = p
    dg = (I - k*dg) / V
end

p = [0.2, 0.18, 10.0]
g₀ = 1.0
t = (0.0, 10.0)
problem = ODEProblem(glucose_volume, g₀, t, p)
sol = solve(problem)
plot(
    sol, xlabel="Time", ylabel="Blood Volume", title="Volume of Blood",
    lw=2.0, label="Numeric Analysis"
)
savefig("ex_112_volume_of_blood.png")

# Example 1.13
function dose_of_medicine(da, k, t)
    da = -k * da
end

k = 0.2
a₀ = 10
t = (0.0, 8.0)

p = plot()
for i = 1:9    
    problem = ODEProblem(dose_of_medicine, a₀, t, k)
    sol = solve(problem)
    plot!(
        p, sol, xlabel="Time", ylabel="Concentration of medicine",
        label=false, title="Drug dosage", color="blue"
    )
    t = t .+ 8.0
    a₀ += sol[end]
end

plot(p, xlim=(0, 72))
savefig("ex_113_drug_dosage.png")

# Example 1.14
const g = 9.81

function rocket_launch(du, u, r, t)
    du[1] = -0.6
    du[2] = (500 * 0.6 / u[1]) - g
end

u0 = [10, 500]
tspan = (0.0, 15.0)
problem = ODEProblem(rocket_launch, u0, tspan)
sol = solve(problem)
plot(sol, idxs=(2))

V(t) = -500 * log((10 - 0.6*t) / 10) - g*t
plot!(
    sol.t,
    t -> V(t),
)

# Problem 1.1
function cable_suspension(dy, p, x)
    W, H = p
    dy = W / H * x
end

y₀ = [0.0]
p = [10, 20]
x = (-10.0, 10.0)
problem = ODEProblem(cable_suspension, y₀, x, p)
sol = solve(problem)
plot(sol)

y(x) = (p[1]/(2p[2])) * x^2
t = collect(-10:0.1:10)
plot!(
    t,
    y.(t),
    label="Analytic Solution"
)


# Problem 1.2
function piston_with_orifices(dv, k, t)
    dv = -k * dv
end

k = 0.5
v₀ = 10.0
t = (0.0, 10.0)
problem = ODEProblem(piston_with_orifices, v₀, t, k)
sol = solve(problem)
plot(sol)

# Problem 1.3
function body_of_snake(u, p, t)
    a, n = p
    return (a/3n) * u
end

p = [1.0, 1.078]
l₀ = [5.0]
t = (0.0, 10.0)
problem = ODEProblem(body_of_snake, l₀, t, p)
sol = solve(problem)
plot(sol)

a, n = p[1], p[2]
c = l₀[1]
L(t) = c*exp(a*t/(3n)) 
t = collect(0.0:0.1:10.0)
plot!(t, L.(t), label="Real Solution", color=:red)

# Problem 1.4
function cancer_research(dl, p, t)
    g, k, ξ, σ = p
    dl = -k * ξ * σ * dl
end 

p = [1.0, 1.2, 2.7, 3.1]
l₀ = 10.0
t = (0.0, 1.0)
problem = ODEProblem(cancer_research, l₀, t, p)
sol = solve(problem)
plot(sol, label="ODE", color=:blue)

L(t) = 10.0 * exp.(-1.2*2.7*3.1*t)
t = collect(0.0:0.01:1.0)
plot!(t, L(t), color=:red, label="Real Solution")

# Problem 1.5

# Problem 1.6

# Problem 1.7
moisture(m, r, t) = -r * m
m₀ = 100.0
r = -log(0.5)/15.0
t = (0.0, 100.0)
problem = ODEProblem(moisture, m₀, t, r)
sol = solve(problem)
plot(sol)

# Problem 1.8


# Problem 1.9

# Problem 1.10
function seasonal_population(u, p, t)
    return u*p*cos(t)
end

p₀ = 1_000_000
p = 0.3
t = (0.0, 10.0)
problem = ODEProblem(seasonal_population, p₀, t, p)
sol = solve(problem)
plot(sol)

P(t) = p₀*exp(p*sin(t))
plot!(sol.t, t->P.(t), label="Analytic solution")

# Problem 1.11
function f(u, p, x)
    α, β, γ, δ = p
    return (α - x) / (β + γ*x + δ*x^2) * u
end

# α = β = δ = 0, γ > 0 (exponential distribution)
p1 = [0.0, 0.0, 0.2, 0.0]

# γ = δ = 0, β > 0 (normal distribution)
p2 = [0.4, 0.7, 0.0, 0.0]

# β = δ = 0, γ > 0, α > −γ (gamma distribution)
p3 = [1.0, 0.0, 0.5, 0.0]

# β = 0, γ = −δ, (α − 1)/γ < 1, α/γ > −1 (beta distribution)
p4 = [0.5 , 0.0, -1.2, 5/-1.2]


y₀ = 2.0
t = (0.1, 100.0)

prob1 = ODEProblem(f, y₀, t, p1)
prob2 = ODEProblem(f, y₀, t, p2)
prob3 = ODEProblem(f, y₀, t, p3)
prob4 = ODEProblem(f, y₀, t, p4)

sol1 = solve(prob1, Tsit5())
sol2 = solve(prob2)
sol3 = solve(prob3)
sol4 = solve(prob4)

plot(sol1, color=:blue, label="Exponential Dist")
plot(sol2, color=:red, label="Normal Dist")
plot(sol3, color=:green, label="Gamma Dist")
plot(sol4, color=:purple, label="Betta Dist")

# Problem 1.12

# Problem 1.13
function company_capital(dy, p, t)
    N, r, s = p
    dy = (1 - N)*r*dy + s
end

p = [0.4, 0.5, 0.7]
y₀ = 100
t = (0.0, 12.0)
problem = ODEProblem(company_capital, y₀, t, p)
sol = solve(problem)
plot(sol)

# Problem 1.14

# Problem 1.15

# Problem 1.16
function deep_diving(dy, u, p, x)
    α, β, a, b = p
    dy[1] = α * u[1] + β + b * exp(-a * x)
end

p = [1.0, 0.0, 0.4, 1.2]
y₀ = [1.0]
x = (0.0, 10.0)
problem = ODEProblem(deep_diving, y₀, x, p)
sol = solve(problem)
plot(sol)

α, β, a, b = p
c = y₀[1]
y(x) = -β/α + (-b/(a + α) * exp(-a*x)) + c*exp(α*x)
plot!(
    sol.t,
    t -> y.(t),
    label="Analytic Solution"
)

# Problem 1.17
function education_investment(dy, p, t)
    dy = 1 - p * dy
end 

p = 1.0
y₀ = 5.0
t = (0.0, 10.0)
problem = ODEProblem(education_investment, y₀, t, p)
sol = solve(problem)
plot(sol)

y(t) = y₀*exp(-p*t) + (1/p)
plot!(
    sol.t,
    t -> y.(t),
    label="Analytic Solution"
)

# Problem 1.18
function memorization_process(dy, p, t)
    a, b, M = p
    dy = a*(M - dy) - b * dy
end

p = [0.2, 0.1, 0.5]
y₀ = 2.0
t = (0.0, 100.0)
problem = ODEProblem(memorization_process, y₀, t, p)
sol = solve(problem)
plot(sol)

a, b, M = p
y(t) = a*M/(a+b) * (1 - exp(-(a+b)*t))
plot!(
    sol.t,
    t -> y.(t),
    label="Analytic Solution"
)


# Problem 1.19

# Problem 1.20
function torque(dw, p, t)
    J, B, T = p
    dw = (T - B*dw) / J
end

p = [1.2, 0.1, 0.4]
w₀ = 10.0
t = (0.0, 10.0)
problem = ODEProblem(torque, w₀, t, p)
sol = solve(problem)
plot(sol)

J, B, T = p
w(t) = w₀*exp(-B*t/J) + ((1 - exp(-B*t/J)) * T/B)
plot!(
    sol.t,
    t -> w.(t),
    label="Analytic Solution"
)


# Problem 1.21
function drug(dy, p, t)
    k, c = p
    dy = c - k*dy
end

p = [0.4, 0.2]
y₀ = 5.0
t = (0.0, 12.0)
problem = ODEProblem(drug, y₀, t, p)
sol = solve(problem)
plot(sol)

# Problem 1.22

# Problem 1.23
function cell_growth(dv, p, t)
    a, b, c, r, M, W = p
    dv = W + a*M + (2*c/r) * (M - b*dv)
end

p = [0.2, 0.4, 1.0, 2.0, 0.3, 0.1]
v₀ = 2.0
t = (0.0, 12.0)
problem = ODEProblem(cell_growth, v₀, t, p)
sol = solve(problem)
plot(sol)


a, b, c, r₀, M, W = p
k = 2*b*c/r₀
A = W/k + M/k*(a + 2*c/r₀)
v(t) = A + (v₀ - A)*exp(-k*t)
plot!(
    sol.t,
    t -> v.(t),
    label = "Analytic Solution"
)

# Problem 1.24

# Problem 1.25
function cell_density(du, u, p, t)
    μ, A = p
    du[1] = μ*(u[1] + A)
end

p = [0.2, 2.0]
u0 = [0.0]
tspan = (0.0, 10.0)
problem = ODEProblem(cell_density, u0, tspan, p)
sol = solve(problem)
plot(sol)

μ, A = p
N(t) = A*(exp(μ*t) - 1)
plot!(
    sol.t,
    t -> N.(t),
    label="Analytic Solution"
)

# Problem 1.26
function nerve(du, u, p, t)
    x₀, y₀ = u
    A, α, B, b, x, y = p 
    du[1] = dx = A - α*(x₀ - x)
    du[2] = dy = B - b*(y₀ - y)
end

p = [0.2, 0.3, 0.5, 0.1, 15, 13]
u₀ = [2.0, 3.0]
tspan = (0.0, 20.0)
problem = ODEProblem(nerve, u₀, tspan, p)
sol = solve(problem)
plot(sol, idxs=(1))

x₀, y₀ = u₀
A, α, B, b, x, y1 = p
x1(t) = x₀ + ((1 - exp(-α*t))*A/α)
y2(t) = y₀ + ((1 + exp(-b*t))*B/b)

plot!(
    sol.t,
    t -> x1.(t),
)

plot!(
    sol.t,
    t -> y2.(t),
)

# Problem 1.27

# Problem 1.28
function secretion_of_hormone(dy, p, t)
    a, b, k = p
    dy = a - b*cos.(π*t/12) - k * dy
end

p = [0.6, 0.3, 0.4]
y₀ = 3.0
t = (0.0, 48.0)
problem = ODEProblem(secretion_of_hormone, y₀, t, p)
sol = solve(problem)
plot(sol)

a, b, k = p
b1(t) = (y₀ - a/k + 144*b*k/(144*k^2 + π^2))*exp(-k*t)
y(t) = a/k - (144*b*k/(144*k^2 + π^2) * (cos(π*t/12) + π/(12*k) * sin(π*t/12))) + b1(t)
plot!(
    sol.t,
    t ->  y.(t)
)


# Problem 1.29
function electric_circuit(di, p, t)
    L, R, E, k = p
    di = (E * sin.(k * t) - R*di)/L
end 

p = [0.4, 0.2, 0.5, 0.1]
i₀ = 3.0
t = (0.0, 10.0)
problem = ODEProblem(electric_circuit, i₀, t, p)
sol = solve(problem)
plot(sol)

L, R, E, k = p
A1(t) = (R*E*sin(k*t) - k*L*E*cos(k*t))/(R^2 + k^2*L^2)
i(t) = A1(t) + i₀*exp(-(R/L)*t)
plot!(
    sol.t,
    t -> i.(t),
    label="Analytic Solution"
)


# Problem 1.30
