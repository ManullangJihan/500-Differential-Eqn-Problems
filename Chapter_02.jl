using DifferentialEquations
using LinearAlgebra
using Plots
using Plots: plot, plot!

# Example 2.1
function population(dp, u, p, t)
    a, b = p
    P = u[1]
    dp[1] = P*(a - b*P)
end

p = [0.02, 0.03]
p₀ = [0.1]
t = (0.0, 1000.0)
problem = ODEProblem(population, p₀, t, p)
sol = solve(problem)
plot(sol, label="Numeric Solution")

a, b = p
p0 = p₀[1]
P(t) = a*p0 / (b*p0 + (a - b*p0) * exp(-a*t))
plot!(
    sol.t,
    t -> P.(t),
    label="Analytic Solution"
)

# Example 2.2
function rumor_spread(dp, p, t)
    k, N = p
    dp = k * dp * (N - dp)
end

p = [0.002, 1000]
p₀ = 10
t = (0.0, 12.0)
problem = ODEProblem(rumor_spread, p₀, t, p)
sol = solve(problem)
plot(sol)

# Example 2.3
function technology_adoption(dn, p, t)
    k, C = p
    dn = k * dn * (C - dn)
end

u0 = 1
p = [0.3, 10]
tspan = (0.0, 10.0)
problem = ODEProblem(technology_adoption, u0, tspan, p)
sol = solve(problem)
plot(sol, label="Numeric Solution")

# Example 2.4
function undissolved_solute(dy, p, t)
    k, v, c = p
    dy = k*dy*()
end

p = [0.01, 20.0, 10.0]
y₀ = 1.0
t = (0.0, 12.0)
problem = ODEProblem(undissolved_solute, y₀, t, p)
sol = solve(problem)
plot(sol)

# Example 2.5
function epidemiology(du, u, p, t)
    k, m = p
    s = u[1]
    du[1] = k*(s) - k/m
end

p = [0.1, 10]
u0 = [1000/800]
tspan = (0.0, 24.0)
problem = ODEProblem(epidemiology, u0, tspan, p)
sol = solve(problem)
plot(sol)

N, k, m = p
s = u0[1]
S(t) = m/(1 + (m-1)*exp(k*t))
plot!(
    sol.t,
    t -> S.(t) ./ N,
    label="Analytic Solution"
)


# Example 2.6
function hmb_temperature(dT, u, p, r)
    Q, k = p
    T₀, T1 = u
    dT[1] = -Q/(4π*k*(T₀ - T1))
end
    
p = [30.0, 1.5]
T₀ = [10.0, 15.0]
r = (0.0, 6.0)
problem = ODEProblem(hmb_temperature, T₀, r, p)
sol = solve(problem)
plot(sol)

# Example 2.7
function hydrogen(dx, p, t)
    a, b, K = p
    dx = K*(a - 2dx)*(b- dx)
end

p = [20.0, 12.0, 0.2]
x₀ = 1.0
t = (0.0, 12.0)
problem = ODEProblem(hydrogen, x₀, t, p)
sol = solve(problem)
plot(sol)

# Example 2.8
function evangelista(dv, p, t)
    A, h, k = p
    dv = k*A*sqrt(2*9.8*h)
end 

p = [0.001, 100.0, -0.6]
v₀ = 1000
t = (0.0, 100.0)
problem = ODEProblem(evangelista, v₀, t, p)
sol = solve(problem)
plot(sol)

# Example 2.9
function mirror_reflection(du, u, p, x)
    y = p[1]
    du[1] = y/(x + sqrt(x^2 + y^2))
end

p = [1.0]
u0 = [0.0]
tspan = (0.0, 10.0)
problem = ODEProblem(mirror_reflection, u0, tspan, p)
sol = solve(problem)
plot(sol)

# Example 2.10
function escape_velocity(du, u, p , t)
    y, R = p
    du[1] = dv = -9.8*R^2/(y + R)^2 * u[1]
end

p = [0, 6371]
v0 = [11.2654]
tspan = (0.0, 1.0)
problem = ODEProblem(escape_velocity, v0, tspan, p)
sol = solve(problem)
plot(sol, label="Analytic Solution")

# Problem 2.1
function infection(di, k, t)
    di = k * di * (100 - di)
end

k = 0.1
i₀ = 1.0
t = (0.0, 10.0)
problem = ODEProblem(infection, i₀, t, k)
sol = solve(problem)
plot(sol)

I(t) = 100exp(100k*t)/(99 + exp(100k*t))
plot!(
    sol.t,
    t -> I.(t),
    label="Analytic Solution"
)

# Problem 2.2
function walruses_population(dm, p, t)
    a, b = p
    dm = a * dm - b * dm * (dm -1)
end

p = [0.4, 0.01]
m₀ = 1000
t = (0.0, 20.0)
problem = ODEProblem(walruses_population, m₀, t, p)
sol = solve(problem)
plot(sol)

a, b = p
male_end = (a + b)/b

println(male_end)
println(sol[end])


# Problem 2.3
g = 9.8
B = 50.0
y = 20.0
s = 0.4
v₀ = 1.0
y₀ = 5.0
c₀ = 0.5 * (g * y₀)
a = g * s / (v₀ + c₀)^2
b = 3 / (v₀ + c₀)^2

#τ(t, x) = (v₀ + c₀) * t - x
#c₁(t, x) = (sqrt(g*y) - c₀) / (τ(t, x))

function flooding(dc, p, x)
    dc = b*dc^2 - a*dc
end

c₀ = 17.9
x = (0.0, 100.0)
problem = ODEProblem(flooding, c₀, x)
sol = solve(problem)
plot(sol)


# Problem 2.4
C = 2.0
k = 1.0

function ionized_gas(dn, p, t)
    dn = C - k*dn^2
end

n₀ = 0.0
t = (0.0, 6.0)
problem = ODEProblem(ionized_gas, n₀, t)
sol = solve(problem)
plot(sol)

# Problem 2.5
function logistic_law(dp, p, t)
    a, b, H = p
    dp = a * dp - b * dp^2 - H
end

p = [0.8, 0.1, 1.7]
p₀ = 10
t = (0.0, 12.0)
problem = ODEProblem(logistic_law, p₀, t, p)
sol = solve(problem)
plot(sol)

a, b, H = p
c = p₀
P1(t) = 1/(2b) * (a + sqrt(4b*H - a^2) * tan(1/2*(c- t) * sqrt(4b*H - a^2)))

plot!(
    sol.t,
    t -> P1.(t),
    label="Analytic solution"
)

# Problem 2.6

function benjamin_population(dp, p, t)
    a, b = p
    dp = dp * (a - b*log(dp))
end

a = 0.6
b = 0.1
p = [a, b]
p₀ = 100
t = (0.0, 1000.0)
problem = ODEProblem(benjamin_population, p₀, t, p)
sol = solve(problem)
plot(sol)

c = a/b - log(p₀)
Benjamin(t) = exp(a/b) * exp(-c*exp(-b*t))
plot!(
    sol.t,
    t -> Benjamin.(t),
    label="Analytic Solution"
)

# Problem 2.7
function critical_population(du, u, p, t)
    k, Pc, Pm = p
    du[1] = k*u[1]*(u[1] - Pc)*(Pm - u[1])
end

k = 0.1
Pc = 100
Pm = 300
p = [k, Pc, Pm]
u0_1 = [80]
u0_2 = [150]
u0_3 = [400]
tspan = (0.0, 5.0)
problem = ODEProblem(critical_population, u0_3, tspan, p)
sol = solve(problem)
plot(sol, label="Numeric Solution", xlabel="Time", ylabel="Population", title="Critical Population")

# Problem 2.8
function population3(du, u, p, t)
    k, n = p
    du[1] = k*u[1]^(n+1)
end

k = 0.1
n = 0.2
p = [k, n]
P₀ = [10]
tspan = (0.0, 1.0)
problem = ODEProblem(population3, P₀, tspan, p)
sol = solve(problem)
plot(sol)

P0 = 10
P3(t) = (P0^(-n) - n*k*t)^(-1/n)
plot!(
    sol.t,
    t -> P3.(t),
    label="Analytic Solution"
)

# Problem 2.9
function allometric_equation1(du, u, p, t)
    k, α = p
    du[1] = dw = k*u[1]^α
end

k = 0.4
α = 0.2
p = [k, α]
w₀ = [1.0]
tspan = (0.0, 6.0)
problem = ODEProblem(allometric_equation1, w₀, tspan, p)
sol = solve(problem)
plot(sol, label="Numeric solution")

function allometric_equation2(du, u, p, t)
    k, α, Wm = p
    μ = 1 - α
    h = 1 - (u[1] / Wm)^μ
    du[1] = dw = k*u[1]^α * h
end

Wm = 10.0
p = [k, α, Wm]
tspan = (0.0, 100.0)
problem2 = ODEProblem(allometric_equation2, w₀, tspan, p)
sol2 = solve(problem2)
plot(sol2)

function raising_fishing_cost(du, u, p, t)
    k, α, Wm, a1, a2 = p
    μ = 1 - α
    h = 1 - (u[1] / Wm)^μ
    du[1] = dw = k*u[1]^α * h
    du[2] = dc = a1 + a2*du[1]
end

a1 = 1
a2 = 0.1
p = [k, α, Wm, a1, a2]
u0 = [1.0, 1.0]
tspan = (0.0, 3.0)
problem3 = ODEProblem(raising_fishing_cost, u0, tspan, p)
sol = solve(problem3)
plot(sol, label=["Fish Weight" "Raising Cost"])

# Problem 2.10
function rhesus_monkey(du, u, p, s)
    k, n = p
    du[1] = k*u[1]^n/s
end

k = 0.2
n = 1
p = [k, n]
r0 = [0.1]
tspan = (0.0, 10.0)
problem = ODEProblem(rhesus_monkey, r0, tspan, p)
sol = solve(problem)

# Problem 2.11
function isometric_contraction(du, u, p, t)
end

# Problem 2.12
function autocatalytic_reaction(du, u, p, t)
    k, β = p
    du[1] = -k*u[1]*(β + α - u[1])
end

k = 0.2
β = 0.3
p = [k, β]
u0 = [1.0]
tspan = (0.0, 6.0)
problem = ODEProblem(autocatalytic_reaction, u0, tspan, p)
sol = solve(problem)
plot(sol)

# Problem 2.13
function chemical_reaction(du, u, p, t)
    k, a, b = p
    du[1] = dx = k*(a- u[1])^2 *(b - u[1]/2)
end

k = 0.2
a = 10.0
b = 20.0
p = [k, a, b]
u0 = [1.0]
tspan = (0.0, 10.0)
problem = ODEProblem(chemical_reaction, u0, tspan, p)
sol = solve(problem)
plot(sol)

# Problem 2.17
k = 0.2

function dna(u, p, t)
    return -k*u[1]^2
end

u0 = 10.0
tspan = (0.0, 10.0)
problem = ODEProblem(dna, u0, tspan)
sol = solve(problem)
plot(sol, label="Numeric Solution")

# Problem 2.18
a = 0.2
b = 0.3

function selection_two_genes(du, u, p, x)
    y = u[1]
    du[1] = dy = y^2*(1-y)*(x^2 - a^2) / (x^2 * (1-x) * (y^2 - b^2))
end

u0 = [100.0]
tspan = (1.0, 10.0)
problem = ODEProblem(selection_two_genes, u0, tspan)
sol = solve(problem)
plot(sol)

# Problem 2.19
a = 0.2
b = 0.25
c = 0.6
d = 0.4

function ownership_consumer_durable(u, p, t)
    return a*u -  a/b * (1 + exp(c-d*t)) * u^2
end

u0 = 10
tspan = (0.0, 6.0)
problem = ODEProblem(ownership_consumer_durable, u0, tspan)
sol = solve(problem)
plot(sol)

# Problem 2.20
a = 0.3
b = 0.2
α = 0.25

function cobb_douglas(u, p, t)
    return a*exp(b*t)*u^α
end

u0 = 10
tspan = (0.0, 6.0)
problem = ODEProblem(cobb_douglas, u0, tspan)
sol = solve(problem)
plot(sol)

analytic_eqn(t) = (a/b * (1 - α) * (exp(b*t) - 1) + u0^(1 - α))^(1/(1-α))
plot!(
    sol.t,
    t-> analytic_eqn.(t),
    label="Analytic Solution"
)

# Problem 2.21
k = 10
m = 40

function state_of_learner(u, p, t)
    return (2k/sqrt(m)) * (u * (1-u))^(3/2)
end

u0 = 0.1
tspan = (1.0, 10.0)
problem = ODEProblem(state_of_learner, u0, tspan)
sol = solve(problem)
plot(sol)

# Domain error = sensitive to initial value

# Problem 2.22
k = 0.001
T = 35

function josef_stefan(u, p, t)
    return k * (u^4 - T^4)
end

u0 = 1.0
tspan = (0.0, 1.0)
problem = ODEProblem(josef_stefan, u0, tspan)
sol = solve(problem)
plot(sol)

# Problem 2.23
μ = 0.4
h = 0.5

function two_body_problems(u, p, t)
    r = u
    return sqrt(2μ/r + 2h)
end

u0 = 1_000_000
tspan = (0.0, 100.0)
problem = ODEProblem(two_body_problems, u0, tspan)
sol = solve(problem)
plot(sol)

# Problem 2.24
