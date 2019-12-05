module GaussLaw1D

export Problem, solve, Constant

using OrdinaryDiffEq: Vern7
using ArgCheck: @argcheck
using BoundaryValueDiffEq: BoundaryValueDiffEq
const BV = BoundaryValueDiffEq
using Unitful: Unitful
const UF = Unitful

# export units for convenience
using Unitful: V, kV
export V, kV
using Unitful: pC, nC, μC, mC, C
export pC, nC, μC, mC, C
using Unitful: F
using Unitful: nm, μm, mm, cm, m
export nm, μm, mm, cm, m
using Unitful: ns, μs, ms, s
export ns, μs, ms, s
using Unitful: nGy, μGy, mGy, Gy
export nGy, μGy, mGy, Gy
using UnitfulRecipes # enable unitful recipes for convenience

const VACUUM_PERMITIVITY = 8.8541878128e-12*F * m^-1
const U_LENGTH = mm
const U_POT = V
const U_E = U_POT / U_LENGTH

const Length = typeof(1.0U_LENGTH)
const Potential = typeof(1.0U_POT)
const EFieldStrength = typeof(1.0U_E)

struct Problem{Fn}
    charge_density::Fn
    limits::NTuple{2, Length}
    shape::Symbol
    potential::NTuple{2, Potential}
    electric_permitivity::typeof(VACUUM_PERMITIVITY)
    radii::Vector{Length}
end

struct Constant{X}
    value::X
end
(c::Constant)(r) = c.value

function Problem(charge_density;
                 potential,
                 limits,
                 shape,
        electric_permitivity=VACUUM_PERMITIVITY,
        radii=collect(range(limits[1], stop=limits[2], length=1000)),
    )
    @argcheck issorted(limits)
    @argcheck limits[1] > 0.0mm
    @argcheck shape in (:parallel, :cylindrical, :spherical)
    r = limits[1]
    cd = charge_density(r)
    actual_cd_dim = UF.dimension(cd)
    expected_cd_dim = UF.dimension(1.0C/(mm^3))
    if  actual_cd_dim != expected_cd_dim
        msg = """Charge density $cd has the wrong unit. Expected unit is C*mm^-3.
        """
        throw(ArgumentError(msg))
    end
    Problem(charge_density, limits, shape, potential, electric_permitivity, radii)
end

function to_compute_unit(pot::Potential)
    UF.ustrip(U_POT, pot)
end

@noinline function bc!(residual, u, w, r)
    # @show residual
    # @show u
    # @show w
    # @show r
    # error()
    wanted_phi_left, wanted_phi_right = w.potential
    actual_phi_left = u[1][1]
    actual_phi_right = u[end][1]
    residual[1] = actual_phi_left  - UF.ustrip(U_POT, wanted_phi_left)
    residual[2] = actual_phi_right - UF.ustrip(U_POT, wanted_phi_right)
    residual
end

@noinline function step!(du, u, w, r0)
    # @show du
    # @show u
    # @show w
    # @show r0
    # error()
    r = r0*U_LENGTH
    phi0, E0 = u
    E = E0 * U_E
    rho = w.charge_density(r)
    dphi = -E0
    dEdr = compute_dEdr(E, rho, r, w)
    du[1] = dphi
    du[2] = UF.ustrip(U_E/U_LENGTH, dEdr)
    du
end

function compute_dEdr(E, rho, r, w)
    shape = w.shape
    ε = w.electric_permitivity
    ret = rho/ε
    if shape == :parallel
        ret -= 0E/r
    elseif shape == :cylindrical
        ret -= 1E/r
    elseif shape == :spherical
        ret -= 2E/r
    else
        error("This error should be unreachable. Unknown shape $shape")
    end
    ret
end

# count = 0
function initial_guess(w, r_strip)
    # @show r_strip
    r = U_LENGTH * r_strip
    phi_left, phi_right = w.potential
    r_left, r_right = w.limits
    E = (phi_right - phi_left) / (r_right - r_left)
    w_left = (r_right - r) / (r_right - r_left)
    w_right = (r - r_left) / (r_right - r_left)
    phi = w_left * phi_left + w_right * phi_right
    # global count += 1
    # @show count
    # if mod(count, 10000) == 0
    #     error()
    # end
    [UF.ustrip(U_POT, phi), UF.ustrip(U_E,E)]
end

function build_problem(w::Problem)
    r1, r2 = w.limits
    rspan = (UF.ustrip(U_LENGTH, r1)::Float64, UF.ustrip(U_LENGTH,r2)::Float64)
    ret = BV.BVProblem(step!, bc!, initial_guess, rspan, w)
    # ret = BV.TwoPointBVProblem(step!, bc!, initial_guess(w, w.limits[1]/U_LENGTH), rspan, w)
end

function solve(w::Problem,
        solver=BV.Shooting(Vern7()),
        ; kw...)
    prob = build_problem(w)
    sol = BV.solve(prob, solver; kw...)

    itp = sol.(UF.ustrip.(U_LENGTH, w.radii))
    potential = first.(itp) .* U_POT
    E = last.(itp) .* U_E
    return (E=E,
        potential=potential,
        rho=w.charge_density.(w.radii),
        radii=w.radii,
        raw=sol)
end

end # module
