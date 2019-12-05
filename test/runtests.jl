@info "@time using GaussLaw1D"
@time using GaussLaw1D
using Test

@testset "Smoketest $shape" for shape in [:parallel, :cylindrical, :spherical]
    prob = Problem(Constant(1.0nC/(mm^3)),
        limits=(1.0mm, 2.0mm),
        potential=(0.0V, 100.0V),
        shape=shape,
    )
    @info "@time solve $shape problem"
    sol = @time solve(prob)
end

@testset "Manufactured solution $shape" for shape in [:parallel, :cylindrical, :spherical]
    phi_true = r-> 1/2*(r/mm)^2*V
    E_true = r -> -(r/mm)*(V/mm)
    ε = GaussLaw1D.VACUUM_PERMITIVITY
    n_div_dims = if shape == :parallel
        1
    elseif shape == :cylindrical
        2
    elseif shape == :spherical
        3
    end
    rho_true = r -> -n_div_dims*ε * (V*mm^-2)
    
    r1 = 1.0mm
    r2 = 2.0mm
    limits = (1.0mm, 2.0mm)
    p = Problem(rho_true,
        limits=(r1,r2),
        potential=phi_true.(limits),
        shape=shape,
    )
    
    sol = solve(p)
    @test E_true.(p.radii) ≈ sol.E
    @test phi_true.(p.radii) ≈ sol.potential
    @test rho_true.(p.radii) ≈ sol.rho
end
