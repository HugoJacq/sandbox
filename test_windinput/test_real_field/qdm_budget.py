"""
We compute the global momentum budget for a forced simulation
"""

# watch out for the pressure. Phi in the code is the NH pressure, I have to add back the hydrostatic pressure.

# Along X

# du/dt + grad(uu) = -1/rho grad(p) - 1/rho grad(ps) + nu lapla(vect(u))
