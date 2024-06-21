using Pkg

Pkg.update()
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("LaTeXStrings")
Pkg.add("LinearAlgebra")
Pkg.add("Plots")
Pkg.add("Printf")
Pkg.add("Symbolics")
Pkg.add(url="https://github.com/jmanthony3/NumericalMethods.jl.git")


# if running WSL, then add PyCall on localhost NOT from within WSL
# Pkg.add("PyCall")
# Pkg.build()
