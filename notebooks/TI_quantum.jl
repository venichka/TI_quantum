module TI_quantum

using BigCombinatorics
using Nemo

export path, heaviside, u, v, c_l_hyper, Factorial_1, Binomial_1

# function definitions go here

"""
    TI_quantum.path()

Returns the paths for data and figures, depending on the home directory of the user. If the home directory is "C:\\Users\\nnefedkin", it returns the paths for Windows, if it's "/home/nikita", it returns the paths for Linux, and if it's "/Users/jimi", it returns the paths for MacOS.

# Returns:
- PATH_FIGS: The path for the figures.
- PATH_DATA: The path for the data.
"""
function path()
    home = homedir()
    if home == "C:\\Users\\nnefedkin"
        PATH_FIGS = "D:/nnefedkin/Google_Drive/Work/In process/Projects/Time_interfaces/Figs/"
        PATH_DATA = "D:/nnefedkin/Google_Drive/Work/In process/Projects/Time_interfaces/Data/"
    elseif home == "/home/nikita"
        PATH_FIGS = "/home/nikita/Documents/figs/time_interfaces/"
        PATH_DATA = "/home/nikita/Documents/data/time_interfaces/"
    elseif home == "/Users/jimi"
        PATH_FIGS = "/Users/jimi/Google Drive/Work/In process/Projects/Time_interfaces/Figs/"
        PATH_DATA = "/Users/jimi/Google Drive/Work/In process/Projects/Time_interfaces/Figs/"
    end
    return PATH_FIGS, PATH_DATA
end


"""
    TI_quantum.heaviside(x)

Computes the Heaviside step function of the input argument `x`, defined as 0.5*(sign(x) + 1).

# Arguments:
- x: The input argument.

# Returns:
- The value of the Heaviside step function of `x`.
"""
function heaviside(x)
    0.5*(sign.(x) .+ 1)
end

"""
    TI_quantum.u(r)

Computes the coefficient `u` in the Bogoliubov transformation for a single time interface, given the input argument `r`.

# Arguments:
- r: The input argument.

# Returns:
- The value of the coefficient `u`.
"""
function u(r::Real)
    0.5*(1.0+r)/sqrt(r)
end


"""
    TI_quantum.v(r)

Computes the coefficient `v` in the Bogoliubov transformation for a single time interface, given the input argument `r`.

# Arguments:
- r: The input argument.

# Returns:
- The value of the coefficient `v`.
"""
function v(r::Real)
    0.5*(-1.0+r)/sqrt(r)
end


"""
    TI_quantum.Factorial_1(n::Integer)::BigInt

# Description:
This function calculates the factorial of a non-negative integer `n`, returning the result as a `BigInt`. It uses memoization to save previously calculated factorials, making subsequent calls to the function faster. It also checks that the argument `n` is non-negative and throws a `DomainError` if it is not.

# See also:

    https://oeis.org/A000142
    Julia's factorial function
    FallingFactorial and RisingFactorial functions in Julia

# Important
- **Initialize before use:** `Factorial_1()`

# Arguments:
- n

# Returns
- n!
"""
function Factorial_1(n::Integer)::BigInt
    if n < 0
        throw(DomainError(n, "argument must be nonngative"))
    end
    if BigCombinatorics._has(Factorial_1, n)
        return BigCombinatorics._get(Factorial_1, n)
    end

    start = BigCombinatorics._max_arg(Factorial_1) + 1
    for m = start:n
        val = m * BigCombinatorics._get(Factorial_1, m - 1)
        BigCombinatorics._save(Factorial_1, m, val)
    end

    return BigCombinatorics._get(Factorial_1, n)
end


"""
    TI_quantum.Factorial_1()

# Description:
This function initializes the memoization table used by `Factorial_1`. It saves the values of `0!` and `1!` as `BigInts` to the table.
"""
function Factorial_1()
    BigCombinatorics._make(Factorial_1, Integer)
    BigCombinatorics._save(Factorial_1, 0, big(1))
    BigCombinatorics._save(Factorial_1, 1, big(1))
end
Factorial_1()


"""
    TI_quantum.Binomial_1(n::Integer, k::Integer)::BigInt

# Description:
This function calculates the binomial coefficient "n choose k", returning the result as a `BigInt`. It uses the memoized `Factorial_1` function to calculate the factorials necessary for the calculation. It requires that `0 <= k <= n`.
"""
function Binomial_1(n::Integer, k::Integer)::BigInt
    return (Factorial_1(n) / (Factorial_1(n-k) * Factorial_1(k)))
end

"""
    TI_quantum.c_l_hyper(n_1, n_2, l, r)

Computes the coefficients in the expansion of the output wavefunction over the Fock states after the time interface, using the hypergeometric function.

# Arguments:
- n_1: The number of particles in the first mode.
- n_2: The number of particles in the second mode.
- l: The difference between the number of particles in the first and second modes.
- r: The input argument.

# Returns:
- The coefficients in the expansion of the output wavefunction over the Fock states after the time interface, using the hypergeometric function.
"""
function c_l_hyper(n_1::Int, n_2::Int, l::Int, r::Real)
    u_0 = big(u(r))
    v_0 = big(v(r))
    (
    sign(u_0) / u_0 * conj(u_0)^n_1 * (- (v_0 / u_0))^l * 
    sqrt((Factorial_1(n_1 + l) / Factorial_1(n_1)) * 
         (Factorial_1(n_2 + l) / Factorial_1(n_2))) *
    u_0^(-n_2) *  
        hypergeometric_2f1(CC(-n_1), CC(n_2 + l + 1), CC(l + 1), 
            CC(abs(v_0)^2 / abs(u_0)^2); flags=1)
    )
end

end # module
