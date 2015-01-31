module Polynomials

export Polynomial, deg, X, Y

abstract Poly

immutable Polynomial{T <: Union(Number, Poly)} <: Poly
  coeffs :: Vector{T}
end

Base.zero{T}(::Polynomial{T}) = zero(Polynomial{T})
Base.zero{T}(::Type{Polynomial{T}}) = Polynomial{T}([])

Base.one{T}(::Polynomial{T}) = one(Polynomial{T})
Base.one{T}(::Type{Polynomial{T}}) = Polynomial{T}([one(T)])

Base.copy{T}(p::Polynomial{T}) = Polynomial(copy(p.coeffs))

Base.convert(::Type{Polynomial}, x) = if isa(x,Polynomial) x else Polynomial([x]) end

Base.convert{R}(::Type{Polynomial{R}}, x) = Polynomial{R}([convert(R,a) for a in convert(Polynomial,x).coeffs])

Base.promote_rule{T,R}(::Type{Polynomial{T}}, ::Type{Polynomial{R}}) = Polynomial{promote_type(T,R)}

Base.promote_rule{T,R}(::Type{T}, ::Type{Polynomial{R}}) = Polynomial{promote_type(T,R)}

function +{T}(p::Polynomial{T}, q::Polynomial{T})
  n = min(length(p.coeffs), length(q.coeffs))
  Polynomial{T}([p.coeffs[1:n] .+ q.coeffs[1:n], p.coeffs[n+1:end], q.coeffs[n+1:end]])
end

function *{T}(p::Polynomial{T}, q::Polynomial{T})
  if isempty(p.coeffs)
    p
  else
    Polynomial(T[p.coeffs[1] * x for x in q.coeffs]) + Polynomial([zero(T), (Polynomial{T}(p.coeffs[2:end]) * q).coeffs])
  end
end

-(p::Polynomial) = Polynomial(map(-, p.coeffs))
-(p::Polynomial, q::Polynomial) = +(p,-q)

X = Polynomial([0,1])
Y = Polynomial([X])

+(x::Union(Number,Poly), y::Union(Number,Poly)) = +(promote(x,y)...)
-(x::Union(Number,Poly), y::Union(Number,Poly)) = -(promote(x,y)...)
*(x::Union(Number,Poly), y::Union(Number,Poly)) = *(promote(x,y)...)

iszero(n::Number) = n == zero(n)

iszero(p::Polynomial) = deg(p) == -1

==(p::Polynomial, q::Polynomial) = iszero(p-q)

==(x::Union(Number,Poly), y::Union(Number,Poly)) = ==(promote(x,y)...)

function deg(p::Polynomial)
  d = length(p.coeffs)
  while d > 0 && iszero(p.coeffs[d])
    d -= 1
  end
  d - 1
end

function Base.divrem{T}(p::Polynomial{T}, q::Polynomial{T})
  if iszero(q)
    error("Divide by zero")
  elseif deg(q) > deg(p)
    (zero(p), p)
  else
    (d,r) = divrem(p.coeffs[1+deg(p)], q.coeffs[1+deg(q)])
    (d,r) = map(x -> Polynomial([x]), (d,r))
    d = d * X^(deg(p)-deg(q))
    r = r * X^deg(p)
    (dd,rr) = divrem(p-(d*q+r), q)
    (d+dd,r+rr)
  end
end

Base.divrem(x::Union(Number,Poly), y::Union(Number,Poly)) = Base.divrem(promote(x,y)...)
Base.divrem(x::Number, y::Number) = (div(x,y), rem(x,y))

Base.div(p::Polynomial, q::Polynomial) = divrem(p,q)[1]
Base.rem(p::Polynomial, q::Polynomial) = divrem(p,q)[2]

end
