begin
  multithreading() = true

  macro threading(test)
      esc(:(if $(@__MODULE__).multithreading()
        Threads.@threads($test)
      else
        $test
      end))
  end
end


@inline function staticselectdim(A, ::Val{D}, i) where {D}
  I = ntuple(k->k==D ? i : (:), Val(max(ndims(A), D)))
  @boundscheck checkbounds(A, I...)
  return @inbounds @view A[I...]
end

# SECTION: Allocator
struct Allocator{T}
  pools::Vector{Dict{UInt, Vector{Vector{T}}}}
end

Allocator(::Type{T}) where {T} = Allocator{T}([Dict() for _ in 1:Threads.nthreads()])

function alloc!(allocator::Allocator{T}, dims::Vararg{<:Integer, N})::Array{T, N} where {T, N}
  pool = @inbounds allocator.pools[Threads.threadid()]
  len = prod(dims)

  haskey(pool, len) || return Array{T, N}(undef, dims)

  buffers = @inbounds pool[len]

  isempty(buffers) && return Array{T, N}(undef, dims)

  return reshape(pop!(buffers), dims...)::Array{T, N}
end


function free!(allocator::Allocator, X::Array{UInt})
  pool = @inbounds allocator.pools[Threads.threadid()]
  len = length(X)

  if !haskey(pool, len)
    pool[len] = []
  end

  buffers = @inbounds pool[len]
  push!(buffers, reshape(X, :))

  return
end
# !SECTION: Allocator