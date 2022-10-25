struct Allocator{T}
  pools::Vector{Dict{Integer, Vector{Vector{T}}}}
end


Allocator(T::DataType) = Allocator{T}([Dict() for _ in 1:Threads.nthreads()])


alloc!(a::Allocator{T}, dims...) where {T} = @inbounds begin
  pool = a.pools[Threads.threadid()]
  len = prod(dims)
  if !haskey(pool, len)
    pool[len] = []
  end

  buffers = pool[len]
  if isempty(buffers)
    push!(buffers, Vector{T}(undef, len))
  end

  return reshape(pop!(buffers), dims...)
end


free!(a::Allocator{T}, X::Array{T}) where {T} = @inbounds begin
  pool = a.pools[Threads.threadid()]
  dims = size(X)
  len = prod(dims)

  if !haskey(pool, len)
    pool[len] = []
  end

  buffers = pool[len]
  push!(buffers, reshape(X, :))

  return
end
