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
