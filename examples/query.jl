using AdaptiveHierarchicalRegularBinning

d=2
nc = 1000
nq = 100

dims=2

C = rand(d, nc)
C[1,:] .= 1:nc
Q = rand(d, nq)
Q[1,:] .= 0 .- (1:nq)

Q[:, 10:30] .= C[:, 30:50]


if dims == 1
  C = C |> transpose |> collect
  Q = Q |> transpose |> collect
end

X = cat(C, Q, dims=dims)

labels = zeros(Bool, nc+nq) # 0 -> C, 1 -> Q
labels[nc+1:end] .= true

tree = ahrb!(X, 5, 1000; dims=dims, QT=UInt128)

labels = labels[tree.info.perm]


function query(tree, labels)
  r = zeros(Bool, length(labels))

  for leaf in Leaves(tree)
    lr      = @view r[range(leaf)]
    llabels = @view labels[range(leaf)]
    lpoints = points(leaf)
    lcorpus = selectdim(lpoints, leaddim(leaf), .!llabels)

    for (i, (label, lpoint)) in enumerate(zip(llabels, eachslice(lpoints; dims=leaddim(leaf))))
      label || continue

      lr[i] = lpoint in eachslice(lcorpus; dims=leaddim(leaf))
    end

  end

  return r
end


el = query(tree, labels)
