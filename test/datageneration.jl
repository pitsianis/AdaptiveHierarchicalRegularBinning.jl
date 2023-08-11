"""
# generate an n-point cloud in d-dimensions that is tall and thin
"""
function lineartreedata(n, d, maxL)
X = zeros(d, n)
r = [zeros(d) ones(d)]

  function helper(r, l, k)
    if k == 0
      return
    end
    s = diff(r, dims=2)
    if l == maxL || k == n
      X[:, k:n] .= rand(d, n-k+1) .* s .+ r[:,1]
    else
      X[:, k] = rand(d) .* s .+ r[:, 1]
      msk = X[:, k] .<= sum(r, dims=2) ./ 2
      r .= [r[:, 1] .+ msk .* s ./ 2    r[:, 2] .- (.!msk) .* s ./ 2]
      helper(r, l+1, k+1)
    end
  end

helper(r, 1, 1)

return 10 .* X .- 5
end

second((x, y)) = y
third((x, y, z)) = z
"""
# generate 2^(maxL) per side regular grid of points in 2 and 3 dimensions that form a full tree
"""
function fulltree(maxL, d)
  r = collect(range(0,1,2^maxL+1))[1:end-1]
  if d == 1
    X = collect(r')
  elseif d == 2
    P = [(x, y) for x in r, y in r]
    X = [first.(P[:])'; second.(P[:])']
  elseif d == 3
    P = [(x, y, z) for x in r, y in r, z in r]
    X = [first.(P[:])'; second.(P[:])'; third.(P[:])']
  else
    error("d must be in 1, 2 or 3")
  end
  return X
end
