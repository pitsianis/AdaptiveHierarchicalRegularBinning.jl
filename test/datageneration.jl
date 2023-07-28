"""
# generate an n-point cloud in d-dimensions that cannot be partitioned in L-levels
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

return X
end