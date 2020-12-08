"""
Subsample a dataframe or matrix.

`df`: A dataframe of matrix (rows is number of observations).
`n::Integer`: Size of subset. If `n >= size(df, 1)`, then the original dataset
              is returned. Otherwise, a dataset of size `(n, size(df, 2)` is
              returned.

**NOTE**: requires `DataFrames.jl`
"""
function subsample(df, n::Integer; seed::Integer=nothing)
  seed == nothing || Random.seed!(seed)

  N = size(df, 1)
  if n >= size(df, 1)
    return df
  else
    return df[sample(1:N, n), :]
  end
end


"""
Partition `data::DataFrame` into control and treatment groups according to
`treatmentsym::Symbol` (whether missing(->control) or not(->treatment)),
for a given marker (`markersym::Symbol`).

**NOTE**: requires `DataFrames.jl`
"""
function partition(data, markersym::Symbol; treatmentsym::Symbol=:treatment,
                   log_response=false)
  marker = data[!, markersym]
  yC = marker[ismissing.(data[!, treatmentsym])]
  yT = marker[.!ismissing.(data[!, treatmentsym])]
  if log_response
    yC .= log.(yC)
    yT .= log.(yT)
  end
  return (yC=yC, yT=yT)
end
