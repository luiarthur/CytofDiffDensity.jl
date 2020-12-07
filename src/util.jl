s3sync(; from, to, tags=``) = run(`aws s3 sync $(from) $(to) $(tags)`)

"""
Redirects stdout to a file.
"""
function redirect_stdout_to_file(f::Function, path::String)
  open(path, "w") do io
    redirect_stdout(io) do
      f()
    end
  end
end
