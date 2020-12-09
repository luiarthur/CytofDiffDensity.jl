module Info

project_name = "cytof-diff-dens"

awsbucket = "s3://$(project_name)"
resultsdir = joinpath(ENV["SCRATCH_DIR"], project_name)

awsbucket_datastudy = joinpath(awsbucket, "datastudy")
resultsdir_datastudy = joinpath(resultsdir, "datastudy")

awsbucket_simstudy = joinpath(awsbucket, "simstudy")
resultsdir_simstudy = joinpath(resultsdir, "simstudy")

awsbucket_compare = joinpath(awsbucket, "compare")
resultsdir_compare = joinpath(resultsdir, "compare")

end  # module
