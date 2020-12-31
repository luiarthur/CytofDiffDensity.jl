SHELL = /bin/bash

awsbucket = `julia -e 'include("../../info.jl"); println(Info.awsbucket_simstudy)'`
resultsdir = `julia -e 'include("../../info.jl"); println(Info.resultsdir_simstudy)'`
simname = "v4"

getresults:
	aws s3 sync $(awsbucket)/$(simname) $(resultsdir)/$(simname) --exclude '*.bson'

getallresults:
	aws s3 sync $(awsbucket)/$(simname) $(resultsdir)/$(simname)

sendresults:
	aws s3 sync $(resultsdir)/$(simname) $(awsbucket)/$(simname)
