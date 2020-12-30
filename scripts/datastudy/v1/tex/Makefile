SHELL = /bin/bash

awsbucket = `julia -e 'include("../../../info.jl"); println(Info.awsbucket_datastudy)'`
resultsdir = `julia -e 'include("../../../info.jl"); println(Info.resultsdir_datastudy)'`
simname = "v1"

link: metrics
	ln -s $(resultsdir)/$(simname) results

metrics:
	python3 scripts/elbo.py

