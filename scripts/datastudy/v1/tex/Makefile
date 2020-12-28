SHELL = /bin/bash

awsbucket = `julia -e 'include("../../../info.jl"); println(Info.awsbucket_datastudy)'`
resultsdir = `julia -e 'include("../../../info.jl"); println(Info.resultsdir_datastudy)'`
simname = "v1"

link:
	ln -s $(resultsdir)/$(simname) results

