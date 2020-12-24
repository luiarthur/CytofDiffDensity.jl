SHELL = /bin/bash

awsbucket = `julia -e 'include("../../../info.jl"); println(Info.awsbucket_simstudy)'`
resultsdir = `julia -e 'include("../../../info.jl"); println(Info.resultsdir_simstudy)'`
simname = "v3"

link:
	ln -s $(resultsdir)/$(simname) results

