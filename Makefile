SHELL = /bin/bash

.PHONY: getdata all

version = 1.0.0
data-url = https://github.com/luiarthur/CytofDiffDensityData/archive/v$(version).zip

all: unzip-data

getdata: data/v$(version).zip

data/v$(version).zip:
	cd data && wget $(data-url)

unzip-data: getdata
	cd data && unzip v$(version).zip && mv CytofDiffDensityData-$(version)/data/* .
	rm -rf data/CytofDiffDensityData* data/*.zip
