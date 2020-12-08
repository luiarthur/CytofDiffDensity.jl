SHELL = /bin/bash

.PHONY: getdata

version = 1.0.0
data-url = https://github.com/luiarthur/CytofDiffDensityData/archive/v$(version).zip

getdata: unzip-data

unzip-data: download-data
	cd data && unzip v$(version).zip && mv CytofDiffDensityData-$(version)/data/* .
	rm -rf data/CytofDiffDensityData* data/*.zip

download-data: data/v$(version).zip

data/v$(version).zip:
	cd data && wget $(data-url)
