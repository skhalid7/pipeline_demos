ROLE=$(/usr/share/google/get_metadata_value attributes/dataproc-role)
if [[ "${ROLE}" == 'Master' ]]; then 
	#install glow.py python package
	gsutil cp ./glow.py-0.6.0-py3-none-any.whl .
	pip3 install glow.py-0.6.0-py3-none-any.whl

	#install custom functions and move them to site packages folder
	gsutil cp ./cncd_pipelines.zip
	unzip cncd_pipelines.zip
	cp -r cncd_pipelines /opt/conda/anaconda/lib/python3.6/site-packages/ ##IF PYTHON VERSION HAS CHANGED UPDATE PATH HERE
fi

#Install dependencies on all nodes
#reference file
gsutil cp ./hg38.nochr.fa .
gsutil cp ./hg38.nochr.fa.fai .

#lift over chain file
gsutil cp ./b37ToHg38.over.chain .
