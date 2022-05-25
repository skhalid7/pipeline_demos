apt-get update
apt-get install build essential
apt-get -y install apt-transport-https ca-certificates curl gnupg2 software-properties-common tabix
gsutil cp gs://shareef-dev/annotation_dependencies/homo_sapiens_ensembl_101.tar.gz .
tar -xvf homo_sapiens_ensembl_101.tar.gz &
apt-get install -y unzip
gsutil cp gs://shareef-dev/init_scripts/vep101.zip . #Get Ensembl101
unzip vep101.zip
#VEP dependencies
apt-get install -y cpanminus
cpanm Archive::Zip
cpanm DBD::mysql
cpanm DBI
#Install VEP
cd /ensembl-vep-release-101/
perl INSTALL.pl -a a -n #Install only api and don't check for updates
cd ..
#Download fasta files
gsutil cp gs://shareef-dev/annotation_dependencies/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz* .
