wget -c http://118.31.70.55/ProGro-neo/data/dbsnp_146.tar.gz
wget -c http://118.31.70.55/ProGro-neo/data/rna.tar.gz
wget -c http://119.3.70.71/ProGeo-neo2.0/data/reference_files.tar.gz
wget -c http://119.3.70.71/ProGeo-neo2.0/data/software.tar.gz
wget -c http://119.3.70.71/ProGeo-neo2.0/data/test.tar.gz
tar -zxvf dbsnp_146.tar.gz
tar -zxvf rna.tar.gz
tar -zxvf reference_files.tar.gz
tar -zxvf software.tar.gz
tar -zxvf test.tar.gz
rm -f *.gz
mv rna test/
mv dbsnp_146 reference_files/
