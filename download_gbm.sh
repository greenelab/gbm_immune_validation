# Gregory Way 2016 - GBM Immune Profiles
# Raw data downloaded on 16 Aug 2016 and deposited into zenodo for versioning

# All files will be downloaded to the `data/` folder
mkdir --parents 'data/' 

# Download data from UCSC Xena (http://xena.ucsc.edu/)

# Download GBM clinical matrix
# url=https://tcga.xenahubs.net/download/TCGA.GBM.sampleMap/GBM_clinicalMatrix
# wget --directory-prefix 'data/' $url

# Download Affy U133A microarray data
# url=https://tcga.xenahubs.net/download/TCGA.GBM.sampleMap/HT_HG-U133A
# wget --directory-prefix 'data/' $url

# Download zenodo tar.gz file
# Note - the zenodo file includes the files downloaded above 
url=https://zenodo.org/record/60304/files/gbm_immune_data.tar.gz
curl $url | tar --extract --gzip --strip-components=1 --directory=data

