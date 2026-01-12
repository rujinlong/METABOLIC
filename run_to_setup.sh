tar zxvf Accessory_scripts.tgz; rm Accessory_scripts.tgz
tar zxvf METABOLIC_hmm_db.tgz;  rm METABOLIC_hmm_db.tgz
tar zxvf METABOLIC_template_and_database.tgz; rm  METABOLIC_template_and_database.tgz
tar zxvf Motif.tgz; rm Motif.tgz
mkdir kofam_database  
cd kofam_database  
wget -c -O ko_list.gz "https://www.genome.jp/ftp/db/kofam/ko_list.gz"
wget -c -O profiles.tar.gz "https://www.genome.jp/ftp/db/kofam/profiles.tar.gz"
gzip --quiet -d ko_list.gz  
tar xzf profiles.tar.gz; rm profiles.tar.gz  
mv ../All_Module_KO_ids.txt profiles

# Create merged HMM databases (significantly reduces file count for HPC compatibility)
# Full database: all prokaryotic KOs
perl ../Accessory_scripts/merge_hmm_profiles.pl profiles kofam_all.hmm profiles/prokaryote.hal
hmmpress kofam_all.hmm

# Small database: only module-related KOs (for -kofam-db small option)
perl ../Accessory_scripts/merge_hmm_profiles.pl profiles kofam_small.hmm profiles/All_Module_KO_ids.txt
hmmpress kofam_small.hmm

cd ../
mkdir dbCAN2
cd dbCAN2
wget https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-old@UGA/dbCAN-fam-HMMs.txt.v10  -O dbCAN-fam-HMMs.txt
perl ../Accessory_scripts/batch_hmmpress_for_dbCAN2_HMMdb.pl
cd ../
mkdir MEROPS
cd MEROPS
wget https://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib
perl ../Accessory_scripts/make_pepunit_db.pl
cd ../
wget -c https://figshare.com/ndownloader/files/43500597 -O METABOLIC_test_files.tgz
tar zxvf METABOLIC_test_files.tgz
rm *.tgz
