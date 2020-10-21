groups $1
mkdir /global/project/projectdirs/metatlas/raw_data/$1
chgrp metatlas /global/project/projectdirs/metatlas/raw_data/$1
chmod 770 -R /global/project/projectdirs/metatlas/raw_data/$1
chmod g+s /global/project/projectdirs/metatlas/raw_data/$1

