-- ##### SECOND SETUP OF TEST DATA ######

/*
This is run against a copy of the meta_atlas mysql database.
Usually this is done by spinning up a mysql container
and loading an all database dump:

mkdir -p dumps
mysqldump -h nerscdb04.nersc.gov -u meta_atlas_admin -p --all-databases > dumps/meta_atlas.sql
docker run -it --rm -e MYSQL_ROOT_PASSWORD=mypw -v $(pwd)/dumps:/docker-entrypoint-initdb.d -v $(pwd):/script mysql:5.7
MYSQL_ID=$(docker ps | grep mysql | cut -f1 -d' ')
docker exec -it $MYSQL_ID /bin/bash

# then within the mysql container's shell:
mysql --password=mypw meta_atlas < /script/extract_test_case_from_db.sql

From within the atlasdb container
apt-get update
apt-get install -y python3 python3-pip
pip3 install mysql-to-sqlite3
mysql2sqlite -f /meta_atlas.sqlite3 -d meta_atlas -p -u root
exit

# back on the host computer you can now copy out the sqlite3 db file
docker cp $MYSQL_ID:/meta_atlas.sqlite3 .
*/

-- results in a database with a 1 atlas, 4 lcmsruns, and 6 compounds
-- the following tables are not modified:
--    mzintensityhpairs
--    fragmentationreferences_mz_intensities
--    compoundidentifications_frag_references
--    fragmentationreferences


-- remove tables that are not used
DROP TABLE IF EXISTS `Compounds`;
DROP TABLE IF EXISTS `Group`;
DROP TABLE IF EXISTS `group`;

-- clean out tables we don't need pre-populated values in
DELETE FROM groups;
DELETE FROM groups_items;
DELETE FROM methods;
DELETE FROM samples;
DELETE FROM mzintensitypairs;
DELETE FROM identificationgrades;
DELETE FROM functionalsets;
DELETE FROM fragmentationreferences_mz_intensities;
DELETE FROM compoundidentifications_frag_references;
DELETE FROM fragmentationreferences;

-- Modify atlas with name 'MSMLS_HILICz150mm_Annot20190824_Predicted_EMA_Unlab_POS_Shih_OakGall_505892_final_20210210.csv'
DELETE FROM atlases
WHERE unique_id!='4b05837a53494dd8b680e6b5059e1934';

UPDATE atlases
SET name='HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_root0'
WHERE unique_id='4b05837a53494dd8b680e6b5059e1934';

UPDATE atlases
SET username='root'
WHERE unique_id='4b05837a53494dd8b680e6b5059e1934';

DELETE l
FROM lcmsruns AS l
LEFT JOIN (
	SELECT unique_id
	FROM lcmsruns AS l1
	JOIN (
		SELECT MAX(creation_time) AS ctime, hdf5_file
		FROM lcmsruns
		WHERE (name LIKE '20201106\_JGI-AK\_PS-KM\_505892\_OakGall\_final\_QE-HF\_HILICZ\_USHXG01583\_POS\_MSMS%') AND
		      (name LIKE '%Cone-S%\_1\_%')
		GROUP BY hdf5_file
	) AS early
	ON l1.creation_time=early.ctime AND l1.hdf5_file=early.hdf5_file
	LIMIT 4
) AS j
ON l.unique_id=j.unique_id
WHERE j.unique_id is NULL;

DELETE FROM compounds
WHERE chebi_id NOT IN ('CHEBI:17256', 'CHEBI:16708', 'CHEBI:16708', 'CHEBI:48517///CHEBI:17712', 'CHEBI:30959///CHEBI:17405', 'CHEBI:16335');

-- work from compounds up to atlases_compound_identifications
DELETE cic
FROM compoundidentifications_compound AS cic
LEFT JOIN compounds AS c
ON cic.target_id=c.unique_id
WHERE c.unique_id is null;

DELETE ci
FROM compoundidentifications AS ci
LEFT JOIN compoundidentifications_compound AS cic
ON ci.unique_id=cic.source_id
WHERE cic.source_id is null;

DELETE aci
FROM atlases_compound_identifications AS aci
LEFT JOIN compoundidentifications AS ci
ON aci.target_id=ci.unique_id
WHERE ci.unique_id is null;

-- work from atlases_compound_identifications down to everything else
DELETE atlases_compound_identifications
FROM atlases_compound_identifications
LEFT JOIN atlases
ON atlases.unique_id=atlases_compound_identifications.source_id
WHERE atlases.unique_id is null;

DELETE compoundidentifications
FROM compoundidentifications
LEFT JOIN atlases_compound_identifications AS aci
ON aci.target_id=compoundidentifications.unique_id
WHERE aci.target_id is null;

DELETE compoundidentifications_compound
FROM compoundidentifications_compound
LEFT JOIN compoundidentifications AS ci
ON ci.unique_id=compoundidentifications_compound.head_id
WHERE ci.unique_id is null;

DELETE compoundidentifications_rt_references
FROM compoundidentifications_rt_references
LEFT JOIN compoundidentifications AS ci
ON ci.unique_id=compoundidentifications_rt_references.head_id
WHERE ci.unique_id is null;

DELETE compoundidentifications_mz_references
FROM compoundidentifications_mz_references
LEFT JOIN compoundidentifications AS ci
ON ci.unique_id=compoundidentifications_mz_references.head_id
WHERE ci.unique_id is null;

DELETE compounds
FROM compounds
LEFT JOIN compoundidentifications_compound AS cic
ON compounds.head_id=cic.target_id
WHERE cic.target_id is null;

DELETE rtreferences
FROM rtreferences
LEFT JOIN compoundidentifications_rt_references AS cirr
ON rtreferences.head_id=cirr.target_id
WHERE cirr.target_id is null;

DELETE mzreferences
FROM mzreferences
LEFT JOIN compoundidentifications_mz_references AS cimr
ON mzreferences.head_id=cimr.target_id
WHERE cimr.target_id is null;

-- SELECT table_name, table_rows FROM information_schema.TABLES WHERE table_schema='meta_atlas' ORDER BY table_rows DESC;
