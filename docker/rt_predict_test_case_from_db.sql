-- ##### SETUP RT PREDICT TEST DATA ######

/*
This is run against a copy of the meta_atlas mysql database.
Usually this is done by spinning up a mysql container
and loading an all database dump:

mysql password is on NERSC at
/global/cfs/cdirs/metatlas/mysql_user.txt

conda install -c ostrokach mysql-client=5.7.10
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
mysql2sqlite -f /docker-entrypoint-initdb.d/meta_atlas.sqlite3 -d meta_atlas -p -u root
exit


*/

SET GLOBAL sql_mode=(SELECT REPLACE(@@sql_mode,'ONLY_FULL_GROUP_BY',''));

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

DELETE l
FROM lcmsruns AS l
LEFT JOIN (
        SELECT unique_id
        FROM lcmsruns AS l1
        JOIN (
                SELECT MAX(creation_time) AS ctime, hdf5_file
                FROM lcmsruns
                WHERE (name LIKE '20201106\_JGI-AK\_PS-KM\_505892\_OakGall\_final\_QE-HF\_HILICZ\_USHXG01583\_POS\_MSMS\_0\_QC\_P%')
                GROUP BY hdf5_file
        ) AS early
        ON l1.creation_time=early.ctime AND l1.hdf5_file=early.hdf5_file
) AS j
ON l.unique_id=j.unique_id
WHERE j.unique_id is NULL;


-- delete atlases that are not the most recent version of the atlases listed
DELETE FROM atlases
WHERE unique_id NOT IN (
	SELECT error_fix.unique_id FROM (  -- extra SELECT required - https://stackoverflow.com/questions/45494/
		SELECT first.unique_id
		FROM atlases AS first
		JOIN (
			SELECT unique_id, max(last_modified) as last_modified
			FROM atlases
			WHERE name IN (
				'HILICz150_ANT20190824_TPL_EMA_Unlab_POS',
				'HILICz150_ANT20190824_TPL_QCv3_Unlab_POS',
				'HILICz150_ANT20190824_TPL_ISv5_Unlab_POS',
				'HILICz150_ANT20190824_TPL_ISv5_13C15N_POS',
				'HILICz150_ANT20190824_TPL_IS_LabUnlab2_POS',
				'HILICz150_ANT20190824_TPL_EMA_Unlab_NEG',
				'HILICz150_ANT20190824_TPL_QCv3_Unlab_NEG',
				'HILICz150_ANT20190824_TPL_ISv5_Unlab_NEG',
				'HILICz150_ANT20190824_TPL_ISv5_13C15N_NEG',
				'HILICz150_ANT20190824_TPL_IS_LabUnlab2_NEG') AND
				username = 'vrsingan'
			GROUP BY name
		) AS newest
		ON first.unique_id = newest.unique_id
	) as error_fix
);


-- delete compounds not in the atlas list
DELETE FROM compounds
WHERE unique_id NOT IN (
	SELECT error_fix.unique_id FROM (  -- extra SELECT required - https://stackoverflow.com/questions/45494/
		SELECT c.unique_id
		FROM atlases_compound_identifications AS aci
		JOIN compoundidentifications AS ci
		ON aci.target_id = ci.unique_id
		JOIN compoundidentifications_compound AS cic
		ON ci.unique_id = cic.source_id
		JOIN compounds as c
		ON cic.target_id = c.unique_id
		WHERE aci.source_id IN (
			SELECT first.unique_id
			FROM atlases AS first
			JOIN (
				SELECT unique_id, max(last_modified) as last_modified
				FROM atlases
				WHERE name IN (
					'HILICz150_ANT20190824_TPL_EMA_Unlab_POS',
					'HILICz150_ANT20190824_TPL_QCv3_Unlab_POS',
					'HILICz150_ANT20190824_TPL_ISv5_Unlab_POS',
					'HILICz150_ANT20190824_TPL_ISv5_13C15N_POS',
					'HILICz150_ANT20190824_TPL_IS_LabUnlab2_POS',
					'HILICz150_ANT20190824_TPL_EMA_Unlab_NEG',
					'HILICz150_ANT20190824_TPL_QCv3_Unlab_NEG',
					'HILICz150_ANT20190824_TPL_ISv5_Unlab_NEG',
					'HILICz150_ANT20190824_TPL_ISv5_13C15N_NEG',
					'HILICz150_ANT20190824_TPL_IS_LabUnlab2_NEG') AND
					username = 'vrsingan'
				GROUP BY name
			) AS newest
			ON first.unique_id = newest.unique_id
		)
	) AS error_fix
);

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
