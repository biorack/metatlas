-- Use this file to generate a near-minimal database for testing
-- to generate an sqlite3 database run:
-- mysql_to_sqlite_filtered.sh c18_test_case_from_db.sql <mysql_password>

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
		WHERE name LIKE '20210915\_JGI-AK\_MK\_506588\_SoilWaterRep\_final\_QE-HF\_C18\_USDAY63680\_NEG\_MSMS%Run2__.mzML'
		      OR
                      name in (
                        '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_FPS_MS1_0_QC_Post_Rg70to1050-CE102040--QC_Run307.mzML',
                        '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_FPS_MS1_0_QC_Pre_Rg70to1050-CE102040--QC_Run6.mzML',
                        '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_NEG_MSMS_0_QC_Post_Rg70to1050-CE102040--QC_Run309.mzML',
                        '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_NEG_MSMS_0_QC_Pre_Rg70to1050-CE102040--QC_Run8.mzML',
                        '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_NEG_MSMS_53_Cone-S1_5_Rg70to1050-CE102040-QlobataAkingi-S1_Run188.mzML',
                        '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_NEG_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run41.mzML',
                        '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_NEG_MSMS_58_Cone-S2_2_Rg70to1050-CE102040-QlobataAkingi-S1_Run56.mzML',
                        '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_NEG_MSMS_59_Cone-S2_3_Rg70to1050-CE102040-QlobataAkingi-S1_Run87.mzML',
                        '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_0_QC_Post_Rg70to1050-CE102040--QC_Run308.mzML',
                        '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_0_QC_Pre_Rg70to1050-CE102040--QC_Run7.mzML',
                        '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.mzML',
                        '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_53_Cone-S1_5_Rg70to1050-CE102040-QlobataAkingi-S1_Run187.mzML',
                        '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_54_Cone-S1_6_Rg70to1050-CE102040-QlobataAkingi-S1_Run221.mzML',
                        '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.mzML',
                        '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_65_Cone-S3_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run16.mzML',
                        '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.mzML'
		      )
		GROUP BY hdf5_file
	) AS early
	ON l1.creation_time=early.ctime AND l1.hdf5_file=early.hdf5_file
) AS j
ON l.unique_id=j.unique_id
WHERE j.unique_id is NULL;

DELETE FROM atlases
WHERE unique_id NOT IN (
	'f74a731c590544aba5c3720b346e508e',
	'19b4c10e304246cbbbe5fd3574770e5d',
	'322ed4c5fabe49349bcbc2857fbcd0dc',
	'4b05837a53494dd8b680e6b5059e1934',
	'5b77242ad9c04e76a745e51e9d6a3c4b',
	'669b750765634159a7f16645e6cf7758',
	'89694aa326cd46958d38d8e9066de16c',
	'a5f7bc81caa94853bbd6ee4b44e09187',
	'c6db576b879043768125c4e03e6a8f6e',
	'db58154082824230be4f7fee93e4ebd9',
	'e299c951fc8b48ea82524a6c9615f418',
	'e7fba1813272439498405436a28b90b2'
);

UPDATE atlases
SET username='root'
WHERE unique_id='f74a731c590544aba5c3720b346e508e';

UPDATE atlases
SET name='HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_root0'
WHERE unique_id='4b05837a53494dd8b680e6b5059e1934';

UPDATE atlases
SET username='root'
WHERE unique_id='4b05837a53494dd8b680e6b5059e1934';

DELETE FROM compounds
WHERE inchi_key NOT IN (
	'OIRDTQYFTABQOQ-KQYNXXCUSA-N',
	'OLXZPDWKRNYJJZ-RRKCRQDMSA-N',
	'LRFVTYWOQMYALW-UHFFFAOYSA-N',
	'GFFGJBXGBJISGV-UHFFFAOYSA-N',
	'HXACOUQIXZGNBF-UHFFFAOYSA-N',
	'CZMRCDWAGMRECN-UGDNZRGBSA-N',
	'ISAKRJDGNUQOIC-UHFFFAOYSA-N',
	'OPTASPLRGRRNAP-UHFFFAOYSA-N',
	'BDJRBEYXGGNYIS-UHFFFAOYSA-N',
	'FBBCSYADXYILEH-UHFFFAOYSA-N',
	'OMZRMXULWNMRAE-BMRADRMJSA-N'
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
