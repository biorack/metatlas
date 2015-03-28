
/*
A web-based atlas to liquid chromatography–mass spectrometry (LCMS) data
*/

module MetaboliteAtlas2 {


    /* Meta data associated with an object stored in a workspace.

        object_id id - ID of the object assigned by the user or retreived from the IDserver (e.g. kb|g.0)
        object_type type - type of the object (e.g. Genome)
        timestamp moddate - date when the object was modified by the user (e.g. 2012-12-17T23:24:06)
        int instance - instance of the object, which is equal to the number of times the user has overwritten the object
        timestamp date_created - time at which the alignment was built/loaded in seconds since the epoch
        string command - name of the command last used to modify or create the object
        username lastmodifier - name of the user who last modified the object
        username owner - name of the user who owns (who created) this object
        workspace_id workspace - ID of the workspace in which the object is currently stored
        workspace_ref ref - a 36 character ID that provides permanent undeniable access to this specific instance of this object
        string chsum - checksum of the associated data object
        string metadata - custom metadata entered for data object during save operation

    */

    /* A string identifier for a workspace. Any string consisting of alphanumeric characters and "-" is acceptable. */
    typedef string workspace_id;

    /* A string indicating the type of an object stored in a workspace. */
    typedef string object_type;

    /* ID of an object stored in the workspace. */
    typedef string object_id;

    /* Login name of KBase user account to which permissions for workspaces are mapped */
    typedef string username;

    /* Exact time for workspace operations. e.g. 2012-12-17T23:24:06 */
    typedef string timestamp;

    /* A permanent reference to an object in a workspace. */
    typedef string workspace_ref;

    /* A permanent reference to an Metabolite Atlas Dictionary.
    @id ws MetaboliteAtlas.MetAtlasDictionary
    */
    typedef string atlas_id;

    typedef structure {
        atlas_id atlas_id;
    } AtlasRef;


    typedef tuple<object_id id, object_type type, timestamp moddate, int instance, string command, username lastmodifier, username owner, workspace_id workspace, workspace_ref ref, string chsum, string metadata> object_metadata;

    /* ************************************************************************************* */
    /* Metabolite Atlas Data Types */
    /* ************************************************************************************* */

    /* Metabolite dictionary object
     *
     *  name - common name of the compound
     *  formula - chemical formula
     *  adducts - adduct ions
     *  mz - mass-to-charge ratio
     *  mz_threshold - threshold in ppm
     *  rt_min - min retention time
     *  rt_max - max retention time
     *  rt_peak - peak retention time
     *
     */

    /* MetAtlas_FileInfo dictionary object
     *
     *  index - arbitrary order
     *  status - 0 means loaded or finished converting; 1 means pending or just uploaded; 2 means loading or in process of converting
     *  name - This is the name and location of the file where it was uploaded to
     *  fid - This is the index in the SciDB array that identifies all the datapoints in your file
     *  polarity - User gets to Specify whether the file was acquired in positive or negative ionization mode
     *  group - A grouping descriptor to be used for analysis. Things like "control" and "treatment" for each file will allow those files to be grouped and processed by later code blocks
     *  inclusion_order - numerical values here to configure your plots to be in the order you want
     *  normalization_factor - A factor to divide the intensity of observed molecules by. This can partially correct for artifacts introduced by loading and extraction inconsistencies. Use with care and avoid by careful experimental design.
     *  retention_correction - A number (in minutes) to shift the retention time by
     *
     */

    typedef structure {
        string name;
        string formula;
        string adducts;
        string mz;
        string mz_threshold;
        string rt_min;
        string rt_max;
        string rt_peak;
        string neutral_mass;
        string pubchem_id;
        string dict_id;
        string creator;
        string creation_date;
    } Compound;


    typedef list<Compound> CompoundList;

    /* note: we might want to use the @optional annotation to flag optional fields
        --Mike Sneddon
    */
    typedef structure {
        string name;
        CompoundList compounds;
        string creator;
        string creation_date;
        string sample;
        string method;
    } MetAtlasDictionary;

    typedef structure {
        string index;
        string status;
        string name;
        string fid;
        list<Atlas> atl;
        string polarity;
        string group;
        string inclusion_order;
        string normalization_factor;
        string retention_correction;
    } MetAtlasFileInfo;


    typedef list<FileInfo> MetAtasFileInfoList;

    typedef structure {
        string name;
        FileInfoList file_info_list;
        string creator;
        string creation_date;
    } MetAtlasFileInfoCollection;


    /* ************************************************************************************* */
    /* FUNCTION DEFINITIONS */
    /* ************************************************************************************* */

    /* All functions require authentication. */
    authentication required;

    /* Input parameters for the "load dictionary" function. */
    typedef structure {
        string workspace;
        string username;
        string dict_id;
    } LoadDictionary;

    /*
        Load a metabolite dictionary given a dictionary id and credentials.  Results are stored in a Dictionary object.  Returns the
        metadata for the metabolite dictionary object.
    */
    funcdef loadDictionary(LoadDictionary params) returns(object_metadata output);


};
