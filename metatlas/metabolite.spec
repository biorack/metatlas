
/*
A web-based atlas to liquid chromatography–mass spectrometry (LCMS) data
*/

module MetaboliteAtlas2 {
{

    /* id ws MetaboliteAtlas2.Run */
    typedef string Run_ref;

    /* @id handle */
    typedef string Run_data_ref;

    typedef structure {
        int fileid;
        float mz;
        float scan_time;
        float i;
        int polarity;
        int ms_level;
        float precursor_MZ;
        float precursor_intensity;
        float collision_energy;
        Run_data_ref data;
    } Run;

    typedef structure {
        list<Run_ref> StudyRuns;
        mapping<string, string> optional_metadata;
        string organisms;
        string study_design_description;
        string publications;
        string experimental_factors;
    } Study;

};
