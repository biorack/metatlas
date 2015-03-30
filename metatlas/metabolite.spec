
/*
A KBase Atlas to liquid chromatography–mass spectrometry (LCMS) data
*/

module MetaboliteAtlas2 {

    /* Metabolite Compound object
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
    } Compound;


    typedef list<Compound> CompoundList;

    typedef structure {
        string name;
        CompoundList compounds;
        string sample;
        string method;
    } MetAtlasDictionary;


    /* @id handle */
    typedef string Run_data_ref;

    typedef structure {
        string name;
        Run_data_ref data;
    } MetAtlasFileInfo;

};
