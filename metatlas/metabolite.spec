
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
        float mz;
        float mz_threshold;
        float rt_min;
        float rt_max;
        float rt_peak;
        float neutral_mass;
        string pubchem_id;
    } MACompound;

    typedef structure {
        string name;
        list<MACompound> compounds;
        string sample;
        string method;
    } MADictionary;


    /* @id handle */
    typedef string Run_data_ref;

    /*
    @optional atlases polarity group inclusion_order normalization_factor retention_correction
    */
    typedef structure {
           string name;
           list<MADictionary> atlases;
           string polarity;
           string group;
           string inclusion_order;
           string normalization_factor;
           string retention_correction;
           Run_data_ref run_file_id;
       } MAFileInfo;


   typedef structure {
        string name;
        string description;
        string reconstitution_method;
        string quenching_method;
        string extraction_method;
        string chromatography_method;
        list<MADictionary> atlas_ids;
    } MAExperiment;

};
