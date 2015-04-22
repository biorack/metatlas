<h1>Readable Labels and Programic IDs for Features in a Metabolite Atlas</h1>
<ul>
	<li>IDs should be flexible enough to permit capturing the quick thoughts about an oberservation</li>
	<li>IDs should be specific enough to permit the comparison of like-for-like metabolites across Atlases</li>
</ul>
Identification of metabolites in untargeted LC/MS and MS/MS experiments is analagous to describing the ingredients in a high-school lunch cafeteria. You might describe the dish as "Steak and Potatoes" and someone else might describe it as "grey-goo".  Both could be correct.  Perhaps "Steak" is an observation that a component of the dish has the characteristics of steak and perhaps its based on actual knowledge of the dish containing one or more components of a cow.  

In an LC/MS experiment, there are observations which can be assigned to a specific chemical entitiy with very high confidence. For these cases, an investigator might specify a precise chemical structure they believe to have detected.  There is no more accurate way to specify an observed metabolite than a complete structural description (as a SMILES or InChi string) that contains the charged form of the detected molecule.  These structures can be matched to chemical databases unambiguously and are the unique description of a compound.  In comparison, there are often times when only partial knowedge about the identity is known.  Difficulty often is due to isomers (chemicals with the same chemical formula but different structure). Many sugars exist with the forumula C6H12O6, but knowing if its glucose or galactose might require extensive investigation. Consequently, an investigator might choose to identify a feature resembling glucose as "A Hexose" to constrain the specification appropriately.  

What is needed is a label-xref service.  There are millions of known metabolite structures.  Tens of thousands of these are captured by biochemical databases.  Exact XREF with boolean logic to specifically show which chemical structures comprise a assertion is the only way to precisely capture the meaning of the assertion. 

The XREF table would have name, structure, and comment as fields

Users can label the peaks they are annotating with free text descriptions.  This description must correspond to an entry in the XREF table as column "name".  

Within the XREF service IDs within metatlas will be
<ol>
<li>m/z with uncertainty</li>
<li>chemical formula & adduct (one or more formulae with boolean logic). formulae are given as neutral formula with adduct (ie: proton, sodium, etc) as two terms</li>
<li>chemical structure & adduct (one of more structures with boolean logic)</li>
<li>Isotopologues and Isotopomers.  This is essential and can be handled adhoc at the present time.</li>
</ol>

A layer that unwraps the xref service to link out will use the following IDs
<ol>
<li>neutral mass with uncertainty</li>
<li>neutral chemical formula (one or more formulae with boolean logic)</li>
<li>chemical structure (one of more neutral structures with boolean logic)</li>
</ol>

Given that these molecules are detected in a mass spectrometer, they will be seen as ions. Consideration of how to maintain the detected ion and also handling the neutral form for identification. The same consieration must be given to isotopologues and isotopomers.

This is the best example of an ontology file for chemical information.
http://semanticscience.org/ontology/cheminf.owl

Adding a few terms specific to electrospray-ionization and liquid chromatography should enable us to relate features to one another.  These terms would include:
<ul>
	<li>Ion Cluster</li>
	<li>Adduct</li>
	<li>Degradation Product</li>
	<li>Pure spectra member</li>
</ul>

Here is the hub that has all their ontologies:
https://code.google.com/p/semanticscience/

Its worth noting that the EC IDs used to build knowledge about enzyme function is a generica and ambiguos way to capture the compounds that participate in many reactions:
ftp://ftp.expasy.org/databases/enzyme/enzyme.dat
Things like, an alcohol, are listed under the enzyme name, alcohol dehydrogenase.  This likely didn't happen overnight.

Some vaguely useful things here:
http://www.semantic-science.org/

Maybe useful talking to this guy:
http://dumontierlab.com/

