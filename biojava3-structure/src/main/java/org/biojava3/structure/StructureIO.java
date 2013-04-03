package org.biojava3.structure;

import java.io.IOException;
import java.util.List;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.quaternary.ModelTransformationMatrix;
import org.biojava.bio.structure.quaternary.BiologicalAssemblyBuilder;
import org.biojava.bio.structure.quaternary.io.BioUnitDataProvider;
import org.biojava.bio.structure.quaternary.io.BioUnitDataProviderFactory;

/** A class that provides static access methods for easy lookup of protein structure related components
 * 
 * @author Andreas Prlic
 *
 * @since 3.0.5
 */
public class StructureIO {


	private static AtomCache cache ;


	/** Loads a structure based on a name. Supported naming conventions are:
	 * 
	 *  * <pre>
		Formal specification for how to specify the <i>name</i>:

		name     := pdbID
		               | pdbID '.' chainID
		               | pdbID '.' range
		               | scopID 
		               | biol
		               | pdp
		range         := '('? range (',' range)? ')'?
		               | chainID
		               | chainID '_' resNum '-' resNum
		pdbID         := [0-9][a-zA-Z0-9]{3}
		chainID       := [a-zA-Z0-9]
		scopID        := 'd' pdbID [a-z_][0-9_]
		biol		  := 'BIOL:' pdbID [:]? [0-9]+
		pdp			  := 'PDP:' pdbID[A-Za-z0-9_]+ 
		resNum        := [-+]?[0-9]+[A-Za-z]?


		Example structures:
		1TIM     	#whole structure - asym unit
		4HHB.C     	#single chain
		4GCR.A_1-83 #one domain, by residue number
		3AA0.A,B    #two chains treated as one structure
		d2bq6a1     #scop domain
		BIOL:1fah   #biological assembly nr 1 for 1fah
		BIOL:1fah:0 #asym unit for 1fah
		BIOL:1fah:1 #biological assembly nr 1 for 1fah
		BIOL:1fah:2 #biological assembly nr 2 for 1fah
		
		</pre>
	 * 
	 * With the additional set of rules:
	 *  
	 *  <ul>
	 *  <li>If only a PDB code is provided, the whole structure will be return including ligands, but the first model only (for NMR).
	 *	<li>Chain IDs are case sensitive, PDB ids are not. To specify a particular chain write as: 4hhb.A or 4HHB.A </li>
	 *  <li>To specify a SCOP domain write a scopId e.g. d2bq6a1. Some flexibility can be allowed in SCOP domain names, see {@link #setStrictSCOP(boolean)}</li>
	 *  <li>URLs are accepted as well</li>
	 *  </ul>
	 *  
	 * @param name
	 * @return a Structure object, or null if name appears improperly formated (eg too short, etc)
	 * @throws IOException The PDB file cannot be cached due to IO errors
	 * @throws StructureException The name appeared valid but did not correspond to a structure.
	 * 	Also thrown by some submethods upon errors, eg for poorly formatted subranges.
	 */
	public static Structure getStructure(String name) throws IOException, StructureException{

		checkInitAtomCache();

		// delegate this functionality to AtomCache...

		return cache.getStructure(name);

	}


	private static void checkInitAtomCache() {
		if ( cache == null){
			cache = new AtomCache();
		}

	}

	public static void setAtomCache(AtomCache c){
		cache = c;
	}

	/** Returns the first biologicalAssembly that is available for a protein structure. For more documentation on quaternary structures see:
	 * {@link http://www.pdb.org/pdb/101/static101.do?p=education_discussion/Looking-at-Structures/bioassembly_tutorial.html}
	 * 
	 * 
	 * @param pdbId
	 * @return a Structure object or null if that assembly is not available
	 * @throws StructureException 
	 * @throws IOException 
	 */
	public static Structure getBiologicalAssembly(String pdbId) throws IOException, StructureException{
				
		return getBiologicalAssembly(pdbId,1);
	}

	/** By default the getStructure method loads asym units. This access method allows to recreate the quaternary structure for a protein if it is available.
	 * 
	 * @param pdbId
	 * @param biolAssemblyNr - the ith biological assembly that is available for a PDB ID (we start counting at 1, 0 represents the asym unit).
	 * @return a Structure object or null if that assembly is not available
	 * @throws StructureException 
	 * @throws IOException 
	 */
	public static Structure getBiologicalAssembly(String pdbId, int biolAssemblyNr) throws IOException, StructureException {

		pdbId = pdbId.toLowerCase();
		
		BioUnitDataProvider provider = BioUnitDataProviderFactory.getBioUnitDataProvider();

		if ( cache != null)
			provider.setAtomCache(cache);
		
		Structure asymUnit = null;

		asymUnit = provider.getAsymUnit(pdbId);
		
		//cleanup to avoid memory leaks
		provider.setAsymUnit(null);
		provider.setAtomCache(null);
		
		// 0 ... asym unit
		if ( biolAssemblyNr == 0) {
			return asymUnit;
		}
		
		//List<ModelTransformationMatrix> transformations = provider.getBioUnitTransformationList(pdbId, biolAssemblyNr -1);
		List<ModelTransformationMatrix> transformations = asymUnit.getPDBHeader().getBioUnitTranformationMap().get(biolAssemblyNr  );
		if ( transformations == null || transformations.size() == 0){
			
			throw new StructureException("Could not load transformations to recreate biological assembly nr " + biolAssemblyNr + " of " + pdbId);
		}
		//System.out.println("StructureIO:" + transformations.get(0).toXML());
		BiologicalAssemblyBuilder builder = new BiologicalAssemblyBuilder();

		return builder.rebuildQuaternaryStructure(asymUnit, transformations);


	}

	/** Does the provider PDB ID have a biological assembly?
	 * 
	 * @param pdbId
	 * @return flag if one or more biological assemblies are available
	 */
	public static boolean hasBiologicalAssembly(String pdbId){

		pdbId = pdbId.toLowerCase();
		
		BioUnitDataProvider provider = BioUnitDataProviderFactory.getBioUnitDataProvider();


		return provider.hasBiolAssembly(pdbId);

	}

	public static int getNrBiologicalAssemblies(String pdbId){

		pdbId = pdbId.toLowerCase();
		
		BioUnitDataProvider provider = BioUnitDataProviderFactory.getBioUnitDataProvider();

		return provider.getNrBiolAssemblies(pdbId);
	}

	private static final String FILE_SEPARATOR = System.getProperty("file.separator");

	/** Utility method to set the location where PDB files can be found
	 * 
	 * @param pathToPDBFiles
	 */
	public static void setPdbPath(String pathToPDBFiles){

		if ( ! pathToPDBFiles.endsWith(FILE_SEPARATOR))
			pathToPDBFiles += FILE_SEPARATOR;
		System.setProperty(AbstractUserArgumentProcessor.PDB_DIR,pathToPDBFiles);
	}
}
