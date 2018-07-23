/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.structure;

import java.io.IOException;
import java.util.Collections;
import java.util.List;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.MMCIFFileReader;
import org.biojava.nbio.structure.io.PDBFileReader;

/**
 * A class that provides static access methods for easy lookup of protein structure related components
 *
 * @author Andreas Prlic
 *
 * @since 3.0.5
 */
public class StructureIO {

	//private static final Logger logger = LoggerFactory.getLogger(StructureIO.class);

	private static AtomCache cache ;


	/** 
	 * Loads a structure based on a name. Supported naming conventions are:
	 *
	 *  <pre>
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
		biol		  := 'BIO:' pdbID [:]? [0-9]+
		pdp			  := 'PDP:' pdbID[A-Za-z0-9_]+
		resNum        := [-+]?[0-9]+[A-Za-z]?


		Example structures:
		1TIM     	#whole structure - asym unit
		4HHB.C     	#single chain
		4GCR.A_1-83 #one domain, by residue number
		3AA0.A,B    #two chains treated as one structure
		d2bq6a1     #scop domain
		BIO:1fah   #biological assembly nr 1 for 1fah
		BIO:1fah:0 #asym unit for 1fah
		BIO:1fah:1 #biological assembly nr 1 for 1fah
		BIO:1fah:2 #biological assembly nr 2 for 1fah

     * </pre>
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

	public static AtomCache getAtomCache() {
		return cache;
	}


	/**
	 * Returns the first biological assembly that is available for the given PDB id.
	 * <p>
	 * The output Structure will be different depending on the multiModel parameter:
	 * <li>
	 * the symmetry-expanded chains are added as new models, one per transformId. All original models but 
	 * the first one are discarded.
	 * </li>
	 * <li>
	 * as original with symmetry-expanded chains added with renamed chain ids and names (in the form 
	 * originalAsymId_transformId and originalAuthId_transformId)
	 * </li> 
	 * <p> 
	 * For more documentation on quaternary structures see:
	 * {@link http://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/biological-assemblies}
	 *
	 *
	 * @param pdbId
	 * @param multiModel if true the output Structure will be a multi-model one with one transformId per model, 
	 * if false the outputStructure will be as the original with added chains with renamed asymIds (in the form originalAsymId_transformId and originalAuthId_transformId).              
	 * @return a Structure object or null if that assembly is not available
	 * @throws StructureException
	 * @throws IOException
	 */
	public static Structure getBiologicalAssembly(String pdbId, boolean multiModel) throws IOException, StructureException{

		checkInitAtomCache();		

		pdbId = pdbId.toLowerCase();
		
		Structure s = cache.getBiologicalAssembly(pdbId, multiModel); 
		
		return s;
	}
	
	/**
	 * Returns the first biological assembly that is available for the given PDB id, 
	 * using multiModel={@value AtomCache#DEFAULT_BIOASSEMBLY_STYLE}
	 * <p> 
	 * For more documentation on quaternary structures see:
	 * {@link http://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/biological-assemblies}
	 *
	 *
	 * @param pdbId
	 * @return a Structure object or null if that assembly is not available
	 * @throws StructureException
	 * @throws IOException
	 */
	public static Structure getBiologicalAssembly(String pdbId) throws IOException, StructureException{
		return getBiologicalAssembly(pdbId, AtomCache.DEFAULT_BIOASSEMBLY_STYLE);
	}	

	/**
	 * Returns the biological assembly for the given PDB id and bioassembly identifier.
	 * <p>
	 * The output Structure will be different depending on the multiModel parameter:
	 * <li>
	 * the symmetry-expanded chains are added as new models, one per transformId. All original models but 
	 * the first one are discarded.
	 * </li>
	 * <li>
	 * as original with symmetry-expanded chains added with renamed chain ids and names (in the form 
	 * originalAsymId_transformId and originalAuthId_transformId)
	 * </li>  
	 * @param pdbId
	 * @param biolAssemblyNr - the ith biological assembly that is available for a PDB ID (we start counting at 1, 0 represents the asym unit).
	 * @param multiModel if true the output Structure will be a multi-model one with one transformId per model, 
	 * if false the outputStructure will be as the original with added chains with renamed asymIds (in the form originalAsymId_transformId and originalAuthId_transformId).              
	 * @return a Structure object or null if that assembly is not available
	 * @throws StructureException if there is no bioassembly available for given biolAssemblyNr or some other problems encountered while loading it
	 * @throws IOException
	 */
	public static Structure getBiologicalAssembly(String pdbId, int biolAssemblyNr, boolean multiModel) throws IOException, StructureException {
		
		checkInitAtomCache();		

		pdbId = pdbId.toLowerCase();
		
		Structure s = cache.getBiologicalAssembly(pdbId, biolAssemblyNr, multiModel); 
		
		return s;
	}
	
	/**
	 * Returns the biological assembly for the given PDB id and bioassembly identifier,
	 * using multiModel={@value AtomCache#DEFAULT_BIOASSEMBLY_STYLE}
	 * @param pdbId
	 * @param biolAssemblyNr - the ith biological assembly that is available for a PDB ID (we start counting at 1, 0 represents the asym unit).
	 * @return a Structure object or null if that assembly is not available
	 * @throws StructureException if there is no bioassembly available for given biolAssemblyNr or some other problems encountered while loading it
	 * @throws IOException
	 */
	public static Structure getBiologicalAssembly(String pdbId, int biolAssemblyNr) throws IOException, StructureException {
		return getBiologicalAssembly(pdbId, biolAssemblyNr, AtomCache.DEFAULT_BIOASSEMBLY_STYLE);
	}
		
	
	/**
	 * Returns all biological assemblies for the given PDB id.
	 * <p>
	 * The output Structure will be different depending on the multiModel parameter:
	 * <li>
	 * the symmetry-expanded chains are added as new models, one per transformId. All original models but 
	 * the first one are discarded.
	 * </li>
	 * <li>
	 * as original with symmetry-expanded chains added with renamed chain ids and names (in the form 
	 * originalAsymId_transformId and originalAuthId_transformId)
	 * </li>  
	 * If only one biological assembly is required use {@link #getBiologicalAssembly(String)} or {@link #getBiologicalAssembly(String, int)} instead.
	 * @param pdbId
	 * @param multiModel if true the output Structure will be a multi-model one with one transformId per model, 
	 * if false the outputStructure will be as the original with added chains with renamed asymIds (in the form originalAsymId_transformId and originalAuthId_transformId).              
	 * @return
	 * @throws IOException
	 * @throws StructureException
	 * @since 5.0
	 */
	public static List<Structure> getBiologicalAssemblies(String pdbId, boolean multiModel) throws IOException, StructureException {

		checkInitAtomCache();		

		pdbId = pdbId.toLowerCase();		
		
		List<Structure> s = cache.getBiologicalAssemblies(pdbId, multiModel); 
		
		return s;
		
	}

	/**
	 * Returns all biological assemblies for the given PDB id,
	 * using multiModel={@value AtomCache#DEFAULT_BIOASSEMBLY_STYLE}
	 * <p>
	 * If only one biological assembly is required use {@link #getBiologicalAssembly(String)} or {@link #getBiologicalAssembly(String, int)} instead.
	 * @param pdbId
	 * @return
	 * @throws IOException
	 * @throws StructureException
	 * @since 5.0
	 */
	public static List<Structure> getBiologicalAssemblies(String pdbId) throws IOException, StructureException {
		return getBiologicalAssemblies(pdbId, AtomCache.DEFAULT_BIOASSEMBLY_STYLE);
	}
	

	private static final String FILE_SEPARATOR = System.getProperty("file.separator");

	/**
	 * Utility method to set the location where PDB files can be found
	 *
	 * @param pathToPDBFiles
	 */
	public static void setPdbPath(String pathToPDBFiles){

		if ( ! pathToPDBFiles.endsWith(FILE_SEPARATOR))
			pathToPDBFiles += FILE_SEPARATOR;
	}


	public static enum StructureFiletype {
		PDB( (new PDBFileReader()).getExtensions()),
		CIF( new MMCIFFileReader().getExtensions()),
		UNKNOWN(Collections.<String>emptyList());

		private List<String> extensions;
		/**
		 * @param extensions List of supported extensions, including leading period
		 */
		private StructureFiletype(List<String> extensions) {
			this.extensions = extensions;
		}
		/**
		 * @return a list of file extensions associated with this type
		 */
		public List<String> getExtensions() {
			return extensions;
		}
	}

	/**
	 * Attempts to guess the type of a structure file based on the extension
	 * @param filename
	 * @return
	 */
	public static StructureFiletype guessFiletype(String filename) {
		String lower = filename.toLowerCase();
		for(StructureFiletype type : StructureFiletype.values()) {
			for(String ext : type.getExtensions()) {
				if(lower.endsWith(ext.toLowerCase())) {
					return type;
				}
			}
		}
		return StructureFiletype.UNKNOWN;
	}
}
