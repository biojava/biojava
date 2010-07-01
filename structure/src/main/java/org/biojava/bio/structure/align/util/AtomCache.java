package org.biojava.bio.structure.align.util;

import java.io.IOException;

import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;


import org.biojava.bio.structure.Atom;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;

import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopInstallation;

import org.biojava.utils.io.InputStreamProvider;

public class AtomCache {
	
	public static final String CHAIN_NR_SYMBOL = ":";
	public static final String CHAIN_SPLIT_SYMBOL = ".";
		
	String path;

	
	// this is a static flag across all instances of AtomCaches, to make sure 
	// only one file is accessed in the file system at a time.
	// this is to avoid reading of partial files that are being
	// automatically downloaded by other instances / threads.	
	private static AtomicBoolean loading  = new AtomicBoolean();
	
	private static ScopInstallation scopInstallation ;
	boolean autoFetch;
	boolean isSplit;
	FileParsingParameters params;
	public AtomCache(String pdbFilePath, boolean isSplit){
		
		// we are caching the binary files that contain the PDBs gzipped
		// that is the most memory efficient way of caching...
		// set the input stream provider to caching mode
		System.setProperty(InputStreamProvider.CACHE_PROPERTY, "true");

		path = pdbFilePath;
		
		
		//this.cache = cache;
		this.isSplit = isSplit;
		
		autoFetch = true;
		loading.set(false);
		params = new FileParsingParameters();
		scopInstallation = null;
	}

	public AtomCache(UserConfiguration config){
		this(config.getPdbFilePath(),config.isSplit());
		autoFetch = config.getAutoFetch();
	}
	

	
	public String getPath() {
		return path;
	}

	public void setPath(String path) {
		this.path = path;
	}

	public boolean isSplit() {
		return isSplit;
	}

	public void setSplit(boolean isSplit) {
		this.isSplit = isSplit;
	}

	public boolean isAutoFetch() {
		return autoFetch;
	}

	public void setAutoFetch(boolean autoFetch) {
		this.autoFetch = autoFetch;
	}

	/** Returns the representation of a ScopDomain as a BioJava Structure object
	 * 
	 * @param domain a scop domain
	 * @return a Structure object.
	 * @throws IOException
	 * @throws StructureException
	 */

	public Structure getStructureForDomain(ScopDomain domain) throws IOException, StructureException{
		
		
		Structure s = null;
		
		String pdbId = domain.getPdbId();
		try {
			s = getStructure(pdbId);
	
		} catch (StructureException ex){
			System.err.println("error getting Structure for " + pdbId);
			throw new StructureException(ex);
		}
	
		String range = "(";
		int rangePos = 0;
		for ( String r : domain.getRanges()) {
			rangePos++;
			range+= r;
			if ( ( domain.getRanges().size()> 1) && (rangePos < domain.getRanges().size())){
				range+=",";
			}
			
		}
		range+=")";
		//System.out.println("getting range for "+ pdbId + " " + range);
		
		Structure n = StructureTools.getSubRanges(s, range);
		
		// get free ligands of first chain...
		if ( n.getChains().size()> 0) {
		   Chain c1 = n.getChains().get(0);
		   for ( Chain c : s.getChains()) {
		      if ( c1.getName().equals(c.getName())) {
		         List<Group> ligands = c.getAtomLigands();
		     
		         for(Group g: ligands){
		            if ( ! c1.getAtomGroups().contains(g)) {
		               c1.addGroup(g);
		            }
		         }
		      }
		   
		   }
		}
		n.setName(domain.getScopId());
		n.setPDBCode(domain.getScopId());
		
		return n;
	}
	

	/** Returns the CA atoms for the provided name. See {@link #getStructure()} for supported naming conventions.
	 * 
	 * @param name
	 * @return
	 * @throws IOException
	 * @throws StructureException
	 */
	public  Atom[] getAtoms(String name) throws IOException,StructureException{
		// synchronizing the whole method now to prevent the same PDB file to be loaded multiple times
	
		Atom[] atoms = null;
		
		//System.out.println("loading " + name);
		Structure s = null;
		try {
			s = getStructure(name);
	
		} catch (StructureException ex){
			System.err.println("error getting Structure for " + name);
			throw new StructureException(ex);
		}
		
		atoms =  StructureTools.getAtomCAArray(s);
		
		/*synchronized (cache){
			cache.put(name, atoms);
		}*/


		return atoms;
	}

	public  Atom[] getAtoms(String name, boolean clone)throws IOException,StructureException{
		Atom[] atoms =  getAtoms(name);

		if ( clone)
			return StructureTools.cloneCAArray(atoms);
		return atoms; 

	}

	

	
	/** Request a Structure based on a <i>name</i>.
	 * The following rules are applied to this name:
	 *  If only a PDB code is provided, the whole structure will be used for the alignment.
	 *	To specify a particular chain write as: 4hhb.A (chain IDs are case sensitive, PDB ids are not)
	 * 	To specify that the 1st chain in a structure should be used write: 4hhb:0 .
	 * 
	 * @param name
	 * @return
	 * @throws IOException
	 * @throws StructureException
	 */
	
   public Structure getStructure(String name) throws IOException, StructureException{
	   
	   if ( name.length() < 4)
	      throw new IllegalArgumentException("Can't interpred IDs that are shorter than 4 residues!");
	   
		while ( loading.get() ){
			// waiting for loading to be finished...
			
		}
		loading.set(true);
				
		
		Structure n = null;
		
		boolean useChainNr = false;
		boolean useDomainInfo = false;
		String range = null;
		int chainNr = -1;
	
		try {
			PDBFileReader reader = new PDBFileReader();
			reader.setPath(path);
			reader.setPdbDirectorySplit(isSplit);
			reader.setAutoFetch(autoFetch);
			
			reader.setFileParsingParameters(params);

			String pdbId   = null;
			String chainId = null;
			
			if ( name.length() == 4){
			
			   pdbId = name; 
		
			} else if ( name.startsWith("d")){
			
			
			   // looks like a SCOP domain!
			   ScopDomain domain = getScopDomain(name);
			   if ( domain == null){
			      System.err.println("Warning, could not find SCOP domain: " + name);
			   }
			   loading.set(false);

			   Structure s = getStructureForDomain(domain);
			   
			   
			   return s;
			
			} else if (name.length() == 6){
				pdbId = name.substring(0,4);
				if ( name.substring(4,5).equals(CHAIN_SPLIT_SYMBOL)) {
					chainId = name.substring(5,6);
				} else if ( name.substring(4,5).equals(CHAIN_NR_SYMBOL)) {

					useChainNr = true;	
					chainNr = Integer.parseInt(name.substring(5,6));
				}
			} else if ( (name.length() > 6) &&  
			      (name.contains(CHAIN_NR_SYMBOL))) {
			   pdbId = name.substring(0,4);
			   // this ID has domain split information...
			   useDomainInfo = true;
			   range = name.substring(5);
			}
			
			//System.out.println("got: " + name + " " + pdbId + " " + chainId + " useChainNr:" + useChainNr + " " +chainNr + " useDomainInfo:" + useDomainInfo + " " + range);

			if (pdbId == null) {
				loading.set(false);
				return null;
			}

			
			//long start  = System.currentTimeMillis();
			Structure s= reader.getStructureById(pdbId.toLowerCase());
			//long end  = System.currentTimeMillis();
			//System.out.println("time to load " + pdbId + " " + (end-start) + "\t  size :" + StructureTools.getNrAtoms(s) + "\t cached: " + cache.size());
			if ( chainId == null && chainNr < 0 && range == null) {
				loading.set(false);
				
				// we only want the 1st model in this case
				return StructureTools.getReducedStructure(s,-1);

			}
			
			
			if ( useChainNr) {
			   //System.out.println("using ChainNr");
			   n = StructureTools.getReducedStructure(s, chainNr);
			} else if ( useDomainInfo) {
			   //System.out.println("calling getSubRanges");
			   n = StructureTools.getSubRanges(s, range);
			} else  {
			   //System.out.println("reducing Chain Id " + chainId);
			   n = StructureTools.getReducedStructure(s, chainId);
			}
			
						
		} catch (Exception e){
			e.printStackTrace();
			loading.set(false);
			throw new StructureException(e.getMessage(),e);

		}
		loading.set(false);
		return n;


	}

   private ScopDomain getScopDomain(String scopId)
   {
    
      if ( scopInstallation == null) {
         scopInstallation = new ScopInstallation(path);
      }
      
      return scopInstallation.getDomainByScopID(scopId);
   }

   public FileParsingParameters getFileParsingParams()
   {
      return params;
   }

   public void setFileParsingParams(FileParsingParameters params)
   {
      this.params = params;
   }
	
	

	

	
	
}
