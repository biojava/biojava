package org.biojava.bio.structure.quaternary.io;

import java.util.List;

import org.biojava.bio.structure.PDBHeader;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.quaternary.ModelTransformationMatrix;
import org.biojava3.core.util.SoftHashMap;


/** A BioUnitDataProvider that extracts the necessary info from PDB files
 * 
 * @author Andreas Prlic
 *
 */
public class PDBBioUnitDataProvider implements BioUnitDataProvider{

	
	SoftHashMap<String, PDBHeader> headerCache = new SoftHashMap<String, PDBHeader>(0);
	
	Structure s = null;
	
	AtomCache cache = new AtomCache();
	
	public PDBHeader loadPDB(String pdbId){
				
		FileParsingParameters params = cache.getFileParsingParams();

		params.setParseBioAssembly(true);		
		params.setAlignSeqRes(true);
		
		PDBHeader header = null;
		try {
			s =  cache.getStructure(pdbId);
			
			header = s.getPDBHeader();
			headerCache.put(s.getPDBCode(),header);
		} catch (Exception e){
			e.printStackTrace();
		}	
		
		
		return header ;
	}
	
	public Structure getAsymUnit(String pdbId){
		
		if (s == null ||( ! s.getPDBCode().equalsIgnoreCase(pdbId))) {
			loadPDB(pdbId);
		}
		
		if ( s.nrModels() > 1)  {
			// temporarily overwrite the NMR setting
			// so we can get rid of multiple modles
			// eg for 2F03 which is an xray with multi-models...
			boolean isNMR = s.isNmr();
			s.setNmr(true);
			s = StructureTools.removeModels(s);
			s.setNmr(isNMR);
		}
		
		
		return s;
	}
	public void setAsymUnit(Structure s){
		this.s = s;
	}
	
	@Override
	public List<ModelTransformationMatrix> getBioUnitTransformationList(
			String pdbId, int biolAssemblyNr) {

	
		PDBHeader header = headerCache.get(pdbId);
		
		if ( header == null) {
			header = loadPDB(pdbId);
		}
		
		return header.getBioUnitTranformationMap().get(biolAssemblyNr);


	}

	@Override
	public int getNrBiolAssemblies(String pdbId) {
		PDBHeader header = headerCache.get(pdbId);
		
		if ( header == null) {
			header = loadPDB(pdbId);
		}
		
		return header.getNrBioAssemblies();
	}

	@Override
	public boolean hasBiolAssembly(String pdbId) {
		PDBHeader header = headerCache.get(pdbId);
		
		if ( header == null) {
			header = loadPDB(pdbId);
		}
		
		if ( header.getNrBioAssemblies() > 0) {
			return true;
		}
		
		return false;
		
	}

	@Override
	public void setAtomCache(AtomCache cache) {
		this.cache = cache;
		
	}

	@Override
	public AtomCache getAtomCache() {
		return cache;
	}

}
