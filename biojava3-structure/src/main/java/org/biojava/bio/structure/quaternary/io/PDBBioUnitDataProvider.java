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
	
	public PDBHeader loadPDB(String pdbId){
		AtomCache cache = new AtomCache();

		FileParsingParameters params = cache.getFileParsingParams();

		params.setParseBioAssembly(true);
		params.setUpdateRemediatedFiles(true);
		
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
	
	public Structure getAsymUnit(){
		
		if ( s.nrModels() > 1) 
			s = StructureTools.removeModels(s);
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

}
