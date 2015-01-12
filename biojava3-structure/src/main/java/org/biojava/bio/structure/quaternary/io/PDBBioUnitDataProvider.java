package org.biojava.bio.structure.quaternary.io;

import java.io.IOException;
import java.util.Collections;
import java.util.List;

import org.biojava.bio.structure.PDBHeader;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssembly;
import org.biojava.bio.structure.quaternary.BiologicalAssemblyTransformation;
import org.biojava3.core.util.SoftHashMap;


/** A BioUnitDataProvider that extracts the necessary info from PDB files
 * 
 * @author Andreas Prlic
 *
 */
public class PDBBioUnitDataProvider implements BioUnitDataProvider{

	
	private SoftHashMap<String, PDBHeader> headerCache = new SoftHashMap<String, PDBHeader>(0);
	
	private Structure s;
	
	// no initialisation here, this gives an opportunity to setAtomCache to initialise it
	private AtomCache cache;
	
	public PDBHeader loadPDB(String pdbId){
			
		
		FileParsingParameters params = null;
		
		if ( cache == null)
			cache = new AtomCache();
		
		params = cache.getFileParsingParams();
				
		if ( params == null)
			params = new FileParsingParameters();

		params.setParseBioAssembly(true);		
		params.setAlignSeqRes(true);
		
		PDBHeader header = null;
		try {
			s =  cache.getStructure(pdbId);
			
			header = s.getPDBHeader();
			headerCache.put(s.getPDBCode(),header);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (StructureException e) {
			e.printStackTrace();
		}
		
		
		return header ;
	}
	
	@Override
	public Structure getAsymUnit(String pdbId){
		
		if (s == null ||( ! s.getPDBCode().equalsIgnoreCase(pdbId))) {
			loadPDB(pdbId);
		}
		
		if ( s.nrModels() > 1)  {
			s = StructureTools.removeModels(s);
		}
		
		
		return s;
	}
	@Override
	public void setAsymUnit(Structure s){
		this.s = s;
	}
	
	@Override
	public List<BiologicalAssemblyTransformation> getBioUnitTransformationList(
			String pdbId, int biolAssemblyNr) {

	
		PDBHeader header = headerCache.get(pdbId);
		
		if ( header == null) {
			header = loadPDB(pdbId);
		}
		
		return header.getBioAssemblies().get(biolAssemblyNr).getTransforms();


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

	@Override
	public List<PdbxStructAssembly> getPdbxStructAssemblies() {
		return Collections.emptyList();		
	}

}
