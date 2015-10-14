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
package org.biojava.nbio.structure.quaternary.io;

import org.biojava.nbio.structure.PDBHeader;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.biojava.nbio.core.util.SoftHashMap;

import java.io.IOException;
import java.util.List;


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

}
