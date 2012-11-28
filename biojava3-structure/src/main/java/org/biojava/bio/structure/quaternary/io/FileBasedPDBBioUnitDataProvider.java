/**
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
 * Created on Oct 22, 2012
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.bio.structure.quaternary.io;

import java.io.IOException;
import java.util.List;

import org.biojava.bio.structure.PDBHeader;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.quaternary.ModelTransformationMatrix;
import org.biojava3.core.util.SoftHashMap;


/** A biounit provider that loads the biol assembly from thepublic PDB file, rather than re-creating it on the fly out of the default PDB file for the asym unit
 * 
 * @author Andreas Prlic
 *
 */
public class FileBasedPDBBioUnitDataProvider implements BioUnitDataProvider{

	SoftHashMap<String, PDBHeader> headerCache = new SoftHashMap<String, PDBHeader>(0);

	Structure s = null;
	
	AtomCache cache = new AtomCache();
	
	public Structure getBioUnit(String pdbId, int bioUnit) throws IOException{
		//System.out.println("load PDB + bioUnit " + bioUnit + " " );
		
		if ( bioUnit == 0){
			PDBBioUnitDataProvider fprov = new PDBBioUnitDataProvider();
			Structure s = fprov.getAsymUnit(pdbId);
			return s;
		}
		
		PDBFileReader reader = new PDBFileReader();
		reader.setAutoFetch(cache.isAutoFetch());
		reader.setFetchCurrent(cache.isFetchCurrent());
		reader.setFetchFileEvenIfObsolete(cache.isFetchFileEvenIfObsolete());
		
		reader.setPath(cache.getPath());
		
		FileParsingParameters params = cache.getFileParsingParams();
		
		if ( bioUnit > 0 ) {
			params.setParseBioAssembly(true);
			reader.setBioAssemblyId(bioUnit);
			reader.setBioAssemblyFallback(false);
		}
			
		params.setHeaderOnly(false);
		
		reader.setFileParsingParameters(params);
		
		return reader.getStructureById(pdbId);
	}
	

	public PDBHeader loadPDB(String pdbId, int bioUnit){
		
		if ( cache == null) {
			cache = new AtomCache();
		}
		
		
		PDBHeader header = null;
		try {
			s = getBioUnit(pdbId,bioUnit);
			
			header = s.getPDBHeader();
			//System.out.println("got header: " + bioUnit + " " + header + " from parsing");
			headerCache.put(s.getPDBCode(),header);
		} catch (Exception e){
			e.printStackTrace();
		}	
		return header ;
	}

	public Structure getAsymUnit(String pdbId){

		if (s == null ||( ! s.getPDBCode().equalsIgnoreCase(pdbId))) {
			loadPDB(pdbId,0);
		}

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
			header = loadPDB(pdbId,biolAssemblyNr);
		}

		return header.getBioUnitTranformationMap().get(biolAssemblyNr);


	}

	@Override
	public int getNrBiolAssemblies(String pdbId) {
		PDBHeader header = headerCache.get(pdbId);

		if ( header == null) {
			header = loadPDB(pdbId,0);
		}

		return header.getNrBioAssemblies();
	}

	@Override
	public boolean hasBiolAssembly(String pdbId) {
		PDBHeader header = headerCache.get(pdbId);

		if ( header == null) {
			
			header = loadPDB(pdbId,0);
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
