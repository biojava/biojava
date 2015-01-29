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
package org.biojava.bio.structure.quaternary.io;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.MMCIFFileReader;
import org.biojava.bio.structure.io.mmcif.SimpleMMcifConsumer;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssembly;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssemblyGen;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructOperList;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;



/** A provider for information about biological units for PDB files that is based on reading local MMcif files.
 * 
 * @author Andreas Prlic
 *
 */
public class MmCifPDBBiolAssemblyProvider implements RawBioUnitDataProvider{

	private String pdbId;
	private List<PdbxStructAssembly> pdbxStructAssemblies;
	private List<PdbxStructAssemblyGen> pdbxStructAssemblyGens;
	private List<PdbxStructOperList> pdbxStructOperList;
	private Structure asymUnit;
	
	private AtomCache cache ;
	
	public MmCifPDBBiolAssemblyProvider(){
		//reset();
	}
	
	@Override
	public void setPdbId(String pdbId) {
		if ( cache == null)
			cache =new AtomCache();
		
		if (
				this.pdbId != null && 
				(this.pdbId.equals(pdbId))
				) {
			// we already have all the data we need...
			return;
		}
		
		this.pdbId= pdbId;
		
		reset();
		
		MMCIFFileReader reader = new MMCIFFileReader(cache.getPath());
		FileParsingParameters params = cache.getFileParsingParams();
		params.setAlignSeqRes(true);
		params.setParseBioAssembly(true);
		reader.setFileParsingParameters(params);
		reader.setFetchBehavior(cache.getFetchBehavior());
		reader.setObsoleteBehavior(cache.getObsoleteBehavior());
		
		try{
			asymUnit = reader.getStructureById(pdbId);
			if ( asymUnit.nrModels() > 1) {
				// why do some NMR structures have bio units???
				asymUnit = StructureTools.removeModels(asymUnit);
			}
				
			SimpleMMcifConsumer consumer = reader.getMMcifConsumer();
					
			
			pdbxStructOperList 		= consumer.getStructOpers();
			pdbxStructAssemblies 	= consumer.getStructAssemblies();
			pdbxStructAssemblyGens 	= consumer.getStructAssemblyGens();
			
			//System.out.println(asymUnit.getPDBHeader());
			
			//System.out.println("OPER:" + pdbxStructOperList);
			//System.out.println("ASSEMBLIES:" + pdbxStructAssemblies);
			//System.out.println("ASSEMBLYGENS:" + pdbxStructAssemblyGens);
			
			// reset the consumer data to avoid memory leaks
			consumer.documentStart();
		} catch (IOException e){
			// TODO shouldn't this be thrown?
			e.printStackTrace();
			
		}
	}

	private void reset() {
		
		pdbxStructOperList   	= new ArrayList<PdbxStructOperList>();
		pdbxStructAssemblies	= new ArrayList<PdbxStructAssembly>();
		pdbxStructAssemblyGens 	= new ArrayList<PdbxStructAssemblyGen>();
		asymUnit = null;
	}

	public String getPdbId(){
		return pdbId;
	}
	
	@Override
	public List<PdbxStructAssembly> getPdbxStructAssemblies() {
		return pdbxStructAssemblies;
	}

	@Override
	public List<PdbxStructAssemblyGen> getPdbxStructAssemblyGens() {
		return pdbxStructAssemblyGens;
	}

	@Override
	public List<PdbxStructOperList> getPdbxStructOperList() {
		return pdbxStructOperList;
	}

	@Override
	public int getNrBiolAssemblies() {
		return pdbxStructAssemblies.size();
	}

	@Override
	public boolean hasBiolAssembly() {
		int nrAssemblies = getNrBiolAssemblies();
		if ( nrAssemblies > 0) 
			return true;
		
		return false;
	}

	@Override
	public PdbxStructAssembly getPdbxStructAssembly(int biolAssemblyNr) {
		if ( biolAssemblyNr < getNrBiolAssemblies())
			return pdbxStructAssemblies.get(biolAssemblyNr);
		return null;
	}

	@Override
	public List<PdbxStructAssemblyGen> getPdbxStructAssemblyGen(int biolAssemblyNr) {
		if ( biolAssemblyNr > getNrBiolAssemblies())
			return null;
		
		List<PdbxStructAssemblyGen> psags = new ArrayList<PdbxStructAssemblyGen>();
		for (PdbxStructAssemblyGen psag : pdbxStructAssemblyGens){
			if ( psag.getAssembly_id().equals((biolAssemblyNr +1)+ ""))
				psags.add(psag);
		}
		return psags;
	}
	
	/** get the asym unit for this PDB ID
	 * 
	 * @return
	 */
	public Structure getAsymUnit(){
		return asymUnit;
	}
	
	public void setAsymUnit(Structure s){
		this.asymUnit = s;
	}

	public AtomCache getAtomCache() {
		return cache;
	}

	public void setAtomCache(AtomCache cache) {
		this.cache = cache;
	}
	
	

}
