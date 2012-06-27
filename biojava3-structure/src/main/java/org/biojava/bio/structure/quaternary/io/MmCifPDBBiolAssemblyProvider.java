package org.biojava.bio.structure.quaternary.io;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.MMCIFFileReader;
import org.biojava.bio.structure.io.mmcif.SimpleMMcifConsumer;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssembly;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssemblyGen;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructOperList;


/** A provider for information about biological units for PDB files that is based on reading local MMcif files.
 * 
 * @author Andreas Prlic
 *
 */
public class MmCifPDBBiolAssemblyProvider implements PDBBioUnitDataProvider{

	String pdbId;
	List<PdbxStructAssembly> pdbxStructAssemblies;
	List<PdbxStructAssemblyGen> pdbxStructAssemblyGens;
	List<PdbxStructOperList> pdbxStructOperList;
	Structure asymUnit;
	
	public MmCifPDBBiolAssemblyProvider(){
		//reset();
	}
	
	@Override
	public void setPdbId(String pdbId) {
		
		if (
				this.pdbId != null && 
				(this.pdbId.equals(pdbId))
				) {
			// we already have all the data we need...
			return;
		}
		
		this.pdbId= pdbId;
		
		reset();
		
		MMCIFFileReader reader = new MMCIFFileReader();
		FileParsingParameters params = new FileParsingParameters();
		reader.setFileParsingParameters(params);
		
		try{
			asymUnit = reader.getStructureById(pdbId);
		
			SimpleMMcifConsumer consumer = reader.getMMcifConsumer();
					
			
			pdbxStructOperList 		= consumer.getStructOpers();
			pdbxStructAssemblies 	= consumer.getStructAssemblies();
			pdbxStructAssemblyGens 	= consumer.getStructAssemblyGens();
			
			//System.out.println(pdbxStructOperList);
			//System.out.println(pdbxStructAssemblies);
			//System.out.println(pdbxStructAssemblyGens);
			
			// reset the consumer data to avoid memory leaks
			consumer.documentStart();
		} catch (Exception e){
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
	public PdbxStructAssemblyGen getPdbxStructAssemblyGen(int biolAssemblyNr) {
		if ( biolAssemblyNr < getNrBiolAssemblies())
			return pdbxStructAssemblyGens.get(biolAssemblyNr);
		return null;
	}
	
	/** get the asym unit for this PDB ID
	 * 
	 * @return
	 */
	public Structure getAsymUnit(){
		return asymUnit;
	}

}
