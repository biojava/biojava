package org.biojava.bio.structure.quaternary.io;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssembly;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssemblyGen;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructOperList;
import org.biojava.bio.structure.quaternary.BiologicalAssemblyBuilder;
import org.biojava.bio.structure.quaternary.ModelTransformationMatrix;
import org.biojava3.structure.StructureIO;

public class MmCifBiolAssemblyProvider implements BioUnitDataProvider {

	public MmCifBiolAssemblyProvider(){
		
	}
	
	public Structure getAsymUnit(String pdbId){
		MmCifPDBBiolAssemblyProvider provider = new MmCifPDBBiolAssemblyProvider();

		provider.setPdbId(pdbId);
		
		return provider.getAsymUnit();
	}
	
	@Override
	public List<ModelTransformationMatrix> getBioUnitTransformationList(
			String pdbId, int biolAssemblyNr) {
		
		MmCifPDBBiolAssemblyProvider provider = new MmCifPDBBiolAssemblyProvider();

		provider.setPdbId(pdbId);
		
		PdbxStructAssembly psa = provider.getPdbxStructAssembly(biolAssemblyNr) ;
		
		PdbxStructAssemblyGen psag = provider.getPdbxStructAssemblyGen(biolAssemblyNr);
		
		if ( psa == null || psag == null) {
			return null;
		}
		//System.out.println(psa);
		//System.out.println(psag);
		
		List<PdbxStructOperList> operators = provider.getPdbxStructOperList();
		//System.out.println(operators);
		
		
		/** now we start to rebuild the quaternary structure
		 * 
		 */
		
		BiologicalAssemblyBuilder builder = new BiologicalAssemblyBuilder();
		
		// these are the transformations that need to be applied to our model
		List<ModelTransformationMatrix> transformations = builder.getBioUnitTransformationList(psa, psag, operators);
		
		return transformations;
	}

	@Override
	public int getNrBiolAssemblies(String pdbId) {
		MmCifPDBBiolAssemblyProvider provider = new MmCifPDBBiolAssemblyProvider();

		provider.setPdbId(pdbId);
		
		return provider.getNrBiolAssemblies();
	}

	@Override
	public boolean hasBiolAssembly(String pdbId) {
		MmCifPDBBiolAssemblyProvider provider = new MmCifPDBBiolAssemblyProvider();

		provider.setPdbId(pdbId);
		
		return provider.hasBiolAssembly();
	}
	
	public Structure getBiolAssembly(String pdbId, int biolAssemblyNr) throws IOException, StructureException{
		PDBBioUnitDataProvider provider = new MmCifPDBBiolAssemblyProvider();
		
		provider.setPdbId(pdbId);
		
		if ( ! provider.hasBiolAssembly()){
			return null;
		}
				
		if (  provider.getNrBiolAssemblies() <= biolAssemblyNr){
			return null;
		}
		/** First we read the required data from wherever we get it from (configured in the factory)
		 * 
		 */
		PdbxStructAssembly psa = provider.getPdbxStructAssembly(biolAssemblyNr) ;
		
		PdbxStructAssemblyGen psag = provider.getPdbxStructAssemblyGen(biolAssemblyNr);
		
		//System.out.println(psa);
		//System.out.println(psag);
		
		List<PdbxStructOperList> operators = provider.getPdbxStructOperList();
		//System.out.println(operators);
		
		
		/** now we start to rebuild the quaternary structure
		 * 
		 */
		
		BiologicalAssemblyBuilder builder = new BiologicalAssemblyBuilder();
		
		// these are the transformations that need to be applied to our model
		ArrayList<ModelTransformationMatrix> transformations = builder.getBioUnitTransformationList(psa, psag, operators);
		
		
		
		Structure asymUnit = null;
		
		if ( provider instanceof MmCifPDBBiolAssemblyProvider){
			MmCifPDBBiolAssemblyProvider mmcifprov = (MmCifPDBBiolAssemblyProvider) provider;
			asymUnit = mmcifprov.getAsymUnit();
		} else {
			
			// how to request internal chain IDs?
			
			asymUnit = StructureIO.getStructure(pdbId);
			
		}
						
		return builder.rebuildQuaternaryStructure(asymUnit, transformations);
	}

}
