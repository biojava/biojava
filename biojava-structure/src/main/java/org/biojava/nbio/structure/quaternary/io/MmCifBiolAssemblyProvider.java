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

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssembly;
import org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssemblyGen;
import org.biojava.nbio.structure.io.mmcif.model.PdbxStructOperList;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyBuilder;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class MmCifBiolAssemblyProvider implements BioUnitDataProvider {

	private MmCifPDBBiolAssemblyProvider provider; 
	
	// no initialisation here, this gives an opportunity to setAtomCache to initialise it
	private AtomCache cache;
	
	public MmCifBiolAssemblyProvider(){
		provider  = new MmCifPDBBiolAssemblyProvider();
	}
	
	@Override
	public Structure getAsymUnit(String pdbId){
	 
		provider.setPdbId(pdbId);
		
		Structure s1 = provider.getAsymUnit();
		return s1;
	}
	
	@Override
	public void setAsymUnit(Structure s){
		provider.setAsymUnit(s);
	}
	
	@Override
	public List<BiologicalAssemblyTransformation> getBioUnitTransformationList(
			String pdbId, int biolAssemblyNr) {
		

		provider.setPdbId(pdbId);

		List<PdbxStructAssembly> psas= provider.getPdbxStructAssemblies();
		
		PdbxStructAssembly psa = null;
		if ( psas.size() >= biolAssemblyNr ){
			psa = psas.get(biolAssemblyNr -1);
		} else {
			throw new IllegalArgumentException("Requested not existing biolAssemblyNr");
			
		}

		// we start counting at 1!
		//PdbxStructAssembly psa = provider.getPdbxStructAssembly(biolAssemblyNr-1) ;
		
		List<PdbxStructAssemblyGen> psags = provider.getPdbxStructAssemblyGen(biolAssemblyNr-1);
				
		if ( psa == null || psags == null) {
			return null;
		}
		//System.out.println(psa);
		//System.out.println(psags);
		
		List<PdbxStructOperList> operators = provider.getPdbxStructOperList();
		//System.out.println(operators);
		
		
		/** 
		 * Now we start to rebuild the quaternary structure
		 */
		
		BiologicalAssemblyBuilder builder = new BiologicalAssemblyBuilder();
		
		// these are the transformations that need to be applied to our model
		List<BiologicalAssemblyTransformation> transformations = builder.getBioUnitTransformationList(psa, psags, operators);
		//System.out.println(transformations);
		return transformations;
	}

	@Override
	public int getNrBiolAssemblies(String pdbId) {
		
		provider.setPdbId(pdbId);
		
		return provider.getNrBiolAssemblies();
	}

	@Override
	public boolean hasBiolAssembly(String pdbId) {


		provider.setPdbId(pdbId);
		
		return provider.hasBiolAssembly();
	}
	
	public Structure getBiolAssembly(String pdbId, int biolAssemblyNr) throws IOException, StructureException{

		
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
		
		List<PdbxStructAssemblyGen> psags = provider.getPdbxStructAssemblyGen(biolAssemblyNr);
		
		
		List<PdbxStructOperList> operators = provider.getPdbxStructOperList();
		
		
		/** now we start to rebuild the quaternary structure
		 * 
		 */
		
		BiologicalAssemblyBuilder builder = new BiologicalAssemblyBuilder();
		
		// these are the transformations that need to be applied to our model
		ArrayList<BiologicalAssemblyTransformation> transformations = builder.getBioUnitTransformationList(psa, psags, operators);
		
		
		
		Structure asymUnit = null;
		
		if ( provider instanceof MmCifPDBBiolAssemblyProvider){
			MmCifPDBBiolAssemblyProvider mmcifprov = provider;
			asymUnit = mmcifprov.getAsymUnit();
		} else {
			
			// how to request internal chain IDs?
			
			asymUnit = StructureIO.getStructure(pdbId);
			
		}
						
		return builder.rebuildQuaternaryStructure(asymUnit, transformations);
	}

	@Override
	public void setAtomCache(AtomCache cache) {

		this.cache =cache;
		provider.setAtomCache(cache);
	}

	@Override
	public AtomCache getAtomCache() {
		return cache;
	}

}
