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
package org.biojava.nbio.structure.io.mmtf;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Bond;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava.nbio.structure.quaternary.BioAssemblyInfo;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.junit.Test;
import org.rcsb.mmtf.decoder.StructureDataToAdapter;
import org.rcsb.mmtf.encoder.AdapterToStructureData;

/**
 * Tests to see if roundtripping of MMTF can be done.
 * @author Anthony Bradley
 *
 */
public class TestMmtfRoundTrip {

	/**
	 * Test that we can round trip a simple structure.
	 * @throws IOException an error reading the file
	 * @throws StructureException an error parsing the structure
	 */
	@Test
	public void testRoundTrip() throws IOException, StructureException {
		AtomCache cache = new AtomCache();
		FileParsingParameters params = new FileParsingParameters();
        params.setParseBioAssembly(true);
        cache.setFileParsingParams(params);
		cache.setUseMmCif(true);
	    StructureIO.setAtomCache(cache);
		
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());
		
		// test case for biojava issue #770, order of subunits
	    Structure structure1 = StructureIO.getStructure("3BW1");
		AdapterToStructureData writerToEncoder = new AdapterToStructureData();
		new MmtfStructureWriter(structure1, writerToEncoder);
		MmtfStructureReader mmtfStructureReader = new MmtfStructureReader();
		new StructureDataToAdapter(writerToEncoder, mmtfStructureReader);
		Structure structure2 = mmtfStructureReader.getStructure();
		
		assertTrue(checkIfAtomsSame(structure1, structure2));
		
		checkBioAssemblies1(structure1, structure2);
	}

	/**
	 * Broad test of atom similarity
	 * @param structOne the first input structure
	 * @param structTwo the second input structure
	 * @param mmtfParams
	 * @return
	 */
	private boolean checkIfAtomsSame(Structure structOne, Structure structTwo) {
		int numModels = structOne.nrModels();
		if(numModels!=structTwo.nrModels()){
			System.out.println("Error - diff number models: "+structOne.getPDBCode());
			return false;
		}
		for(int i=0;i<numModels;i++){
			List<Chain> chainsOne = structOne.getChains(i);
			List<Chain> chainsTwo = structTwo.getChains(i);
			if(chainsOne.size()!=chainsTwo.size()){
				System.out.println("Error - diff number chains: "+structOne.getPDBCode());
				return false;
			}
			// Now make sure they're sorted in the right order
			sortChains(chainsOne, chainsTwo);
			// Check that each one has the same number of poly, non-poly and water chains
			checkDiffChains(structOne, structTwo, i);
			// Now loop over
			for(int j=0; j<chainsOne.size();j++){
				Chain chainOne = chainsOne.get(j);
				Chain chainTwo = chainsTwo.get(j);
				// Check they have the same chain id
				assertEquals(chainOne.getId(), chainTwo.getId());
				checkSeqresGroups(chainOne, chainTwo);
				List<Group> groupsOne = chainOne.getAtomGroups();
				List<Group> groupsTwo = chainTwo.getAtomGroups();
				if(groupsOne.size()!=groupsTwo.size()){
					System.out.println("Error - diff number groups: "+structOne.getPDBCode());
					System.out.println(chainOne.getId()+":"+groupsOne.size()+" "+groupsTwo.size());
					return false;
				}
				for(int k=0; k<groupsOne.size();k++){
					Group groupOne = groupsOne.get(k);
					Group groupTwo = groupsTwo.get(k);
					// Check if the groups are of the same type
					if(!groupOne.getType().equals(groupTwo.getType())){
						System.out.println("Error - diff group type: "+structOne.getPDBCode());
						System.out.println(groupOne.getPDBName() + " and type: "+groupOne.getType());
						System.out.println(groupTwo.getPDBName() + " and type: "+groupTwo.getType());;
					}
					// Check the single letter amino acid is correct
					if(groupOne.getChemComp().getOne_letter_code().length()==1 && groupTwo.getChemComp().getOne_letter_code().length()==1){
						if(!groupOne.getChemComp().getOne_letter_code().equals(groupTwo.getChemComp().getOne_letter_code())){
							System.out.println(groupOne.getPDBName());
						}
						assertEquals(groupOne.getChemComp().getOne_letter_code(), groupTwo.getChemComp().getOne_letter_code());
					}
					assertEquals(groupOne.getType(), groupTwo.getType());
					assertEquals(groupOne.getPDBName(), groupTwo.getPDBName());
					assertEquals(groupOne.getResidueNumber().getSeqNum(), groupTwo.getResidueNumber().getSeqNum());
					assertEquals(groupOne.getResidueNumber().getInsCode(), groupTwo.getResidueNumber().getInsCode());
					assertEquals(groupOne.getResidueNumber().getChainName(), groupTwo.getResidueNumber().getChainName());
					if(groupTwo.getAltLocs().size()!=groupOne.getAltLocs().size()){
						System.out.println("Error - diff number alt locs: "+structOne.getPDBCode()+" "+groupOne.getPDBName()+" "+groupOne.getResidueNumber().getSeqNum());
						System.out.println(groupOne.getAltLocs().size());
						System.out.println(groupTwo.getAltLocs().size());

					}
					// Get the first conf
					List<Atom> atomsOne = new ArrayList<>(groupOne.getAtoms());
					List<Atom> atomsTwo = new ArrayList<>(groupTwo.getAtoms());

					for(Group altLocOne: groupOne.getAltLocs()){
						for(Atom atomAltLocOne: altLocOne.getAtoms()){
							atomsOne.add(atomAltLocOne);
						}
					}
					for(Group altLocTwo: groupTwo.getAltLocs()){
						for(Atom atomAltLocTwo: altLocTwo.getAtoms()){
							atomsTwo.add(atomAltLocTwo);
						}
					}
					if(atomsOne.size()!=atomsTwo.size()){
						System.out.println("Error - diff number atoms: "+structOne.getPDBCode());
						System.out.println(groupOne.getResidueNumber());
						System.out.println(groupOne.getPDBName()+" vs "+groupTwo.getPDBName());
						System.out.println(atomsOne.size()+" vs "+atomsTwo.size());
						return false;           
					}
					// Now sort the atoms 
					sortAtoms(atomsOne, atomsTwo);
					// Now loop through the atoms
					for(int l=0;l<atomsOne.size();l++){
						Atom atomOne = atomsOne.get(l);
						Atom atomTwo = atomsTwo.get(l);
						assertEquals(atomOne.getGroup().getPDBName(), atomTwo.getGroup().getPDBName());
						assertEquals(atomOne.getCharge(), atomTwo.getCharge());
						// Check the coords are the same to three db
						assertArrayEquals(atomOne.getCoords(), atomTwo.getCoords(), 0.0009999999);
						assertEquals(atomOne.getTempFactor(), atomTwo.getTempFactor(), 0.009999999);
						assertEquals(atomOne.getOccupancy(), atomTwo.getOccupancy(), 0.009999999);
						assertEquals(atomOne.getElement(), atomTwo.getElement());
						assertEquals(atomOne.getName(),atomTwo.getName());
						assertEquals(atomOne.getAltLoc(), atomTwo.getAltLoc());
						if(i==0){
							if(atomOne.getBonds()==null){
								if(atomTwo.getBonds()!=null){
									System.out.println("Null bonds in one and not the other");
									return false;
								}
							}
							else if(atomTwo.getBonds()==null){
								System.out.println("Null bonds in one and not the other");
								return false;
							}
							else if(atomOne.getBonds().size()!=atomTwo.getBonds().size()){
								System.out.println("Error different number of bonds: "+structOne.getPDBCode());
								System.out.println(atomOne.getBonds().size()+" vs. "+atomTwo.getBonds().size());
								System.out.println(atomOne);
								System.out.println(atomTwo);
								for(Bond bond : atomOne.getBonds()) {
									System.out.println(bond);
								}
								for(Bond bond : atomTwo.getBonds()) {
									System.out.println(bond);
								}
								return false;
							}
						}
					}
				}
			}
		} 
		return true;
	}
	/**
	 * Check both structures have the same number of poly,non-poly and water chains
	 * @param structOne the first structure
	 * @param structTwo the second structure
	 * @param i the model index
	 */
	private void checkDiffChains(Structure structOne, Structure structTwo, int i) {
		assertEquals(structOne.getPolyChains(i).size(), structTwo.getPolyChains(i).size());
		assertEquals(structOne.getNonPolyChains(i).size(), structTwo.getNonPolyChains(i).size());
		assertEquals(structOne.getWaterChains(i).size(), structTwo.getWaterChains(i).size());
	}
	/**
	 * Sort the atom based on PDB serial id
	 * @param atomsOne the first list
	 * @param atomsTwo  the second list
	 */
	private void sortAtoms(List<Atom> atomsOne, List<Atom> atomsTwo) {
		atomsOne.sort(new Comparator<Atom>() {
			@Override
			public int compare(Atom o1, Atom o2) {
				//  
				if (o1.getPDBserial()<o2.getPDBserial()){
					return -1;
				}
				else{
					return 1;
				}
			}
		});
		atomsTwo.sort(new Comparator<Atom>() {
			@Override
			public int compare(Atom o1, Atom o2) {
				//  
				if (o1.getPDBserial()<o2.getPDBserial()){
					return -1;
				}
				else{
					return 1;
				}
			}
		});  		
	}

	/**
	 * Sort the chains based on chain id.
	 * @param chainsOne the first list of chains
	 * @param chainsTwo the second list of chains
	 */
	private void sortChains(List<Chain> chainsOne, List<Chain> chainsTwo) {
		Collections.sort(chainsOne, new Comparator<Chain>() {
			@Override
			public int compare(Chain o1, Chain o2) {
				return o1.getId().compareTo(o2.getId());
			}
		});
		Collections.sort(chainsTwo, new Comparator<Chain>() {
			@Override
			public int compare(Chain o1, Chain o2) {
				return o1.getId().compareTo(o2.getId());
			}
		});

	}

	private void checkSeqresGroups(Chain chainOne, Chain chainTwo) {
		
		assertEquals(chainOne.getSeqResGroups().size(), chainTwo.getSeqResGroups().size());
		
		for (int i=0; i<chainOne.getSeqResGroups().size(); i++) {
			Group gOne = chainOne.getSeqResGroup(i);
			Group gTwo = chainTwo.getSeqResGroup(i);
			assertNotNull(gOne.getChemComp());
			assertNotNull(gTwo.getChemComp());
			assertEquals(gOne.getChemComp().getOne_letter_code(), gTwo.getChemComp().getOne_letter_code());
			
			assertEquals(gOne.getResidueNumber(), gTwo.getResidueNumber());
			//assertEquals(gOne.getPDBName(), gTwo.getPDBName());
		}
		
		// this is to test issue #517: a null pointer happens if the group hasn't been cloned (including the chem comp associated to it)
		Chain clonedChain = (Chain)chainTwo.clone();
		assertEquals(chainTwo.getSeqResGroups().size(), clonedChain.getSeqResGroups().size());
		for (int i=0; i<clonedChain.getSeqResGroups().size(); i++) {
			Group g = clonedChain.getSeqResGroup(i);
			assertNotNull(g.getChemComp());
			
		}
	}
	
	/**
     * Checks consistency of bioassemblies
     * @param structOne the first input structure
     * @param structTwo the second input structure
     */
    private void checkBioAssemblies1(Structure structOne, Structure structTwo) throws IOException {
        
        Map<Integer, BioAssemblyInfo> expecteds = structOne.getPDBHeader().getBioAssemblies();
        Map<Integer, BioAssemblyInfo> actuals = structTwo.getPDBHeader().getBioAssemblies();
        assertEquals(expecteds.size(), actuals.size());
        
        assertEquals(new ArrayList<>(expecteds.keySet()), new ArrayList<>(actuals.keySet()));
       
        List<BioAssemblyInfo> assemblies1 = new ArrayList<>(expecteds.values());
        List<BioAssemblyInfo> assemblies2 = new ArrayList<>(actuals.values());
        
        for (int i = 0; i < assemblies1.size(); i++) {
            BioAssemblyInfo info1 = assemblies1.get(i);
            BioAssemblyInfo info2 = assemblies2.get(i);
            assertEquals(info1.getId(), info2.getId());
            assertEquals(info1.getTransforms().size(), info2.getTransforms().size());

            for (int j = 0; j < info1.getTransforms().size(); j++) {
                BiologicalAssemblyTransformation trans1 = info1.getTransforms().get(j);
                BiologicalAssemblyTransformation trans2 = info2.getTransforms().get(j);
                
                assertEquals(trans1.getChainId(), trans2.getChainId());
                assertTrue(trans1.getTransformationMatrix().epsilonEquals(trans2.getTransformationMatrix(), 0.000001));
            }
        }
    }
}
