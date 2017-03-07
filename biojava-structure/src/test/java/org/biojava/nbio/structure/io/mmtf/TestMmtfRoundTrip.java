package org.biojava.nbio.structure.io.mmtf;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Bond;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
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
		cache.setUseMmCif(true);
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());
		
		
		StructureIO.setAtomCache(cache);
		Structure structure = StructureIO.getStructure("4CUP");
		AdapterToStructureData writerToEncoder = new AdapterToStructureData();
		new MmtfStructureWriter(structure, writerToEncoder);
		MmtfStructureReader mmtfStructureReader = new MmtfStructureReader();
		new StructureDataToAdapter(writerToEncoder, mmtfStructureReader);
		assertTrue(checkIfAtomsSame(structure,mmtfStructureReader.getStructure()));
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
}
