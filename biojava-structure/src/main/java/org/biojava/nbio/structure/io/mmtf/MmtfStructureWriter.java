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

import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Bond;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.EntityInfo;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.PDBCrystallographicInfo;
import org.biojava.nbio.structure.PDBHeader;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.biojava.nbio.structure.quaternary.BioAssemblyInfo;
import org.rcsb.mmtf.api.StructureAdapterInterface;
import org.rcsb.mmtf.dataholders.MmtfStructure;

/**
 * Class to take Biojava structure data and covert to the DataApi for encoding. 
 * Must implement all the functions in {@link StructureAdapterInterface}.
 * @author Anthony Bradley
 *
 */
public class MmtfStructureWriter {


	private StructureAdapterInterface mmtfDecoderInterface;

	/**
	 * Pass data from Biojava structure  to another generic output type. Loops through the data 
	 * structure and calls all the set functions.
	 * @param structure the input {@link Structure} to write
	 * @param dataTransferInterface the generic interface that 
	 * implements all the set methods.
	 */
	public MmtfStructureWriter(Structure structure, StructureAdapterInterface dataTransferInterface) {
		this.mmtfDecoderInterface = dataTransferInterface;
		// Reset structure to consider altloc groups with the same residue number but different group names as seperate groups
		MmtfUtils.fixMicroheterogenity(structure);
		// Get the chain name to index map
		MmtfSummaryDataBean mmtfSummaryDataBean = MmtfUtils.getStructureInfo(structure);
		Map<String, Integer> chainIdToIndexMap = mmtfSummaryDataBean.getChainIdToIndexMap();
		List<Atom> allAtoms = mmtfSummaryDataBean.getAllAtoms();
		int numBonds = mmtfSummaryDataBean.getNumBonds();
		List<Chain> allChains = mmtfSummaryDataBean.getAllChains();
		mmtfDecoderInterface.initStructure(numBonds, allAtoms.size(), MmtfUtils.getNumGroups(structure), allChains.size(), structure.nrModels(), structure.getPDBCode());
		// Generate the secondary structure
		MmtfUtils.calculateDsspSecondaryStructure(structure);
		// Get the header and the xtal info.
		PDBHeader pdbHeader = structure.getPDBHeader();
		PDBCrystallographicInfo xtalInfo = pdbHeader.getCrystallographicInfo();
		mmtfDecoderInterface.setHeaderInfo(pdbHeader.getRfree(), pdbHeader.getRwork(), pdbHeader.getResolution(), pdbHeader.getTitle(), MmtfUtils.dateToIsoString(pdbHeader.getDepDate()), 
				MmtfUtils.dateToIsoString(pdbHeader.getRelDate()), MmtfUtils.techniquesToStringArray(pdbHeader.getExperimentalTechniques()));
		mmtfDecoderInterface.setXtalInfo(MmtfUtils.getSpaceGroupAsString(xtalInfo.getSpaceGroup()), MmtfUtils.getUnitCellAsArray(xtalInfo), MmtfUtils.getNcsAsArray(xtalInfo.getNcsOperators()));
		// Store the bioassembly data
		storeBioassemblyInformation(chainIdToIndexMap, pdbHeader.getBioAssemblies());
		// Store the entity data
		storeEntityInformation(allChains, structure.getEntityInfos());
		// Now loop through the data structure
		for (int modelIndex=0; modelIndex<structure.nrModels(); modelIndex++) {
			List<Chain> modelChains = structure.getChains(modelIndex);
			// Set this model
			mmtfDecoderInterface.setModelInfo(modelIndex, modelChains.size());
			for(int chainInModelIndex=0; chainInModelIndex<modelChains.size(); chainInModelIndex++) {
				Chain chain = modelChains.get(chainInModelIndex);
				List<Group> groups = chain.getAtomGroups();
				List<Group> sequenceGroups = chain.getSeqResGroups();
				mmtfDecoderInterface.setChainInfo(chain.getId(), chain.getName(), groups.size());
				for(int groupInChainIndex=0; groupInChainIndex<groups.size(); groupInChainIndex++){
					Group group = groups.get(groupInChainIndex);
					List<Atom> atomsInGroup = MmtfUtils.getAtomsForGroup(group);
					ChemComp chemComp = group.getChemComp();
					Character insCode = group.getResidueNumber().getInsCode();
					if(insCode==null || insCode.equals(' ')){
						insCode=MmtfStructure.UNAVAILABLE_CHAR_VALUE;
					}
					char singleLetterCode = 'X';
					if (chemComp.getOne_letter_code().length()==1){
						singleLetterCode = chemComp.getOne_letter_code().charAt(0);
					}
					mmtfDecoderInterface.setGroupInfo(group.getPDBName(), group.getResidueNumber().getSeqNum(), insCode.charValue(), 
							chemComp.getType().toUpperCase(), atomsInGroup.size(), MmtfUtils.getNumBondsInGroup(atomsInGroup), singleLetterCode,
							sequenceGroups.indexOf(group), MmtfUtils.getSecStructType(group));
					for (Atom atom : atomsInGroup){
						char altLoc = MmtfStructure.UNAVAILABLE_CHAR_VALUE;
						if(atom.getAltLoc()!=null){
							if(atom.getAltLoc().charValue()!=' '){
								altLoc=atom.getAltLoc().charValue();
							}
						}
						mmtfDecoderInterface.setAtomInfo(atom.getName(), atom.getPDBserial(), altLoc, (float) atom.getX(), 
								(float) atom.getY(), (float) atom.getZ(), atom.getOccupancy(), 
								atom.getTempFactor(), atom.getElement().toString(), atom.getCharge());
						addBonds(atom, atomsInGroup, allAtoms);
					}
				}
			}
		}
		mmtfDecoderInterface.finalizeStructure();

	}

	/**
	 * Add the bonds for a given atom.
	 * @param atom the atom for which bonds are to be formed
	 * @param atomsInGroup the list of atoms in the group
	 * @param allAtoms the list of atoms in the whole structure
	 */
	private void addBonds(Atom atom, List<Atom> atomsInGroup, List<Atom> allAtoms) {
		if(atom.getBonds()==null){
			return;
		}
		for(Bond bond : atom.getBonds()) {
			// Now set the bonding information.
			Atom other = bond.getOther(atom);
			// If both atoms are in the group
			if (atomsInGroup.indexOf(other)!=-1){
				Integer firstBondIndex = atomsInGroup.indexOf(atom);
				Integer secondBondIndex = atomsInGroup.indexOf(other);
				// Don't add the same bond twice
				if(firstBondIndex>secondBondIndex){
					int bondOrder = bond.getBondOrder();
					mmtfDecoderInterface.setGroupBond(firstBondIndex, secondBondIndex, bondOrder);
				}
			}
			// Otherwise it's an inter group bond - so add it here
			else {
				Integer firstBondIndex = allAtoms.indexOf(atom);
				Integer secondBondIndex = allAtoms.indexOf(other);
				if(firstBondIndex>secondBondIndex){
					// Don't add the same bond twice
					int bondOrder = bond.getBondOrder();							
					mmtfDecoderInterface.setInterGroupBond(firstBondIndex, secondBondIndex, bondOrder);
				}
			}
		}		
	}


	/**
	 * Store the entity information for a given structure.
	 * @param allChains a list of all the chains in a structure
	 * @param entityInfos a list of the entity information
	 */
	private void storeEntityInformation(List<Chain> allChains, List<EntityInfo> entityInfos) {
		for (EntityInfo entityInfo : entityInfos) {
			String description = entityInfo.getDescription();
			String type;
			if (entityInfo.getType()==null){
				type = null;
			}
			else{
				type = entityInfo.getType().getEntityType();
			}
			List<Chain> entityChains = entityInfo.getChains();
			if (entityChains.isEmpty()){
				// Error mapping chain to entity
				System.err.println("ERROR MAPPING CHAIN TO ENTITY: "+description);
				mmtfDecoderInterface.setEntityInfo(new int[0], "", description, type);
				continue;
			}
			else{
				int[] chainIndices = new int[entityChains.size()];
				for (int i=0; i<entityChains.size(); i++) {
					chainIndices[i] = allChains.indexOf(entityChains.get(i));
				}
				Chain chain = entityChains.get(0);
				ChainImpl chainImpl;
				if (chain instanceof ChainImpl){
					chainImpl = (ChainImpl) entityChains.get(0);
				}
				else{
					throw new RuntimeException();
				}
				String sequence = chainImpl.getSeqResOneLetterSeq();
				mmtfDecoderInterface.setEntityInfo(chainIndices, sequence, description, type);
			}
		}		
	}


	/**
	 * Generate the bioassembly information on in the desired form.
	 * @param bioJavaStruct the Biojava structure
	 * @param header the header
	 */
	private void storeBioassemblyInformation(Map<String, Integer> chainIdToIndexMap, Map<Integer, BioAssemblyInfo> inputBioAss) {
		int bioAssemblyIndex = 0;
		for (Entry<Integer, BioAssemblyInfo> entry : inputBioAss.entrySet()) {
			Map<double[], int[]> transformMap = MmtfUtils.getTransformMap(entry.getValue(), chainIdToIndexMap);
			for(Entry<double[], int[]> transformEntry : transformMap.entrySet()) {
				mmtfDecoderInterface.setBioAssemblyTrans(bioAssemblyIndex, transformEntry.getValue(), transformEntry.getKey(), entry.getKey().toString());
			}
			bioAssemblyIndex++;
		}
	}

}
