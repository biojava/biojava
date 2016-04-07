package org.biojava.nbio.structure.io.mmtf;

import org.rcsb.mmtf.api.MmtfDecodedDataInterface;
import org.rcsb.mmtf.api.MmtfInputDataInterface;
import org.rcsb.mmtf.api.ObjectToByteArrayConverterInterface;
import org.rcsb.mmtf.dataholders.MmtfBean;

/**
 * Class to take Biojava structure data and covert to the DataApi for encoding.
 * @author Anthony Bradley
 *
 */
public class MmtfStructureEncoder implements MmtfInputDataInterface {

	@Override
	public void setxCoords(float[] xCoords) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setyCoords(float[] yCoords) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setzCoords(float[] zCoords) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setbFactors(float[] bFactors) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setOccupancies(float[] occupancies) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setAtomIds(int[] atomIds) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setAltLocIds(char[] altLocIds) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setInsCodes(char[] insCodes) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setGroupIds(int[] groupIds) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setGroupName(int groupInd, String groupName) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setNumAtomsInGroup(int groupInd, int numAtomsInGroup) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setGroupAtomNames(int groupInd, String[] groupAtomNames) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setGroupElementNames(int groupInd, String[] groupElements) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setGroupBondOrders(int groupInd, int[] groupBondOrders) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setGroupBondIndices(int groupInd, int[] groupBondIndices) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setGroupAtomCharges(int groupInd, int[] groupAtomCharges) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setGroupSingleLetterCode(int groupInd, char groupSingleLetterCode) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setGroupChemCompType(int groupInd, String groupChemCompType) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setGroupTypeIndices(int[] groupTypeIndices) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setGroupSequenceIndices(int[] groupSequenceIndices) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setChainIds(String[] chainIds) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setChainNames(String[] chainNames) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setChainsPerModel(int[] chainsPerModel) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setGroupsPerChain(int[] groupsPerChain) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setSpaceGroup(String spaceGroup) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setUnitCell(float[] unitCell) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setNumBioassemblies(int numBioassemblies) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setNumTransInBioassembly(int bioassemblyIndex, int numTransInBioassembly) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setChainIndexListForTransform(int bioassemblyIndex, int transformationIndex,
			int[] transChainIndexList) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setMatrixForTransform(int bioassemblyIndex, int transformationIndex, double[] transformationMatrix) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setInterGroupBondIndices(int[] interGroupBondIndices) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setInterGroupBondOrders(int[] interGroupBondOrders) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setMmtfVersion(String mmtfVersion) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setMmtfProducer(String mmtfProducer) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setNumEntities(int numEntities) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setEntityDescription(int entityInd, String entityDescription) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setEntityType(int entityInd, String entityType) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setEntityChainIndexList(int entityInd, int[] chainIndexList) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setEntitySequence(int entityInd, String sequence) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setStructureId(String structureId) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setNumModels(int numModels) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setNumChains(int numChains) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setNumGroups(int numGroups) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setNumAtoms(int numAtoms) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setRfree(float rFree) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setRwork(float rWork) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setResolution(float resolution) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setTitle(String title) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setExperimentalMethods(String[] experimentalMethods) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setDepositionDate(String depositionDate) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public MmtfDecodedDataInterface getDataAsDecodedDataInterface() {
		// TODO Auto-generated method stub
		return null;
	}

}
