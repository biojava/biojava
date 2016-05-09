package org.biojava.nbio.structure.io.mmtf;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.AminoAcid;
import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.BondImpl;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.EntityInfo;
import org.biojava.nbio.structure.EntityType;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.HetatomImpl;
import org.biojava.nbio.structure.NucleotideImpl;
import org.biojava.nbio.structure.PDBCrystallographicInfo;
import org.biojava.nbio.structure.PDBHeader;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureImpl;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.biojava.nbio.structure.quaternary.BioAssemblyInfo;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.biojava.nbio.structure.xtal.CrystalCell;
import org.biojava.nbio.structure.xtal.SpaceGroup;
import org.rcsb.mmtf.api.StructureAdapterInterface;
import org.rcsb.mmtf.dataholders.MmtfStructure;


/**
 * A biojava specific structure inflator for MMTF.
 * Should be ported to biojava code.
 *
 * @author Anthony Bradley
 */
public class MmtfStructureReader implements StructureAdapterInterface, Serializable {

	/** The Constant serialVersionUID. */
	private static final long serialVersionUID = 6772030485225130853L;

	/** The structure. */
	private Structure structure;

	/** The model number. */
	private int modelNumber;

	/** The chain. */
	private Chain chain;

	/** The group. */
	private Group group;

	/** The atoms in a group. */
	private List<Atom> atomsInGroup;

	/** All the atoms. */
	private Atom[] allAtoms;
	private int atomCounter;

	/** The list of EntityInformation */
	private List<EntityInfo> entityInfoList;

	/** All the chains */
	private List<Chain> chainList; 

	/**
	 * Instantiates a new bio java structure decoder.
	 */
	public MmtfStructureReader() {
		structure = new StructureImpl();
		modelNumber = 0;
		entityInfoList = new ArrayList<>();
		chainList = new ArrayList<>();
	}

	/**
	 * Gets the structure.
	 *
	 * @return the structure
	 */
	public Structure getStructure() {
		return structure;
	}

	@Override
	public void finalizeStructure() {
		// Ensure all altlocs have all atoms
		StructureTools.cleanUpAltLocs(structure);
		// Number the remaining ones
		int counter =0;
		for (EntityInfo entityInfo : entityInfoList) {
			counter++;
			entityInfo.setMolId(counter);
		}
		structure.setEntityInfos(entityInfoList);
	}

	@Override
	public void initStructure(int totalNumBonds, int totalNumAtoms, int totalNumGroups, 
			int totalNumChains, int totalNumModels, String modelId) {
		structure.setPDBCode(modelId);
		allAtoms = new Atom[totalNumAtoms];
	}


	/* (non-Javadoc)
	 * @see org.rcsb.mmtf.decoder.StructureDecoderInterface#setModelInfo(int, int)
	 */
	@Override
	public void setModelInfo(int inputModelNumber,
			int chainCount) {
		modelNumber = inputModelNumber;
		structure.addModel(new ArrayList<>(chainCount));
	}

	/* (non-Javadoc)
	 * @see org.rcsb.mmtf.decoder.StructureDecoderInterface
	 * #setChainInfo(java.lang.String, int)
	 */
	@Override
	public void setChainInfo(String chainId, String chainName, int groupCount) {
		// First check to see if the chain exists
		boolean newChain = true;
		for (Chain c: structure.getChains(modelNumber)) {
			if (c.getChainID().equals(chainId)) {
				newChain = false;
				chain = c;
				break;
			}
		}
		// If we need to set a new chain do this
		if (newChain){
			chain = new ChainImpl();
			chain.setChainID(chainId.trim());
			structure.addChain(chain, modelNumber);
			chainList.add(chain);
		}
	}


	/* (non-Javadoc)
	 * @see org.rcsb.mmtf.decoder.StructureDecoderInterface
	 * #setGroupInfo(java.lang.String, int, char, int, int)
	 */
	@Override
	public void setGroupInfo(String groupName, int groupNumber,
			char insertionCode, String chemCompType, int atomCount, int bondCount, 
			char singleLetterCode, int sequenceIndexId, int secStructType) {
		// Get the polymer type
		int polymerType = getGroupTypIndicator(chemCompType);
		switch (polymerType) {
		case 1:
			AminoAcid aa = new AminoAcidImpl();
			// Now set the one letter code
			aa.setAminoType(StructureTools.get1LetterCodeAmino(groupName));
			group = aa;
			break;
		case 2:
			group = new NucleotideImpl();
			break;
		default:
			group = new HetatomImpl();
		}
		atomsInGroup = new ArrayList<>();
		// Set the CC -> empty but not null
		ChemComp chemComp = new ChemComp();
		chemComp.setOne_letter_code("" + singleLetterCode);
		group.setChemComp(chemComp);
		group.setPDBName(groupName);
		if (insertionCode == MmtfStructure.UNAVAILABLE_CHAR_VALUE) {
			group.setResidueNumber(chain.getChainID().trim(), groupNumber, null);
		} else {
			group.setResidueNumber(chain.getChainID().trim(),
					groupNumber, insertionCode);
		}
		group.setAtoms(new ArrayList<>(atomCount));
		if (polymerType != 0) {
			chain.getSeqResGroups().add(group);
		}
		if (atomCount > 0) {
			chain.addGroup(group);
		}
		MmtfUtils.setSecStructType(group, secStructType);
	}

	/**
	 * 
	 * @return
	 */
	private Group getGroupWithSameResNumButDiffPDBName() {
		// If this chain already has this group number
		for (Group g : chain.getAtomGroups() ) {
			if (g.getResidueNumber().getSeqNum()==group.getResidueNumber().getSeqNum()) {
				if( ! g.getPDBName().equals(group.getPDBName() )){
					return g;
				}
			}
		}
		return null;
	}

	/* (non-Javadoc)
	 * @see org.rcsb.mmtf.decoder.StructureDecoderInterface#
	 * setAtomInfo(java.lang.String, int, char, float, float,
	 * float, float, float, java.lang.String, int)
	 */
	@Override
	public void setAtomInfo(String atomName,
			int serialNumber, char alternativeLocationId, float x,
			float y, float z, float occupancy,
			float temperatureFactor,
			String element, int charge) {
		Atom atom = new AtomImpl();
		Group altGroup = null;
		atom.setPDBserial(serialNumber);
		atom.setName(atomName.trim());
		atom.setElement(Element.valueOfIgnoreCase(element));
		if (alternativeLocationId != ' ') {
			// Get the altGroup
			altGroup = getCorrectAltLocGroup(alternativeLocationId);
			atom.setAltLoc(alternativeLocationId);
		} else {
			atom.setAltLoc(Character.valueOf(' '));
		}
		atom.setX(x);
		atom.setY(y);
		atom.setZ(z);
		atom.setOccupancy(occupancy);
		atom.setTempFactor(temperatureFactor);
		atom.setCharge((short) charge);
		if (altGroup == null) {
			group.addAtom(atom);
		} else {
			altGroup.setChain(chain);
			altGroup.addAtom(atom);
		}

		// IF the main group doesn't have this atom
		if (!group.hasAtom(atom.getName())) {
			// If it's not a microheterogenity case
			if (group.getPDBName().equals(atom.getGroup().getPDBName())) {
				group.addAtom(atom);
			}
		}
		atomsInGroup.add(atom);
		allAtoms[atomCounter] = atom;
		atomCounter++;
	}

	/* (non-Javadoc)
	 * @see org.rcsb.mmtf.decoder.StructureDecoderInter
	 * face#setGroupBonds(int, int, int)
	 */
	@Override
	public void setGroupBond(int indOne,
			int indTwo, int bondOrder) {
		// Get the atom
		Atom atomOne = atomsInGroup.get(indOne);
		Atom atomTwo = atomsInGroup.get(indTwo);
		// set the new bond
		@SuppressWarnings("unused")
		BondImpl bond = new BondImpl(atomOne, atomTwo, bondOrder);
	}

	/* (non-Javadoc)
	 * @see org.rcsb.mmtf.decoder.StructureDecoder
	 * Interface#setInterGroupBonds(int, int, int)
	 */
	@Override
	public void setInterGroupBond(int indOne,
			int indTwo, int bondOrder) {
		// Get the atom
		Atom atomOne = allAtoms[indOne];
		Atom atomTwo = allAtoms[indTwo];
		// set the new bond
		@SuppressWarnings("unused")
		BondImpl bond = new BondImpl(atomOne, atomTwo, bondOrder);

	}


	/**
	 * Generates Alternate location groups.
	 *
	 * @param altLoc the alt loc
	 * @return the correct alt loc group
	 */
	private Group getCorrectAltLocGroup(Character altLoc) {
		// see if we know this altLoc already;
		List<Atom> atoms = group.getAtoms();
		if (atoms.size() > 0) {
			Atom a1 = atoms.get(0);
			// we are just adding atoms to the current group
			// probably there is a second group following later...
			if (a1.getAltLoc().equals(altLoc)) {
				return group;	}
		}

		// Get the altLocGroup
		Group altLocgroup = group.getAltLocGroup(altLoc);
		if (altLocgroup != null) {
			return altLocgroup;
		}
		// If the group already exists (microheterogenity).
		Group oldGroup = getGroupWithSameResNumButDiffPDBName();
		if (oldGroup!= null){
			Group altLocG = group;
			group = oldGroup;
			group.addAltLoc(altLocG);
			chain.getAtomGroups().remove(altLocG);
			return altLocG;
		}
		// no matching altLoc group found.
		// build it up.
		if (group.getAtoms().size() == 0) {
			return group;
		}
		Group altLocG = (Group) group.clone();
		// drop atoms from cloned group...
		// https://redmine.open-bio.org/issues/3307
		altLocG.setAtoms(new ArrayList<>());
		altLocG.getAltLocs().clear();
		group.addAltLoc(altLocG);
		return altLocG;

	}


	/* (non-Javadoc)
	 * @see org.rcsb.mmtf.decoder.StructureDecoderInterface#
	 * setXtalInfo(java.lang.String, java.util.List)
	 */
	@Override
	public void setXtalInfo(String spaceGroupString,
			float[] unitCell) {
		// Now set the xtalographic information
		PDBCrystallographicInfo pci = new PDBCrystallographicInfo();
		SpaceGroup spaceGroup = SpaceGroup.parseSpaceGroup(spaceGroupString);
		pci.setSpaceGroup(spaceGroup);
		if (unitCell.length > 0) {
			CrystalCell cell = new CrystalCell(unitCell[0], unitCell[1],
					unitCell[2], unitCell[3], unitCell[4], unitCell[5]);
			pci.setCrystalCell(cell);
			structure.setCrystallographicInfo(pci);
		}
	}


	/**
	 * Get the type of group (0,1 or 2) depending on whether it is an amino aicd (1), nucleic acid (2) or ligand (0)
	 * @param currentGroup
	 * @return The type of group. (0,1 or 2) depending on whether it is an amino aicd (1), nucleic acid (2) or ligand (0)
	 */
	private int getGroupTypIndicator(String currentGroupType) {
		// At the moment - peptide like is a HETATM group (consistent with biojava)
		if(currentGroupType.toUpperCase().equals("PEPTIDE-LIKE")){
			return 0;
		}
		// Again to correspond with Biojava - but I suspect we really want this to be 1
		if(currentGroupType.toUpperCase().equals("D-PEPTIDE LINKING")){
			return 0;
		}
		if(currentGroupType.toUpperCase().contains("PEPTIDE")){
			return 1;
		}
		if(currentGroupType.toUpperCase().contains("DNA") || currentGroupType.toUpperCase().contains("RNA")){
			return 2;
		}
		else{
			return 0;
		}
	}


	@Override
	public void setBioAssemblyTrans(int bioAssemblyId, int[] inputChainIndices, double[] inputTransform) {
		PDBHeader pdbHeader = structure.getPDBHeader();
		List<Chain> totChainList = new ArrayList<>(); 
		for (int i=0; i<structure.nrModels(); i++) { 
			totChainList.addAll(structure.getChains(i));
		}
		// Get the bioassembly data
		Map<Integer, BioAssemblyInfo> bioAssemblies = pdbHeader.getBioAssemblies();
		// Get the bioassembly itself (if it exists
		BioAssemblyInfo bioAssInfo;
		if (bioAssemblies.containsKey(bioAssemblyId)){
			bioAssInfo = bioAssemblies.get(bioAssemblyId);
		}
		else{
			bioAssInfo = new  BioAssemblyInfo();
			bioAssInfo.setTransforms(new ArrayList<>());
			bioAssemblies.put(bioAssemblyId, bioAssInfo);
			bioAssInfo.setId(bioAssemblyId);
		}

		for(int currChainIndex : inputChainIndices){
			BiologicalAssemblyTransformation bioAssTrans = new BiologicalAssemblyTransformation();
			Integer transId = bioAssInfo.getTransforms().size()+1;
			bioAssTrans.setId(transId.toString());
			// If it actually has an index - if it doesn't it is because the chain has no density.
			if (currChainIndex!=-1){
				bioAssTrans.setChainId(totChainList.get(currChainIndex).getChainID());
			}
			else {
				continue;
			}
			// Now set matrix
			Matrix4d mat4d = new Matrix4d(inputTransform);
			bioAssTrans.setTransformationMatrix(mat4d);
			// Now add this
			bioAssInfo.getTransforms().add(bioAssTrans);
		}
	}

	@Override
	public void setEntityInfo(int[] chainIndices, String sequence, String description, String type) {
		// First get the chains
		EntityInfo entityInfo = new EntityInfo();
		entityInfo.setDescription(description);
		entityInfo.setType(EntityType.entityTypeFromString(type));
		List<Chain> chains = new ArrayList<>(); 
		// Now loop through the chain ids and make a list of them
		for( int index : chainIndices) {
			chains.add(chainList.get(index));
			chainList.get(index).setEntityInfo(entityInfo);
		}
		entityInfo.setChains(chains);
		entityInfoList.add(entityInfo);
	}

	@Override
	public void setHeaderInfo(float rFree, float rWork, float resolution, String title, String depositionDate,
			String releaseDate, String[] experimnetalMethods) {
		// Get the pdb header
		PDBHeader pdbHeader = structure.getPDBHeader();
		pdbHeader.setTitle(title);
		pdbHeader.setResolution(resolution);
		pdbHeader.setRfree(rFree);
		// Now loop through the techniques and add them in
		for (String techniqueStr : experimnetalMethods) {
			pdbHeader.setExperimentalTechnique(techniqueStr);
		}
	}


}
