package org.biojava.nbio.structure.io.mmtf;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.AminoAcid;
import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.BondImpl;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.Element;
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
import org.rcsb.mmtf.api.StructureDecoderInterface;
import org.rcsb.mmtf.dataholders.BioAssemblyData;
import org.rcsb.mmtf.dataholders.BioAssemblyTrans;


/**
 * A biojava specific structure inflator for MMTF.
 * Should be ported to biojava code.
 *
 * @author Anthony Bradley
 */
public class BioJavaStructureDecoder implements StructureDecoderInterface, Serializable {

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

	/** The chemical component group. */
	private ChemComp chemicalComponentGroup;

	/** All the atoms. */
	private List<Atom> allAtoms;

	/**
	 * Instantiates a new bio java structure decoder.
	 */
	public BioJavaStructureDecoder() {
		structure = new StructureImpl();
		modelNumber = 0;
		chemicalComponentGroup = new ChemComp();
		allAtoms  = new ArrayList<Atom>();
	}

	/**
	 * Gets the structure.
	 *
	 * @return the structure
	 */
	public final Structure getStructure() {
		return structure;
	}

	/* (non-Javadoc)
	 * @see org.rcsb.mmtf.decoder.StructureDecoderInterface#setModelCount(int)
	 */
	@Override
	public void setModelCount(final int inputModelCount) {
	}

	/* (non-Javadoc)
	 * @see org.rcsb.mmtf.decoder.StructureDecoderInterface#setModelInfo(int, int)
	 */
	@Override
	public final void setModelInfo(final int inputModelNumber,
			final int chainCount) {
		modelNumber = inputModelNumber;
		structure.addModel(new ArrayList<Chain>(chainCount));
	}

	/* (non-Javadoc)
	 * @see org.rcsb.mmtf.decoder.StructureDecoderInterface
	 * #setChainInfo(java.lang.String, int)
	 */
	@Override
	public final void setChainInfo(final String chainId, final int groupCount) {
		chain = new ChainImpl();
		chain.setChainID(chainId.trim());
		structure.addChain(chain, modelNumber);
	}


	/* (non-Javadoc)
	 * @see org.rcsb.mmtf.decoder.StructureDecoderInterface
	 * #setGroupInfo(java.lang.String, int, char, int, int)
	 */
	@Override
	public final void setGroupInfo(final String groupName, final int groupNumber,
			final char insertionCode, final String chemCompType, final int atomCount) {
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
		atomsInGroup = new ArrayList<Atom>();
		// Set the CC -> empty but not null
		group.setChemComp(chemicalComponentGroup);
		group.setPDBName(groupName);
		if (insertionCode == '?') {
			group.setResidueNumber(chain.getChainID().trim(), groupNumber, null);
		} else {
			group.setResidueNumber(chain.getChainID().trim(),
					groupNumber, insertionCode);
		}
		group.setAtoms(new ArrayList<Atom>(atomCount));
		if (polymerType != 0) {
			chain.getSeqResGroups().add(group);
		}
		if (atomCount > 0) {
			chain.addGroup(group);
		}
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
	public final void setAtomInfo(final String atomName,
			final int serialNumber, final char alternativeLocationId, final float x,
			final float y, final float z, final float occupancy,
			final float temperatureFactor,
			final String element, final int charge) {
		Atom atom = new AtomImpl();
		Group altGroup = null;
		atom.setPDBserial(serialNumber);
		atom.setName(atomName.trim());
		atom.setElement(Element.valueOfIgnoreCase(element));
		if (alternativeLocationId != '?') {
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
		allAtoms.add(atom);
	}

	/* (non-Javadoc)
	 * @see org.rcsb.mmtf.decoder.StructureDecoderInter
	 * face#setGroupBonds(int, int, int)
	 */
	@Override
	public final void setGroupBonds(final int indOne,
			final int indTwo, final int bondOrder) {
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
	public final void setInterGroupBonds(final int indOne,
			final int indTwo, final int bondOrder) {
		// Get the atom
		Atom atomOne = allAtoms.get(indOne);
		Atom atomTwo = allAtoms.get(indTwo);
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
	private Group getCorrectAltLocGroup(final Character altLoc) {

		// see if we know this altLoc already;
		List<Atom> atoms = group.getAtoms();
		if (atoms.size() > 0) {
			Atom a1 = atoms.get(0);
			// we are just adding atoms to the current group
			// probably there is a second group following later...
			if (a1.getAltLoc().equals(altLoc)) {
				return group;
			}
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
		altLocG.setAtoms(new ArrayList<Atom>());
		altLocG.getAltLocs().clear();
		group.addAltLoc(altLocG);
		return altLocG;

	}


	/* (non-Javadoc)
	 * @see org.rcsb.mmtf.decoder.StructureDecoderInterface#
	 * setXtalInfo(java.lang.String, java.util.List)
	 */
	@Override
	public final void setXtalInfo(final String spaceGroupString,
			final List<Float> unitCell) {
		// Now set the xtalographic information
		PDBCrystallographicInfo pci = new PDBCrystallographicInfo();
		SpaceGroup spaceGroup = SpaceGroup.parseSpaceGroup(spaceGroupString);
		pci.setSpaceGroup(spaceGroup);
		if (unitCell.size() > 0) {
			CrystalCell cell = new CrystalCell(unitCell.get(0), unitCell.get(1),
					unitCell.get(2), unitCell.get(3), unitCell.get(4), unitCell.get(5));
			pci.setCrystalCell(cell);
			structure.setCrystallographicInfo(pci);
		}
	}

	/* (non-Javadoc)
	 * @see org.rcsb.mmtf.decoder.StructureDecoder
	 * Interface#updateChainInfo(java.lang.String, int)
	 */
	@Override
	public final void updateChainInfo(final String chainId,
			final int groupCount) {
		for (Chain c: structure.getChains(modelNumber)) {
			if (c.getChainID().equals(chainId)) {
				chain = c;
				break;
			}
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
	public void setBioAssembly(Map<Integer, BioAssemblyData> inputBioassemblies) {

		PDBHeader pdbHeader = structure.getPDBHeader();
		// Get the bioassebly data
		Map<Integer, BioAssemblyInfo> bioAssemblies = new HashMap<>();
		int bioassemlyCounter = 0;
		for (Entry<Integer, BioAssemblyData> entry: inputBioassemblies.entrySet()) {
			bioassemlyCounter++;
			// Get the key and the value
			Integer key = entry.getKey();
			BioAssemblyData value = entry.getValue();
			// Make a biojava bioassembly object
			BioAssemblyInfo bioAssInfo = new  BioAssemblyInfo();
			bioAssInfo.setId(bioassemlyCounter);
			// Now get the new sizes
			List<BiologicalAssemblyTransformation> newList = new ArrayList<>();
			Integer transIdCounter =0;
			for (BioAssemblyTrans transform : value.getTransforms()) {
				transIdCounter++;
				// Now loop over the chains
				for(String currChainId : transform.getChainIdList()){
					BiologicalAssemblyTransformation bioAssTrans =
							new BiologicalAssemblyTransformation();
					bioAssTrans.setId(transIdCounter.toString());
					bioAssTrans.setChainId(currChainId);
					// Now set matrix
					Matrix4d mat4d = new Matrix4d(transform.getTransformation());
					bioAssTrans.setTransformationMatrix(mat4d);
					// Now add this
					newList.add(bioAssTrans);
				}
			}
			// Now set the transform list
			bioAssInfo.setTransforms(newList);
			// Now set this
			bioAssemblies.put(key, bioAssInfo);
		}
		// Now actually set them
		pdbHeader.setBioAssemblies(bioAssemblies);
		structure.setPDBHeader(pdbHeader);		
	}

	@Override
	public void cleanUpStructure() {
		// Ensure all altlocs have all atoms
		for (int i =0; i< structure.nrModels() ; i++){
			for (Chain chain : structure.getModel(i)) {
				for (Group group : chain.getAtomGroups()) {
					for (Group altLocGroup : group.getAltLocs()) { 
						for ( Atom groupAtom : group.getAtoms()) {
							// If this alt loc doesn't have this atom
							if (! altLocGroup.hasAtom(groupAtom.getName())) {
								// Make sure it's not a microheterogenity
								if (altLocGroup.getPDBName().equals(group.getPDBName())) {
									altLocGroup.addAtom(groupAtom);
								}
							}
						}
					}
				}
			}
		}
	}

	@Override
	public void prepareStructure(int totalNumAtoms, int totalNumGroups, int totalNumChains, int totalNumModels) {

	}


}
