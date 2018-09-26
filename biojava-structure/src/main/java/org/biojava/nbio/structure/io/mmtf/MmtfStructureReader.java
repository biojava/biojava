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

import java.io.Serializable;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
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
import org.biojava.nbio.structure.io.mmcif.chem.PolymerType;
import org.biojava.nbio.structure.io.mmcif.chem.ResidueType;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.biojava.nbio.structure.quaternary.BioAssemblyInfo;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.biojava.nbio.structure.xtal.CrystalCell;
import org.biojava.nbio.structure.xtal.SpaceGroup;
import org.rcsb.mmtf.api.StructureAdapterInterface;
import org.rcsb.mmtf.dataholders.MmtfStructure;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * A biojava specific structure inflator for MMTF.
 * Should be ported to biojava code.
 *
 * @author Anthony Bradley
 */
public class MmtfStructureReader implements StructureAdapterInterface, Serializable {

	/** The Constant serialVersionUID. */
	private static final long serialVersionUID = 6772030485225130853L;

	/** The logger */
	private static final Logger logger = LoggerFactory.getLogger(MmtfStructureReader.class);

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

	/** All the chains as a list of maps */
	private List<Map<String,Chain>> chainMap;

	private List<double[]> transformList;

	private int bioassIndex;

	private Map<String,String> chainSequenceMap;

	/**
	 * Instantiates a new bio java structure decoder.
	 */
	public MmtfStructureReader() {
		structure = new StructureImpl();
		modelNumber = 0;
		entityInfoList = new ArrayList<>();
		chainList = new ArrayList<>();
		chainMap = new ArrayList<>();
		transformList = new ArrayList<>();
		chainSequenceMap = new HashMap<>();
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
		// Number the remaining ones
		int counter =0;
		// Add the entity info
		for (EntityInfo entityInfo : entityInfoList) {
			counter++;
			entityInfo.setMolId(counter);
		}
		structure.setEntityInfos(entityInfoList);
		// Add the actual chains
		for(int i=0; i<chainMap.size(); i++) {
			// Now add the chain information
			Map<String, Chain> modelChainMap = chainMap.get(i);
			for(Chain modelChain : modelChainMap.values()){
				structure.addChain(modelChain, i);
				String sequence = chainSequenceMap.get(modelChain.getId());
				if (sequence == null) {
					logger.warn("Sequence is null for chain with asym_id {}. Most likely the chain is non-polymeric. Will not add seqres groups for it.", modelChain.getId());
					continue;
				}
				MmtfUtils.addSeqRes(modelChain, sequence);
			}
		}
		StructureTools.cleanUpAltLocs(structure);
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
		structure.addModel(new ArrayList<Chain>(chainCount));
		chainMap.add(new HashMap<>());
	}

	/* (non-Javadoc)
	 * @see org.rcsb.mmtf.decoder.StructureDecoderInterface
	 * #setChainInfo(java.lang.String, int)
	 */
	@Override
	public void setChainInfo(String chainId, String chainName, int groupCount) {
		// First check to see if the chain exists
		Map<String, Chain> modelChainMap = chainMap.get(modelNumber);
		if(modelChainMap.containsKey(chainId)){
			chain = modelChainMap.get(chainId);
		}
		// If we need to set a new chain do this
		else{
			chain = new ChainImpl();
			chain.setId(chainId.trim());
			chain.setName(chainName);
			chain.setAtomGroups(new ArrayList<>(groupCount));
			modelChainMap.put(chainId, chain);
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
		ResidueType residueType = ResidueType.getResidueTypeFromString(chemCompType);
		int polymerType = getGroupTypIndicator(residueType.polymerType);
		switch (polymerType) {
		case 1:
			AminoAcid aa = new AminoAcidImpl();
			// Now set the one letter code
			aa.setAminoType(singleLetterCode);
			group = aa;
			break;
		case 2:
			group = new NucleotideImpl();
			break;
		default:
			group = new HetatomImpl();
			break;
		}
		atomsInGroup = new ArrayList<Atom>();
		ChemComp chemComp = new ChemComp();
		chemComp.setOne_letter_code(String.valueOf(singleLetterCode));
		chemComp.setType(chemCompType.toUpperCase());
		chemComp.setResidueType(residueType);
		chemComp.setPolymerType(residueType.polymerType);
		group.setChemComp(chemComp);
		group.setPDBName(groupName);
		if (insertionCode == MmtfStructure.UNAVAILABLE_CHAR_VALUE) {
			group.setResidueNumber(chain.getName().trim(), groupNumber, null);
		} else {
			group.setResidueNumber(chain.getName().trim(),
					groupNumber, insertionCode);
		}
		group.setAtoms(new ArrayList<Atom>(atomCount));
		if (polymerType==1 || polymerType==2) {
			MmtfUtils.insertSeqResGroup(chain, group, sequenceIndexId);
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
			if (g.getResidueNumber().equals(group.getResidueNumber())) {
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
		if(alternativeLocationId==MmtfStructure.UNAVAILABLE_CHAR_VALUE){
			alternativeLocationId = ' ';
		}
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
				// And it's not a deuterated case.  'nanoheterogenity'
				if(!StructureTools.hasNonDeuteratedEquiv(atom,group)){
					group.addAtom(atom);
				}
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
	public void setXtalInfo(String spaceGroupString, float[] unitCell, double[][] ncsOperMatrixList) {
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

		pci.setNcsOperators(MmtfUtils.getNcsAsMatrix4d(ncsOperMatrixList));
	}


	/**
	 * Get the type of group (0,1 or 2) depending on whether it is an amino aicd (1), nucleic acid (2) or ligand (0)
	 * @param polymerType
	 * @return The type of group. (0,1 or 2) depending on whether it is an amino aicd (1), nucleic acid (2) or ligand (0)
	 */
	private int getGroupTypIndicator(PolymerType polymerType) {
		if(PolymerType.PROTEIN_ONLY.contains(polymerType)){
			return 1;
		}
		if(PolymerType.POLYNUCLEOTIDE_ONLY.contains(polymerType)){
			return 2;
		}
		return 0;
	}


	@Override
	public void setBioAssemblyTrans(int bioAssemblyId, int[] inputChainIndices, double[] inputTransform, String name) {
		// Biojava uses this as a one indexed id.
		bioAssemblyId++;
		if(bioassIndex!=bioAssemblyId){
			transformList = new ArrayList<>();
			bioassIndex = bioAssemblyId;
		}
		PDBHeader pdbHeader = structure.getPDBHeader();
		// Get the bioassembly data
		Map<Integer, BioAssemblyInfo> bioAssemblies = pdbHeader.getBioAssemblies();
		// Get the bioassembly itself (if it exists
		BioAssemblyInfo bioAssInfo;
		if (bioAssemblies.containsKey(bioAssemblyId)){
			bioAssInfo = bioAssemblies.get(bioAssemblyId);
		}
		else{
			bioAssInfo = new  BioAssemblyInfo();
			bioAssInfo.setTransforms(new ArrayList<BiologicalAssemblyTransformation>());
			bioAssemblies.put(bioAssemblyId, bioAssInfo);
			bioAssInfo.setId(bioAssemblyId);
		}

		for(int currChainIndex : inputChainIndices){
			BiologicalAssemblyTransformation bioAssTrans = new BiologicalAssemblyTransformation();
			Integer transId = transformList.indexOf(inputTransform)+1;
			if(transId==0){
				transformList.add(inputTransform);
				transId = transformList.indexOf(inputTransform)+1;
			}
			bioAssTrans.setId(transId.toString());
			// If it actually has an index - if it doesn't it is because the chain has no density.
			if (currChainIndex!=-1){
				bioAssTrans.setChainId(chainList.get(currChainIndex).getId());
			}
			else {
				continue;
			}
			// Now set matrix
			Matrix4d mat4d = new Matrix4d(inputTransform);
			bioAssTrans.setTransformationMatrix(mat4d);
			// Now add this
			bioAssInfo.getTransforms().add(bioAssTrans);
			// sort transformations into a unique order
			Collections.sort(bioAssInfo.getTransforms());
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
			chainSequenceMap.put(chainList.get(index).getId(), sequence);
		}
		entityInfo.setChains(chains);
		entityInfoList.add(entityInfo);
	}

	@Override
	public void setHeaderInfo(float rFree, float rWork, float resolution, String title, String depositionDate,
			String releaseDate, String[] experimentalMethods) {
		SimpleDateFormat formatter = new SimpleDateFormat("yyyy-MM-dd");
		// Get the pdb header
		PDBHeader pdbHeader = structure.getPDBHeader();
		pdbHeader.setTitle(title);
		pdbHeader.setResolution(resolution);
		pdbHeader.setRfree(rFree);
		pdbHeader.setRwork(rWork);
		// Now loop through the techniques and add them in
		if (experimentalMethods!=null) {
			for (String techniqueStr : experimentalMethods) {
				pdbHeader.setExperimentalTechnique(techniqueStr);
			}
		}
		// Set the dates
		if(depositionDate!=null){
			try {
				Date depDate = formatter.parse(depositionDate);
				pdbHeader.setDepDate(depDate);
			} catch (ParseException e) {
				logger.warn("Could not parse date string '{}', depositon date will be unavailable", depositionDate);
			}
		}
		else{
			pdbHeader.setDepDate(new Date(0));
		}
		if(releaseDate!=null){
			try {
				Date relDate = formatter.parse(releaseDate);
				pdbHeader.setRelDate(relDate);
			} catch (ParseException e) {
				logger.warn("Could not parse date string '{}', release date will be unavailable", releaseDate);
			}
		}
		else{
			pdbHeader.setRelDate(new Date(0));
		}
	}
}
