package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.AminoAcid;
import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.DBRef;
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.EntityInfo;
import org.biojava.nbio.structure.EntityType;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.HetatomImpl;
import org.biojava.nbio.structure.NucleotideImpl;
import org.biojava.nbio.structure.PDBCrystallographicInfo;
import org.biojava.nbio.structure.PDBHeader;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.SeqMisMatch;
import org.biojava.nbio.structure.SeqMisMatchImpl;
import org.biojava.nbio.structure.Site;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureImpl;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.BondMaker;
import org.biojava.nbio.structure.io.ChargeAdder;
import org.biojava.nbio.structure.io.EntityFinder;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.SeqRes2AtomAligner;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.model.DatabasePdbrevRecord;
import org.biojava.nbio.structure.quaternary.BioAssemblyInfo;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyBuilder;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.biojava.nbio.structure.xtal.CrystalCell;
import org.biojava.nbio.structure.xtal.SpaceGroup;
import org.biojava.nbio.structure.xtal.SymoplibParser;
import org.rcsb.cif.model.FloatColumn;
import org.rcsb.cif.model.IntColumn;
import org.rcsb.cif.model.StrColumn;
import org.rcsb.cif.model.generated.AtomSite;
import org.rcsb.cif.model.generated.AtomSites;
import org.rcsb.cif.model.generated.AuditAuthor;
import org.rcsb.cif.model.generated.Cell;
import org.rcsb.cif.model.generated.ChemComp;
import org.rcsb.cif.model.generated.ChemCompBond;
import org.rcsb.cif.model.generated.DatabasePDBRemark;
import org.rcsb.cif.model.generated.DatabasePDBRev;
import org.rcsb.cif.model.generated.DatabasePDBRevRecord;
import org.rcsb.cif.model.generated.Entity;
import org.rcsb.cif.model.generated.EntityPoly;
import org.rcsb.cif.model.generated.EntityPolySeq;
import org.rcsb.cif.model.generated.EntitySrcGen;
import org.rcsb.cif.model.generated.EntitySrcNat;
import org.rcsb.cif.model.generated.Exptl;
import org.rcsb.cif.model.generated.PdbxAuditRevisionHistory;
import org.rcsb.cif.model.generated.PdbxChemCompIdentifier;
import org.rcsb.cif.model.generated.PdbxDatabaseStatus;
import org.rcsb.cif.model.generated.PdbxEntityDescriptor;
import org.rcsb.cif.model.generated.PdbxEntitySrcSyn;
import org.rcsb.cif.model.generated.PdbxMolecule;
import org.rcsb.cif.model.generated.PdbxMoleculeFeatures;
import org.rcsb.cif.model.generated.PdbxNonpolyScheme;
import org.rcsb.cif.model.generated.PdbxReferenceEntityLink;
import org.rcsb.cif.model.generated.PdbxReferenceEntityList;
import org.rcsb.cif.model.generated.PdbxReferenceEntityPolyLink;
import org.rcsb.cif.model.generated.PdbxStructAssembly;
import org.rcsb.cif.model.generated.PdbxStructAssemblyGen;
import org.rcsb.cif.model.generated.PdbxStructModResidue;
import org.rcsb.cif.model.generated.PdbxStructOperList;
import org.rcsb.cif.model.generated.Refine;
import org.rcsb.cif.model.generated.Struct;
import org.rcsb.cif.model.generated.StructAsym;
import org.rcsb.cif.model.generated.StructConf;
import org.rcsb.cif.model.generated.StructConn;
import org.rcsb.cif.model.generated.StructConnType;
import org.rcsb.cif.model.generated.StructKeywords;
import org.rcsb.cif.model.generated.StructNcsOper;
import org.rcsb.cif.model.generated.StructRef;
import org.rcsb.cif.model.generated.StructRefSeq;
import org.rcsb.cif.model.generated.StructRefSeqDif;
import org.rcsb.cif.model.generated.StructSheetRange;
import org.rcsb.cif.model.generated.StructSite;
import org.rcsb.cif.model.generated.StructSiteGen;
import org.rcsb.cif.model.generated.Symmetry;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.Matrix4d;
import java.time.LocalDate;
import java.time.ZoneId;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeFormatterBuilder;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

// TODO detach the impl from the redundant mmCIF impl

/**
 * An implementation of a CifFileConsumer for BioJava. Will process the information provided by a CifFile instance and
 * use it to build up a {@link Structure} object. The implementation is for the most part really close to that in
 * {@link org.biojava.nbio.structure.io.mmcif.SimpleMMcifConsumer} and associated classes. The main difference is that
 * all internally used model classes are generated from the MMCIF schema and a standardized interface to CifFile data is
 * provided. This allows to readily parse files in CIF format as well as binary CIF (BCIF) format.
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @since 5.3.0
 */
class CifFileConsumerImpl implements CifFileConsumer<Structure> {
    private static final Logger logger = LoggerFactory.getLogger(CifFileConsumerImpl.class);
    private static final DateTimeFormatter DATE_FORMAT = new DateTimeFormatterBuilder()
            .parseCaseInsensitive()
            .appendPattern("yyyy-MM-dd")
            .toFormatter(Locale.US);

    private Structure structure;
    private Chain currentChain;
    private Group currentGroup;
    private List<List<Chain>> allModels;
    private List<Chain> currentModel;
    private PDBHeader pdbHeader;
    private String currentNmrModelNumber;
    private List<Chain> entityChains;

    private Entity entity;
    private EntityPoly entityPoly;
    private EntitySrcGen entitySrcGen;
    private EntitySrcNat entitySrcNat;
    private PdbxEntitySrcSyn entitySrcSyn;
    private List<Chain> seqResChains;
    private PdbxStructAssembly structAssembly;
    private PdbxStructAssemblyGen structAssemblyGen;
    private StructAsym structAsym;
    private StructConn structConn;
    private StructNcsOper structNcsOper;
    private PdbxStructOperList structOpers;
    private StructRef structRef;
    private StructRefSeqDif structRefSeqDif;
    private StructSiteGen structSiteGen;

    private Map<String, String> asymId2entityId;
    private Map<String, String> asymId2authorId;
    private Matrix4d parsedScaleMatrix;

    private FileParsingParameters params;

    public CifFileConsumerImpl(FileParsingParameters params) {
        this.params = params;
    }

    @Override
    public void prepare() {
        this.structure = new StructureImpl();
        this.pdbHeader = new PDBHeader();
        structure.setPDBHeader(pdbHeader);

        this.allModels = new ArrayList<>();
        this.currentModel = new ArrayList<>();

        this.seqResChains  = new ArrayList<>();
        this.asymId2entityId = new HashMap<>();
        this.asymId2authorId = new HashMap<>();

        this.entityChains = new ArrayList<>();
    }

    @Override
    public void consumeAtomSite(AtomSite atomSite) {
        if (params.isHeaderOnly()) {
            return;
        }

        StrColumn labelAsymId = atomSite.getLabelAsymId();
        StrColumn authAsymId = atomSite.getAuthAsymId();

        StrColumn groupPDB = atomSite.getGroupPDB();
        IntColumn authSeqId = atomSite.getAuthSeqId();

        StrColumn labelCompId = atomSite.getLabelCompId();

        IntColumn id = atomSite.getId();
        StrColumn labelAtomId = atomSite.getLabelAtomId();

        FloatColumn cartnX = atomSite.getCartnX();
        FloatColumn cartnY = atomSite.getCartnY();
        FloatColumn cartnZ = atomSite.getCartnZ();

        FloatColumn occupancy = atomSite.getOccupancy();
        FloatColumn bIsoOrEquiv = atomSite.getBIsoOrEquiv();

        StrColumn labelAltId = atomSite.getLabelAltId();
        StrColumn typeSymbol = atomSite.getTypeSymbol();

        StrColumn pdbxPDBInsCode = atomSite.getPdbxPDBInsCode();
        IntColumn labelSeqId = atomSite.getLabelSeqId();
        IntColumn pdbx_pdb_model_num = atomSite.getPdbxPDBModelNum();

        for (int atomIndex = 0; atomIndex < atomSite.getRowCount(); atomIndex++) {
            boolean startOfNewChain = false;
            Character oneLetterCode = StructureTools.get1LetterCodeAmino(labelCompId.get(atomIndex));

            boolean isHetAtmInFile = false;
            if (!"ATOM".equals(groupPDB.get(atomIndex))) {
                if (oneLetterCode != null && oneLetterCode.equals(StructureTools.UNKNOWN_GROUP_LABEL)) {
                    oneLetterCode = null;
                }

                isHetAtmInFile = true;
            }

            String insCodeString = pdbxPDBInsCode.get(atomIndex);
            Character insCode = null;
            if (insCodeString != null && !insCodeString.isEmpty() && !"?".equals(insCodeString)) {
                insCode = insCodeString.charAt(0);
            }

            // non polymer chains (ligands and small molecules) will have a label_seq_id set to '.'
            long seqId = labelSeqId.get(atomIndex);

            String nmrModelNumber = pdbx_pdb_model_num.getStringData(atomIndex);

            if (currentNmrModelNumber == null) {
                currentNmrModelNumber = nmrModelNumber;
            }
            if (!currentNmrModelNumber.equals(nmrModelNumber)) {
                currentNmrModelNumber = nmrModelNumber;

                if (currentChain != null) {
                    currentChain.addGroup(currentGroup);
                    currentGroup.trimToSize();
                }

                allModels.add(currentModel);
                currentModel = new ArrayList<>();
                currentChain = null;
                currentGroup = null;
            }

            String asymId = labelAsymId.get(atomIndex);
            String authId = authAsymId.get(atomIndex);
            if (currentChain == null) {
                currentChain = new ChainImpl();
                currentChain.setName(authId);
                currentChain.setId(asymId);
                currentModel.add(currentChain);
                startOfNewChain = true;
            }

            if (!asymId.equals(currentChain.getId())) {
                startOfNewChain = true;

                currentChain.addGroup(currentGroup);

                Optional<Chain> testChain = currentModel.stream()
                        .filter(chain -> chain.getId().equals(asymId))
                        .findFirst();

                if (testChain.isPresent()) {
                    currentChain = testChain.get();
                } else {
                    currentChain = new ChainImpl();
                    currentChain.setName(authId);
                    currentChain.setId(asymId);
                }

                if (!currentModel.contains(currentChain)) {
                    currentModel.add(currentChain);
                }
            }

            ResidueNumber residueNumber = new ResidueNumber(authId, authSeqId.get(atomIndex), insCode);

            String recordName = groupPDB.get(atomIndex);
            String compId = labelCompId.get(atomIndex);
            if (currentGroup == null) {
                currentGroup = createGroup(recordName, oneLetterCode, compId, seqId);
                currentGroup.setResidueNumber(residueNumber);
                currentGroup.setPDBName(compId);
                currentGroup.setHetAtomInFile(isHetAtmInFile);
            }

            Group altGroup = null;
            String altLocation = labelAltId.get(atomIndex);

            if (startOfNewChain) {
                currentGroup = createGroup(recordName, oneLetterCode, compId, seqId);
                currentGroup.setResidueNumber(residueNumber);
                currentGroup.setPDBName(compId);
                currentGroup.setHetAtomInFile(isHetAtmInFile);
            } else {
                if (!residueNumber.equals(currentGroup.getResidueNumber())) {
                    currentChain.addGroup(currentGroup);
                    currentGroup.trimToSize();
                    currentGroup = createGroup(recordName, oneLetterCode, compId, seqId);
                    currentGroup.setPDBName(compId);
                    currentGroup.setResidueNumber(residueNumber);
                    currentGroup.setHetAtomInFile(isHetAtmInFile);
                } else {
                    if (altLocation != null && !altLocation.isEmpty() && !altLocation.equals(".")) {
                        altGroup = getAltLocGroup(recordName, altLocation.charAt(0), oneLetterCode, compId, seqId);
                        if (altGroup.getChain() == null) {
                            altGroup.setChain(currentChain);
                        }
                    }
                }
            }

            if (params.isParseCAOnly()) {
                if (!labelAtomId.get(atomIndex).equals(StructureTools.CA_ATOM_NAME) && "C".equals(typeSymbol.get(atomIndex))) {
                    continue;
                }
            }

            Atom atom = new AtomImpl();

            atom.setPDBserial(id.get(atomIndex));
            atom.setName(labelAtomId.get(atomIndex));

            atom.setX(cartnX.get(atomIndex));
            atom.setY(cartnY.get(atomIndex));
            atom.setZ(cartnZ.get(atomIndex));

            atom.setOccupancy((float) occupancy.get(atomIndex));
            atom.setTempFactor((float) bIsoOrEquiv.get(atomIndex));

            if (altLocation == null || altLocation.isEmpty() || altLocation.equals(".")) {
                atom.setAltLoc(' ');
            } else {
                atom.setAltLoc(altLocation.charAt(0));
            }

            String ts = typeSymbol.get(atomIndex);
            try {
                Element element = Element.valueOfIgnoreCase(ts);
                atom.setElement(element);
            }  catch (IllegalArgumentException e) {
                logger.info("Element {} was not recognised as a BioJava-known element, the element will be " +
                        "represented as the generic element {}", typeSymbol, Element.R.name());
                atom.setElement(Element.R);
            }

            if (altGroup != null) {
                altGroup.addAtom(atom);
            } else {
                currentGroup.addAtom(atom);
            }

            String atomName = atom.getName();
            if (!currentGroup.hasAtom(atomName)) {
                if (currentGroup.getPDBName().equals(atom.getGroup().getPDBName())) {
                    if (!StructureTools.hasNonDeuteratedEquiv(atom, currentGroup)) {
                        currentGroup.addAtom(atom);
                    }
                }
            }
        }
    }

    private Group getAltLocGroup(String recordName, Character altLoc, Character oneLetterCode, String threeLetterCode,
                                 long seqId) {
        List<Atom> atoms = currentGroup.getAtoms();
        if (atoms.size() > 0) {
            if (atoms.get(0).getAltLoc().equals(altLoc)) {
                return currentGroup;
            }
        }

        List<Group> altLocs = currentGroup.getAltLocs();
        for (Group altLocGroup : altLocs) {
            atoms = altLocGroup.getAtoms();
            if (atoms.size() > 0) {
                for (Atom a1 : atoms) {
                    if (a1.getAltLoc().equals(altLoc)) {
                        return altLocGroup;
                    }
                }
            }
        }

        if (threeLetterCode.equals(currentGroup.getPDBName())) {
            if (currentGroup.getAtoms().isEmpty()) {
                return currentGroup;
            }

            Group altLocGroup = (Group) currentGroup.clone();
            altLocGroup.setAtoms(new ArrayList<>());
            altLocGroup.getAltLocs().clear();
            currentGroup.addAltLoc(altLocGroup);
            return altLocGroup;
        }

        Group altLocGroup = createGroup(recordName, oneLetterCode, threeLetterCode, seqId);
        altLocGroup.setPDBName(threeLetterCode);
        altLocGroup.setResidueNumber(currentGroup.getResidueNumber());
        currentGroup.addAltLoc(altLocGroup);
        return altLocGroup;
    }

    private Group createGroup(String record, Character oneLetterCode, String threeLetterCode, long seqId) {
        Group group = ChemCompGroupFactory.getGroupFromChemCompDictionary(threeLetterCode);
        if (group != null && !group.getChemComp().isEmpty()) {
            if (group instanceof AminoAcidImpl) {
                AminoAcidImpl aminoAcid = (AminoAcidImpl) group;
                aminoAcid.setId(seqId);
            } else if (group instanceof NucleotideImpl) {
                NucleotideImpl nucleotide = (NucleotideImpl) group;
                nucleotide.setId(seqId);
            } else if (group instanceof HetatomImpl) {
                HetatomImpl hetatom = (HetatomImpl) group;
                hetatom.setId(seqId);
            }
            return group;
        }

        if ("ATOM".equals(record)) {
            if (StructureTools.isNucleotide(threeLetterCode)) {
                NucleotideImpl nucleotide = new NucleotideImpl();
                group = nucleotide;
                nucleotide.setId(seqId);
            } else if (oneLetterCode == null || oneLetterCode == StructureTools.UNKNOWN_GROUP_LABEL) {
                HetatomImpl hetatom = new HetatomImpl();
                group = hetatom;
                hetatom.setId(seqId);
            } else {
                AminoAcidImpl aminoAcid = new AminoAcidImpl();
                group = aminoAcid;
                aminoAcid.setAminoType(oneLetterCode);
                aminoAcid.setId(seqId);
            }
        } else {
            if (StructureTools.isNucleotide(threeLetterCode)) {
                NucleotideImpl nucleotide = new NucleotideImpl();
                group = nucleotide;
                nucleotide.setId(seqId);
            } else if (oneLetterCode != null) {
                AminoAcidImpl aminoAcid = new AminoAcidImpl();
                group = aminoAcid;
                aminoAcid.setAminoType(oneLetterCode);
                aminoAcid.setId(seqId);
            } else {
                HetatomImpl hetatom = new HetatomImpl();
                hetatom.setId(seqId);
                group = hetatom;
            }
        }
        return group;
    }

    @Override
    public void consumeAtomSites(AtomSites atomSites) {
        // no atom sites present
        if (!atomSites.isDefined() || atomSites.getRowCount() == 0) {
            return;
        }

        try {
            parsedScaleMatrix = new Matrix4d(
                    atomSites.getFractTransfMatrix11().get(0),
                    atomSites.getFractTransfMatrix12().get(0),
                    atomSites.getFractTransfMatrix13().get(0),
                    atomSites.getFractTransfVector1().get(0),

                    atomSites.getFractTransfMatrix21().get(0),
                    atomSites.getFractTransfMatrix22().get(0),
                    atomSites.getFractTransfMatrix23().get(0),
                    atomSites.getFractTransfVector2().get(0),

                    atomSites.getFractTransfMatrix31().get(0),
                    atomSites.getFractTransfMatrix32().get(0),
                    atomSites.getFractTransfMatrix33().get(0),
                    atomSites.getFractTransfVector3().get(0),

                    0,
                    0,
                    0,
                    1
            );
        } catch (NumberFormatException e) {
            logger.warn("Some values in _atom_sites.fract_transf_matrix or _atom_sites.fract_transf_vector could not " +
                    "be parsed as numbers. Can't check whether coordinate frame convention is correct! Error: {}",
                    e.getMessage());
            structure.getPDBHeader().getCrystallographicInfo().setNonStandardCoordFrameConvention(false);
        }
    }

    @Override
    public void consumeAuditAuthor(AuditAuthor auditAuthor) {
        for (int rowIndex = 0; rowIndex < auditAuthor.getRowCount(); rowIndex++) {
            String name = auditAuthor.getName().get(rowIndex);

            StringBuilder last = new StringBuilder();
            StringBuilder initials = new StringBuilder();
            boolean afterComma = false;
            for (char c : name.toCharArray()) {
                if (c == ' ') {
                    continue;
                }
                if (c == ',') {
                    afterComma = true;
                    continue;
                }

                if (afterComma) {
                    initials.append(c);
                } else {
                    last.append(c);
                }
            }

            StringBuilder newaa = new StringBuilder();
            newaa.append(initials);
            newaa.append(last);

            String auth = pdbHeader.getAuthors();
            if (auth == null) {
                pdbHeader.setAuthors(newaa.toString());
            } else {
                auth += "," + newaa.toString();
                pdbHeader.setAuthors(auth);
            }
        }
    }

    @Override
    public void consumeCell(Cell cell) {
        if (!cell.isDefined() || cell.getRowCount() == 0) {
            return;
        }

        try {
            float a = (float) cell.getLengthA().get(0);
            float b = (float) cell.getLengthB().get(0);
            float c = (float) cell.getLengthC().get(0);
            float alpha = (float) cell.getAngleAlpha().get(0);
            float beta = (float) cell.getAngleBeta().get(0);
            float gamma = (float) cell.getAngleGamma().get(0);

            CrystalCell crystalCell = new CrystalCell();
            crystalCell.setA(a);
            crystalCell.setB(b);
            crystalCell.setC(c);
            crystalCell.setAlpha(alpha);
            crystalCell.setBeta(beta);
            crystalCell.setGamma(gamma);

            if (!crystalCell.isCellReasonable()) {
                // If the entry describes a structure determined by a technique other than X-ray crystallography,
                // cell is (sometimes!) a = b = c = 1.0, alpha = beta = gamma = 90 degrees
                // if so we don't add and CrystalCell will be null
                logger.debug("The crystal cell read from file does not have reasonable dimensions (at least one " +
                        "dimension is below {}), discarding it.", CrystalCell.MIN_VALID_CELL_SIZE);
                return;
            }

            structure.getPDBHeader()
                    .getCrystallographicInfo()
                    .setCrystalCell(crystalCell);

        } catch (NumberFormatException e){
            structure.getPDBHeader()
                    .getCrystallographicInfo()
                    .setCrystalCell(null);
            logger.info("could not parse some cell parameters ({}), ignoring _cell", e.getMessage());
        }
    }

    @Override
    public void consumeChemComp(ChemComp chemComp) {
        // TODO not impled in ref
    }

    @Override
    public void consumeChemCompBond(ChemCompBond chemCompBond) {
        // TODO not impled in ref
    }

    @Override
    public void consumeDatabasePDBremark(DatabasePDBRemark databasePDBremark) {
        for (int rowIndex = 0; rowIndex < databasePDBremark.getRowCount(); rowIndex++) {
            int id = databasePDBremark.getId().get(rowIndex);
            if (id == 2) {
                String line = databasePDBremark.getText().get(rowIndex);
                int i = line.indexOf("ANGSTROM");

                if (i > 5) {
                    // line contains ANGSTROM info...
                    String resolution = line.substring(i - 5, i).trim();
                    // convert string to float
                    try {
                        float res = Float.parseFloat(resolution);
                        pdbHeader.setResolution(res);
                    } catch (NumberFormatException e) {
                        logger.info("could not parse resolution from line and ignoring it {}", line);
                        return;
                    }
                }
            }
        }
    }

    private Date convert(LocalDate localDate) {
        return Date.from(localDate.atStartOfDay().atZone(ZoneId.systemDefault()).toInstant());
    }

    @Override
    public void consumeDatabasePDBrev(DatabasePDBRev databasePDBrev) {
        logger.debug("got a database revision:" + databasePDBrev);

        for (int rowIndex = 0; rowIndex < databasePDBrev.getRowCount(); rowIndex++) {
            if (databasePDBrev.getNum().get(rowIndex) == 1) {
                String dateOriginal = databasePDBrev.getDateOriginal().get(rowIndex);
                pdbHeader.setDepDate(convert(LocalDate.parse(dateOriginal, DATE_FORMAT)));

                String date = databasePDBrev.getDate().get(rowIndex);
                pdbHeader.setRelDate(convert(LocalDate.parse(date, DATE_FORMAT)));
            } else {
                String dbrev = databasePDBrev.getDate().get(rowIndex);
                pdbHeader.setModDate(convert(LocalDate.parse(dbrev, DATE_FORMAT)));
            }
        }
    }

    @Override
    public void consumeDatabasePDBrevRecord(DatabasePDBRevRecord databasePDBrevRecord) {
        List<DatabasePdbrevRecord> revRecords = pdbHeader.getRevisionRecords();
        if (revRecords == null) {
            revRecords = new ArrayList<>();
            pdbHeader.setRevisionRecords(revRecords);
        }

        revRecords.addAll(convert(databasePDBrevRecord));
    }

    private List<DatabasePdbrevRecord> convert(DatabasePDBRevRecord databasePDBrevRecord) {
        List<DatabasePdbrevRecord> revRecords = new ArrayList<>();
        for (int rowIndex = 0; rowIndex < databasePDBrevRecord.getRowCount(); rowIndex++) {
            DatabasePdbrevRecord revRecord = new DatabasePdbrevRecord();
            revRecord.setDetails(databasePDBrevRecord.getDetails().get(rowIndex));
            revRecord.setRev_num(databasePDBrevRecord.getRevNum().getStringData(rowIndex));
            revRecord.setType(databasePDBrevRecord.getType().get(rowIndex));
            revRecords.add(revRecord);
        }
        return revRecords;
    }

    @Override
    public void consumeEntity(Entity entity) {
        this.entity = entity;
    }

    @Override
    public void consumeEntityPoly(EntityPoly entityPoly) {
        this.entityPoly = entityPoly;
    }

    @Override
    public void consumeEntitySrcGen(EntitySrcGen entitySrcGen) {
        this.entitySrcGen = entitySrcGen;
    }

    @Override
    public void consumeEntitySrcNat(EntitySrcNat entitySrcNat) {
        this.entitySrcNat = entitySrcNat;
    }

    @Override
    public void consumeEntitySrcSyn(PdbxEntitySrcSyn entitySrcSyn) {
        this.entitySrcSyn = entitySrcSyn;
    }

    @Override
    public void consumeEntityPolySeq(EntityPolySeq entityPolySeq) {
        for (int rowIndex = 0; rowIndex < entityPolySeq.getRowCount(); rowIndex++) {
            Chain entityChain = getEntityChain(entityPolySeq.getEntityId().get(rowIndex));

            // first we check through the chemcomp provider, if it fails we do some heuristics to guess the type of group
            // TODO some of this code is analogous to getNewGroup() and we should try to unify them - JD 2016-03-08

            Group g = ChemCompGroupFactory.getGroupFromChemCompDictionary(entityPolySeq.getMonId().get(rowIndex));
            //int seqId = Integer.parseInt(entityPolySeq.getNum());
            if (g != null && !g.getChemComp().isEmpty()) {
                if (g instanceof AminoAcidImpl) {
                    AminoAcidImpl aa = (AminoAcidImpl) g;
                    aa.setRecordType(AminoAcid.SEQRESRECORD);
                }
            } else {
                if (entityPolySeq.getMonId().get(rowIndex).length() == 3 &&
                        StructureTools.get1LetterCodeAmino(entityPolySeq.getMonId().get(rowIndex)) != null) {
                    AminoAcidImpl a = new AminoAcidImpl();
                    a.setRecordType(AminoAcid.SEQRESRECORD);
                    Character code1 = StructureTools.get1LetterCodeAmino(entityPolySeq.getMonId().get(rowIndex));
                    a.setAminoType(code1);
                    g = a;

                } else if (StructureTools.isNucleotide(entityPolySeq.getMonId().get(rowIndex))) {
                    // the group is actually a nucleotide group...
                    g = new NucleotideImpl();
                } else {
                    logger.debug("Residue {} {} is not a standard aminoacid or nucleotide, will create a het group " +
                            "for it", entityPolySeq.getNum().get(rowIndex), entityPolySeq.getMonId().get(rowIndex));
                    g = new HetatomImpl();
                }
            }
            // at this stage we don't know about author residue numbers (insertion codes)
            // we abuse now the ResidueNumber field setting the internal residue numbers (label_seq_id, strictly
            // sequential and follow the seqres sequence 1 to n)
            // later the actual ResidueNumbers (author residue numbers) have to be corrected in alignSeqRes()
            g.setResidueNumber(ResidueNumber.fromString(entityPolySeq.getNum().getStringData(rowIndex)));
            g.setPDBName(entityPolySeq.getMonId().get(rowIndex));
            entityChain.addGroup(g);
        }
    }

    private Chain getEntityChain(String entityId) {
        for (Chain chain : entityChains) {
            if (chain.getId().equals(entityId)) {
                return chain;
            }
        }

        // does not exist yet, so create...
        Chain chain = new ChainImpl();
        chain.setId(entityId);
        entityChains.add(chain);

        return chain;
    }

    @Override
    public void consumeExptl(Exptl exptl) {
        for (int rowIndex = 0; rowIndex < exptl.getRowCount(); rowIndex++) {
            pdbHeader.setExperimentalTechnique(exptl.getMethod().get(rowIndex));
        }
    }

    @Override
    public void consumePdbxAuditRevisionHistory(PdbxAuditRevisionHistory pdbxAuditRevisionHistory) {
        for (int rowIndex = 0; rowIndex < pdbxAuditRevisionHistory.getRowCount(); rowIndex++) {
            // first entry in revision history is the release date
            if (pdbxAuditRevisionHistory.getOrdinal().get(rowIndex) == 1) {
                String release = pdbxAuditRevisionHistory.getRevisionDate().get(rowIndex);
                pdbHeader.setRelDate(convert(LocalDate.parse(release, DATE_FORMAT)));
            } else {
                // all other dates are revision dates;
                // since this method may be called multiple times,
                // the last revision date will "stick"
                String revision = pdbxAuditRevisionHistory.getRevisionDate().get(rowIndex);
                pdbHeader.setModDate(convert(LocalDate.parse(revision, DATE_FORMAT)));
            }
        }
    }

    @Override
    public void consumePdbxChemCompIdentifier(PdbxChemCompIdentifier pdbxChemCompIdentifier) {
        // TODO not impled in ref
    }

    @Override
    public void consumePdbxDatabaseStatus(PdbxDatabaseStatus pdbxDatabaseStatus) {
        for (int rowIndex = 0; rowIndex < pdbxDatabaseStatus.getRowCount(); rowIndex++) {
            // the deposition date field is only available in mmCIF 5.0
            StrColumn recvdInitialDepositionDate = pdbxDatabaseStatus.getRecvdInitialDepositionDate();
            if (recvdInitialDepositionDate.isDefined()) {
                String deposition = recvdInitialDepositionDate.get(rowIndex);
                pdbHeader.setDepDate(convert(LocalDate.parse(deposition, DATE_FORMAT)));
            }
        }
    }

    @Override
    public void consumePdbxEntityDescriptor(PdbxEntityDescriptor pdbxEntityDescriptor) {
        // TODO not considered in ref
    }

    @Override
    public void consumePdbxMolecule(PdbxMolecule pdbxMolecule) {
        // TODO not considered in ref
    }

    @Override
    public void consumePdbxMoleculeFeatures(PdbxMoleculeFeatures pdbxMoleculeFeatures) {
        // TODO not considered in ref
    }

    @Override
    public void consumePdbxNonpolyScheme(PdbxNonpolyScheme pdbxNonpolyScheme) {
        // TODO not impled in ref
    }

    @Override
    public void consumePdbxReferenceEntityLink(PdbxReferenceEntityLink pdbxReferenceEntityLink) {
        // TODO not considered in ref
    }

    @Override
    public void consumePdbxReferenceEntityList(PdbxReferenceEntityList pdbxReferenceEntityList) {
        // TODO not considered in ref
    }

    @Override
    public void consumePdbxReferenceEntityPolyLink(PdbxReferenceEntityPolyLink pdbxReferenceEntityPolyLink) {
        // TODO not considered in ref
    }

    @Override
    public void consumePdbxStructAssembly(PdbxStructAssembly pdbxStructAssembly) {
        this.structAssembly = pdbxStructAssembly;
    }

    @Override
    public void consumePdbxStructAssemblyGen(PdbxStructAssemblyGen pdbxStructAssemblyGen) {
        this.structAssemblyGen = pdbxStructAssemblyGen;
    }

    @Override
    public void consumePdbxStructModResidue(PdbxStructModResidue pdbxStructModResidue) {
        // TODO not considered in ref
    }

    @Override
    public void consumePdbxStructOperList(PdbxStructOperList pdbxStructOperList) {
        this.structOpers = pdbxStructOperList;
    }

    @Override
    public void consumeRefine(Refine refine) {
        for (int rowIndex = 0; rowIndex < refine.getRowCount(); rowIndex++) {
            // RESOLUTION
            // in very rare cases (for instance hybrid methods x-ray + neutron diffraction, e.g. 3ins, 4n9m)
            // there are 2 resolution values, one for each method
            // we take the last one found so that behaviour is like in PDB file parsing
            double lsDResHigh = refine.getLsDResHigh().get(rowIndex);
            // TODO this could use a check to keep reasonable values - 1.5 may be overwritten by 0.0
            if (pdbHeader.getResolution() != PDBHeader.DEFAULT_RESOLUTION) {
                logger.warn("More than 1 resolution value present, will use last one {} and discard previous {}",
                        lsDResHigh, String.format("%4.2f",pdbHeader.getResolution()));
            }
            pdbHeader.setResolution((float) lsDResHigh);

            FloatColumn lsRFactorRFree = refine.getLsRFactorRFree();
            // RFREE
            if (pdbHeader.getRfree() != PDBHeader.DEFAULT_RFREE) {
                logger.warn("More than 1 Rfree value present, will use last one {} and discard previous {}",
                        lsRFactorRFree, String.format("%4.2f",pdbHeader.getRfree()));
            }
            if (lsRFactorRFree.isDefined()) {
                pdbHeader.setRfree((float) lsRFactorRFree.get(rowIndex));
            } else {
                // some entries like 2ifo haven't got this field at all
                logger.info("_refine.ls_R_factor_R_free not present, not parsing Rfree value");
            }

            // RWORK
            FloatColumn lsRFactorRWork = refine.getLsRFactorRWork();
            if(pdbHeader.getRwork() != PDBHeader.DEFAULT_RFREE) {
                logger.warn("More than 1 R work value present, will use last one {} and discard previous {} ",
                        lsRFactorRWork, String.format("%4.2f",pdbHeader.getRwork()));
            }
            if (lsRFactorRWork.isDefined()) {
                pdbHeader.setRwork((float) lsRFactorRWork.get(rowIndex));
            } else {
                logger.info("_refine.ls_R_factor_R_work not present, not parsing R-work value");
            }
        }
    }

    @Override
    public void consumeStruct(Struct struct) {
        if (struct.isDefined() && struct.getTitle().isDefined()) {
            pdbHeader.setTitle(struct.getTitle().get(0));
        }

        if (struct.isDefined() && struct.getEntryId().isDefined()) {
            pdbHeader.setIdCode(struct.getEntryId().get(0));
            structure.setPDBCode(struct.getEntryId().get(0));
        }
    }

    @Override
    public void consumeStructAsym(StructAsym structAsym) {
        this.structAsym = structAsym;
    }

    @Override
    public void consumeStructConf(StructConf structConf) {
        // TODO not considered in ref
    }

    @Override
    public void consumeStructConn(StructConn structConn) {
        this.structConn = structConn;
    }

    @Override
    public void consumeStructConnType(StructConnType structConnType) {
        // TODO not considered in ref
    }

    @Override
    public void consumeStructKeywords(StructKeywords structKeywords) {
        StrColumn pdbxKeywords = structKeywords.getPdbxKeywords();
        // TODO what is the correct format for these?
        pdbHeader.setDescription(pdbxKeywords.values().collect(Collectors.joining(", ")));
        pdbHeader.setClassification(pdbxKeywords.values().collect(Collectors.joining(", ")));
    }

    @Override
    public void consumeStructNcsOper(StructNcsOper structNcsOper) {
        this.structNcsOper = structNcsOper;
    }

    @Override
    public void consumeStructRef(StructRef structRef) {
        this.structRef = structRef;
    }

    @Override
    public void consumeStructRefSeq(StructRefSeq structRefSeq) {
        for (int rowIndex = 0; rowIndex < structRefSeq.getRowCount(); rowIndex++) {
            String refId = structRefSeq.getRefId().get(rowIndex);

            DBRef dbRef = new DBRef();

            dbRef.setIdCode(structRefSeq.getPdbxPDBIdCode().get(rowIndex));
            dbRef.setDbAccession(structRefSeq.getPdbxDbAccession().get(rowIndex));
            dbRef.setDbIdCode(structRefSeq.getPdbxDbAccession().get(rowIndex));
            dbRef.setChainName(structRefSeq.getPdbxStrandId().get(rowIndex));

            OptionalInt structRefRowIndex = IntStream.range(0, structRef.getRowCount())
                    .filter(i -> structRef.getId().get(i).equals(refId))
                    .findFirst();

            if (structRefRowIndex.isPresent()) {
                dbRef.setDatabase(structRef.getDbName().get(structRefRowIndex.getAsInt()));
                dbRef.setDbIdCode(structRef.getDbCode().get(structRefRowIndex.getAsInt()));
            } else {
                logger.info("could not find StructRef `{} for StructRefSeq {}", refId, rowIndex);
            }

            int seqBegin;
            int seqEnd;

            try {
                seqBegin = Integer.parseInt(structRefSeq.getPdbxAuthSeqAlignBeg().get(rowIndex));
                seqEnd = Integer.parseInt(structRefSeq.getPdbxAuthSeqAlignEnd().get(rowIndex));
            } catch (NumberFormatException e) {
                // this happens in a few entries, annotation error? e.g. 6eoj
                logger.warn("Couldn't parse pdbx_auth_seq_align_beg/end in _struct_ref_seq. Will not store dbref " +
                        "alignment info for accession {}. Error: {}", dbRef.getDbAccession(), e.getMessage());
                return;
            }

            char beginInsCode = ' ';
            String pdbxSeqAlignBegInsCode = structRefSeq.getPdbxSeqAlignBegInsCode().get(rowIndex);
            if (pdbxSeqAlignBegInsCode.length() > 0) {
                beginInsCode = pdbxSeqAlignBegInsCode.charAt(0);
            }

            char endInsCode = ' ';
            String pdbxSeqAlignEndInsCode = structRefSeq.getPdbxSeqAlignEndInsCode().get(rowIndex);
            if (pdbxSeqAlignEndInsCode.length() > 0) {
                endInsCode = pdbxSeqAlignEndInsCode.charAt(0);
            }

            if (beginInsCode == '?') {
                beginInsCode = ' ';
            }
            if (endInsCode == '?') {
                endInsCode = ' ';
            }

            dbRef.setSeqBegin(seqBegin);
            dbRef.setInsertBegin(beginInsCode);
            dbRef.setSeqEnd(seqEnd);
            dbRef.setInsertEnd(endInsCode);

            int dbSeqBegin = structRefSeq.getDbAlignBeg().get(rowIndex);
            int dbSeqEnd = structRefSeq.getDbAlignEnd().get(rowIndex);

            char dbBeginInsCode = ' ';
            StrColumn pdbxDbAlignBegInsCodeCol = structRefSeq.getPdbxDbAlignBegInsCode();
            if (pdbxDbAlignBegInsCodeCol.isDefined()) {
                String pdbxDbAlignBegInsCode = pdbxDbAlignBegInsCodeCol.get(rowIndex);
                if (pdbxDbAlignBegInsCode.length() > 0) {
                    dbBeginInsCode = pdbxDbAlignBegInsCode.charAt(0);
                }
            }

            char dbEndInsCode = ' ';
            StrColumn pdbxDbAlignEndInsCodeCol = structRefSeq.getPdbxDbAlignEndInsCode();
            if (pdbxDbAlignEndInsCodeCol.isDefined()) {
                String pdbxDbAlignEndInsCode = pdbxDbAlignEndInsCodeCol.get(rowIndex);
                if (pdbxDbAlignEndInsCode.length() > 0) {
                    dbEndInsCode = pdbxDbAlignEndInsCode.charAt(0);
                }
            }

            if (dbBeginInsCode == '?') {
                dbBeginInsCode = ' ';
            }
            if (dbEndInsCode == '?') {
                dbEndInsCode = ' ';
            }

            dbRef.setDbSeqBegin(dbSeqBegin);
            dbRef.setIdbnsBegin(dbBeginInsCode);
            dbRef.setDbSeqEnd(dbSeqEnd);
            dbRef.setIdbnsEnd(dbEndInsCode);

            List<DBRef> dbrefs = structure.getDBRefs();
            if (dbrefs == null) {
                dbrefs = new ArrayList<>();
            }
            dbrefs.add(dbRef);

            logger.debug(dbRef.toPDB());

            structure.setDBRefs(dbrefs);
        }
    }

    @Override
    public void consumeStructRefSeqDif(StructRefSeqDif structRefSeqDif) {
        this.structRefSeqDif = structRefSeqDif;
    }

    @Override
    public void consumeStructSheetRange(StructSheetRange structSheetRange) {
        // TODO not considered in ref
    }

    @Override
    public void consumeStructSite(StructSite structSite) {
        if (params.isHeaderOnly()) {
            return;
        }

        List<Site> sites = structure.getSites();
        if (sites == null) {
            sites = new ArrayList<>();
        }

        for (int rowIndex = 0; rowIndex < structSite.getRowCount(); rowIndex++) {
            Site site = null;
            for (Site asite : sites) {
                if (asite.getSiteID().equals(structSite.getId().get(rowIndex))) {
                    site = asite; // prevent duplicate siteIds
                }
            }

            boolean addSite = false;
            if (site == null) {
                site = new Site();
                addSite = true;
            }

            site.setSiteID(structSite.getId().get(rowIndex));
            site.setDescription(structSite.getDetails().get(rowIndex));
            site.setEvCode(structSite.getPdbxEvidenceCode().get(rowIndex));

            if (addSite) {
                sites.add(site);
            }
        }

        structure.setSites(sites);
    }

    @Override
    public void consumeStructSiteGen(StructSiteGen structSiteGen) {
        this.structSiteGen = structSiteGen;
    }

    @Override
    public void consumeSymmetry(Symmetry symmetry) {
        for (int rowIndex = 0; rowIndex < symmetry.getRowCount(); rowIndex++) {
            String spaceGroupString = symmetry.getSpaceGroupNameH_M().get(rowIndex);
            SpaceGroup spaceGroup = SymoplibParser.getSpaceGroup(spaceGroupString);
            if (spaceGroup == null) {
                logger.warn("Space group '{}' not recognised as a standard space group", spaceGroupString);
                structure.getPDBHeader()
                        .getCrystallographicInfo()
                        .setNonStandardSg(true);
            } else {
                structure.getPDBHeader()
                        .getCrystallographicInfo()
                        .setSpaceGroup(spaceGroup);
                structure.getPDBHeader()
                        .getCrystallographicInfo()
                        .setNonStandardSg(false);
            }
        }
    }

    @Override
    public void finish() {
        if (currentChain != null) {
            currentChain.addGroup(currentGroup);

            Optional<Chain> testChain = currentModel.stream()
                    .filter(chain -> chain.getId().equals(currentChain.getId()))
                    .findFirst();

            if (!testChain.isPresent()) {
                currentModel.add(currentChain);
            }
        } else if (!params.isHeaderOnly()) {
            logger.warn("current chain is null at end of document.");
        }

        allModels.add(currentModel);

        initMaps();

        for (int rowIndex = 0; rowIndex < structAsym.getRowCount(); rowIndex++) {
            String id = structAsym.getId().get(rowIndex);
            String entityId = structAsym.getEntityId().get(rowIndex);
            logger.debug("Entity {} matches asym_id: {}", entityId, id);

            Chain chain = getEntityChain(entityId);
            Chain seqRes = (Chain) chain.clone();
            // to solve issue #160 (e.g. 3u7t)
            seqRes = removeSeqResHeterogeneity(seqRes);
            seqRes.setId(id);
            seqRes.setName(asymId2authorId.getOrDefault(id, id));

            EntityType type = EntityType.entityTypeFromString(getEntityType(entityId));
            if (type == null || type == EntityType.POLYMER) {
                seqResChains.add(seqRes);
            }

            logger.debug(" seqres: {} {}<", id, seqRes);
            addEntity(rowIndex, entityId, getEntityDescription(entityId), getEntityType(entityId));
        }

        if (!structAsym.isDefined() || structAsym.getRowCount() == 0) {
            logger.warn("No _struct_asym category in file, no SEQRES groups will be added.");
        }

        // entities
        // In addEntities above we created the entities if they were present in the file
        // Now we need to make sure that they are linked to chains and also that if they are not present in the file we
        // need to add them now
        linkEntities();

        // now that we know the entities, we can add all chains to structure so that they are stored
        // properly as polymer/nonpolymer/water chains inside structure
        allModels.forEach(structure::addModel);

        // Only align if requested (default) and not when headerOnly mode with no Atoms.
        // Otherwise, we store the empty SeqRes Groups unchanged in the right chains.
        if (params.isAlignSeqRes() && !params.isHeaderOnly()){
            logger.debug("Parsing mode align_seqres, will parse SEQRES and align to ATOM sequence");
            alignSeqRes();
        } else {
            logger.debug("Parsing mode unalign_seqres, will parse SEQRES but not align it to ATOM sequence");
            SeqRes2AtomAligner.storeUnAlignedSeqRes(structure, seqResChains, params.isHeaderOnly());
        }

        // Now make sure all altlocgroups have all the atoms in all the groups
        StructureTools.cleanUpAltLocs(structure);

        // NOTE bonds and charges can only be done at this point that the chain id mapping is properly sorted out
        if (!params.isHeaderOnly()) {
            if (params.shouldCreateAtomBonds()) {
                addBonds();
            }

            if (params.shouldCreateAtomCharges()) {
                addCharges();
            }
        }

        if (!params.isHeaderOnly()) {
            addSites();
        }

        // set the oligomeric state info in the header...
        if (params.isParseBioAssembly()) {
            // the more detailed mapping of chains to rotation operations happens in StructureIO...

            Map<Integer, BioAssemblyInfo> bioAssemblies = new LinkedHashMap<>();
            List<org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssembly> structAssemblies = convert(structAssembly);
            List<org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssemblyGen> structAssemblyGens = convert(structAssemblyGen);

            for (org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssembly pdbxStructAssembly : structAssemblies) {
                List<org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssemblyGen> pdbxStructAssemblyGens = structAssemblyGens.stream()
                        .filter(sag -> sag.getAssembly_id().equals(pdbxStructAssembly.getId()))
                        .collect(Collectors.toList());

                BiologicalAssemblyBuilder builder = new BiologicalAssemblyBuilder();

                // these are the transformations that need to be applied to our model
                List<BiologicalAssemblyTransformation> transformations = builder.getBioUnitTransformationList(pdbxStructAssembly,
                        pdbxStructAssemblyGens, convert(structOpers));

                int bioAssemblyId = -1;
                try {
                    bioAssemblyId = Integer.parseInt(pdbxStructAssembly.getId());
                } catch (NumberFormatException e) {
                    logger.info("Could not parse a numerical bio assembly id from '{}'", pdbxStructAssembly.getId());
                }

                // if bioassembly id is not numerical we throw it away
                // this happens usually for viral capsid entries, like 1ei7
                // see issue #230 in github
                if (bioAssemblyId != -1) {
                    int mmSize = 0;
                    // note that the transforms contain asym ids of both polymers and non-polymers
                    // For the mmsize, we are only interested in the polymers
                    for (BiologicalAssemblyTransformation transf : transformations) {
                        Chain c = structure.getChain(transf.getChainId());
                        if (c == null) {
                            logger.info("Could not find asym id {} specified in struct_assembly_gen", transf.getChainId());
                            continue;
                        }
                        if (c.getEntityType() == EntityType.POLYMER &&
                                // for entries like 4kro, sugars are annotated as polymers but we
                                // don't want them in the macromolecularSize count
                                !c.getEntityInfo().getDescription().contains("SUGAR")) {
                            mmSize++;
                        }
                    }

                    BioAssemblyInfo bioAssembly = new BioAssemblyInfo();
                    bioAssembly.setId(bioAssemblyId);
                    bioAssembly.setMacromolecularSize(mmSize);
                    bioAssembly.setTransforms(transformations);
                    bioAssemblies.put(bioAssemblyId, bioAssembly);
                }

            }
            structure.getPDBHeader()
                    .setBioAssemblies(bioAssemblies);
        }

        setStructNcsOps();
        setCrystallographicInfoMetadata();

        Map<String, List<SeqMisMatch>> misMatchMap = new HashMap<>();
        for (int rowIndex = 0; rowIndex < structRefSeqDif.getRowCount(); rowIndex++) {
            SeqMisMatch seqMisMatch = new SeqMisMatchImpl();
            seqMisMatch.setDetails(structRefSeqDif.getDetails().get(rowIndex));

            String insCode = structRefSeqDif.getPdbxPdbInsCode().get(rowIndex);
                if (insCode != null && insCode.equals("?")) {
                insCode = null;
            }
            seqMisMatch.setInsCode(insCode);
            seqMisMatch.setOrigGroup(structRefSeqDif.getDbMonId().get(rowIndex));
            seqMisMatch.setPdbGroup(structRefSeqDif.getMonId().get(rowIndex));
            seqMisMatch.setPdbResNum(structRefSeqDif.getPdbxAuthSeqNum().get(rowIndex));
            seqMisMatch.setUniProtId(structRefSeqDif.getPdbxSeqDbAccessionCode().get(rowIndex));
            seqMisMatch.setSeqNum(structRefSeqDif.getSeqNum().get(rowIndex));

            String strandId = structRefSeqDif.getPdbxPdbStrandId().get(rowIndex);
            List<SeqMisMatch> seqMisMatches = misMatchMap.computeIfAbsent(strandId, k -> new ArrayList<>());
            seqMisMatches.add(seqMisMatch);
        }

        for (String chainId : misMatchMap.keySet()){
            Chain chain = structure.getPolyChainByPDB(chainId);
            if (chain == null) {
                logger.warn("Could not set mismatches for chain with author id {}", chainId);
                continue;
            }

            chain.setSeqMisMatches(misMatchMap.get(chainId));
        }
    }

    private String getEntityType(String entityId) {
        return IntStream.range(0, entity.getRowCount())
                .filter(i -> entity.getId().get(i).equals(entityId))
                .mapToObj(i -> entity.getType().get(i))
                .findFirst()
                .orElseThrow(() -> new NoSuchElementException("could not find entity with id " + entityId));
    }

    private String getEntityDescription(String entityId) {
        return IntStream.range(0, entity.getRowCount())
                .filter(i -> entity.getId().get(i).equals(entityId))
                .mapToObj(i -> entity.getPdbxDescription().get(i))
                .findFirst()
                .orElseThrow(() -> new NoSuchElementException("could not find entity with id " + entityId));
    }

    private void addEntity(int asymRowIndex, String entityId, String pdbxDescription, String type) {
        int eId = 0;
        try {
            eId = Integer.parseInt(entityId);
        } catch (NumberFormatException e) {
            logger.warn("Could not parse mol_id from string {}. Will use 0 for creating Entity", entityId);
        }
        
        int entityRowIndex = IntStream.range(0, entity.getRowCount())
                .filter(i -> entity.getId().get(i).equals(entityId))
                .findFirst()
                .orElse(-1);
        
        EntityInfo entityInfo = structure.getEntityById(eId);
        
        if (entityInfo == null) {
            entityInfo = new EntityInfo();
            entityInfo.setMolId(eId);
            // we only add the compound if a polymeric one (to match what the PDB parser does)
            if (entityRowIndex != -1) {
                entityInfo.setDescription(pdbxDescription);

                EntityType eType = EntityType.entityTypeFromString(type);
                if (eType != null) {
                    entityInfo.setType(eType);
                } else {
                    logger.warn("Type '{}' is not recognised as a valid entity type for entity {}", type, eId);
                }
                addAncilliaryEntityData(asymRowIndex, entityInfo);
                structure.addEntityInfo(entityInfo);
                logger.debug("Adding Entity with entity id {} from _entity, with name: {}", eId, 
                        entityInfo.getDescription());
            }
        }
    }

    private void addAncilliaryEntityData(int asymRowIndex, EntityInfo entityInfo) {
        // Loop through each of the entity types and add the corresponding data
        // We're assuming if data is duplicated between sources it is consistent
        // This is a potentially huge assumption...

        for (int rowIndex = 0; rowIndex < entitySrcGen.getRowCount(); rowIndex++) {
            if (entitySrcGen.getEntityId().get(rowIndex).equals(structAsym.getEntityId().get(asymRowIndex))) {
                continue;
            }

            addInformationFromEntitySrcGen(rowIndex, entityInfo);
        }

        for (int rowIndex = 0; rowIndex < entitySrcNat.getRowCount(); rowIndex++) {
            if (entitySrcNat.getEntityId().get(rowIndex).equals(structAsym.getEntityId().get(asymRowIndex))) {
                continue;
            }

            addInformationFromEntitySrcNat(rowIndex, entityInfo);
        }

        for (int rowIndex = 0; rowIndex < entitySrcSyn.getRowCount(); rowIndex++) {
            if (entitySrcSyn.getEntityId().get(rowIndex).equals(structAsym.getEntityId().get(asymRowIndex))) {
                continue;
            }

            addInformationFromEntitySrcSyn(rowIndex, entityInfo);
        }
    }

    private void addInformationFromEntitySrcSyn(int rowIndex, EntityInfo entityInfo) {
        entityInfo.setOrganismCommon(entitySrcSyn.getOrganismCommonName().get(rowIndex));
        entityInfo.setOrganismScientific(entitySrcSyn.getOrganismScientific().get(rowIndex));
        entityInfo.setOrganismTaxId(entitySrcSyn.getNcbiTaxonomyId().get(rowIndex));
    }

    private void addInformationFromEntitySrcNat(int rowIndex, EntityInfo entityInfo) {
        entityInfo.setAtcc(entitySrcNat.getPdbxAtcc().get(rowIndex));
        entityInfo.setCell(entitySrcNat.getPdbxCell().get(rowIndex));
        entityInfo.setOrganismCommon(entitySrcNat.getCommonName().get(rowIndex));
        entityInfo.setOrganismScientific(entitySrcNat.getPdbxOrganismScientific().get(rowIndex));
        entityInfo.setOrganismTaxId(entitySrcNat.getPdbxNcbiTaxonomyId().get(rowIndex));
    }

    private void addInformationFromEntitySrcGen(int rowIndex, EntityInfo entityInfo) {
        entityInfo.setAtcc(entitySrcGen.getPdbxGeneSrcAtcc().get(rowIndex));
        entityInfo.setCell(entitySrcGen.getPdbxGeneSrcCell().get(rowIndex));
        entityInfo.setOrganismCommon(entitySrcGen.getGeneSrcCommonName().get(rowIndex));
        entityInfo.setOrganismScientific(entitySrcGen.getPdbxGeneSrcScientificName().get(rowIndex));
        entityInfo.setOrganismTaxId(entitySrcGen.getPdbxGeneSrcNcbiTaxonomyId().get(rowIndex));
        entityInfo.setExpressionSystemTaxId(entitySrcGen.getPdbxHostOrgNcbiTaxonomyId().get(rowIndex));
        entityInfo.setExpressionSystem(entitySrcGen.getPdbxHostOrgScientificName().get(rowIndex));
    }

    private void setStructNcsOps() {
        List<Matrix4d> ncsOperators = new ArrayList<>();

        for (int rowIndex = 0; rowIndex < structNcsOper.getRowCount(); rowIndex++) {
            if (!"generate".equals(structNcsOper.getCode().get(rowIndex))) {
                continue;
            }

            try {
                Matrix4d operator = new Matrix4d();

                operator.setElement(0, 0, structNcsOper.getMatrix11().get(rowIndex));
                operator.setElement(0, 1, structNcsOper.getMatrix12().get(rowIndex));
                operator.setElement(0, 2, structNcsOper.getMatrix13().get(rowIndex));

                operator.setElement(1, 0, structNcsOper.getMatrix21().get(rowIndex));
                operator.setElement(1, 1, structNcsOper.getMatrix22().get(rowIndex));
                operator.setElement(1, 2, structNcsOper.getMatrix23().get(rowIndex));

                operator.setElement(2, 0, structNcsOper.getMatrix31().get(rowIndex));
                operator.setElement(2, 1, structNcsOper.getMatrix32().get(rowIndex));
                operator.setElement(2, 2, structNcsOper.getMatrix33().get(rowIndex));

                operator.setElement(3, 0, 0);
                operator.setElement(3, 1, 0);
                operator.setElement(3, 2, 0);
                operator.setElement(3, 3, 1);

                ncsOperators.add(operator);
            } catch (NumberFormatException e) {
                logger.warn("Error parsing doubles in NCS operator list, skipping operator {}", rowIndex + 1);
            }
        }

        if (ncsOperators.size() > 0) {
            structure.getCrystallographicInfo()
                    .setNcsOperators(ncsOperators.toArray(new Matrix4d[0]));
        }
    }

    private List<org.biojava.nbio.structure.io.mmcif.model.PdbxStructOperList> convert(PdbxStructOperList structOpers) {
        List<org.biojava.nbio.structure.io.mmcif.model.PdbxStructOperList> re = new ArrayList<>();
        for (int rowIndex = 0; rowIndex < structOpers.getRowCount(); rowIndex++) {
            org.biojava.nbio.structure.io.mmcif.model.PdbxStructOperList pdbxStructOperList =
                    new org.biojava.nbio.structure.io.mmcif.model.PdbxStructOperList();

            pdbxStructOperList.setId(structOpers.getId().get(rowIndex));
            pdbxStructOperList.setName(structOpers.getName().get(rowIndex));
            pdbxStructOperList.setSymmetry_operation(structOpers.getSymmetryOperation().get(rowIndex));
            pdbxStructOperList.setType(structOpers.getType().get(rowIndex));

            pdbxStructOperList.setMatrix11(String.valueOf(structOpers.getMatrix11().get(rowIndex)));
            pdbxStructOperList.setMatrix12(String.valueOf(structOpers.getMatrix12().get(rowIndex)));
            pdbxStructOperList.setMatrix13(String.valueOf(structOpers.getMatrix13().get(rowIndex)));
            pdbxStructOperList.setMatrix21(String.valueOf(structOpers.getMatrix21().get(rowIndex)));
            pdbxStructOperList.setMatrix22(String.valueOf(structOpers.getMatrix22().get(rowIndex)));
            pdbxStructOperList.setMatrix23(String.valueOf(structOpers.getMatrix23().get(rowIndex)));
            pdbxStructOperList.setMatrix31(String.valueOf(structOpers.getMatrix31().get(rowIndex)));
            pdbxStructOperList.setMatrix32(String.valueOf(structOpers.getMatrix32().get(rowIndex)));
            pdbxStructOperList.setMatrix33(String.valueOf(structOpers.getMatrix33().get(rowIndex)));

            pdbxStructOperList.setVector1(String.valueOf(structOpers.getVector1().get(rowIndex)));
            pdbxStructOperList.setVector2(String.valueOf(structOpers.getVector2().get(rowIndex)));
            pdbxStructOperList.setVector3(String.valueOf(structOpers.getVector3().get(rowIndex)));

            re.add(pdbxStructOperList);
        }
        return re;
    }

    private List<org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssemblyGen> convert(PdbxStructAssemblyGen structAssemblyGen) {
        List<org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssemblyGen> re = new ArrayList<>();
        for (int rowIndex = 0; rowIndex < structAssemblyGen.getRowCount(); rowIndex++) {
            org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssemblyGen pdbxStructAssemblyGen =
                    new org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssemblyGen();

            pdbxStructAssemblyGen.setAssembly_id(structAssemblyGen.getAssemblyId().get(rowIndex));
            pdbxStructAssemblyGen.setAsym_id_list(structAssemblyGen.getAsymIdList().get(rowIndex));
            pdbxStructAssemblyGen.setOper_expression(structAssemblyGen.getOperExpression().get(rowIndex));

            re.add(pdbxStructAssemblyGen);
        }
        return re;
    }

    private List<org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssembly> convert(PdbxStructAssembly structAssembly) {
        List<org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssembly> re = new ArrayList<>();
        for (int rowIndex = 0; rowIndex < structAssembly.getRowCount(); rowIndex++) {
            org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssembly pdbxStructAssembly =
                    new org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssembly();

            pdbxStructAssembly.setDetails(structAssembly.getDetails().get(rowIndex));
            pdbxStructAssembly.setId(structAssembly.getId().get(rowIndex));
            pdbxStructAssembly.setMethod_details(structAssembly.getMethodDetails().get(rowIndex));
            pdbxStructAssembly.setOligomeric_count(structAssembly.getOligomericCount().getStringData(rowIndex));
            pdbxStructAssembly.setOligomeric_details(structAssembly.getOligomericDetails().get(rowIndex));

            re.add(pdbxStructAssembly);
        }
        return re;
    }

    private void setCrystallographicInfoMetadata() {
        if (parsedScaleMatrix != null) {
            PDBCrystallographicInfo crystalInfo = structure.getCrystallographicInfo();
            boolean nonStd = false;
            if (crystalInfo.getCrystalCell() != null && !crystalInfo.getCrystalCell().checkScaleMatrix(parsedScaleMatrix)) {
                nonStd = true;
            }

            crystalInfo.setNonStandardCoordFrameConvention(nonStd);
        }
    }

    private void addSites() {
        List<Site> sites = structure.getSites();
        if (sites == null) sites = new ArrayList<>();

        for (int rowIndex = 0; rowIndex < structSiteGen.getRowCount(); rowIndex++) {
            // For each StructSiteGen, find the residues involved, if they exist then
            String site_id = structSiteGen.getSiteId().get(rowIndex); // multiple could be in same site.
            if (site_id == null) {
                site_id = "";
            }
            String comp_id = structSiteGen.getLabelCompId().get(rowIndex);  // PDBName

            // Assumption: the author chain ID and residue number for the site is consistent with the original
            // author chain id and residue numbers.

            String asymId = structSiteGen.getLabelAsymId().get(rowIndex); // chain name
            String authId = structSiteGen.getAuthAsymId().get(rowIndex); // chain Id
            String auth_seq_id = structSiteGen.getAuthSeqId().get(rowIndex); // Res num

            String insCode = structSiteGen.getPdbxAuthInsCode().get(rowIndex);
            if (insCode != null && insCode.equals("?")) {
                insCode = null;
            }

            // Look for asymID = chainID and seqID = seq_ID.  Check that comp_id matches the resname.
            Group g = null;
            try {
                Chain chain = structure.getChain(asymId);

                if (null != chain) {
                    try {
                        Character insChar = null;
                        if (null != insCode && insCode.length() > 0) {
                            insChar = insCode.charAt(0);
                        }
                        g = chain.getGroupByPDB(new ResidueNumber(null, Integer.parseInt(auth_seq_id), insChar));
                    } catch (NumberFormatException e) {
                        logger.warn("Could not lookup residue : {}{}", authId, auth_seq_id);
                    }
                }
            } catch (StructureException e) {
                logger.warn("Problem finding residue in site entry {} - {}",
                        structSiteGen.getSiteId().get(rowIndex), e.getMessage());
            }

            if (g != null) {
                // 2. find the site_id, if not existing, create anew.
                Site site = null;
                for (Site asite : sites) {
                    if (site_id.equals(asite.getSiteID())) {
                        site = asite;
                    }
                }

                boolean addSite = false;

                // 3. add this residue to the site.
                if (site == null) {
                    addSite = true;
                    site = new Site();
                    site.setSiteID(site_id);
                }

                List<Group> groups = site.getGroups();
                if (groups == null) {
                    groups = new ArrayList<>();
                }

                // Check the self-consistency of the residue reference from auth_seq_id and chain_id
                if (!comp_id.equals(g.getPDBName())) {
                    logger.warn("comp_id doesn't match the residue at {} {} - skipping", authId, auth_seq_id);
                } else {
                    groups.add(g);
                    site.setGroups(groups);
                }
                if (addSite) {
                    sites.add(site);
                }
            }
        }
        structure.setSites(sites);
    }

    private void addCharges() {
        ChargeAdder.addCharges(structure);
    }

    /**
     * The method will return a new reference to a Chain with any consecutive groups
     * having same residue numbers removed.
     * This is necessary to solve the microheterogeneity issue in entries like 3u7t (see github issue #160)
     */
    private static Chain removeSeqResHeterogeneity(Chain c) {
        Chain trimmedChain = new ChainImpl();
        ResidueNumber lastResNum = null;

        for (Group g : c.getAtomGroups()) {
            // note we have to deep copy this, otherwise they stay linked and would get altered in addGroup(g)
            ResidueNumber currentResNum = new ResidueNumber(
                    g.getResidueNumber().getChainName(),
                    g.getResidueNumber().getSeqNum(),
                    g.getResidueNumber().getInsCode());

            if (lastResNum == null || !lastResNum.equals(currentResNum)) {
                trimmedChain.addGroup(g);
            } else {
                logger.debug("Removing seqres group because it seems to be repeated in entity_poly_seq, most likely " +
                        "has hetero='y': {}", g);
            }
            lastResNum = currentResNum;

        }
        return trimmedChain;
    }

    private void addBonds() {
        BondMaker maker = new BondMaker(structure, params);
        maker.makeBonds();
        maker.formBondsFromStructConn(convert(structConn));
    }

    private List<org.biojava.nbio.structure.io.mmcif.model.StructConn> convert(StructConn structConn) {
        return IntStream.range(0, structConn.getRowCount())
                .mapToObj(rowIndex -> {
                    org.biojava.nbio.structure.io.mmcif.model.StructConn sc =
                            new org.biojava.nbio.structure.io.mmcif.model.StructConn();

                    sc.setPdbx_ptnr1_PDB_ins_code(structConn.getPdbxPtnr1PDBInsCode().get(rowIndex));
                    sc.setPdbx_ptnr2_PDB_ins_code(structConn.getPdbxPtnr2PDBInsCode().get(rowIndex));
                    sc.setPtnr1_auth_seq_id(structConn.getPtnr1AuthSeqId().getStringData(rowIndex));
                    sc.setPtnr2_auth_seq_id(structConn.getPtnr2AuthSeqId().getStringData(rowIndex));
                    sc.setPtnr1_label_comp_id(structConn.getPtnr1LabelCompId().get(rowIndex));
                    sc.setPtnr2_label_comp_id(structConn.getPtnr2LabelCompId().get(rowIndex));
                    sc.setPtnr1_label_atom_id(structConn.getPtnr1LabelAtomId().get(rowIndex));
                    sc.setPtnr2_label_atom_id(structConn.getPtnr2LabelAtomId().get(rowIndex));
                    sc.setPdbx_ptnr1_label_alt_id(structConn.getPdbxPtnr1LabelAltId().get(rowIndex));
                    sc.setPdbx_ptnr2_label_alt_id(structConn.getPdbxPtnr2LabelAltId().get(rowIndex));
                    sc.setPtnr1_symmetry(structConn.getPtnr1Symmetry().get(rowIndex));
                    sc.setPtnr2_symmetry(structConn.getPtnr2Symmetry().get(rowIndex));
                    sc.setConn_type_id(structConn.getConnTypeId().get(rowIndex));

                    return sc;
                })
                .collect(Collectors.toList());
    }

    private void alignSeqRes() {
        logger.debug("Parsing mode align_seqres, will align to ATOM to SEQRES sequence");

        // fix SEQRES residue numbering for all models

        for (int model = 0; model < structure.nrModels(); model++) {
            List<Chain> atomList   = structure.getModel(model);

            for (Chain seqResChain : seqResChains){

                // this extracts the matching atom chain from atomList
                Chain atomChain = SeqRes2AtomAligner.getMatchingAtomRes(seqResChain, atomList, true);

                if (atomChain == null) {
                    // most likely there's no observed residues at all for the seqres chain: can't map
                    // e.g. 3zyb: chains with asym_id L,M,N,O,P have no observed residues
                    logger.info("Could not map SEQRES chain with asym_id={} to any ATOM chain. Most likely there's " +
                            "no observed residues in the chain.", seqResChain.getId());
                    continue;
                }

                //map the atoms to the seqres...

                // we need to first clone the seqres so that they stay independent for different models
                List<Group> seqResGroups = new ArrayList<>();
                for (int i = 0; i < seqResChain.getAtomGroups().size(); i++) {
                    seqResGroups.add((Group)seqResChain.getAtomGroups().get(i).clone());
                }

                for (int seqResPos = 0 ; seqResPos < seqResGroups.size(); seqResPos++) {
                    Group seqresG = seqResGroups.get(seqResPos);
                    boolean found = false;
                    for (Group atomG : atomChain.getAtomGroups()) {

                        int internalNr = getInternalNr(atomG);

                        if (seqresG.getResidueNumber().getSeqNum() == internalNr) {
                            seqResGroups.set(seqResPos, atomG);
                            found = true;
                            break;
                        }
                    }

                    if (!found)
                        // so far the residue number has tracked internal numbering.
                        // however there are no atom records, as such this can't be a PDB residue number...
                        seqresG.setResidueNumber(null);
                }
                atomChain.setSeqResGroups(seqResGroups);
            }
        }
    }

    private int getInternalNr(Group atomG) {
        if (atomG.getType().equals(GroupType.AMINOACID)) {
            AminoAcidImpl aa = (AminoAcidImpl) atomG;
            return (int) aa.getId();
        } else if (atomG.getType().equals(GroupType.NUCLEOTIDE)) {
            NucleotideImpl nu = (NucleotideImpl) atomG;
            return (int) nu.getId();
        } else {
            HetatomImpl he = (HetatomImpl) atomG;
            return (int) he.getId();
        }
    }

    private void linkEntities() {
        for (List<Chain> allModel : allModels) {
            for (Chain chain : allModel) {
                //logger.info("linking entities for " + chain.getId() + " "  + chain.getName());
                String entityId = asymId2entityId.get(chain.getId());

                if (entityId == null) {
                    // this can happen for instance if the cif file didn't have _struct_asym category at all
                    // and thus we have no asymId2entityId mapping at all
                    logger.info("No entity id could be found for chain {}", chain.getId());
                    continue;
                }

                int eId = Integer.parseInt(entityId);

                // Entities are not added for non-polymeric entities, if a chain is non-polymeric its entity won't be found.
                // TODO: add all entities and unique compounds and add methods to directly get polymer or non-polymer
                // asyms (chains).  Either create a unique StructureImpl or modify existing for a better representation of the
                // mmCIF internal data structures but is compatible with Structure interface.
                // Some examples of PDB entries with this kind of problem:
                //   - 2uub: asym_id X, chainName Z, entity_id 24: fully non-polymeric but still with its own chainName
                //   - 3o6j: asym_id K, chainName Z, entity_id 6 : a single water molecule
                //   - 1dz9: asym_id K, chainName K, entity_id 6 : a potassium ion alone

                EntityInfo entityInfo = structure.getEntityById(eId);
                if (entityInfo == null) {
                    // Supports the case where the only chain members were from non-polymeric entity that is missing.
                    // Solved by creating a new Compound(entity) to which this chain will belong.
                    logger.info("Could not find an Entity for entity_id {}, for chain id {}, creating a new Entity.",
                            eId, chain.getId());
                    entityInfo = new EntityInfo();
                    entityInfo.setMolId(eId);
                    entityInfo.addChain(chain);
                    if (chain.isWaterOnly()) {
                        entityInfo.setType(EntityType.WATER);
                    } else {
                        entityInfo.setType(EntityType.NONPOLYMER);
                    }
                    chain.setEntityInfo(entityInfo);
                    structure.addEntityInfo(entityInfo);
                } else {
                    logger.debug("Adding chain with chain id {} (auth id {}) to Entity with entity_id {}",
                            chain.getId(), chain.getName(), eId);
                    entityInfo.addChain(chain);
                    chain.setEntityInfo(entityInfo);
                }

            }

        }

        // if no entity information was present in file we then go and find the entities heuristically with EntityFinder
        List<EntityInfo> entityInfos = structure.getEntityInfos();
        if (entityInfos == null || entityInfos.isEmpty()) {
            List<List<Chain>> polyModels = new ArrayList<>();
            List<List<Chain>> nonPolyModels = new ArrayList<>();
            List<List<Chain>> waterModels = new ArrayList<>();

            for (List<Chain> model : allModels) {
                List<Chain> polyChains = new ArrayList<>();
                List<Chain> nonPolyChains = new ArrayList<>();
                List<Chain> waterChains = new ArrayList<>();

                polyModels.add(polyChains);
                nonPolyModels.add(nonPolyChains);
                waterModels.add(waterChains);

                for (Chain chain : model) {
                    // we only have entities for polymeric chains, all others are ignored for assigning entities
                    if (chain.isWaterOnly()) {
                        waterChains.add(chain);
                    } else if (chain.isPureNonPolymer()) {
                        nonPolyChains.add(chain);
                    } else {
                        polyChains.add(chain);
                    }
                }
            }

            entityInfos = EntityFinder.findPolyEntities(polyModels);
            EntityFinder.createPurelyNonPolyEntities(nonPolyModels, waterModels, entityInfos);

            structure.setEntityInfos(entityInfos);
        }

        // final sanity check: it can happen that from the annotated entities some are not linked to any chains
        // e.g. 3s26: a sugar entity does not have any chains associated to it (it seems to be happening with many sugar compounds)
        // we simply log it, this can sign some other problems if the entities are used down the line
        for (EntityInfo e : entityInfos) {
            if (e.getChains().isEmpty()) {
                logger.info("Entity {} '{}' has no chains associated to it",
                        e.getMolId() < 0 ? "with no entity id" : e.getMolId(), e.getDescription());
            }
        }
    }

    private void initMaps() {
        if (structAsym == null || !structAsym.isDefined() || structAsym.getRowCount() == 0) {
            logger.info("No _struct_asym category found in file. No asym id to entity_id mapping will be available");
            return;
        }

        Map<String, List<String>> entityId2asymId = new HashMap<>();
        for (int rowIndex = 0; rowIndex < structAsym.getRowCount(); rowIndex++) {
            String id = structAsym.getId().get(rowIndex);
            String entityId = structAsym.getEntityId().get(rowIndex);

            logger.debug("Entity {} matches asym_id: {}", entityId, id);

            asymId2entityId.put(id, entityId);

            if (entityId2asymId.containsKey(entityId)) {
                List<String> asymIds = entityId2asymId.get(entityId);
                asymIds.add(id);
            } else {
                List<String> asymIds = new ArrayList<>();
                asymIds.add(id);
                entityId2asymId.put(entityId, asymIds);
            }
        }

        if (entityPoly == null || !entityPoly.isDefined() || entityPoly.getRowCount() == 0) {
            logger.info("No _entity_poly category found in file. No asym id to author id mapping will be available " +
                    "for header only parsing");
            return;
        }

        for (int rowIndex = 0; rowIndex < entityPoly.getRowCount(); rowIndex++) {
            if (!entityPoly.getPdbxStrandId().isDefined()) {
                logger.info("_entity_poly.pdbx_strand_id is null for entity {}. Won't be able to map asym ids to " +
                        "author ids for this entity.", entityPoly.getEntityId().get(rowIndex));
                break;
            }

            String[] chainNames = entityPoly.getPdbxStrandId().get(rowIndex).split(",");
            List<String> asymIds = entityId2asymId.get(entityPoly.getEntityId().get(rowIndex));
            if (chainNames.length != asymIds.size()) {
                logger.warn("The list of asym ids (from _struct_asym) and the list of author ids (from _entity_poly) " +
                        "for entity {} have different lengths! Can't provide a mapping from asym ids to author chain " +
                        "ids", entityPoly.getEntityId().get(rowIndex));
                break;
            }

            for (int i = 0; i < chainNames.length; i++) {
                asymId2authorId.put(asymIds.get(i), chainNames[i]);
            }
        }
    }

    @Override
    public Structure getContainer() {
        return structure;
    }
}
