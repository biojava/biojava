package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.EntityType;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.xtal.CrystalCell;
import org.biojava.nbio.structure.xtal.SpaceGroup;
import org.rcsb.cif.CifBuilder;
import org.rcsb.cif.model.Category;
import org.rcsb.cif.model.CifFile;
import org.rcsb.cif.model.FloatColumnBuilder;
import org.rcsb.cif.model.IntColumnBuilder;
import org.rcsb.cif.model.StrColumnBuilder;
import org.rcsb.cif.schema.StandardSchemata;
import org.rcsb.cif.schema.mm.MmCifBlockBuilder;
import org.rcsb.cif.schema.mm.MmCifCategoryBuilder;
import org.rcsb.cif.schema.mm.MmCifFileBuilder;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;
import java.util.stream.Collector;

/**
 * Convert a BioJava {@link Structure} to a CifFile.
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @since 5.3.0
 */
class CifFileSupplierImpl implements CifFileSupplier<Structure> {
    @Override
    public CifFile get(Structure structure) {
        // for now BioJava only considered 3 categories for create a Cif representation of a structure

        // cell
        CrystalCell crystalCell = structure.getPDBHeader().getCrystallographicInfo().getCrystalCell();
        // symmetry
        SpaceGroup spaceGroup = structure.getPDBHeader().getCrystallographicInfo().getSpaceGroup();
        // atom_site
        List<WrappedAtom> wrappedAtoms = collectWrappedAtoms(structure);
        Category atomSite = wrappedAtoms.stream().collect(toAtomSite());

        MmCifBlockBuilder blockBuilder = CifBuilder.enterFile(StandardSchemata.MMCIF)
                .enterBlock(structure.getPDBCode());

        if (atomSite.isDefined() && atomSite.getRowCount() > 0) {
            // set atom site
            blockBuilder.addCategory(atomSite);
        }

        if (crystalCell != null) {
            // set cell category
            blockBuilder.enterCell()
                    .enterLengthA()
                    .add(crystalCell.getA())
                    .leaveColumn()

                    .enterLengthB()
                    .add(crystalCell.getB())
                    .leaveColumn()

                    .enterLengthC()
                    .add(crystalCell.getC())
                    .leaveColumn()

                    .enterAngleAlpha()
                    .add(crystalCell.getAlpha())
                    .leaveColumn()

                    .enterAngleBeta()
                    .add(crystalCell.getBeta())
                    .leaveColumn()

                    .enterAngleGamma()
                    .add(crystalCell.getGamma())
                    .leaveColumn()
                    .leaveCategory();
        }

        if (spaceGroup != null) {
            // set symmetry category
            blockBuilder.enterSymmetry()
                    .enterSpaceGroupNameH_M()
                    .add(spaceGroup.getShortSymbol())
                    .leaveColumn()
                    .leaveCategory();
        }

        return blockBuilder.leaveBlock().leaveFile();
    }

    private static List<WrappedAtom> collectWrappedAtoms(Structure structure) {
        List<WrappedAtom> wrappedAtoms = new ArrayList<>();

        for (int modelIndex = 0; modelIndex < structure.nrModels(); modelIndex++) {
            final int model = modelIndex + 1;
            for (Chain chain : structure.getChains(modelIndex)) {
                final String chainName = chain.getName();
                final String chainId = chain.getId();
                for (Group group : chain.getAtomGroups()) {
                    // The alt locs can have duplicates, since at parsing time we make sure that all alt loc groups have
                    // all atoms (see StructureTools#cleanUpAltLocs)
                    // Thus we have to remove duplicates here by using the atom id
                    // See issue https://github.com/biojava/biojava/issues/778 and
                    // TestAltLocs.testMmcifWritingAllAltlocs/testMmcifWritingPartialAltlocs
                    Map<Integer, WrappedAtom> uniqueAtoms = new LinkedHashMap<>();
                    for (int atomIndex = 0; atomIndex < group.size(); atomIndex++) {
                        Atom atom = group.getAtom(atomIndex);
                        if (atom == null) {
                            continue;
                        }

                        uniqueAtoms.put(atom.getPDBserial(), new WrappedAtom(chain, model, chainName, chainId, atom, atom.getPDBserial()));
                    }

                    if (group.hasAltLoc()) {
                        for (Group alt : group.getAltLocs()) {
                            for (int atomIndex = 0; atomIndex < alt.size(); atomIndex++) {
                                Atom atom = alt.getAtom(atomIndex);
                                if (atom == null) {
                                    continue;
                                }

                                uniqueAtoms.put(atom.getPDBserial(), new WrappedAtom(chain, model, chainName, chainId, atom, atom.getPDBserial()));
                            }
                        }
                    }

                    wrappedAtoms.addAll(uniqueAtoms.values());
                }
            }
        }

        return wrappedAtoms;
    }

    static class WrappedAtom {
        private final Chain chain;
        private final int model;
        private final String chainName;
        private final String chainId;
        private final Atom atom;
        private final int atomId;

        WrappedAtom(Chain chain, int model, String chainName, String chainId, Atom atom, int atomId) {
            this.chain = chain;
            this.model = model;
            this.chainName = chainName;
            this.chainId = chainId;
            this.atom = atom;
            this.atomId = atomId;
        }

        Chain getChain() {
            return chain;
        }

        int getModel() {
            return model;
        }

        String getChainName() {
            return chainName;
        }

        String getChainId() {
            return chainId;
        }

        Atom getAtom() {
            return atom;
        }

        int getAtomId() {
            return atomId;
        }
    }

    private static Collector<WrappedAtom, ?, Category> toAtomSite() {
        return Collector.of(AtomSiteCollector::new,
                AtomSiteCollector::accept,
                AtomSiteCollector::combine,
                AtomSiteCollector::get);
    }

    static class AtomSiteCollector implements Consumer<WrappedAtom> {
        private final MmCifCategoryBuilder.AtomSiteBuilder atomSiteBuilder;
        private final StrColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> groupPDB;
        private final IntColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> id;
        private final StrColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> typeSymbol;
        private final StrColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> labelAtomId;
        private final StrColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> labelAltId;
        private final StrColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> labelCompId;
        private final StrColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> labelAsymId;
        private final StrColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> labelEntityId;
        private final IntColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> labelSeqId;
        private final StrColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> pdbxPDBInsCode;
        private final FloatColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> cartnX;
        private final FloatColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> cartnY;
        private final FloatColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> cartnZ;
        private final FloatColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> occupancy;
        private final FloatColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> bIsoOrEquiv;
        private final IntColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> authSeqId;
        private final StrColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> authCompId;
        private final StrColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> authAsymId;
        private final StrColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> authAtomId;
        private final IntColumnBuilder<MmCifCategoryBuilder.AtomSiteBuilder, MmCifBlockBuilder, MmCifFileBuilder> pdbxPDBModelNum;

        AtomSiteCollector() {
            this.atomSiteBuilder = new MmCifCategoryBuilder.AtomSiteBuilder(null);
            this.groupPDB = atomSiteBuilder.enterGroupPDB();
            this.id = atomSiteBuilder.enterId();
            this.typeSymbol = atomSiteBuilder.enterTypeSymbol();
            this.labelAtomId = atomSiteBuilder.enterLabelAtomId();
            this.labelAltId = atomSiteBuilder.enterLabelAltId();
            this.labelCompId = atomSiteBuilder.enterLabelCompId();
            this.labelAsymId = atomSiteBuilder.enterLabelAsymId();
            this.labelEntityId = atomSiteBuilder.enterLabelEntityId();
            this.labelSeqId = atomSiteBuilder.enterLabelSeqId();
            this.pdbxPDBInsCode = atomSiteBuilder.enterPdbxPDBInsCode();
            this.cartnX = atomSiteBuilder.enterCartnX();
            this.cartnY = atomSiteBuilder.enterCartnY();
            this.cartnZ = atomSiteBuilder.enterCartnZ();
            this.occupancy = atomSiteBuilder.enterOccupancy();
            this.bIsoOrEquiv = atomSiteBuilder.enterBIsoOrEquiv();
            this.authSeqId = atomSiteBuilder.enterAuthSeqId();
            this.authCompId = atomSiteBuilder.enterAuthCompId();
            this.authAsymId = atomSiteBuilder.enterAuthAsymId();
            this.authAtomId = atomSiteBuilder.enterAuthAtomId();
            this.pdbxPDBModelNum = atomSiteBuilder.enterPdbxPDBModelNum();
        }

        @Override
        public void accept(WrappedAtom wrappedAtom) {
            Atom atom = wrappedAtom.getAtom();
            Group group = atom.getGroup();
            Chain chain = group.getChain();

            groupPDB.add(group.getType().equals(GroupType.HETATM) ? "HETATM" : "ATOM");
            id.add(wrappedAtom.getAtomId());
            Element element = atom.getElement();
            typeSymbol.add(element.equals(Element.R) ? "X" : element.toString().toUpperCase());
            labelAtomId.add(atom.getName());
            Character altLoc = atom.getAltLoc();
            if (altLoc == null || altLoc == ' ') {
                labelAltId.markNextNotPresent();
            } else {
                labelAltId.add(String.valueOf(altLoc));
            }
            labelCompId.add(group.getPDBName());
            labelAsymId.add(wrappedAtom.getChainId());
            String entityId = "0";
            int seqId = group.getResidueNumber().getSeqNum();
            if (chain.getEntityInfo() != null) {
                entityId = Integer.toString(chain.getEntityInfo().getMolId());
                if (chain.getEntityInfo().getType() == EntityType.POLYMER) {
                    // this only makes sense for polymeric chains, non-polymer chains will never have seqres groups and
                    // there's no point in calling getAlignedResIndex
                    seqId = chain.getEntityInfo().getAlignedResIndex(group, chain);
                }
            }
            labelEntityId.add(entityId);
            labelSeqId.add(seqId);
            String insCode = "";
            if (group.getResidueNumber().getInsCode() != null ) {
                insCode = Character.toString(group.getResidueNumber().getInsCode());
            }
            if (insCode.isEmpty()) {
                pdbxPDBInsCode.markNextUnknown();
            } else {
                pdbxPDBInsCode.add(insCode);
            }
            cartnX.add(atom.getX());
            cartnY.add(atom.getY());
            cartnZ.add(atom.getZ());
            occupancy.add(atom.getOccupancy());
            bIsoOrEquiv.add(atom.getTempFactor());
            authSeqId.add(group.getResidueNumber().getSeqNum());
            authCompId.add(group.getPDBName());
            authAsymId.add(wrappedAtom.getChainName());
            authAtomId.add(atom.getName());
            pdbxPDBModelNum.add(wrappedAtom.getModel());
        }

        AtomSiteCollector combine(AtomSiteCollector other) {
            throw new UnsupportedOperationException("impl by calling addAll for all collection");
        }

        Category get() {
            groupPDB.leaveColumn();
            id.leaveColumn();
            typeSymbol.leaveColumn();
            labelAtomId.leaveColumn();
            labelAltId.leaveColumn();
            labelCompId.leaveColumn();
            labelAsymId.leaveColumn();
            labelEntityId.leaveColumn();
            labelSeqId.leaveColumn();
            pdbxPDBInsCode.leaveColumn();
            cartnX.leaveColumn();
            cartnY.leaveColumn();
            cartnZ.leaveColumn();
            occupancy.leaveColumn();
            bIsoOrEquiv.leaveColumn();
            authSeqId.leaveColumn();
            authCompId.leaveColumn();
            authAsymId.leaveColumn();
            authAtomId.leaveColumn();
            pdbxPDBModelNum.leaveColumn();
            return atomSiteBuilder.build();
        }
    }
}
