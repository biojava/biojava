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
import org.rcsb.cif.model.Category;
import org.rcsb.cif.model.CifFile;
import org.rcsb.cif.model.builder.BlockBuilder;
import org.rcsb.cif.model.builder.CategoryBuilder;
import org.rcsb.cif.model.builder.CifBuilder;
import org.rcsb.cif.model.builder.FloatColumnBuilder;
import org.rcsb.cif.model.builder.IntColumnBuilder;
import org.rcsb.cif.model.builder.StrColumnBuilder;

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

        BlockBuilder blockBuilder = new CifBuilder()
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

        public WrappedAtom(Chain chain, int model, String chainName, String chainId, Atom atom, int atomId) {
            this.chain = chain;
            this.model = model;
            this.chainName = chainName;
            this.chainId = chainId;
            this.atom = atom;
            this.atomId = atomId;
        }

        public Chain getChain() {
            return chain;
        }

        public int getModel() {
            return model;
        }

        public String getChainName() {
            return chainName;
        }

        public String getChainId() {
            return chainId;
        }

        public Atom getAtom() {
            return atom;
        }

        public int getAtomId() {
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
        private final StrColumnBuilder<CategoryBuilder> groupPDB;
        private final IntColumnBuilder<CategoryBuilder> id;
        private final StrColumnBuilder<CategoryBuilder> typeSymbol;
        private final StrColumnBuilder<CategoryBuilder> labelAtomId;
        private final StrColumnBuilder<CategoryBuilder> labelAltId;
        private final StrColumnBuilder<CategoryBuilder> labelCompId;
        private final StrColumnBuilder<CategoryBuilder> labelAsymId;
        private final StrColumnBuilder<CategoryBuilder> labelEntityId;
        private final IntColumnBuilder<CategoryBuilder> labelSeqId;
        private final StrColumnBuilder<CategoryBuilder> pdbxPDBInsCode;
        private final FloatColumnBuilder<CategoryBuilder> cartnX;
        private final FloatColumnBuilder<CategoryBuilder> cartnY;
        private final FloatColumnBuilder<CategoryBuilder> cartnZ;
        private final FloatColumnBuilder<CategoryBuilder> occupancy;
        private final FloatColumnBuilder<CategoryBuilder> bIsoOrEquiv;
        private final IntColumnBuilder<CategoryBuilder> authSeqId;
        private final StrColumnBuilder<CategoryBuilder> authCompId;
        private final StrColumnBuilder<CategoryBuilder> authAsymId;
        private final StrColumnBuilder<CategoryBuilder> authAtomId;
        private final IntColumnBuilder<CategoryBuilder> pdbxPDBModelNum;

        AtomSiteCollector() {
            // TODO this doesn't really make the case for the builder ;)
            this.groupPDB = new StrColumnBuilder<>("atom_site", "group_PDB", null);
            this.id = new IntColumnBuilder<>("atom_site", "id", null);
            this.typeSymbol = new StrColumnBuilder<>("atom_site", "type_symbol", null);
            this.labelAtomId = new StrColumnBuilder<>("atom_site", "label_atom_id", null);
            this.labelAltId = new StrColumnBuilder<>("atom_site", "label_alt_id", null);
            this.labelCompId = new StrColumnBuilder<>("atom_site", "label_comp_id", null);
            this.labelAsymId = new StrColumnBuilder<>("atom_site", "label_asym_id", null);
            this.labelEntityId = new StrColumnBuilder<>("atom_site", "label_entity_id", null);
            this.labelSeqId = new IntColumnBuilder<>("atom_site", "label_seq_id", null);
            this.pdbxPDBInsCode = new StrColumnBuilder<>("atom_site", "pdbx_PDB_ins_code", null);
            this.cartnX = new FloatColumnBuilder<>("atom_site", "Cartn_x", null);
            this.cartnY = new FloatColumnBuilder<>("atom_site", "Cartn_y", null);
            this.cartnZ = new FloatColumnBuilder<>("atom_site", "Cartn_z", null);
            this.occupancy = new FloatColumnBuilder<>("atom_site", "occupancy", null);
            this.bIsoOrEquiv = new FloatColumnBuilder<>("atom_site", "B_iso_or_equiv", null);
            this.authSeqId = new IntColumnBuilder<>("atom_site", "auth_seq_id", null);
            this.authCompId = new StrColumnBuilder<>("atom_site", "auth_comp_id", null);
            this.authAsymId = new StrColumnBuilder<>("atom_site", "auth_asym_id", null);
            this.authAtomId = new StrColumnBuilder<>("atom_site", "auth_atom_id", null);
            this.pdbxPDBModelNum = new IntColumnBuilder<>("atom_site", "pdbx_PDB_model_num", null);
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
            throw new UnsupportedOperationException("impl by calling addAll for all collection - not feeling like writing that code");
        }

        Category get() {
            return new CategoryBuilder("atom_site", null)
                    .addColumn(groupPDB.build())
                    .addColumn(id.build())
                    .addColumn(typeSymbol.build())
                    .addColumn(labelAtomId.build())
                    .addColumn(labelAltId.build())
                    .addColumn(labelCompId.build())
                    .addColumn(labelAsymId.build())
                    .addColumn(labelEntityId.build())
                    .addColumn(labelSeqId.build())
                    .addColumn(pdbxPDBInsCode.build())
                    .addColumn(cartnX.build())
                    .addColumn(cartnY.build())
                    .addColumn(cartnZ.build())
                    .addColumn(occupancy.build())
                    .addColumn(bIsoOrEquiv.build())
                    .addColumn(authSeqId.build())
                    .addColumn(authCompId.build())
                    .addColumn(authAsymId.build())
                    .addColumn(authAtomId.build())
                    .addColumn(pdbxPDBModelNum.build())
                    .build();
        }
    }
}
