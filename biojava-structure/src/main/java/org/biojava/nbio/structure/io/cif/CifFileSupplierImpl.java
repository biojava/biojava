package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.xtal.CrystalCell;
import org.biojava.nbio.structure.xtal.SpaceGroup;
import org.rcsb.cif.model.Block;
import org.rcsb.cif.model.Category;
import org.rcsb.cif.model.CifFile;
import org.rcsb.cif.model.Column;
import org.rcsb.cif.model.atomsite.AtomSite;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;
import java.util.stream.Collector;

/**
 * Convert a BioJava {@link Structure} to a CifFile.
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @since 5.2.1
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
        AtomSite atomSite = wrappedAtoms.stream().collect(toAtomSite());

        Block.BlockBuilder blockBuilder = CifFile.enterFile()
                .enterBlock(structure.getPDBCode());

        if (atomSite.isDefined() && atomSite.getRowCount() > 0) {
            // set atom site
            blockBuilder.addCategory(atomSite);
        }

        if (crystalCell != null) {
            // set cell category
            blockBuilder.enterCategory("cell")
                    .enterColumn("length_a")
                    .floatValues(crystalCell.getA())
                    .leaveColumn()

                    .enterColumn("length_b")
                    .floatValues(crystalCell.getB())
                    .leaveColumn()

                    .enterColumn("length_c")
                    .floatValues(crystalCell.getC())
                    .leaveColumn()

                    .enterColumn("angle_alpha")
                    .floatValues(crystalCell.getAlpha())
                    .leaveColumn()

                    .enterColumn("angle_beta")
                    .floatValues(crystalCell.getBeta())
                    .leaveColumn()

                    .enterColumn("angle_gamma")
                    .floatValues(crystalCell.getGamma())
                    .leaveColumn()
                    .leaveCategory();
        }

        if (spaceGroup != null) {
            // set symmetry category
            blockBuilder.enterCategory("symmetry")
                    .enterColumn("space_group_name_H-M")
                    .stringValues(spaceGroup.getShortSymbol())
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

    private static Collector<WrappedAtom, ?, AtomSite> toAtomSite() {
        return Collector.of(AtomSiteCollector::new,
                AtomSiteCollector::accept,
                AtomSiteCollector::combine,
                AtomSiteCollector::get);
    }

    static class AtomSiteCollector implements Consumer<WrappedAtom> {
        private final Column.StrColumnBuilder groupPDB;
        private final Column.IntColumnBuilder id;
        private final Column.StrColumnBuilder typeSymbol;
        private final Column.StrColumnBuilder labelAtomId;
        private final Column.StrColumnBuilder labelAltId;
        private final Column.StrColumnBuilder labelCompId;
        private final Column.StrColumnBuilder labelAsymId;
        private final Column.StrColumnBuilder labelEntityId;
        private final Column.IntColumnBuilder labelSeqId;
        private final Column.StrColumnBuilder pdbxPDBInsCode;
        private final Column.FloatColumnBuilder cartnX;
        private final Column.FloatColumnBuilder cartnY;
        private final Column.FloatColumnBuilder cartnZ;
        private final Column.FloatColumnBuilder occupancy;
        private final Column.FloatColumnBuilder bIsoOrEquiv;
        private final Column.IntColumnBuilder authSeqId;
        private final Column.StrColumnBuilder authCompId;
        private final Column.StrColumnBuilder authAsymId;
        private final Column.StrColumnBuilder authAtomId;
        private final Column.IntColumnBuilder pdbxPDBModelNum;
        private int atomId = 1;

        AtomSiteCollector() {
            this.groupPDB = Column.enterStrColumn("group_PDB");
            this.id = Column.enterIntColumn("id");
            this.typeSymbol = Column.enterStrColumn("type_symbol");
            this.labelAtomId = Column.enterStrColumn("label_atom_id");
            this.labelAltId = Column.enterStrColumn("label_alt_id");
            this.labelCompId = Column.enterStrColumn("label_comp_id");
            this.labelAsymId = Column.enterStrColumn("label_asym_id");
            this.labelEntityId = Column.enterStrColumn("label_entity_id");
            this.labelSeqId = Column.enterIntColumn("label_seq_id");
            this.pdbxPDBInsCode = Column.enterStrColumn("pdbx_PDB_ins_code");
            this.cartnX = Column.enterFloatColumn("Cartn_x");
            this.cartnY = Column.enterFloatColumn("Cartn_y");
            this.cartnZ = Column.enterFloatColumn("Cartn_z");
            this.occupancy = Column.enterFloatColumn("occupancy");
            this.bIsoOrEquiv = Column.enterFloatColumn("B_iso_or_equiv");
            this.authSeqId = Column.enterIntColumn("auth_seq_id");
            this.authCompId = Column.enterStrColumn("auth_comp_id");
            this.authAsymId = Column.enterStrColumn("auth_asym_id");
            this.authAtomId = Column.enterStrColumn("auth_atom_id");
            this.pdbxPDBModelNum = Column.enterIntColumn("pdbx_PDB_model_num");
        }

        @Override
        public void accept(WrappedAtom wrappedAtom) {
            Atom atom = wrappedAtom.getAtom();
            Group group = atom.getGroup();
            Chain chain = group.getChain();

            groupPDB.stringValues(group.getType().equals(GroupType.HETATM) ? "HETATM" : "ATOM");
            id.intValues(wrappedAtom.getAtomId());
            Element element = atom.getElement();
            typeSymbol.stringValues(element.equals(Element.R) ? "X" : element.toString().toUpperCase());
            labelAtomId.stringValues(atom.getName());
            Character altLoc = atom.getAltLoc();
            if (altLoc == null || altLoc == ' ') {
                labelAltId.markNextNotPresent();
            } else {
                labelAltId.stringValues(String.valueOf(altLoc));
            }
            labelCompId.stringValues(group.getPDBName());
            labelAsymId.stringValues(wrappedAtom.getChainId());
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
            labelEntityId.stringValues(entityId);
            labelSeqId.intValues(seqId);
            String insCode = "";
            if (group.getResidueNumber().getInsCode() != null ) {
                insCode = Character.toString(group.getResidueNumber().getInsCode());
            }
            if (insCode.isEmpty()) {
                pdbxPDBInsCode.markNextUnknown();
            } else {
                pdbxPDBInsCode.stringValues(insCode);
            }
            cartnX.floatValues(atom.getX());
            cartnY.floatValues(atom.getY());
            cartnZ.floatValues(atom.getZ());
            occupancy.floatValues(atom.getOccupancy());
            bIsoOrEquiv.floatValues(atom.getTempFactor());
            authSeqId.intValues(group.getResidueNumber().getSeqNum());
            authCompId.stringValues(group.getPDBName());
            authAsymId.stringValues(wrappedAtom.getChainName());
            authAtomId.stringValues(atom.getName());
            pdbxPDBModelNum.intValues(wrappedAtom.getModel());

            atomId++;
        }

        AtomSiteCollector combine(AtomSiteCollector other) {
            throw new UnsupportedOperationException("impl by calling addAll for all collection - not feeling like writing that code");
        }
        
        @SuppressWarnings("Duplicates")
        AtomSite get() {
            Map<String, Column> columns = new LinkedHashMap<>();
            put(columns, groupPDB.build());
            put(columns, groupPDB.build());
            put(columns, id.build());
            put(columns, typeSymbol.build());
            put(columns, labelAtomId.build());
            put(columns, labelAltId.build());
            put(columns, labelCompId.build());
            put(columns, labelAsymId.build());
            put(columns, labelEntityId.build());
            put(columns, labelSeqId.build());
            put(columns, pdbxPDBInsCode.build());
            put(columns, cartnX.build());
            put(columns, cartnY.build());
            put(columns, cartnZ.build());
            put(columns, occupancy.build());
            put(columns, bIsoOrEquiv.build());
            put(columns, authSeqId.build());
            put(columns, authCompId.build());
            put(columns, authAsymId.build());
            put(columns, authAtomId.build());
            put(columns, pdbxPDBModelNum.build());

            return (AtomSite) Category.enterCategory("atom_site", columns, null).build();
        }

        private void put(Map<String, Column> columns, Column column) {
            columns.put(column.getColumnName(), column);
        }
    }
}
