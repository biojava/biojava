package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.xtal.CrystalCell;
import org.biojava.nbio.structure.xtal.SpaceGroup;
import org.rcsb.cif.model.Block;
import org.rcsb.cif.model.Category;
import org.rcsb.cif.model.CifFile;
import org.rcsb.cif.model.Column;
import org.rcsb.cif.model.atomsite.AtomSite;

import java.util.Collection;
import java.util.LinkedHashMap;
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
        AtomSite atomSite = structure.getChains()
                .stream()
                .map(Chain::getAtomGroups)
                .flatMap(Collection::stream)
                .map(Group::getAtoms)
                .flatMap(Collection::stream)
                .collect(toAtomSite());

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

    private static Collector<Atom, ?, AtomSite> toAtomSite() {
        return Collector.of(AtomSiteCollector::new,
                AtomSiteCollector::accept,
                AtomSiteCollector::combine,
                AtomSiteCollector::get);
    }

    static class AtomSiteCollector implements Consumer<Atom> {
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
        public void accept(Atom atom) {
            Group group = atom.getGroup();
            Chain chain = group.getChain();

            groupPDB.stringValues(group.getType().equals(GroupType.HETATM) ? "HETATM" : "ATOM");
            id.intValues(atomId);
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
            labelAsymId.stringValues(chain.getId());
            String entityId = "0";
            int seqId = group.getResidueNumber().getSeqNum();
            if (chain.getEntityInfo() != null) {
                entityId = Integer.toString(chain.getEntityInfo().getMolId());
                if (chain.getEntityInfo().getType() == EntityType.POLYMER) {
                    // this only makes sense for polymeric chains, non-polymer chains will never have seqres groups and there's no point in calling getAlignedResIndex
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
            authAsymId.stringValues(chain.getId());
            authAtomId.stringValues(atom.getName());
            pdbxPDBModelNum.intValues(1);

            // TODO respect symmetry
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
