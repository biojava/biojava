package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.chem.ChemCompAtom;
import org.biojava.nbio.structure.chem.ChemCompDescriptor;
import org.rcsb.cif.schema.mm.AtomSite;
import org.rcsb.cif.schema.mm.AtomSites;
import org.rcsb.cif.schema.mm.AuditAuthor;
import org.rcsb.cif.schema.mm.Cell;
import org.rcsb.cif.schema.mm.ChemComp;
import org.rcsb.cif.schema.mm.ChemCompBond;
import org.rcsb.cif.schema.mm.DatabasePDBRemark;
import org.rcsb.cif.schema.mm.DatabasePDBRev;
import org.rcsb.cif.schema.mm.DatabasePDBRevRecord;
import org.rcsb.cif.schema.mm.Entity;
import org.rcsb.cif.schema.mm.EntityPoly;
import org.rcsb.cif.schema.mm.EntityPolySeq;
import org.rcsb.cif.schema.mm.EntitySrcGen;
import org.rcsb.cif.schema.mm.EntitySrcNat;
import org.rcsb.cif.schema.mm.Exptl;
import org.rcsb.cif.schema.mm.PdbxAuditRevisionHistory;
import org.rcsb.cif.schema.mm.PdbxChemCompIdentifier;
import org.rcsb.cif.schema.mm.PdbxDatabaseStatus;
import org.rcsb.cif.schema.mm.PdbxEntityBranchDescriptor;
import org.rcsb.cif.schema.mm.PdbxEntitySrcSyn;
import org.rcsb.cif.schema.mm.PdbxMolecule;
import org.rcsb.cif.schema.mm.PdbxMoleculeFeatures;
import org.rcsb.cif.schema.mm.PdbxNonpolyScheme;
import org.rcsb.cif.schema.mm.PdbxReferenceEntityLink;
import org.rcsb.cif.schema.mm.PdbxReferenceEntityList;
import org.rcsb.cif.schema.mm.PdbxReferenceEntityPolyLink;
import org.rcsb.cif.schema.mm.PdbxStructAssembly;
import org.rcsb.cif.schema.mm.PdbxStructAssemblyGen;
import org.rcsb.cif.schema.mm.PdbxStructModResidue;
import org.rcsb.cif.schema.mm.PdbxStructOperList;
import org.rcsb.cif.schema.mm.Refine;
import org.rcsb.cif.schema.mm.Struct;
import org.rcsb.cif.schema.mm.StructAsym;
import org.rcsb.cif.schema.mm.StructConf;
import org.rcsb.cif.schema.mm.StructConn;
import org.rcsb.cif.schema.mm.StructConnType;
import org.rcsb.cif.schema.mm.StructKeywords;
import org.rcsb.cif.schema.mm.StructNcsOper;
import org.rcsb.cif.schema.mm.StructRef;
import org.rcsb.cif.schema.mm.StructRefSeq;
import org.rcsb.cif.schema.mm.StructRefSeqDif;
import org.rcsb.cif.schema.mm.StructSheetRange;
import org.rcsb.cif.schema.mm.StructSite;
import org.rcsb.cif.schema.mm.StructSiteGen;
import org.rcsb.cif.schema.mm.Symmetry;

/**
 * Defines a rather generic interface which allows to populate some data structure with data parsed from a CIF file.
 *
 * @param <S> the type of container an implementing class will return
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @since 5.3.0
 */
public interface CifFileConsumer<S> {
    /**
     * Setup routine which initializes a new container.
     */
    void prepare();

    /**
     * Ultimate setup which can include steps which require several categories to be available and integrate them into
     * the final container.
     */
    void finish();

    /**
     * Retrieve the created container representing a CIF file.
     * @return all desired information wrapped as object of type <code>S</code>
     */
    S getContainer();
}
