package org.biojava.nbio.structure.io.cif;

import org.rcsb.cif.model.Category;
import org.rcsb.cif.model.atomsite.AtomSite;
import org.rcsb.cif.model.atomsites.AtomSites;
import org.rcsb.cif.model.cell.Cell;
import org.rcsb.cif.model.chemcomp.ChemComp;
import org.rcsb.cif.model.chemcompbond.ChemCompBond;
import org.rcsb.cif.model.entity.Entity;
import org.rcsb.cif.model.entitypoly.EntityPoly;
import org.rcsb.cif.model.entitypolyseq.EntityPolySeq;
import org.rcsb.cif.model.exptl.Exptl;
import org.rcsb.cif.model.pdbxchemcompidentifier.PdbxChemCompIdentifier;
import org.rcsb.cif.model.pdbxentitydescriptor.PdbxEntityDescriptor;
import org.rcsb.cif.model.pdbxmolecule.PdbxMolecule;
import org.rcsb.cif.model.pdbxmoleculefeatures.PdbxMoleculeFeatures;
import org.rcsb.cif.model.pdbxnonpolyscheme.PdbxNonpolyScheme;
import org.rcsb.cif.model.pdbxreferenceentitylink.PdbxReferenceEntityLink;
import org.rcsb.cif.model.pdbxreferenceentitylist.PdbxReferenceEntityList;
import org.rcsb.cif.model.pdbxreferenceentitypolylink.PdbxReferenceEntityPolyLink;
import org.rcsb.cif.model.pdbxstructassembly.PdbxStructAssembly;
import org.rcsb.cif.model.pdbxstructassemblygen.PdbxStructAssemblyGen;
import org.rcsb.cif.model.pdbxstructmodresidue.PdbxStructModResidue;
import org.rcsb.cif.model.pdbxstructoperlist.PdbxStructOperList;
import org.rcsb.cif.model.struct.Struct;
import org.rcsb.cif.model.structasym.StructAsym;
import org.rcsb.cif.model.structconf.StructConf;
import org.rcsb.cif.model.structconn.StructConn;
import org.rcsb.cif.model.structconntype.StructConnType;
import org.rcsb.cif.model.structkeywords.StructKeywords;
import org.rcsb.cif.model.structncsoper.StructNcsOper;
import org.rcsb.cif.model.structsheetrange.StructSheetRange;
import org.rcsb.cif.model.structsite.StructSite;
import org.rcsb.cif.model.structsitegen.StructSiteGen;
import org.rcsb.cif.model.symmetry.Symmetry;

public interface CifFileConsumer<S> {
    void prepare();

    void consumeAtomSite(AtomSite atomSite);

    void consumeAtomSites(AtomSites atomSites);

    void consumeAuditAuthor(Category auditAuthor);

    void consumeCell(Cell cell);

    void consumeChemComp(ChemComp chemComp);

    void consumeChemCompBond(ChemCompBond chemCompBond);

    void consumeDatabasePDBremark(Category databasePDBremark);

    void consumeDatabasePDBrev(Category databasePDBrev);

    void consumeDatabasePDBrevRecord(Category databasePDBrevRecord);

    void consumeEntity(Entity e);

    void consumeEntityPoly(EntityPoly entityPoly);

    void consumeEntitySrcGen(Category entitySrcGen);

    void consumeEntitySrcNat(Category entitySrcNat);

    void consumeEntitySrcSyn(Category entitySrcSyn);

    void consumeEntityPolySeq(EntityPolySeq entityPolySeq);

    void consumeExptl(Exptl exptl);

    void consumePdbxAuditRevisionHistory(Category pdbxAuditRevisionHistory);

    void consumePdbxChemCompIdentifier(PdbxChemCompIdentifier pdbxChemCompIdentifier);

    void consumePdbxDatabaseStatus(Category pdbxDatabaseStatus);

    void consumePdbxEntityDescriptor(PdbxEntityDescriptor pdbxEntityDescriptor);

    void consumePdbxMolecule(PdbxMolecule pdbxMolecule);

    void consumePdbxMoleculeFeatures(PdbxMoleculeFeatures pdbxMoleculeFeatures);

    void consumePdbxNonpolyScheme(PdbxNonpolyScheme pdbxNonpolyScheme);

    void consumePdbxReferenceEntityLink(PdbxReferenceEntityLink pdbxReferenceEntityLink);

    void consumePdbxReferenceEntityList(PdbxReferenceEntityList pdbxReferenceEntityList);

    void consumePdbxReferenceEntityPolyLink(PdbxReferenceEntityPolyLink pdbxReferenceEntityPolyLink);

    void consumePdbxStructAssembly(PdbxStructAssembly pdbxStructAssembly);

    void consumePdbxStructAssemblyGen(PdbxStructAssemblyGen pdbxStructAssemblyGen);

    void consumePdbxStructModResidue(PdbxStructModResidue pdbxStructModResidue);

    void consumePdbxStructOperList(PdbxStructOperList pdbxStructOperList);

    void consumeRefine(Category refine);

    void consumeStruct(Struct struct);

    void consumeStructAsym(StructAsym structAsym);

    void consumeStructConf(StructConf structConf);

    void consumeStructConn(StructConn structConn);

    void consumeStructConnType(StructConnType structConnType);

    void consumeStructKeywords(StructKeywords structKeywords);

    void consumeStructNcsOper(StructNcsOper structNcsOper);

    void consumeStructRef(Category structRef);

    void consumeStructRefSeq(Category structRefSeq);

    void consumeStructRefSeqDif(Category structRefSeqDif);

    void consumeStructSheetRange(StructSheetRange structSheetRange);

    void consumeStructSite(StructSite structSite);

    void consumeStructSiteGen(StructSiteGen structSiteGen);

    void consumeSymmetry(Symmetry symmetry);

    void finish();

    S getContainer();
}
