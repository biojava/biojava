package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.rcsb.cif.model.Block;
import org.rcsb.cif.model.CifFile;

/**
 * Convert BioJava structures to CifFiles and vice versa.
 */
public class CifFileConverter {
    /**
     * @see CifFileConverter#convert(CifFile, FileParsingParameters)
     */
    public static Structure convert(CifFile cifFile) {
        return convert(cifFile, new FileParsingParameters());
    }

    /**
     * Convert CifFile to Structure.
     * @param cifFile the source
     * @param parameters parameters for parsing
     * @return the target
     */
    public static Structure convert(CifFile cifFile, FileParsingParameters parameters) {
        // initialize consumer
        CifFileConsumer<Structure> consumer = new CifFileConsumerImpl(parameters);

        // init structure
        consumer.prepare();

        // feed individual categories to consumer
        Block cifBlock = cifFile.getBlocks().get(0);

        // TODO maybe integrate rogue categories into Cif schema
        consumer.consumeAuditAuthor(cifBlock.getCategory("audit_author"));
        consumer.consumeAtomSite(cifBlock.getAtomSite());
        consumer.consumeAtomSites(cifBlock.getAtomSites());
        consumer.consumeCell(cifBlock.getCell());
        consumer.consumeChemComp(cifBlock.getChemComp());
        consumer.consumeChemCompBond(cifBlock.getChemCompBond());
        consumer.consumeDatabasePDBremark(cifBlock.getCategory("database_PDB_remark"));
        consumer.consumeDatabasePDBrev(cifBlock.getCategory("database_PDB_rev"));
        consumer.consumeDatabasePDBrevRecord(cifBlock.getCategory("database_PDB_rev_record"));
        consumer.consumeEntity(cifBlock.getEntity());
        consumer.consumeEntityPoly(cifBlock.getEntityPoly());
        consumer.consumeEntitySrcGen(cifBlock.getCategory("entity_src_gen"));
        consumer.consumeEntitySrcNat(cifBlock.getCategory("entity_src_nat"));
        consumer.consumeEntitySrcSyn(cifBlock.getCategory("entity_src_syn"));
        consumer.consumeEntityPolySeq(cifBlock.getEntityPolySeq());
        consumer.consumeExptl(cifBlock.getExptl());
        consumer.consumePdbxAuditRevisionHistory(cifBlock.getCategory("pdbx_audit_revision_history"));
        consumer.consumePdbxChemCompIdentifier(cifBlock.getPdbxChemCompIdentifier());
        consumer.consumePdbxDatabaseStatus(cifBlock.getCategory("pdbx_database_status"));
        consumer.consumePdbxEntityDescriptor(cifBlock.getPdbxEntityDescriptor());
        consumer.consumePdbxMolecule(cifBlock.getPdbxMolecule());
        consumer.consumePdbxMoleculeFeatures(cifBlock.getPdbxMoleculeFeatures());
        consumer.consumePdbxNonpolyScheme(cifBlock.getPdbxNonpolyScheme());
        consumer.consumePdbxReferenceEntityLink(cifBlock.getPdbxReferenceEntityLink());
        consumer.consumePdbxReferenceEntityList(cifBlock.getPdbxReferenceEntityList());
        consumer.consumePdbxReferenceEntityPolyLink(cifBlock.getPdbxReferenceEntityPolyLink());
        consumer.consumePdbxStructAssembly(cifBlock.getPdbxStructAssembly());
        consumer.consumePdbxStructAssemblyGen(cifBlock.getPdbxStructAssemblyGen());
        consumer.consumePdbxStructModResidue(cifBlock.getPdbxStructModResidue());
        consumer.consumePdbxStructOperList(cifBlock.getPdbxStructOperList());
        consumer.consumeRefine(cifBlock.getCategory("refine"));
        consumer.consumeStruct(cifBlock.getStruct());
        consumer.consumeStructAsym(cifBlock.getStructAsym());
        consumer.consumeStructConf(cifBlock.getStructConf());
        consumer.consumeStructConn(cifBlock.getStructConn());
        consumer.consumeStructConnType(cifBlock.getStructConnType());
        consumer.consumeStructKeywords(cifBlock.getStructKeywords());
        consumer.consumeStructNcsOper(cifBlock.getStructNcsOper());
        consumer.consumeStructRef(cifBlock.getCategory("struct_ref"));
        consumer.consumeStructRefSeq(cifBlock.getCategory("struct_ref_seq"));
        consumer.consumeStructRefSeqDif(cifBlock.getCategory("struct_ref_seq_dif"));
        consumer.consumeStructSheetRange(cifBlock.getStructSheetRange());
        consumer.consumeStructSite(cifBlock.getStructSite());
        consumer.consumeStructSiteGen(cifBlock.getStructSiteGen());
        consumer.consumeSymmetry(cifBlock.getSymmetry());

        // prepare structure to be retrieved
        consumer.finish();

        return consumer.getContainer();
    }

    /**
     * Convert Structure to CifFile.
     * @param structure the source
     * @return the target
     */
    public static CifFile convert(Structure structure) {
        return new CifFileSupplierImpl().get(structure);
    }
}
