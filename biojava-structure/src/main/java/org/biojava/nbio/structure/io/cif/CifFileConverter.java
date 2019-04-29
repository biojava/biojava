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
        Block cifBlock = cifFile.getFirstBlock();

        consumer.consumeAuditAuthor(cifBlock.getAuditAuthor());
        consumer.consumeAtomSite(cifBlock.getAtomSite());
        consumer.consumeAtomSites(cifBlock.getAtomSites());
        consumer.consumeCell(cifBlock.getCell());
        consumer.consumeChemComp(cifBlock.getChemComp());
        consumer.consumeChemCompBond(cifBlock.getChemCompBond());
        consumer.consumeDatabasePDBremark(cifBlock.getDatabasePDBRemark());
        consumer.consumeDatabasePDBrev(cifBlock.getDatabasePDBRev());
        consumer.consumeDatabasePDBrevRecord(cifBlock.getDatabasePDBRevRecord());
        consumer.consumeEntity(cifBlock.getEntity());
        consumer.consumeEntityPoly(cifBlock.getEntityPoly());
        consumer.consumeEntitySrcGen(cifBlock.getEntitySrcGen());
        consumer.consumeEntitySrcNat(cifBlock.getEntitySrcNat());
        consumer.consumeEntitySrcSyn(cifBlock.getPdbxEntitySrcSyn());
        consumer.consumeEntityPolySeq(cifBlock.getEntityPolySeq());
        consumer.consumeExptl(cifBlock.getExptl());
        consumer.consumePdbxAuditRevisionHistory(cifBlock.getPdbxAuditRevisionHistory());
        consumer.consumePdbxChemCompIdentifier(cifBlock.getPdbxChemCompIdentifier());
        consumer.consumePdbxDatabaseStatus(cifBlock.getPdbxDatabaseStatus());
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
        consumer.consumeRefine(cifBlock.getRefine());
        consumer.consumeStruct(cifBlock.getStruct());
        consumer.consumeStructAsym(cifBlock.getStructAsym());
        consumer.consumeStructConf(cifBlock.getStructConf());
        consumer.consumeStructConn(cifBlock.getStructConn());
        consumer.consumeStructConnType(cifBlock.getStructConnType());
        consumer.consumeStructKeywords(cifBlock.getStructKeywords());
        consumer.consumeStructNcsOper(cifBlock.getStructNcsOper());
        consumer.consumeStructRef(cifBlock.getStructRef());
        consumer.consumeStructRefSeq(cifBlock.getStructRefSeq());
        consumer.consumeStructRefSeqDif(cifBlock.getStructRefSeqDif());
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
