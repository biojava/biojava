package org.biojava.nbio.structure.io;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.cif.CifFileConsumer;
import org.biojava.nbio.structure.io.cif.CifFileConsumerImpl;
import org.rcsb.cif.CifReader;
import org.rcsb.cif.model.Block;
import org.rcsb.cif.model.CifFile;

import java.io.IOException;
import java.io.InputStream;

public class CIFFileReader extends LocalPDBDirectory {
    // TODO paths are meaningless for now
    public static final String[] BCIF_SPLIT_DIR = new String[] { "data", "structures", "divided", "bcif" };
    public static final String[] BCIF_OBSOLETE_DIR = new String[] { "data", "structures", "obsolete", "bcif" };

    public CIFFileReader() {
        this(null);
    }

    public CIFFileReader(String path) {
        super(path);
        // TODO
        addExtension(".cif");
        addExtension(".cif.gz");
        addExtension(".bcif");
        addExtension(".bcif.gz");
    }

    @Override
    public Structure getStructure(InputStream inStream) throws IOException {
        return read(inStream);
    }

    @Override
    protected String getFilename(String pdbId) {
        // TODO
        return pdbId.toLowerCase() + ".bcif.gz";
    }

    @Override
    protected String[] getSplitDirPath() {
        // TODO
        return BCIF_SPLIT_DIR;
    }

    @Override
    protected String[] getObsoleteDirPath() {
        // TODO
        return BCIF_OBSOLETE_DIR;
    }

    private Structure read(InputStream inputStream) throws IOException {
        // TODO allow to switch dynamically between text and binary?
        // obtain the decoded, flat representation of the input
        return read(CifReader.parseBinary(inputStream));
    }

    private Structure read(CifFile cifFile) {
        // initialize consumer
        CifFileConsumer<Structure> consumer = new CifFileConsumerImpl(getFileParsingParameters());

        // init structure
        consumer.prepare();

        // feed individual categories to consumer
        Block cifBlock = cifFile.getBlocks().get(0);

        // TODO maybe integrate rogue categories into Cif schema
        // TODO rethink way to access missing categories
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
}
