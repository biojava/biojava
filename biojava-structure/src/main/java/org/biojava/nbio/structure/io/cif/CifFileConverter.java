package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.rcsb.cif.CifIO;
import org.rcsb.cif.model.Block;
import org.rcsb.cif.model.CifFile;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Convert BioJava structures to CifFiles and vice versa.
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @since 5.3.0
 */
public class CifFileConverter {
    /**
     * Read data from a file and convert to Structure without any FileParsingParameters.
     * @param path the source of information - can be gzipped or binary or text data
     * @return the target
     */
    public static Structure fromPath(Path path) throws IOException {
        return fromInputStream(Files.newInputStream(path), new FileParsingParameters());
    }

    /**
     * Read data from a file and convert to Structure.
     * @param path the source of information - can be gzipped or binary or text data
     * @param parameters parameters for parsing
     * @return the target
     */
    public static Structure fromPath(Path path, FileParsingParameters parameters) throws IOException {
        return fromInputStream(Files.newInputStream(path), parameters);
    }

    /**
     * Get data from a URL and convert to Structure without any FileParsingParameters.
     * @param url the source of information - can be gzipped or binary or text data
     * @return the target
     * @throws IOException thrown when reading fails
     */
    public static Structure fromURL(URL url) throws IOException {
        return fromURL(url, new FileParsingParameters());
    }

    /**
     * Get data from a URL and convert to Structure.
     * @param url the source of information - can be gzipped or binary or text data
     * @param parameters parameters for parsing
     * @return the target
     * @throws IOException thrown when reading fails
     */
    private static Structure fromURL(URL url, FileParsingParameters parameters) throws IOException {
        return fromInputStream(url.openStream(), parameters);
    }

    /**
     * Convert InputStream to Structure without any FileParsingParameters.
     * @param inputStream the InputStream of information - can be gzipped or binary or text data
     * @return the target
     * @throws IOException thrown when reading fails
     * @see CifFileConverter#fromInputStream(InputStream, FileParsingParameters)
     */
    public static Structure fromInputStream(InputStream inputStream) throws IOException {
        return fromInputStream(inputStream, new FileParsingParameters());
    }

    /**
     * Convert InputStream to Structure.
     * @param inputStream the InputStream of information - can be gzipped or binary or text data
     * @param parameters parameters for parsing
     * @return the target
     * @throws IOException thrown when reading fails
     */
    public static Structure fromInputStream(InputStream inputStream, FileParsingParameters parameters) throws IOException {
        return fromCifFile(CifIO.readFromInputStream(inputStream), parameters);
    }

    /**
     * Convert CifFile to Structure without any FileParsingParameters.
     * @param cifFile the source
     * @return the target
     * @see CifFileConverter#fromCifFile(CifFile, FileParsingParameters)
     */
    public static Structure fromCifFile(CifFile cifFile) {
        return fromCifFile(cifFile, new FileParsingParameters());
    }

    /**
     * Convert CifFile to Structure.
     * @param cifFile the source
     * @param parameters parameters for parsing
     * @return the target
     */
    public static Structure fromCifFile(CifFile cifFile, FileParsingParameters parameters) {
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
     * Write a structure to a CIF file.
     * @param structure the source
     * @param path where to write to
     * @throws IOException thrown when writing fails
     */
    public static void toTextFile(Structure structure, Path path) throws IOException {
        CifIO.writeText(toCifFile(structure), path);
    }

    /**
     * Write a structure to a BCIF file.
     * @param structure the source
     * @param path where to write to
     * @throws IOException thrown when writing fails
     */
    public static void toBinaryFile(Structure structure, Path path) throws IOException {
        CifIO.writeBinary(toCifFile(structure), path);
    }

    /**
     * Convert a structure to BCIF format.
     * @param structure the source
     * @return the binary representation of the structure
     * @throws IOException thrown when writing fails
     */
    public static byte[] toBinary(Structure structure) throws IOException {
        return CifIO.writeText(toCifFile(structure));
    }

    /**
     * Convert a structure to mmCIF format.
     * @param structure the source
     * @return the mmCIF String representation of the structure
     * @throws IOException thrown when writing fails
     */
    public static String toText(Structure structure) throws IOException {
        return new String(CifIO.writeText(toCifFile(structure)));
    }

    /**
     * Convert Structure to CifFile.
     * @param structure the source
     * @return the target
     */
    public static CifFile toCifFile(Structure structure) {
        return new CifFileSupplierImpl().get(structure);
    }
}
