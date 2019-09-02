package org.biojava.nbio.structure.io;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.io.cif.CifFileConverter;

import java.io.IOException;
import java.io.InputStream;

/**
 * Parse text Cif files and provide capabilities to store them locally.
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @since 5.3.0
 */
public class CifFileReader extends LocalPDBDirectory {
    public static final String[] CIF_SPLIT_DIR    = new String[] { "data", "structures", "divided", "mmCIF" };
    public static final String[] CIF_OBSOLETE_DIR = new String[] { "data", "structures", "obsolete", "mmCIF" };

    /**
     * Constructs a new CifFileReader, initializing the extensions member variable.
     * The path is initialized in the same way as {@link UserConfiguration},
     * i.e. to system property/environment variable {@link UserConfiguration#PDB_DIR}.
     * Both autoFetch and splitDir are initialized to false
     */
    public CifFileReader() {
        this(null);
    }

    /**
     * Constructs a new CifFileReader, initializing the extensions member variable.
     * The path is initialized to the given path, both autoFetch and splitDir are initialized to false.
     */
    public CifFileReader(String path) {
        super(path);
        addExtension(".cif");
        addExtension(".mmcif");
        addExtension(".cif.gz");
        addExtension(".mmcif.gz");
    }

    @Override
    public Structure getStructure(InputStream inStream) throws IOException{
        return CifFileConverter.fromInputStream(inStream, getFileParsingParameters());
    }

    @Override
    protected String getFilename(String pdbId) {
        return pdbId.toLowerCase() + ".cif.gz";
    }

    @Override
    protected String[] getSplitDirPath() {
        return CIF_SPLIT_DIR;
    }

    @Override
    protected String[] getObsoleteDirPath() {
        return CIF_OBSOLETE_DIR;
    }
}
