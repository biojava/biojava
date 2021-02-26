package org.biojava.nbio.structure.io;

import java.util.Collections;
import java.util.List;

/**
 * An enum of supported file formats.
 * @author Sebastian Bittrich
 * @since 6.0.0
 */
public enum StructureFiletype {
    PDB(new PDBFileReader().getExtensions()),
    CIF(new CifFileReader().getExtensions()),
    BCIF(new BcifFileReader().getExtensions()),
    MMTF(new MMTFFileReader().getExtensions()),
    UNKNOWN(Collections.emptyList());

    private final List<String> extensions;

    /**
     * @param extensions List of supported extensions, including leading period
     */
    StructureFiletype(List<String> extensions) {
        this.extensions = extensions;
    }

    /**
     * @return a list of file extensions associated with this type
     */
    public List<String> getExtensions() {
        return extensions;
    }
}
