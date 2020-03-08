package org.biojava.nbio.structure;

import org.biojava.nbio.structure.io.MMCIFFileReader;
import org.biojava.nbio.structure.io.PDBFileReader;

import java.util.Collections;
import java.util.List;

public enum StructureFiletype {
    PDB( (new PDBFileReader()).getExtensions()),
    CIF( new MMCIFFileReader().getExtensions()),
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
