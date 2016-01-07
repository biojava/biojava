package org.biojava.nbio.structure.io.mmcif.model;

/**
 * Created by Matt on 11/1/2015.
 */
public class StructSite {
    String id;
    String details;
    String pdbx_evidence_code;

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getDetails() {
        return details;
    }

    public void setDetails(String details) {
        this.details = details;
    }

    public String getPdbx_evidence_code() {
        return pdbx_evidence_code;
    }

    public void setPdbx_evidence_code(String pdbx_evidence_code) {
        this.pdbx_evidence_code = pdbx_evidence_code;
    }
}
