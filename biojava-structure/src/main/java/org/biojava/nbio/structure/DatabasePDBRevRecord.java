package org.biojava.nbio.structure;

import org.biojava.nbio.structure.io.cif.CifBean;

/**
 * Represents revision records for use by {@link PDBHeader}.
 * @author Sebastian Bittrich
 * @since 6.0.0
 */
public class DatabasePDBRevRecord implements CifBean {
    private static final long serialVersionUID = 1L;
    private String revNum;
    private String type;
    private String details;

    public DatabasePDBRevRecord() {

    }

    public DatabasePDBRevRecord(String revNum, String type, String details) {
        this.revNum = revNum;
        this.type = type;
        this.details = details;
    }

    public DatabasePDBRevRecord(org.rcsb.cif.schema.mm.DatabasePDBRevRecord cif, int row) {
        this(cif.getDetails().get(row),
                cif.getRevNum().getStringData(row),
                cif.getType().get(row));
    }

    public String getRevNum() {
        return revNum;
    }

    public void setRevNum(String revNum) {
        this.revNum = revNum;
    }

    public String getType() {
        return type;
    }

    public void setType(String type) {
        this.type = type;
    }

    public String getDetails() {
        return details;
    }

    public void setDetails(String details) {
        this.details = details;
    }

    @Override
    public String toString() {
        return "DatabasePDBRevRecord{" +
                "revNum='" + revNum + '\'' +
                ", type='" + type + '\'' +
                ", details='" + details + '\'' +
                '}';
    }
}
