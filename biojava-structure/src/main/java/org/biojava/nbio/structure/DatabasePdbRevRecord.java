package org.biojava.nbio.structure;

import java.io.Serializable;

public class DatabasePdbRevRecord implements Serializable {
    private static final long serialVersionUID = -791924804009516791L;
    private String rev_num;
    private String type;
    private String details;

    public String getRev_num() {
        return rev_num;
    }

    public void setRev_num(String rev_num) {
        this.rev_num = rev_num;
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
        return "DatabasePdbrevRecord{" +
                "rev_num='" + rev_num + '\'' +
                ", type='" + type + '\'' +
                ", details='" + details + '\'' +
                '}';
    }
}
