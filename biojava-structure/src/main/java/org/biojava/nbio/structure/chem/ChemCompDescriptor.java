package org.biojava.nbio.structure.chem;

import org.biojava.nbio.structure.io.cif.CifBean;

import java.util.Objects;

/**
 * Properties of the chemical component descriptor.
 * @author Sebastian Bittrich
 * @since 6.0.0
 */
public class ChemCompDescriptor implements CifBean {
    private static final long serialVersionUID = 1078685833800736278L;
    private String compId;
    private String type;
    private String program;
    private String programVersion;
    private String descriptor;

    public String getCompId() {
        return compId;
    }

    public void setCompId(String compId) {
        this.compId = compId;
    }

    public String getType() {
        return type;
    }

    public void setType(String type) {
        this.type = type;
    }

    public String getProgram() {
        return program;
    }

    public void setProgram(String program) {
        this.program = program;
    }

    public String getProgramVersion() {
        return programVersion;
    }

    public void setProgramVersion(String programVersion) {
        this.programVersion = programVersion;
    }

    public String getDescriptor() {
        return descriptor;
    }

    public void setDescriptor(String descriptor) {
        this.descriptor = descriptor;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        ChemCompDescriptor that = (ChemCompDescriptor) o;
        return Objects.equals(compId, that.compId) &&
                Objects.equals(type, that.type) &&
                Objects.equals(program, that.program) &&
                Objects.equals(programVersion, that.programVersion) &&
                Objects.equals(descriptor, that.descriptor);
    }

    @Override
    public int hashCode() {
        return Objects.hash(compId, type, program, programVersion, descriptor);
    }

    @Override
    public String toString() {
        return "ChemCompDescriptor [comp_id=" + compId +
                ", type=" + type +
                ", program=" + program +
                ", program_version=" + programVersion +
                ", descriptor=" + descriptor + "]";
    }
}
