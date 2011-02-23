package org.biojava.bio.structure.io.mmcif.model;


/** Container object for _pdbx_chem_comp_descriptor
 * 
 * @author Andreas Prlic
 * @since 3.2
 *
 */
public class ChemCompDescriptor {
	String comp_id; 
	String type; 
	String program;
	String program_version;
	String descriptor;
	public String getComp_id() {
		return comp_id;
	}
	public void setComp_id(String comp_id) {
		this.comp_id = comp_id;
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
	public String getProgram_version() {
		return program_version;
	}
	public void setProgram_version(String program_version) {
		this.program_version = program_version;
	}
	public String getDescriptor() {
		return descriptor;
	}
	public void setDescriptor(String descriptor) {
		this.descriptor = descriptor;
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((comp_id == null) ? 0 : comp_id.hashCode());
		result = prime * result
				+ ((descriptor == null) ? 0 : descriptor.hashCode());
		result = prime * result + ((program == null) ? 0 : program.hashCode());
		result = prime * result
				+ ((program_version == null) ? 0 : program_version.hashCode());
		result = prime * result + ((type == null) ? 0 : type.hashCode());
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ChemCompDescriptor other = (ChemCompDescriptor) obj;
		if (comp_id == null) {
			if (other.comp_id != null)
				return false;
		} else if (!comp_id.equals(other.comp_id))
			return false;
		if (descriptor == null) {
			if (other.descriptor != null)
				return false;
		} else if (!descriptor.equals(other.descriptor))
			return false;
		if (program == null) {
			if (other.program != null)
				return false;
		} else if (!program.equals(other.program))
			return false;
		if (program_version == null) {
			if (other.program_version != null)
				return false;
		} else if (!program_version.equals(other.program_version))
			return false;
		if (type == null) {
			if (other.type != null)
				return false;
		} else if (!type.equals(other.type))
			return false;
		return true;
	}
	@Override
	public String toString() {
		return "ChemCompDescriptior [comp_id=" + comp_id + ", type=" + type
				+ ", program=" + program + ", program_version="
				+ program_version + ", descriptor=" + descriptor + "]";
	}
	 
	
	
}
